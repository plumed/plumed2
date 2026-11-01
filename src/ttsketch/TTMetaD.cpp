/* Copyright (c) 2026, Nils E. Strand, Siyao Yang, Yuehaw Khoo, and Aaron R. Dinner
 *
 * Tensor train (TT)-metadynamics.
 * https://arxiv.org/abs/2603.13549
 * See module.md for installation instructions.
 */

#include "TTHelper.h"
#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "tools/File.h"
#include "tools/OpenMP.h"
#include "tools/OFile.h"
#include "tools/IFile.h"
#include <algorithm>
#include <ctime>
#include <filesystem>
#include <random>

namespace PLMD {
namespace ttsketch {

//+PLUMEDOC BIAS TTMETAD
/*
Performs metadynamics with a tensor-train representation of the bias potential.

This action runs metadynamics while periodically compressing the accumulated
Gaussian hills into a low-rank tensor train (TT) using the TT-Sketch algorithm.
The TT representation keeps memory and evaluation cost bounded throughout long
simulations, in contrast to standard metadynamics where cost grows linearly with
the number of deposited hills.

The basis functions of each argument are defined on [INTERVAL_MIN, INTERVAL_MAX] and are
periodic over this interval. For a periodic argument the interval should thus coincide with the
periodicity domain. For a non periodic argument the interval should contain the region that is
explored: outside it the bias is the one at the closest limit and there is no biasing force along
that argument, and hills deposited outside are folded back into the interval when they are
compressed into the TT (a warning is written in the log).

After every compression the TT is appended to the HDF5 file given by TT_FILE, and the hills file
restarts empty. When restarting, the last TT in TT_FILE is read back together with the hills in FILE,
i.e. those deposited after it. With multiple replicas each replica has its own TT_FILE (with the
replica suffix), unless WALKERS_MPI is used, in which case the walkers share one TT_FILE next to
the shared hills file. Two TTMETAD actions must use different TT_FILE.

## Examples

The following input runs TT-MetaD on two distances. Hills are deposited every
500 steps and compressed into a TT every 1e6 steps with target rank 10.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
TTMETAD ...
  ARG=d1,d2
  SIGMA=0.1,0.1
  HEIGHT=1.0
  PACE=500
  SKETCH_PACE=1000000
  SKETCH_RANK=10
  SKETCH_INITRANK=50
  INTERVAL_MIN=0.0,0.0
  INTERVAL_MAX=3.0,3.0
  FILE=HILLS
  LABEL=tt
... TTMETAD
PRINT ARG=tt.bias FILE=COLVAR STRIDE=100
```

*/
//+ENDPLUMEDOC

using namespace PLMD::bias;

#if !defined(__PLUMED_HAS_LIBITENSOR) || !defined(__PLUMED_HAS_LIBHDF5)

class TTMetaD : public Bias {
public:
  static void registerKeywords(Keywords& keys);
  explicit TTMetaD(const ActionOptions& ao) : Action(ao), Bias(ao) {
    error("TTMETAD requires libitensor and libhdf5. "
          "Reconfigure with --enable-libitensor.");
  }
  void calculate() override {}
  void update() override {}
};

#else

using namespace itensor;

class TTMetaD : public Bias {

private:
  struct Gaussian {
    double height;
    std::vector<double> center;
    std::vector<double> sigma;
    std::vector<double> invsigma;
    Gaussian(const double h, const std::vector<double>& c, const std::vector<double>& s)
      : height(h), center(c), sigma(s), invsigma(s) {
      for(unsigned i = 0; i < invsigma.size(); ++i) {
        if(abs(invsigma[i]) > 1.e-20) {
          invsigma[i] = 1.0 / invsigma[i];
        } else {
          invsigma[i] = 0.0;
        }
      }
    }
  };
  // standard metadynamics parameters
  double kbt_;
  int stride_;
  bool welltemp_;
  double biasf_;
  std::string fmt_;
  bool isFirstStep_;
  double height0_;
  std::vector<double> sigma0_;
  std::vector<Gaussian> hills_;   // Gaussians deposited since the last TT update
  OFile hillsOfile_;
  std::string mw_dir_;
  std::string hillsfname_;
  bool walkers_mpi_;
  int mpi_size_;
  int mpi_rank_;
  // TT-sketch parameters
  unsigned d_;               // number of CV dimensions
  int sketch_rc_;            // initial (sketch) rank for the random TT coefficient tensor
  int sketch_r_;             // target TT bond dimension after SVD truncation (0 = use cutoff)
  double sketch_cutoff_;     // relative SVD truncation threshold (used when sketch_r_=0)
  int sketch_stride_;        // steps between TT bias updates (tau in the paper)
  double sketch_alpha_;      // weight for non-constant basis functions in the sketch coefficient TT
  std::vector<BasisFunc> sketch_basis_;  // one BasisFunc per CV dimension
  unsigned sketch_count_;    // number of TT updates performed so far (starts at 1, increments after each update)
  MPS vb_;                   // current TT approximation of the accumulated bias (V_bias^TT)
  DenseTT vbDense_;          // plain-array copy of vb_ used to evaluate the bias (rebuilt by setVb())
  std::string ttfilename_;   // HDF5 file where the TTs are stored (TT_FILE, with replica suffix or in the walkers directory)
  std::vector<double> imin_; // INTERVAL_MIN
  std::vector<double> imax_; // INTERVAL_MAX
  double sketch_until_;      // simulation time after which the bias is frozen (no further updates)
  bool frozen_;              // true once the bias has been frozen (after sketch_until_)
  bool sketch_conv_;         // true if Gaussian smoothing (conv mode) is active for basis evaluation
  bool nonintrusive_;        // if true, uses the intrusive-free variant: accumulates Bemp directly
  // instead of re-projecting the previous vb_ through the Gram matrix
  bool deterministic_;       // if true, use a fixed RNG seed (42) in createTTCoeff() for reproducibility
  MPS B_prev_;               // accumulated tensor moment from all previous sketches (nonintrusive only)
  std::vector<ITensor> A_prev_;   // accumulated cross-gram matrices per bond (nonintrusive only)

  void readGaussians(IFile *ifile);
  void setVb(const MPS& vb);
  // clamp cv to [INTERVAL_MIN, INTERVAL_MAX]; out[i] is true if cv[i] was outside
  void clampToInterval(std::vector<double>& cv, std::vector<char>& out) const;
  bool writerReplica() const {
    return !this->walkers_mpi_ || this->mpi_rank_ == 0;
  }
  void writeGaussian(const Gaussian& hill, OFile& file);
  double getHeight(const std::vector<double>& cv);
  double getBias(const std::vector<double>& cv);
  double getBiasAndDerivatives(const std::vector<double>& cv, std::vector<double>& der);
  double evaluateGaussian(const std::vector<double>& cv, const Gaussian& hill);
  double evaluateGaussianAndDerivatives(const std::vector<double>& cv, const Gaussian& hill, std::vector<double>& der, std::vector<double>& dp);
  bool scanOneHill(IFile* ifile, std::vector<Value>& tmpvalues, std::vector<double>& center, std::vector<double>& sigma, double& height);
  void paraSketch();
  MPS createTTCoeff() const;
  std::vector<Matrix<double>> intBasisSample() const;
  std::pair<MPS, std::vector<Matrix<double>>> formTensorMoment(const std::vector<Matrix<double>>& M, const MPS& coeff) const;
  std::tuple<MPS, std::vector<ITensor>, std::vector<ITensor>> formTensorMomentVb(const MPS& coeff);

public:
  explicit TTMetaD(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
};

TTMetaD::TTMetaD(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao),
  kbt_(0.0),
  stride_(0),
  welltemp_(false),
  biasf_(-1.0),
  isFirstStep_(true),
  height0_(std::numeric_limits<double>::max()),
  mw_dir_(""),
  walkers_mpi_(false),
  mpi_size_(0),
  mpi_rank_(0),
  sketch_r_(0),
  sketch_cutoff_(0.0),
  sketch_count_(1),
  sketch_until_(std::numeric_limits<double>::max()),
  frozen_(false),
  sketch_conv_(false),
  nonintrusive_(false),
  deterministic_(false) {
  bool kernel;
  parseFlag("KERNEL_BASIS", kernel);
  this->d_ = getNumberOfArguments();
  if(this->d_ < 2) {
    error("Number of arguments must be at least 2");
  }
  parse("FMT", this->fmt_);
  parseVector("SIGMA", this->sigma0_);
  if(this->sigma0_.size() != d_) {
    error("number of arguments does not match number of SIGMA parameters");
  }
  parse("HEIGHT", this->height0_);
  parse("PACE", this->stride_);
  if(stride_ <= 0) {
    error("frequency for hill addition is nonsensical");
  }
  this->hillsfname_ = "HILLS";
  parse("FILE", this->hillsfname_);
  parse("BIASFACTOR", this->biasf_);
  if(this->biasf_ < 1.0 && this->biasf_ != -1.0) {
    error("well tempered bias factor is nonsensical");
  }
  this->kbt_ = getkBT();
  if(this->biasf_ >= 1.0) {
    if(this->kbt_ == 0.0) {
      error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    }
    this->welltemp_ = true;
  }

  parseFlag("WALKERS_MPI", this->walkers_mpi_);
  parse("WALKERS_DIR", this->mw_dir_);
  if(this->walkers_mpi_ && this->mw_dir_== "") {
    const std::string ret = std::filesystem::current_path();
    this->mw_dir_ = ret + "/";
    // multi_sim_comm is valid only on the first rank of each replica
    if(comm.Get_rank() == 0) {
      multi_sim_comm.Bcast(this->mw_dir_, 0);
    }
    comm.Bcast(this->mw_dir_, 0);
  }
  if(this->walkers_mpi_ && this->mw_dir_ != "") {
    this->hillsfname_ = this->mw_dir_ + "/" + this->hillsfname_;
  }
  std::string ttfile = "ttsketch.h5";
  parse("TT_FILE", ttfile);
  if(this->walkers_mpi_) {
    // shared by all the walkers, next to the shared hills file
    this->ttfilename_ = this->mw_dir_ + "/" + ttfile;
  } else {
    // one file per replica, as for the other output files
    this->ttfilename_ = FileBase::appendSuffix(ttfile, plumed.getSuffix());
  }
  for(const auto* other : plumed.getActionSet().select<const TTMetaD*>()) {
    if(other->ttfilename_ == this->ttfilename_) {
      error("TT_FILE " + this->ttfilename_ + " is already used by action " + other->getLabel() + ": use a different TT_FILE");
    }
  }

  parse("SKETCH_RANK", this->sketch_r_);
  parse("SKETCH_CUTOFF", this->sketch_cutoff_);
  if(this->sketch_r_ <= 0 && (this->sketch_cutoff_ <= 0.0 || this->sketch_cutoff_ > 1.0)) {
    error("Valid SKETCH_RANK or SKETCH_CUTOFF needs to be specified");
  }
  parse("SKETCH_INITRANK", this->sketch_rc_);
  if(this->sketch_rc_ <= 0) {
    error("SKETCH_INITRANK must be positive");
  }
  parse("SKETCH_PACE", this->sketch_stride_);
  if(this->sketch_stride_ <= 0) {
    error("SKETCH_PACE must be positive");
  }
  std::vector<double> interval_min;
  parseVector("INTERVAL_MIN", interval_min);
  if(interval_min.size() != this->d_) {
    error("Number of arguments does not match number of INTERVAL_MIN parameters");
  }
  std::vector<double> interval_max;
  parseVector("INTERVAL_MAX", interval_max);
  if(interval_max.size() != this->d_) {
    error("Number of arguments does not match number of INTERVAL_MAX parameters");
  }
  int nbasis = 20;
  parse("SKETCH_NBASIS", nbasis);
  if(nbasis <= 1) {
    error("SKETCH_NBASIS must be greater than 1");
  }
  if(!kernel && nbasis % 2 == 0) {
    ++nbasis;
  }
  parse("SKETCH_ALPHA", this->sketch_alpha_);
  if(this->sketch_alpha_ <= 0.0 || this->sketch_alpha_ > 1.0) {
    error("SKETCH_ALPHA must be positive and no greater than 1");
  }
  std::vector<double> w;
  parseVector("SKETCH_WIDTH", w);
  if(w.size() == 0) {
    w.assign(this->d_, 0.0);
  }
  if(w.size() != this->d_) {
    error("Number of arguments does not match number of SKETCH_WIDTH parameters");
  }
  for (double val : w) {
    if (val != 0.0) {
      this->sketch_conv_ = true;
    }
  }
  std::vector<double> dx;
  parseVector("KERNEL_DX", dx);
  if(dx.size() == 0) {
    dx.resize(this->d_, 0.0);
  }
  if(dx.size() != this->d_) {
    error("Number of arguments does not match number of KERNEL_DX parameters");
  }
  for(unsigned i = 0; i < this->d_; ++i) {
    if(this->sketch_conv_ && w[i] <= 0.0) {
      error("Gaussian smoothing requires positive WIDTH");
    }
    if(kernel && dx[i] < 0.0) {
      error("Kernel basis requires positive KERNEL_DX");
    }
    if(interval_max[i] <= interval_min[i]) {
      error("INTERVAL_MAX parameters need to be greater than respective INTERVAL_MIN parameters");
    }
    if(getPntrToArgument(i)->isPeriodic()) {
      // the basis is periodic over the interval: for a periodic CV the two periods should coincide
      double pmin, pmax;
      getPntrToArgument(i)->getDomain(pmin, pmax);
      const double tol = 1e-4 * (pmax - pmin);
      if(std::abs(interval_min[i] - pmin) > tol || std::abs(interval_max[i] - pmax) > tol) {
        warning("the interval of the periodic argument " + getPntrToArgument(i)->getName() +
                " does not coincide with its periodicity domain: the bias will be discontinuous");
      }
    }
    this->sketch_basis_.push_back(BasisFunc(std::make_pair(interval_min[i], interval_max[i]), nbasis, w[i], kernel, dx[i]));
  }
  if(kernel && this->sketch_conv_) {
    error("kernel smoothing incompatible with kernel basis");
  }
  this->imin_ = interval_min;
  this->imax_ = interval_max;
  if(this->walkers_mpi_) {
    this->mpi_size_ = getNumberOfReplicas();
    this->mpi_rank_ = getReplicaIndex();
  }

  parse("SKETCH_UNTIL", this->sketch_until_);

  parseFlag("NONINTRUSIVE", this->nonintrusive_);
  parseFlag("DETERMINISTIC", this->deterministic_);

  if(getRestart()) {
    // the last TT in the file is the bias; the hills file holds the hills deposited after it
    const unsigned ntt = ttCount(this->ttfilename_);
    if(ntt > 0) {
      this->sketch_count_ = ntt + 1;
      setVb(ttRead(this->ttfilename_, this->sketch_count_));
      if(length(this->vb_) != static_cast<int>(this->d_)) {
        error("the TT in " + this->ttfilename_ + " has " + std::to_string(length(this->vb_)) +
              " dimensions, but TTMETAD has " + std::to_string(this->d_) + " arguments");
      }
      for(unsigned i = 0; i < this->d_; ++i) {
        if(this->vbDense_.nbasis(i) != this->sketch_basis_[i].nbasis()) {
          error("the TT in " + this->ttfilename_ + " has " + std::to_string(this->vbDense_.nbasis(i)) +
                " basis functions along argument " + std::to_string(i + 1) + ", but SKETCH_NBASIS gives " +
                std::to_string(this->sketch_basis_[i].nbasis()));
        }
      }
      if(this->nonintrusive_ && !ttMomentsRead(this->ttfilename_, this->sketch_count_, this->d_ - 1, this->B_prev_, this->A_prev_)) {
        error("NONINTRUSIVE restart needs the tensor moments of the last sketch, which are not in " + this->ttfilename_);
      }
    }
    log << "  restarting from TT " << this->sketch_count_ - 1 << " in " << this->ttfilename_ << "\n";
    IFile hills_ifile;
    hills_ifile.link(*this);
    if(this->walkers_mpi_) {
      hills_ifile.enforceSuffix("");
    }
    if(!hills_ifile.FileExist(this->hillsfname_)) {
      error("The hills file cannot be found");
    }
    hills_ifile.open(this->hillsfname_);
    readGaussians(&hills_ifile);
    hills_ifile.close();
  } else if(writerReplica()) {
    // a fresh run: keep a TT file left by a previous run as a backup, as for the other output files,
    // so that it cannot be mistaken for the TTs of this run at a later restart
    OFile backup;
    backup.link(*this);
    backup.backupFile("bck", this->ttfilename_);
  }

  log << "  TTs are stored in " << this->ttfilename_ << "\n";
  log << "  Bibliography " << plumed.cite("Strand, Yang, Khoo, and Dinner, https://arxiv.org/abs/2603.13549 (2026)");
  log << plumed.cite("Fishman, White, and Stoudenmire, https://arxiv.org/abs/2007.14822 (2020)");
  log << "\n";

  if(writerReplica()) {
    this->hillsOfile_.link(*this);
    if(this->walkers_mpi_) {
      // a single hills file shared by all the walkers (written by the first one);
      // otherwise each replica has its own, with the replica suffix
      this->hillsOfile_.enforceSuffix("");
    }
    this->hillsOfile_.open(this->hillsfname_);
    if(this->fmt_.length() > 0) {
      this->hillsOfile_.fmtField(this->fmt_);
    }
    hillsOfile_.setHeavyFlush();
    for(unsigned i = 0; i < this->d_; ++i) {
      hillsOfile_.setupPrintValue(getPntrToArgument(i));
    }
  }
}

void TTMetaD::readGaussians(IFile *ifile) {
  std::vector<double> center(this->d_);
  std::vector<double> sigma(this->d_);
  double height;
  int nhills = 0;

  std::vector<Value> tmpvalues;
  for(unsigned j = 0; j < this->d_; ++j) {
    tmpvalues.push_back(Value(this, getPntrToArgument(j)->getName(), false));
  }

  while(scanOneHill(ifile, tmpvalues, center, sigma, height)) {
    ++nhills;
    if(this->welltemp_ && this->biasf_ > 1.0) {
      height *= (this->biasf_ - 1.0) / this->biasf_;
    }
    this->hills_.push_back(Gaussian(height, center, sigma));
  }
  log << "  " << nhills << " hills retrieved\n";
}

void TTMetaD::setVb(const MPS& vb) {
  this->vb_ = vb;
  this->vbDense_ = length(vb) == 0 ? DenseTT() : DenseTT(vb);
}

void TTMetaD::clampToInterval(std::vector<double>& cv, std::vector<char>& out) const {
  out.assign(this->d_, 0);
  for(unsigned i = 0; i < this->d_; ++i) {
    if(cv[i] < this->imin_[i]) {
      cv[i] = this->imin_[i];
      out[i] = 1;
    } else if(cv[i] > this->imax_[i]) {
      cv[i] = this->imax_[i];
      out[i] = 1;
    }
  }
}

bool TTMetaD::scanOneHill(IFile* ifile, std::vector<Value>& tmpvalues, std::vector<double>& center, std::vector<double>& sigma, double& height) {
  double dummy;
  if(ifile->scanField("time", dummy)) {
    unsigned ncv = tmpvalues.size();
    for(unsigned i = 0; i < ncv; ++i) {
      ifile->scanField(&tmpvalues[i]);
      if(tmpvalues[i].isPeriodic() && !getPntrToArgument(i)->isPeriodic()) {
        error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
      } else if(tmpvalues[i].isPeriodic()) {
        std::string imin, imax;
        tmpvalues[i].getDomain(imin, imax);
        std::string rmin, rmax;
        getPntrToArgument(i)->getDomain(rmin, rmax);
        if(imin != rmin || imax != rmax) {
          error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
        }
      }
      center[i] = tmpvalues[i].get();
    }
    for(unsigned i = 0; i < ncv; ++i) {
      ifile->scanField("sigma_" + getPntrToArgument(i)->getName(), sigma[i]);
    }

    ifile->scanField("height", height);
    ifile->scanField("biasf", dummy);
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}

void TTMetaD::writeGaussian(const Gaussian& hill, OFile&file) {
  file.printField("time", getTimeStep() * getStep());
  for(unsigned i = 0; i < this->d_; ++i) {
    file.printField(getPntrToArgument(i), hill.center[i]);
  }
  for(unsigned i = 0; i < this->d_; ++i) {
    file.printField("sigma_" + getPntrToArgument(i)->getName(), hill.sigma[i]);
  }
  double height = hill.height;
  if(this->welltemp_ && this->biasf_ > 1.0) {
    height *= this->biasf_ / (this->biasf_ - 1.0);
  }
  file.printField("height", height).printField("biasf", this->biasf_);
  file.printField();
}

// Called every MD step: evaluate total bias (TT + remaining Gaussians) and set forces.
void TTMetaD::calculate() {
  std::vector<double> cv(this->d_);
  for(unsigned i = 0; i < this->d_; ++i) {
    cv[i] = getArgument(i);
  }

  // outside the interval the bias is the one at the closest boundary and the force along that CV is zero
  std::vector<char> out;
  clampToInterval(cv, out);

  std::vector<double> der(this->d_, 0.0);

  double ene = getBiasAndDerivatives(cv, der);
  for(unsigned i = 0; i < this->d_; ++i) {
    if(out[i]) {
      der[i] = 0.0;
    }
  }
  setBias(ene);
  for(unsigned i = 0; i < this->d_; ++i) {
    setOutputForce(i, -der[i]);
  }
}

// Called every MD step after forces are applied.
// Two separate triggers:
//   1. Every sketch_stride_ steps: run paraSketch() to compress accumulated hills into vb_,
//      clear the hills list, and rewind the HILLS file for the next interval.
//   2. Every stride_ steps: deposit a new Gaussian hill and write it to the HILLS file.
void TTMetaD::update() {
  // SKETCH_UNTIL: at the first step after this time the pending hills are sketched
  // one last time and the bias is frozen
  const bool freezeNow = !this->frozen_ && !this->isFirstStep_ && getTime() >= this->sketch_until_;
  // all the walkers hold the same list of hills, so they all take the same decision here
  const bool nowAddATT = !this->frozen_ && !this->isFirstStep_ && !this->hills_.empty() &&
                         (getStep() % this->sketch_stride_ == 0 || freezeNow);

  if(nowAddATT) {
    if(writerReplica()) {
      this->hillsOfile_.flush();
      unsigned N = this->hills_.size();
      log << "Sample limits\n";
      for(unsigned i = 0; i < this->d_; ++i) {
        auto [large, small] = this->sketch_basis_[i].dom();
        unsigned nout = 0;
        for(unsigned j = 0; j < N; ++j) {
          const double c = this->hills_[j].center[i];
          if(c < this->imin_[i] || c > this->imax_[i]) {
            ++nout;
          }
          large = std::max(large, c);
          small = std::min(small, c);
        }
        log << small << " " << large << "\n";
        if(nout > 0) {
          log << "  WARNING: " << nout << " hills are centered outside the interval of argument " << getPntrToArgument(i)->getName()
              << ": the periodic basis folds them back into the interval\n";
        }
      }

      log << "\nEmpirical means:\n";
      std::vector<double> muhat(this->d_, 0.0);
      for(unsigned k = 0; k < this->d_; ++k) {
        for(unsigned j = 0; j < N; ++j) {
          muhat[k] += this->hills_[j].center[k] / N;
        }
        log << muhat[k] << " ";
      }
      if(N > 1) {
        log << "\nEmpirical covariance matrix:\n";
        Matrix<double> sigmahat(this->d_, this->d_);
        for(unsigned k = 0; k < this->d_; ++k) {
          for(unsigned l = k; l < this->d_; ++l) {
            double sum = 0.0;
            for(unsigned j = 0; j < N; ++j) {
              sum += (this->hills_[j].center[k] - muhat[k]) * (this->hills_[j].center[l] - muhat[l]);
            }
            sigmahat(k, l) = sigmahat(l, k) = sum / (N - 1);
          }
        }
        matrixOut(log, sigmahat);
      } else {
        log << "\n";
      }

      // record current bias at hill centers before update, for computing relative error after
      std::vector<std::vector<double>> x(N);
      std::vector<double> A0(N);
      for(unsigned i = 0; i < N; ++i) {
        x[i] = this->hills_[i].center;
        A0[i] = getBias(x[i]);
      }

      log << "\nStarting TT-sketch...\n";
      log.flush();
      paraSketch();
      ++this->sketch_count_;

      this->hills_.clear();

      // compute relative L2 approximation error at hill centers: ||V_new - V_old|| / ||V_old||
      std::vector<double> diff(N);
      for(unsigned i = 0; i < N; ++i) {
        diff[i] = getBias(x[i]) - A0[i];
      }
      log << "Relative l2 error = " << sqrt(norm(diff) / norm(A0)) << "\n\n";
      log.flush();

      // all ranks of the replica hold the same vb_; one of them writes it
      if(comm.Get_rank() == 0) {
        ttWrite(this->ttfilename_, this->vb_, this->sketch_count_);
        if(this->nonintrusive_) {
          ttMomentsWrite(this->ttfilename_, this->B_prev_, this->A_prev_, this->sketch_count_);
        }
      }
    }

    if(this->walkers_mpi_) {
      if(comm.Get_rank() == 0) {
        multi_sim_comm.Bcast(this->sketch_count_, 0);
      }
      comm.Bcast(this->sketch_count_, 0);
      if(this->mpi_rank_ != 0) {
        this->hills_.clear();
        setVb(ttRead(this->ttfilename_, this->sketch_count_));
      }
    }
    // the hills are now in the TT: the hills file restarts empty
    if(writerReplica()) {
      this->hillsOfile_.rewind();
      this->hillsOfile_.clearFields();
      if(this->fmt_.length() > 0) {
        this->hillsOfile_.fmtField(this->fmt_);
      }
      hillsOfile_.setHeavyFlush();
      for(unsigned i = 0; i < this->d_; ++i) {
        hillsOfile_.setupPrintValue(getPntrToArgument(i));
      }
    }
  }
  if(freezeNow) {
    this->frozen_ = true;
    log << "TT bias frozen at time " << getTime() << "\n";
  }

  bool nowAddAHill;
  if(getStep() % this->stride_ == 0 && !isFirstStep_ && !this->frozen_) {
    nowAddAHill = true;
  } else {
    nowAddAHill = false;
    this->isFirstStep_ = false;
  }

  std::vector<double> cv(this->d_);
  for(unsigned i = 0; i < this->d_; ++i) {
    cv[i] = getArgument(i);
  }

  if(nowAddAHill) {
    double height = getHeight(cv);

    if(this->walkers_mpi_) {
      std::vector<double> all_cv(this->mpi_size_ * this->d_, 0.0);
      std::vector<double> all_sigma(this->mpi_size_ * this->sigma0_.size(), 0.0);
      std::vector<double> all_height(this->mpi_size_, 0.0);
      if(comm.Get_rank() == 0) {
        multi_sim_comm.Allgather(cv, all_cv);
        multi_sim_comm.Allgather(this->sigma0_, all_sigma);
        multi_sim_comm.Allgather(height * (this->biasf_ > 1.0 ? this->biasf_ / (this->biasf_ - 1.0) : 1.0), all_height);
      }
      comm.Bcast(all_cv, 0);
      comm.Bcast(all_sigma, 0);
      comm.Bcast(all_height, 0);

      for(int i = 0; i < this->mpi_size_; i++) {
        std::vector<double> cv_now(this->d_);
        std::vector<double> sigma_now(this->sigma0_.size());
        for(unsigned j = 0; j < this->d_; j++) {
          cv_now[j] = all_cv[i * this->d_ + j];
        }
        for(unsigned j = 0; j < this->sigma0_.size(); j++) {
          sigma_now[j] = all_sigma[i * this->sigma0_.size() + j];
        }
        double fact = (this->biasf_ > 1.0 ? (this->biasf_ - 1.0) / this->biasf_ : 1.0);
        Gaussian newhill(all_height[i] * fact, cv_now, sigma_now);
        this->hills_.push_back(newhill);
        if(this->mpi_rank_ == 0) {
          writeGaussian(newhill, hillsOfile_);
        }
      }
    } else {
      Gaussian newhill(height, cv, this->sigma0_);
      this->hills_.push_back(newhill);
      writeGaussian(newhill, hillsOfile_);
    }
  }

  if(getStep() % this->sketch_stride_ == 1 && !this->frozen_) {
    log << "Vbias update " << this->sketch_count_ << "...\n\n";
    log.flush();
  }
}

double TTMetaD::getHeight(const std::vector<double>& cv) {
  double height = this->height0_;
  if(this->welltemp_) {
    // the bias felt by the system (see calculate())
    std::vector<double> cvc(cv);
    std::vector<char> out;
    clampToInterval(cvc, out);
    double vbias = getBias(cvc);
    if(this->biasf_ > 1.0) {
      height = this->height0_ * exp(-vbias / (this->kbt_ * (this->biasf_ - 1.0)));
    } else {
      height = this->height0_ * exp(-vbias / this->kbt_);
    }
  }
  return height;
}

// Total bias = TT approximation of accumulated hills + sum of recently deposited Gaussians.
double TTMetaD::getBias(const std::vector<double>& cv) {
  double bias = this->vbDense_.empty() ? 0.0 : this->vbDense_.eval(this->sketch_basis_, cv, this->sketch_conv_, nullptr);
  unsigned nt = OpenMP::getNumThreads();
  #pragma omp parallel num_threads(nt)
  {
    #pragma omp for reduction(+:bias) nowait
    for(unsigned i = 0; i < hills_.size(); ++i) {
      bias += evaluateGaussian(cv, this->hills_[i]);
    }
  }
  return bias;
}

double TTMetaD::getBiasAndDerivatives(const std::vector<double>& cv, std::vector<double>& der) {
  double bias = 0.0;
  if (!this->vbDense_.empty()) {
    bias = this->vbDense_.eval(this->sketch_basis_, cv, this->sketch_conv_, &der);
  } else {
    der.assign(this->d_, 0.0);
  }
  unsigned nt = OpenMP::getNumThreads();
  if(this->hills_.size() < 2 * nt || nt == 1) {
    std::vector<double> dp(this->d_);
    for(unsigned i = 0; i < this->hills_.size(); ++i) {
      bias += evaluateGaussianAndDerivatives(cv, this->hills_[i], der, dp);
    }
  } else {
    #pragma omp parallel num_threads(nt)
    {
      std::vector<double> omp_deriv(this->d_, 0.0);
      std::vector<double> dp(this->d_);
      #pragma omp for reduction(+:bias) nowait
      for(unsigned i = 0; i < this->hills_.size(); ++i) {
        bias += evaluateGaussianAndDerivatives(cv, this->hills_[i], omp_deriv, dp);
      }
      #pragma omp critical
      for(unsigned i = 0; i < this->d_; ++i) {
        der[i] += omp_deriv[i];
      }
    }
  }
  return bias;
}

double TTMetaD::evaluateGaussian(const std::vector<double>& cv, const Gaussian& hill) {
  double dp2 = 0.0;
  for(unsigned i = 0; i < this->d_; i++) {
    double dp = difference(i, hill.center[i], cv[i]) * hill.invsigma[i];
    dp2 += dp * dp;
  }
  dp2 *= 0.5;

  double bias = 0.0;
  if(dp2 < dp2cutoff) {
    bias = hill.height * exp(-dp2);
  }

  return bias;
}

double TTMetaD::evaluateGaussianAndDerivatives(const std::vector<double>& cv, const Gaussian& hill, std::vector<double>& der, std::vector<double>& dp) {
  double dp2 = 0.0;
  double bias = 0.0;
  for(unsigned i = 0; i < this->d_; i++) {
    dp[i] = difference(i, hill.center[i], cv[i]) * hill.invsigma[i];
    dp2 += dp[i] * dp[i];
  }
  dp2 *= 0.5;
  if(dp2 < dp2cutoff) {
    bias = hill.height * exp(-dp2);
    for(unsigned i = 0; i < this->d_; i++) {
      der[i] -= bias * dp[i] * hill.invsigma[i];
    }
  }

  return bias;
}

// TT-Sketch algorithm: compresses the N accumulated Gaussian hills plus the previous
// TT bias vb_ into a new low-rank TT approximation G, which becomes the updated vb_.
//
// Overview (see paper Algorithm 2):
//   1. createTTCoeff(): build a random sketch coefficient TT `coeff` of rank sketch_rc_.
//   2. intBasisSample(): compute M[i](j,k) = <phi_k, g_j>_{dim i}, the inner product of
//      each basis function with each Gaussian along dimension i.
//   3. formTensorMoment(): contract coeff with M to get the tensor moment Bemp and the
//      cross-gram matrices A of the bonds (normal equations).
//   4. For each bond between cores k and k+1:
//      - Add the contribution of the previous bias (or of the previous moments) to A.
//      - SVD A to get the rank-trimmed projection V (truncating to sketch_r_ or sketch_cutoff_).
//   5. Recover TT cores G from Bemp and the pseudo-inverse of U*S.
//   6. Apply Gram matrix correction (kernel basis only) to convert to dual basis coefficients.
void TTMetaD::paraSketch() {
  auto coeff = createTTCoeff();
  auto M = intBasisSample();
  MPS G(this->d_);

  auto [Bemp, AMats] = formTensorMoment(M, coeff);
  MPS Bemp_Vb;
  std::vector<ITensor> envi_L_Vb;
  std::vector<ITensor> envi_R_Vb;
  if(this->sketch_count_ != 1) {
    if(this->nonintrusive_) {
      // Non-intrusive variant: accumulate tensor moments directly across sketches.
      // B_prev_ stores the sum of Bemp from all prior intervals; re-index its bond/site
      // indices to match the current coeff TT's indices before adding.
      auto sites = siteInds(Bemp);
      auto sites_prev = siteInds(this->B_prev_);
      auto links = linkInds(Bemp);
      auto links_prev = linkInds(this->B_prev_);
      for(unsigned i = 1; i <= this->d_; ++i) {
        this->B_prev_.ref(i) *= delta(sites(i), sites_prev(i));
        if(i != 1) {
          this->B_prev_.ref(i) *= delta(links(i - 1), links_prev(i - 1));
        }
        if (i != this->d_) {
          this->B_prev_.ref(i) *= delta(links(i), links_prev(i));
        }
        Bemp.ref(i) += this->B_prev_(i);
      }
    } else {
      // Intrusive variant: project the previous TT bias vb_ onto the current sketch's
      // basis by pre-multiplying each core by the Gram matrix (kernel basis) or the
      // convolution inner-product matrix (Fourier conv mode), then add its tensor moment.
      if(this->sketch_basis_[0].kernel()) {
        // Convert vb_ from dual basis to standard basis by applying Gram matrix per core
        for(unsigned i = 1; i <= this->d_; ++i) {
          auto s = siteIndex(this->vb_, i);
          ITensor gram(s, prime(s));
          for(int j = 1; j <= dim(s); ++j) {
            for(int l = 1; l <= dim(s); ++l) {
              gram.set(s = j, prime(s) = l, this->sketch_basis_[i - 1].gram()(j - 1, l - 1));
            }
          }
          this->vb_.ref(i) *= gram;
          this->vb_.ref(i).noPrime();
        }
      }
      // Note: in SKETCH_WIDTH (conv) mode the coefficients of vb_ are those of the unsmoothed bias,
      // the smoothing being applied when the basis is evaluated, so they are added as they are.
      // Damping them here as well (as done before) smoothed the old bias again at every sketch,
      // so that the intrusive variant did not match NONINTRUSIVE even when the TT is exact.
      auto vbresult = formTensorMomentVb(coeff);
      Bemp_Vb = std::get<0>(vbresult);
      envi_L_Vb = std::get<1>(vbresult);
      envi_R_Vb = std::get<2>(vbresult);
      for(unsigned i = 1; i <= this->d_; ++i) {
        Bemp.ref(i) += Bemp_Vb(i);
      }
    }
  }
  if(this->nonintrusive_) {
    this->B_prev_ = Bemp;
    if(this->sketch_count_ == 1) {
      this->A_prev_.resize(this->d_ - 1);
    }
  }
  // For each bond between cores k-1 and k, form the cross-gram matrix A = L^T * R,
  // where L = envi_L[k] and R = envi_R[k-2] are the left/right environment projections
  // of coeff onto the sample indices. SVD A to find the rank-trimmed subspace V.
  auto links = linkInds(coeff);
  std::vector<ITensor> U(this->d_), S(this->d_), V(this->d_);
  std::vector<Index> links_trimmed;
  for(unsigned core_id = 2; core_id <= this->d_; ++core_id) {
    int rank = dim(links(core_id - 1));
    const Matrix<double>& AMat = AMats[core_id - 2];
    Matrix<double> LMat, RMat, Lt, AMat_Vb;

    if(this->sketch_count_ != 1 && !this->nonintrusive_) {
      auto ivb = linkIndex(this->vb_, core_id - 1);
      int rank_vb = dim(ivb);
      LMat = Matrix<double>(rank_vb, rank);
      RMat = Matrix<double>(rank_vb, rank);
      for(int i = 1; i <= rank_vb; ++i) {
        for(int j = 1; j <= rank; ++j) {
          LMat(i - 1, j - 1) = envi_L_Vb[core_id - 1].elt(ivb = i, links(core_id - 1) = j);
          RMat(i - 1, j - 1) = envi_R_Vb[core_id - 2].elt(ivb = i, links(core_id - 1) = j);
        }
      }
      transpose(LMat, Lt);
      mult(Lt, RMat, AMat_Vb);
    }

    ITensor A(prime(links(core_id - 1)), links(core_id - 1));
    for(int i = 1; i <= rank; ++i) {
      for(int j = 1; j <= rank; ++j) {
        A.set(prime(links(core_id - 1)) = i, links(core_id - 1) = j, AMat(i - 1, j - 1));
      }
    }
    if(this->sketch_count_ != 1) {
      if(this->nonintrusive_) {
        this->A_prev_[core_id - 2] *= delta(this->A_prev_[core_id - 2].index(1), A.index(1));
        this->A_prev_[core_id - 2] *= delta(this->A_prev_[core_id - 2].index(2), A.index(2));
        A += this->A_prev_[core_id - 2];
      } else {
        ITensor A_Vb(prime(links(core_id - 1)), links(core_id - 1));
        for(int i = 1; i <= rank; ++i) {
          for(int j = 1; j <= rank; ++j) {
            A_Vb.set(prime(links(core_id - 1)) = i, links(core_id - 1) = j, AMat_Vb(i - 1, j - 1));
          }
        }
        A += A_Vb;
      }
    }
    if(this->nonintrusive_) {
      this->A_prev_[core_id - 2] = A;
    }
    auto original_link_tags = tags(links(core_id - 1));
    V[core_id - 1] = ITensor(links(core_id - 1));
    if(this->sketch_r_ > 0) {
      svd(A, U[core_id - 1], S[core_id - 1], V[core_id - 1],
      {"Cutoff=", this->sketch_cutoff_, "RightTags=", original_link_tags, "MaxDim=", this->sketch_r_});
    } else {
      svd(A, U[core_id - 1], S[core_id - 1], V[core_id - 1], {"Cutoff=", this->sketch_cutoff_, "RightTags=", original_link_tags});
    }
    links_trimmed.push_back(commonIndex(S[core_id - 1], V[core_id - 1]));
  }

  // Recover TT cores G from Bemp and the trimmed subspaces.
  // G[1] = Bemp[1] * V[1] (project first core onto rank-trimmed subspace)
  // G[k] = pinv(U[k]*S[k]) * Bemp[k] * V[k]  for 2 <= k < d
  // G[d] = pinv(U[d]*S[d]) * Bemp[d]
  G.ref(1) = Bemp(1) * V[1];
  for(unsigned core_id = 2; core_id <= this->d_; ++core_id) {
    int rank = dim(links(core_id - 1)), rank_trimmed = dim(links_trimmed[core_id - 2]);
    ITensor A = U[core_id - 1] * S[core_id - 1];
    ITensor Pinv(links_trimmed[core_id - 2], links(core_id - 1));
    Matrix<double> AMat(rank, rank_trimmed), PMat;
    for(int i = 1; i <= rank; ++i) {
      for(int j = 1; j <= rank_trimmed; ++j) {
        AMat(i - 1, j - 1) = A.elt(prime(links(core_id - 1)) = i, links_trimmed[core_id - 2] = j);
      }
    }
    pseudoInvert(AMat, PMat);

    for(int i = 1; i <= rank_trimmed; ++i) {
      for(int j = 1; j <= rank; ++j) {
        Pinv.set(links_trimmed[core_id - 2] = i, links(core_id - 1) = j, PMat(i - 1, j - 1));
      }
    }
    G.ref(core_id) = Pinv * Bemp(core_id);
    if(core_id != this->d_) {
      G.ref(core_id) *= V[core_id];
    }
  }

  log << "Final ranks ";
  for(unsigned i = 1; i < this->d_; ++i) {
    log << dim(linkIndex(G, i)) << " ";
  }
  log << "\n";
  log.flush();

  if(this->sketch_basis_[0].kernel()) {
    // Apply pseudo-inverse of the Gram matrix to each core to convert from the
    // primal Gram representation back to dual basis coefficients (G^+ * G_core).
    for(unsigned i = 1; i <= this->d_; ++i) {
      auto s = siteIndex(G, i);
      ITensor ginv(s, prime(s));
      for(int j = 1; j <= dim(s); ++j) {
        for(int l = 1; l <= dim(s); ++l) {
          ginv.set(s = j, prime(s) = l, this->sketch_basis_[i - 1].ginv()(j - 1, l - 1));
        }
      }
      G.ref(i) *= ginv;
      G.ref(i).noPrime();
    }
  }

  setVb(G);
}

// Build the random sketch coefficient TT `coeff` of bond dimension sketch_rc_.
// Each element is drawn from N(0,1); then each core is multiplied by a diagonal
// matrix diag(1, alpha, alpha, ...) so the constant basis function (pos=1) has
// full weight while all other harmonics are scaled down by sketch_alpha_. This
// ensures the sketch captures the large constant component of the bias accurately.
MPS TTMetaD::createTTCoeff() const {
  unsigned seed = this->deterministic_ ? 42u : static_cast<unsigned int>(time(nullptr));
  // all ranks of a replica must build the same sketch
  comm.Bcast(seed, 0);
  std::default_random_engine generator(seed);
  std::normal_distribution<double> distribution(0.0, 1.0);
  int n = this->sketch_basis_[0].nbasis();
  auto sites = SiteSet(this->d_, n);
  auto coeff = MPS(sites, this->sketch_rc_);
  for(int j = 1; j <= n; ++j) {
    for(int k = 1; k <= this->sketch_rc_; ++k) {
      coeff.ref(1).set(sites(1) = j, linkIndex(coeff, 1) = k, distribution(generator));
    }
  }
  for(unsigned i = 2; i <= this->d_ - 1; ++i) {
    for(int j = 1; j <= n; ++j) {
      for(int k = 1; k <= this->sketch_rc_; ++k) {
        for(int l = 1; l <= this->sketch_rc_; ++l) {
          coeff.ref(i).set(sites(i) = j, linkIndex(coeff, i - 1) = k, linkIndex(coeff, i) = l, distribution(generator));
        }
      }
    }
  }
  for(int j = 1; j <= n; ++j) {
    for(int k = 1; k <= this->sketch_rc_; ++k) {
      coeff.ref(this->d_).set(sites(this->d_) = j, linkIndex(coeff, this->d_ - 1) = k, distribution(generator));
    }
  }
  for(unsigned i = 1; i <= this->d_; ++i) {
    auto s = sites(i);
    auto sp = prime(s);
    std::vector<double> Avec(n, this->sketch_alpha_);
    Avec[0] = 1.0;
    auto A = diagITensor(Avec, s, sp);
    coeff.ref(i) *= A;
    coeff.ref(i).noPrime();
  }
  return coeff;
}

// Compute the matrix M[dim i][sample j, basis k] = <phi_k, g_j>_{dim i}, the
// analytical inner product of each 1D basis function with the 1D marginal of each
// Gaussian hill along dimension i. The d-dimensional Gaussian factorizes as a
// product of d 1D Gaussians, so we distribute height as h^(1/d) per dimension.
// These integrals are computed in closed form:
//   - Kernel basis: convolution of the kernel function with the 1D Gaussian hill
//   - Fourier basis: Fourier transform of the 1D Gaussian (exponential decay * harmonic)
// Returns M as a vector of ITensors (one per dimension) and a new IndexSet is_new
// where is_new(i) indexes the N sample points for dimension i.
std::vector<Matrix<double>> TTMetaD::intBasisSample() const {
  unsigned N = this->hills_.size();
  int nb = this->sketch_basis_[0].nbasis();
  std::vector<Matrix<double>> M(this->d_, Matrix<double>(N, nb));
  for(unsigned i = 1; i <= this->d_; ++i) {
    double L = (this->sketch_basis_[i - 1].dom().second - this->sketch_basis_[i - 1].dom().first) / 2;
    double a = (this->sketch_basis_[i - 1].dom().second + this->sketch_basis_[i - 1].dom().first) / 2;
    for(unsigned j = 1; j <= N; ++j) {
      double x = this->hills_[j - 1].center[i - 1];
      double w = this->hills_[j - 1].sigma[i - 1];
      double h = pow(this->hills_[j - 1].height, 1.0 / this->d_);  // per-dimension d-th root of height
      for(int pos = 1; pos <= nb; ++pos) {
        double result = 0.0;
        if(this->sketch_basis_[0].kernel()) {
          if(pos == 1) {
            // <1, g_j> = integral of 1D Gaussian = sqrt(2*pi)*sigma
            result = h * sqrt(2 * M_PI) * w;
          } else {
            // <kernel_c, g_j> = conv of two Gaussians (widths dx and w), with periodic images
            double c = this->sketch_basis_[i - 1].center(pos - 1);
            double dx = this->sketch_basis_[i - 1].dx();
            for(int k = -1; k <= 1; ++k) {
              result += exp(-pow(x - c + 2 * k * L, 2) / (2 * (pow(dx, 2) + pow(w, 2)))) * h * sqrt(2 * M_PI) * w /
                        (sqrt(1 / pow(dx, 2) + 1 / pow(w, 2)) * w);
            }
          }
        } else {
          if(pos == 1) {
            // <1/sqrt(2L), g_j> = sqrt(pi/L) * sigma (Gaussian times constant)
            result = h * sqrt(M_PI / L) * w;
          } else if(pos % 2 == 0) {
            // <cos harmonic, g_j>: Fourier transform of Gaussian times exponential damping
            result = exp(-pow(M_PI * w * (pos / 2), 2) / (2 * pow(L, 2))) * h * sqrt(2 * M_PI / L) * w * cos(M_PI * (x - a) * (pos / 2) / L);
          } else {
            result = exp(-pow(M_PI * w * (pos / 2), 2) / (2 * pow(L, 2))) * h * sqrt(2 * M_PI / L) * w * sin(M_PI * (x - a) * (pos / 2) / L);
          }
        }
        M[i - 1](j - 1, pos - 1) = result;
      }
    }
  }
  return M;
}

// Compute the tensor moment B of the new Gaussian hills and the cross-gram matrices A.
//
// With K_i(j) = sum_p coeff_i[:, p, :] M_i(j, p) the core i of coeff contracted with the
// basis integrals of hill j (a r_{i-1} x r_i matrix), the left and right environments are
//   EL_i(j, :) = K_1(j) ... K_i(j)       (i = 1..d-1, size r_i)
//   ER_i(j, :) = K_{i+2}(j) ... K_d(j)   (i = 0..d-2, size r_{i+1})
// The tensor moment is
//   B_1[p, b]    = sum_j ER_0(j, b) M_1(j, p)
//   B_k[a, p, b] = sum_j EL_{k-1}(j, a) ER_{k-1}(j, b) M_k(j, p)   (1 < k < d)
//   B_d[a, p]    = sum_j EL_{d-1}(j, a) M_d(j, p)
// and the cross-gram matrix of bond k = 1..d-1 is A_k = EL_k^T ER_{k-1}.
// Everything is done on plain arrays (the environments are independent across hills),
// only B is converted to ITensors at the end.
std::pair<MPS, std::vector<Matrix<double>>> TTMetaD::formTensorMoment(const std::vector<Matrix<double>>& M, const MPS& coeff) const {
  const unsigned d = this->d_;
  const unsigned N = M[0].nrows();
  const DenseTT C(coeff);
  const auto sites = siteInds(coeff);
  const auto links = linkInds(coeff);
  std::vector<int> r(d + 1), n(d);
  for(unsigned i = 0; i <= d; ++i) {
    r[i] = C.rank(i);
  }
  for(unsigned i = 0; i < d; ++i) {
    n[i] = C.nbasis(i);
  }
  // EL[i] for i = 1..d-1 (N x r_i), ER[i] for i = 0..d-2 (N x r_{i+1}), row-major
  std::vector<std::vector<double>> EL(d), ER(d - 1);
  for(unsigned i = 1; i < d; ++i) {
    EL[i].assign(static_cast<size_t>(N) * r[i], 0.0);
    ER[i - 1].assign(static_cast<size_t>(N) * r[i], 0.0);
  }
  const unsigned nt = OpenMP::getNumThreads();
  #pragma omp parallel num_threads(nt)
  {
    std::vector<std::vector<double>> K(d);
    std::vector<double> tmp;
    #pragma omp for schedule(static)
    for(unsigned j = 0; j < N; ++j) {
      // contracted cores K[i] = K_{i+1}(j), r_i x r_{i+1}
      for(unsigned i = 0; i < d; ++i) {
        const int rl = r[i], rr = r[i + 1], nb = n[i];
        K[i].assign(static_cast<size_t>(rl) * rr, 0.0);
        const double* c = C.core(i).data();
        for(int a = 0; a < rl; ++a) {
          double* k = K[i].data() + static_cast<size_t>(a) * rr;
          for(int p = 0; p < nb; ++p) {
            const double m = M[i](j, p);
            const double* cp = c + (static_cast<size_t>(a) * nb + p) * rr;
            for(int b = 0; b < rr; ++b) {
              k[b] += cp[b] * m;
            }
          }
        }
      }
      // left environments
      for(int b = 0; b < r[1]; ++b) {
        EL[1][static_cast<size_t>(j) * r[1] + b] = K[0][b];
      }
      for(unsigned i = 2; i < d; ++i) {
        const int rl = r[i - 1], rr = r[i];
        const double* prev = EL[i - 1].data() + static_cast<size_t>(j) * rl;
        double* out = EL[i].data() + static_cast<size_t>(j) * rr;
        for(int a = 0; a < rl; ++a) {
          const double* k = K[i - 1].data() + static_cast<size_t>(a) * rr;
          for(int b = 0; b < rr; ++b) {
            out[b] += prev[a] * k[b];
          }
        }
      }
      // right environments
      for(int a = 0; a < r[d - 1]; ++a) {
        ER[d - 2][static_cast<size_t>(j) * r[d - 1] + a] = K[d - 1][a];
      }
      for(int i = static_cast<int>(d) - 3; i >= 0; --i) {
        // ER[i] = K_{i+2} ER[i+1], K_{i+2} = K[i+1] is r_{i+1} x r_{i+2}
        const int rl = r[i + 1], rr = r[i + 2];
        const double* next = ER[i + 1].data() + static_cast<size_t>(j) * rr;
        double* out = ER[i].data() + static_cast<size_t>(j) * rl;
        for(int a = 0; a < rl; ++a) {
          const double* k = K[i + 1].data() + static_cast<size_t>(a) * rr;
          double sum = 0.0;
          for(int b = 0; b < rr; ++b) {
            sum += k[b] * next[b];
          }
          out[a] = sum;
        }
      }
    }
  }

  // cross-gram matrices
  std::vector<Matrix<double>> A(d - 1);
  for(unsigned k = 1; k < d; ++k) {
    const int rk = r[k];
    A[k - 1] = Matrix<double>(rk, rk);
    #pragma omp parallel for num_threads(nt) schedule(static)
    for(int a = 0; a < rk; ++a) {
      for(int b = 0; b < rk; ++b) {
        double sum = 0.0;
        for(unsigned j = 0; j < N; ++j) {
          sum += EL[k][static_cast<size_t>(j) * rk + a] * ER[k - 1][static_cast<size_t>(j) * rk + b];
        }
        A[k - 1](a, b) = sum;
      }
    }
  }

  // tensor moment
  MPS B(d);
  {
    const int nb = n[0], rr = r[1];
    ITensor core(sites(1), links(1));
    for(int p = 0; p < nb; ++p) {
      for(int b = 0; b < rr; ++b) {
        double sum = 0.0;
        for(unsigned j = 0; j < N; ++j) {
          sum += ER[0][static_cast<size_t>(j) * rr + b] * M[0](j, p);
        }
        core.set(sites(1) = p + 1, links(1) = b + 1, sum);
      }
    }
    B.ref(1) = core;
  }
  for(unsigned k = 2; k < d; ++k) {
    const int rl = r[k - 1], rr = r[k], nb = n[k - 1];
    std::vector<double> Bk(static_cast<size_t>(rl) * rr * nb, 0.0); // [a][b][p]
    #pragma omp parallel for num_threads(nt) schedule(static)
    for(int a = 0; a < rl; ++a) {
      for(unsigned j = 0; j < N; ++j) {
        const double la = EL[k - 1][static_cast<size_t>(j) * rl + a];
        const double* er = ER[k - 1].data() + static_cast<size_t>(j) * rr;
        for(int b = 0; b < rr; ++b) {
          const double w = la * er[b];
          double* out = Bk.data() + (static_cast<size_t>(a) * rr + b) * nb;
          for(int p = 0; p < nb; ++p) {
            out[p] += w * M[k - 1](j, p);
          }
        }
      }
    }
    ITensor core(links(k - 1), sites(k), links(k));
    for(int a = 0; a < rl; ++a) {
      for(int b = 0; b < rr; ++b) {
        for(int p = 0; p < nb; ++p) {
          core.set(links(k - 1) = a + 1, sites(k) = p + 1, links(k) = b + 1, Bk[(static_cast<size_t>(a) * rr + b) * nb + p]);
        }
      }
    }
    B.ref(k) = core;
  }
  {
    const int nb = n[d - 1], rl = r[d - 1];
    ITensor core(links(d - 1), sites(d));
    for(int a = 0; a < rl; ++a) {
      for(int p = 0; p < nb; ++p) {
        double sum = 0.0;
        for(unsigned j = 0; j < N; ++j) {
          sum += EL[d - 1][static_cast<size_t>(j) * rl + a] * M[d - 1](j, p);
        }
        core.set(links(d - 1) = a + 1, sites(d) = p + 1, sum);
      }
    }
    B.ref(d) = core;
  }
  return std::make_pair(B, A);
}

// Compute the tensor moment B and environments for the previous TT bias vb_.
// This is the intrusive variant's contribution from the prior accumulated bias.
// Analogous to formTensorMoment but uses vb_ in place of M*coeff: the "samples"
// are the vb_ bond indices, and the environments are formed by contracting coeff
// with vb_ core by core. The result has the same structure as Bemp from the new hills.
std::tuple<MPS, std::vector<ITensor>, std::vector<ITensor>> TTMetaD::formTensorMomentVb(const MPS& coeff) {
  // align vb_ site indices with coeff's site indices for contraction
  for(unsigned i = 1; i <= this->d_; ++i) {
    this->vb_.ref(i) *= delta(siteIndex(this->vb_, i), siteIndex(coeff, i));
  }
  std::vector<ITensor> envi_L(this->d_);
  envi_L[1] = coeff(1) * this->vb_(1);
  for(unsigned i = 2; i < this->d_; ++i) {
    envi_L[i] = envi_L[i - 1] * coeff(i) * this->vb_(i);
  }

  std::vector<ITensor> envi_R(this->d_);
  envi_R[this->d_ - 2] = coeff(this->d_) * this->vb_(this->d_);
  for(int i = this->d_ - 3; i >= 0; --i) {
    envi_R[i] = envi_R[i + 1] * coeff(i + 2) * this->vb_(i + 2);
  }

  MPS B(this->d_);
  B.ref(1) = this->vb_(1) * envi_R[0];
  for(unsigned core_id = 2; core_id < this->d_; ++core_id) {
    B.ref(core_id) = envi_L[core_id - 1] * this->vb_(core_id) * envi_R[core_id - 1];
  }
  B.ref(this->d_) = envi_L[this->d_ - 1] * this->vb_(this->d_);

  return std::make_tuple(B, envi_L, envi_R);
}

#endif // !defined(__PLUMED_HAS_LIBITENSOR) || !defined(__PLUMED_HAS_LIBHDF5)

void TTMetaD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.addFlag("KERNEL_BASIS", false, "Specifies that local kernel basis should be used instead of Fourier basis");
  keys.add("compulsory", "SIGMA", "the widths of the Gaussian hills");
  keys.add("compulsory", "PACE", "the frequency for hill addition");
  keys.add("compulsory", "FILE", "HILLS", "a file in which the list of added hills is stored");
  keys.add("compulsory", "TT_FILE", "ttsketch.h5", "the HDF5 file where the TT bias is stored after every sketch, and read at restart. "
           "With multiple replicas the replica suffix is added, with WALKERS_MPI the file is shared and is in the walkers directory");
  keys.add("compulsory", "HEIGHT", "the heights of the Gaussian hills");
  keys.add("optional", "FMT", "specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional", "BIASFACTOR", "use well tempered metadynamics and use this bias factor. Please note you must also specify temp");
  keys.add("optional", "TEMP", "the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.addFlag("WALKERS_MPI", false, "To be used when gromacs + multiple walkers are used");
  keys.add("optional", "WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.use("RESTART");
  keys.add("optional", "SKETCH_RANK", "Target rank for TTSketch algorithm - compulsory if SKETCH_CUTOFF is not specified");
  keys.add("optional", "SKETCH_CUTOFF", "Truncation error cutoff for singular value decomposition - compulsory if SKETCH_RANK is not specified");
  keys.add("compulsory", "SKETCH_INITRANK", "Initial rank for TTSketch algorithm");
  keys.add("compulsory", "SKETCH_PACE", "1e6", "The frequency for TT Vbias updates");
  keys.add("compulsory", "INTERVAL_MIN", "Lower limits of the domain of the basis functions. Outside the limits the bias is the one at the closest limit and the system does not feel the biasing force along that argument");
  keys.add("compulsory", "INTERVAL_MAX", "Upper limits of the domain of the basis functions. Outside the limits the bias is the one at the closest limit and the system does not feel the biasing force along that argument");
  keys.add("compulsory", "SKETCH_NBASIS", "20", "Number of basis functions per dimension");
  keys.add("compulsory", "SKETCH_ALPHA", "0.05", "Weight coefficient for random tensor train construction");
  keys.add("optional", "SKETCH_UNTIL", "After this time, the hills deposited so far are compressed into the TT one last time and the bias potential freezes");
  keys.add("optional", "SKETCH_WIDTH", "Width of Gaussian kernels for smoothing");
  keys.add("optional", "KERNEL_DX", "Width of basis function kernels");
  keys.addFlag("NONINTRUSIVE", false, "Sketching uses previous exact sum of Gaussians instead of TT approximation");
  keys.addFlag("DETERMINISTIC", false, "Use a fixed random seed for TT sketch construction, ensuring reproducible results. Intended for regression testing.");
}

PLUMED_REGISTER_ACTION(TTMetaD, "TTMETAD")

} // namespace ttsketch
} // namespace PLMD
