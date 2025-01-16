/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "core/FlexibleBin.h"
#include "tools/Exception.h"
#include "tools/Grid.h"
#include "tools/Matrix.h"
#include "tools/OpenMP.h"
#include "tools/Random.h"
#include "tools/File.h"
#include <ctime>
#include <numeric>
#if defined(__PLUMED_HAS_GETCWD)
#include <unistd.h>
#endif

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS PBMETAD
/*
Used to performed Parallel Bias metadynamics.

This action activate Parallel Bias Metadynamics (PBMetaD) \cite pbmetad, a version of metadynamics \cite metad in which
multiple low-dimensional bias potentials are applied in parallel.
In the current implementation, these have the form of mono-dimensional metadynamics bias
potentials:

\f[
{V(s_1,t), ..., V(s_N,t)}
\f]

where:

\f[
V(s_i,t) = \sum_{ k \tau < t} W_i(k \tau)
\exp\left(
- \frac{(s_i-s_i^{(0)}(k \tau))^2}{2\sigma_i^2}
\right).
\f]

To ensure the convergence of each mono-dimensional bias potential to the corresponding free energy,
at each deposition step the Gaussian heights are multiplied by the so-called conditional term:

\f[
W_i(k \tau)=W_0 \frac{\exp\left(
- \frac{V(s_i,k \tau)}{k_B T}
\right)}{\sum_{i=1}^N
\exp\left(
- \frac{V(s_i,k \tau)}{k_B T}
\right)}
\f]

where \f$W_0\f$ is the initial Gaussian height.

The PBMetaD bias potential is defined by:

\f[
V_{PB}(\vec{s},t) = -k_B T \log{\sum_{i=1}^N
\exp\left(
- \frac{V(s_i,t)}{k_B T}
\right)}.
\f]


Information on the Gaussian functions that build each bias potential are printed to
multiple HILLS files, which
are used both to restart the calculation and to reconstruct the mono-dimensional
free energies as a function of the corresponding CVs.
These can be reconstructed using the \ref sum_hills utility because the final bias is given
by:

\f[
V(s_i) = -F(s_i)
\f]

Currently, only a subset of the \ref METAD options are available in PBMetaD.

The bias potentials can be stored on a grid to increase performances of long PBMetaD simulations.
You should
provide either the number of bins for every collective variable (GRID_BIN) or
the desired grid spacing (GRID_SPACING). In case you provide both PLUMED will use
the most conservative choice (highest number of bins) for each dimension.
In case you do not provide any information about bin size (neither GRID_BIN nor GRID_SPACING)
and if Gaussian width is fixed PLUMED will use 1/5 of the Gaussian width as grid spacing.
This default choice should be reasonable for most applications.

Another option that is available is well-tempered metadynamics \cite Barducci:2008. In this
variant of PBMetaD the heights of the Gaussian hills are scaled at each step by the
additional well-tempered metadynamics term.
This  ensures that each bias converges more smoothly. It should be noted that, in the case of well-tempered metadynamics, in
the output printed the Gaussian height is re-scaled using the bias factor.
Also notice that with well-tempered metadynamics the HILLS files do not contain the bias,
but the negative of the free-energy estimate. This choice has the advantage that
one can restart a simulation using a different value for the \f$\Delta T\f$. The applied bias will be scaled accordingly.

Note that you can use here also the flexible Gaussian approach  \cite Branduardi:2012dl
in which you can adapt the Gaussian to the extent of Cartesian space covered by a variable or
to the space in collective variable covered in a given time. In this case the width of the deposited
Gaussian potential is denoted by one value only that is a Cartesian space (ADAPTIVE=GEOM) or a time
(ADAPTIVE=DIFF). Note that a specific integration technique for the deposited Gaussian kernels
should be used in this case. Check the documentation for utility sum_hills.

With the keyword INTERVAL one changes the metadynamics algorithm setting the bias force equal to zero
outside boundary \cite baftizadeh2012protein. If, for example, metadynamics is performed on a CV s and one is interested only
to the free energy for s > boundary, the history dependent potential is still updated according to the above
equations but the metadynamics force is set to zero for s < boundary. Notice that Gaussian kernels are added also
if s < boundary, as the tails of these Gaussian kernels influence VG in the relevant region s > boundary. In this way, the
force on the system in the region s > boundary comes from both metadynamics and the force field, in the region
s < boundary only from the latter. This approach allows obtaining a history-dependent bias potential VG that
fluctuates around a stable estimator, equal to the negative of the free energy far enough from the
boundaries. Note that:
- It works only for one-dimensional biases;
- It works both with and without GRID;
- The interval limit boundary in a region where the free energy derivative is not large;
- If in the region outside the limit boundary the system has a free energy minimum, the INTERVAL keyword should
  be used together with a \ref UPPER_WALLS or \ref LOWER_WALLS at boundary.

For systems with multiple CVs that share identical properties, PBMetaD with partitioned families can be used
to group them under one bias potential that each contributes to \cite Prakash2018PF. This is done with a list
of PF keywords, where each PF* argument contains the list of CVs from ARG to be placed in that family. Once
invoked, each CV in ARG must be placed in exactly one PF, even if it results in families containing only one CV.
Additionally, in cases where each of SIGMA or GRID entry would correspond to each ARG entry, they now correspond to
each PF and must be adjusted accordingly.

Multiple walkers  \cite multiplewalkers can also be used. See below the examples.

\par Examples

The following input is for PBMetaD calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the PBMetaD bias potential are written to the COLVAR file every 100 steps.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 LABEL=pb FILE=HILLS_d1,HILLS_d2
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
\endplumedfile
(See also \ref DISTANCE and \ref PRINT).

\par
If you use well-tempered metadynamics, you should specify a single bias factor and initial
Gaussian height.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
... PBMETAD
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
\endplumedfile

\par
Using partitioned families, each CV in ARG must be in exactly one family. Here,
the distance between atoms 1,2 is degenerate with 2,4, but not with the
distance between 3,5. Note that two SIGMA are provided to match the two PF.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
DISTANCE ATOMS=1,2 LABEL=d3
PBMETAD ...
ARG=d1,d2,d3 SIGMA=0.2,0.2 HEIGHT=0.3
PF0=d1 PF1=d2,d3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
... PBMETAD
PRINT ARG=d1,d2,d3,pb.bias STRIDE=100 FILE=COLVAR
\endplumedfile

\par
The following input enables the MPI version of multiple-walkers.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
WALKERS_MPI
... PBMETAD
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
\endplumedfile

\par
The disk version of multiple-walkers can be
enabled by setting the number of walker used, the id of the
current walker which interprets the input file, the directory where the
hills containing files resides, and the frequency to read the other walkers.
Here is an example
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
WALKERS_N=10
WALKERS_ID=3
WALKERS_DIR=../
WALKERS_RSTRIDE=100
... PBMETAD
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
\endplumedfile
where  WALKERS_N is the total number of walkers, WALKERS_ID is the
id of the present walker (starting from 0 ) and the WALKERS_DIR is the directory
where all the walkers are located. WALKERS_RSTRIDE is the number of step between
one update and the other.

*/
//+ENDPLUMEDOC

class PBMetaD : public Bias {

private:
  struct Gaussian {
    std::vector<double> center;
    std::vector<double> sigma;
    double height;
    bool   multivariate; // this is required to discriminate the one dimensional case
    std::vector<double> invsigma;
    Gaussian(const std::vector<double> & center,const std::vector<double> & sigma, double height, bool multivariate):
      center(center),sigma(sigma),height(height),multivariate(multivariate),invsigma(sigma) {
      // to avoid troubles from zero element in flexible hills
      for(unsigned i=0; i<invsigma.size(); ++i)
        if(std::abs(invsigma[i])>1.e-20) {
          invsigma[i]=1.0/invsigma[i] ;
        } else {
          invsigma[i]=0.0;
        }
    }
  };
  // general setup
  double kbt_;
  int stride_;
  // well-tempered MetaD
  bool welltemp_;
  double biasf_;
  // output files format
  std::string fmt_;
  // first step?
  bool isFirstStep_;
  // Gaussian starting parameters
  double height0_;
  std::vector<double> sigma0_;
  std::vector<double> sigma0min_;
  std::vector<double> sigma0max_;
  // Gaussians
  std::vector<std::vector<Gaussian> > hills_;
  std::vector<FlexibleBin> flexbin_;
  int adaptive_;
  std::vector<std::unique_ptr<OFile>> hillsOfiles_;
  std::vector<std::unique_ptr<IFile>> ifiles_;
  std::vector<std::string> ifilesnames_;
  // Grids
  bool grid_;
  std::vector<std::unique_ptr<GridBase>> BiasGrids_;
  std::vector<std::unique_ptr<OFile>> gridfiles_;
  int wgridstride_;
  // Partitioned Families
  unsigned int pf_n_; // initialize number of partitioned families
  std::vector<int> pfs_; //std::vector length of arguments that holds which pf# each cv belongs in
  std::vector<Value*> pfhold_; // std::vector length of pf_n which stores a pointer to the first argument fed to each family
  bool do_pf_; // if partitioned families are enabled
  // multiple walkers
  int mw_n_;
  std::string mw_dir_;
  int mw_id_;
  int mw_rstride_;
  bool walkers_mpi_;
  size_t mpi_nw_;
  unsigned mpi_id_;
  std::vector<std::string> hillsfname_;
  // intervals
  std::vector<double> uppI_;
  std::vector<double> lowI_;
  std::vector<bool>  doInt_;
  // variable for selector
  std::string selector_;
  bool  do_select_;
  unsigned select_value_;
  unsigned current_value_;

  double stretchA=1.0;
  double stretchB=0.0;

  bool noStretchWarningDone=false;

  void noStretchWarning() {
    if(!noStretchWarningDone) {
      log<<"\nWARNING: you are using a HILLS file with Gaussian kernels, PLUMED 2.8 uses stretched Gaussians by default\n";
    }
    noStretchWarningDone=true;
  }

  void   readGaussians(unsigned iarg, IFile*);
  void   writeGaussian(unsigned iarg, const Gaussian&, OFile*);
  void   addGaussian(unsigned iarg, const Gaussian&);
  double getBiasAndDerivatives(unsigned iarg, const std::vector<double>&, double* der=NULL);
  double evaluateGaussian(unsigned iarg, const std::vector<double>&, const Gaussian&,double* der=NULL);
  std::vector<unsigned> getGaussianSupport(unsigned iarg, const Gaussian&);
  bool   scanOneHill(unsigned iarg, IFile *ifile,  std::vector<Value> &v, std::vector<double> &center, std::vector<double>  &sigma, double &height, bool &multivariate);

public:
  explicit PBMetaD(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
  bool checkNeedsGradients()const override;
};

PLUMED_REGISTER_ACTION(PBMetaD,"PBMETAD")

void PBMetaD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition, one for all biases");
  keys.add("optional","FILE","files in which the lists of added hills are stored, default names are assigned using arguments if FILE is not found");
  keys.add("optional","HEIGHT","the height of the Gaussian hills, one for all biases. Compulsory unless TAU, TEMP and BIASFACTOR are given");
  keys.add("optional","FMT","specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics with this bias factor, one for all biases.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","TAU","in well tempered metadynamics, sets height to (k_B Delta T*pace*timestep)/tau");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","GRID_WSTRIDE", "frequency for dumping the grid");
  keys.add("optional","GRID_WFILES", "dump grid for the bias, default names are used if GRID_WSTRIDE is used without GRID_WFILES.");
  keys.add("optional","GRID_RFILES", "read grid for the bias");
  keys.add("optional","ADAPTIVE","use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or timestep dimensions");
  keys.add("optional","SIGMA_MAX","the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","SIGMA_MIN","the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("numbered","PF", "specify which CVs belong in a partitioned family. Once a PF is specified, all CVs in ARG must be placed in a PF even if there is one CV per PF”");
  keys.add("optional","SELECTOR", "add forces and do update based on the value of SELECTOR");
  keys.add("optional","SELECTOR_ID", "value of SELECTOR");
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers - not compatible with WALKERS_* options other than WALKERS_DIR");
  keys.add("optional","INTERVAL_MIN","one dimensional lower limits, outside the limits the system will not feel the biasing force.");
  keys.add("optional","INTERVAL_MAX","one dimensional upper limits, outside the limits the system will not feel the biasing force.");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

PBMetaD::PBMetaD(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao),
  kbt_(0.0),
  stride_(0),
  welltemp_(false),
  biasf_(1.0),
  isFirstStep_(true),
  height0_(std::numeric_limits<double>::max()),
  adaptive_(FlexibleBin::none),
  grid_(false),
  wgridstride_(0),
  pf_n_(0), do_pf_(false),
  mw_n_(1), mw_dir_(""), mw_id_(0), mw_rstride_(1),
  walkers_mpi_(false), mpi_nw_(0),
  do_select_(false) {

  // parse the flexible hills
  std::string adaptiveoption;
  adaptiveoption="NONE";
  parse("ADAPTIVE",adaptiveoption);
  if(adaptiveoption=="GEOM") {
    log.printf("  Uses Geometry-based hills width: sigma must be in distance units and only one sigma is needed\n");
    adaptive_=FlexibleBin::geometry;
  } else if(adaptiveoption=="DIFF") {
    log.printf("  Uses Diffusion-based hills width: sigma must be in time steps and only one sigma is needed\n");
    adaptive_=FlexibleBin::diffusion;
  } else if(adaptiveoption=="NONE") {
    adaptive_=FlexibleBin::none;
  } else {
    error("I do not know this type of adaptive scheme");
  }

  parse("FMT",fmt_);

  // Partitioned Families - fill with -1 to mark as invalid
  pfs_.assign(getNumberOfArguments(), -1);
  pfhold_.resize(getNumberOfArguments());
  std::vector<Value*> familyargs;
  for(int i = 0;; i++) {
    parseArgumentList("PF", i, familyargs);
    if (familyargs.empty()) {
      break;
    }

    do_pf_ = true;
    log << "  Identified Partitioned Family " << i << ":";
    for (unsigned j = 0; j < familyargs.size(); j++) {
      log << " " << familyargs[j]->getName();
      // loop through the argument list to make sure it exists and assign it
      bool foundArg = false;
      for (unsigned argnum = 0; argnum < getNumberOfArguments(); argnum++) {
        if (familyargs[j]->getName() == getPntrToArgument(argnum)->getName()) {
          foundArg = true;
          if (pfs_[argnum] != -1) {
            error(familyargs[j]->getName() + " already present in PF" + std::to_string(pfs_[argnum]));
          }
          pfs_[argnum] = i;  // store the pf# for each cv
          if (pfhold_[i] == nullptr) {
            // if this is the first argument in the family, store a pointer for it (this is for HILLS & GRID files)
            pfhold_[i] = getPntrToArgument(argnum);
          }
        }
      }
      if (!foundArg) {
        error(familyargs[j]->getName() + " in PF" + std::to_string(i) + " not found in ARG");
      }
    }
    log << "\n";
    pf_n_++;
  }

  // if PF were specified, every argument gets treated as its own PF
  if (!do_pf_) {
    pf_n_ = getNumberOfArguments();
    for(unsigned i=0; i < pf_n_; i++) {
      pfhold_[i] = getPntrToArgument(i);
      pfs_[i] = i;
    }
  } else {
    // If we are doing PF, make sure each argument got assigned to a family.
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      if (pfs_[i] == -1) {
        error(getPntrToArgument(i)->getName() + " was not assigned a PF");
      }
    }
  }

  // parse the sigma
  parseVector("SIGMA",sigma0_);
  if(adaptive_==FlexibleBin::none) {
    // if you use normal sigma you need one sigma per argument
    if( sigma0_.size()!=pf_n_ ) {
      std::string fields = do_pf_ ? "PFs" : "arguments";
      error("number of " + fields + " does not match number of SIGMA parameters");
    }
  } else {
    // if you use flexible hills you need one sigma
    if(sigma0_.size()!=1) {
      error("If you choose ADAPTIVE you need only one sigma according to your choice of type (GEOM/DIFF)");
    }
    // if adaptive then the number must be an integer
    if(adaptive_==FlexibleBin::diffusion) {
      if(int(sigma0_[0])-sigma0_[0]>1.e-9 || int(sigma0_[0])-sigma0_[0] <-1.e-9 || int(sigma0_[0])<1 ) {
        error("In case of adaptive hills with diffusion, the sigma must be an integer which is the number of time steps\n");
      }
    }
    // here evtl parse the sigma min and max values
    parseVector("SIGMA_MIN",sigma0min_);
    if(sigma0min_.size()>0 && sigma0min_.size()!=pf_n_) {
      error("the number of SIGMA_MIN values be the same of the number of the arguments/PF");
    } else if(sigma0min_.size()==0) {
      sigma0min_.resize(pf_n_);
      for(unsigned i=0; i<pf_n_; i++) {
        sigma0min_[i]=-1.;
      }
    }

    parseVector("SIGMA_MAX",sigma0max_);
    if(sigma0max_.size()>0 && sigma0max_.size()!=pf_n_) {
      error("the number of SIGMA_MAX values be the same of the number of the arguments/PF");
    } else if(sigma0max_.size()==0) {
      sigma0max_.resize(pf_n_);
      for(unsigned i=0; i<pf_n_; i++) {
        sigma0max_[i]=-1.;
      }
    }

    for(unsigned i=0; i<pf_n_; i++) {
      std::vector<double> tmp_smin, tmp_smax;
      tmp_smin.resize(1,sigma0min_[i]);
      tmp_smax.resize(1,sigma0max_[i]);
      flexbin_.push_back(FlexibleBin(adaptive_,this,i,sigma0_[0],tmp_smin,tmp_smax));
    }
  }

  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  parse("PACE",stride_);
  if(stride_<=0) {
    error("frequency for hill addition is nonsensical");
  }


  parseVector("FILE",hillsfname_);
  if(hillsfname_.size()==0) {
    for(unsigned i=0; i< pf_n_; i++) {
      std::string name = do_pf_ ? "HILLS.PF"+std::to_string(i) : "HILLS."+getPntrToArgument(i)->getName();
      hillsfname_.push_back(name);
    }
  }
  if( hillsfname_.size()!=pf_n_ ) {
    error("number of FILE arguments does not match number of HILLS files");
  }

  parse("BIASFACTOR",biasf_);
  if( biasf_<1.0 ) {
    error("well tempered bias factor is nonsensical");
  }
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) {
    kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  } else {
    kbt_=plumed.getAtoms().getKbT();
  }
  if(biasf_>1.0) {
    if(kbt_==0.0) {
      error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    }
    welltemp_=true;
  }
  double tau=0.0;
  parse("TAU",tau);
  if(tau==0.0) {
    if(height0_==std::numeric_limits<double>::max()) {
      error("At least one between HEIGHT and TAU should be specified");
    }
    // if tau is not set, we compute it here from the other input parameters
    if(welltemp_) {
      tau=(kbt_*(biasf_-1.0))/height0_*getTimeStep()*stride_;
    }
  } else {
    if(!welltemp_) {
      error("TAU only makes sense in well-tempered metadynamics");
    }
    if(height0_!=std::numeric_limits<double>::max()) {
      error("At most one between HEIGHT and TAU should be specified");
    }
    height0_=(kbt_*(biasf_-1.0))/tau*getTimeStep()*stride_;
  }


  // Multiple walkers
  parse("WALKERS_N",mw_n_);
  parse("WALKERS_ID",mw_id_);
  if(mw_n_<=mw_id_) {
    error("walker ID should be a numerical value less than the total number of walkers");
  }
  parse("WALKERS_DIR",mw_dir_);
  parse("WALKERS_RSTRIDE",mw_rstride_);

  // MPI version
  parseFlag("WALKERS_MPI",walkers_mpi_);

  // Grid file
  parse("GRID_WSTRIDE",wgridstride_);
  std::vector<std::string> gridfilenames_;
  parseVector("GRID_WFILES",gridfilenames_);
  if (wgridstride_ == 0 && gridfilenames_.size() > 0) {
    error("frequency with which to output grid not specified use GRID_WSTRIDE");
  }
  if(gridfilenames_.size()==0 && wgridstride_ > 0) {
    for(unsigned i=0; i<pf_n_; i++) {
      std::string name = do_pf_ ? "GRID.PF"+std::to_string(i) : "GRID."+getPntrToArgument(i)->getName();
      gridfilenames_.push_back(name);
    }
  }
  if(gridfilenames_.size() > 0 && hillsfname_.size() > 0 && gridfilenames_.size() != hillsfname_.size()) {
    error("number of GRID_WFILES arguments does not match number of HILLS files");
  }

  // Read grid
  std::vector<std::string> gridreadfilenames_;
  parseVector("GRID_RFILES",gridreadfilenames_);

  // Grid Stuff
  std::vector<std::string> gmin(pf_n_);
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=pf_n_ && gmin.size()!=0) {
    error("not enough values for GRID_MIN");
  }
  std::vector<std::string> gmax(pf_n_);
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=pf_n_ && gmax.size()!=0) {
    error("not enough values for GRID_MAX");
  }
  std::vector<unsigned> gbin(pf_n_);
  std::vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=pf_n_ && gbin.size()!=0) {
    error("not enough values for GRID_BIN");
  }
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=pf_n_ && gspacing.size()!=0) {
    error("not enough values for GRID_SPACING");
  }
  if(gmin.size()!=gmax.size()) {
    error("GRID_MAX and GRID_MIN should be either present or absent");
  }
  if(gspacing.size()!=0 && gmin.size()==0) {
    error("If GRID_SPACING is present also GRID_MIN and GRID_MAX should be present");
  }
  if(gbin.size()!=0     && gmin.size()==0) {
    error("If GRID_BIN is present also GRID_MIN and GRID_MAX should be present");
  }
  if(gmin.size()!=0) {
    if(gbin.size()==0 && gspacing.size()==0) {
      if(adaptive_==FlexibleBin::none) {
        log<<"  Binsize not specified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0_.size()==pf_n_);
        gspacing.resize(pf_n_);
        for(unsigned i=0; i<gspacing.size(); i++) {
          gspacing[i]=0.2*sigma0_[i];
        }
      } else {
        // with adaptive hills and grid a sigma min must be specified
        for(unsigned i=0; i<sigma0min_.size(); i++)
          if(sigma0min_[i]<=0) {
            error("When using ADAPTIVE Gaussians on a grid SIGMA_MIN must be specified");
          }
        log<<"  Binsize not specified, 1/5 of sigma_min will be be used\n";
        gspacing.resize(pf_n_);
        for(unsigned i=0; i<gspacing.size(); i++) {
          gspacing[i]=0.2*sigma0min_[i];
        }
      }
    } else if(gspacing.size()!=0 && gbin.size()==0) {
      log<<"  The number of bins will be estimated from GRID_SPACING\n";
    } else if(gspacing.size()!=0 && gbin.size()!=0) {
      log<<"  You specified both GRID_BIN and GRID_SPACING\n";
      log<<"  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if(gbin.size()==0) {
      gbin.assign(pf_n_,1);
    }
    if(gspacing.size()!=0)
      for(unsigned i=0; i<pf_n_; i++) {
        double a,b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        unsigned n=std::ceil(((b-a)/gspacing[i]));
        if(gbin[i]<n) {
          gbin[i]=n;
        }
      }
  }
  if(gbin.size()>0) {
    grid_=true;
  }

  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  if(!grid_&&gridfilenames_.size() > 0) {
    error("To write a grid you need first to define it!");
  }
  if(!grid_&&gridreadfilenames_.size() > 0) {
    error("To read a grid you need first to define it!");
  }

  doInt_.resize(pf_n_,false);
  // Interval keyword
  parseVector("INTERVAL_MIN",lowI_);
  parseVector("INTERVAL_MAX",uppI_);
  // various checks
  if(lowI_.size()!=uppI_.size()) {
    error("both a lower and an upper limits must be provided with INTERVAL");
  }
  if(lowI_.size()!=0 && lowI_.size()!=pf_n_) {
    error("check number of argument of INTERVAL");
  }
  for(unsigned i=0; i<lowI_.size(); ++i) {
    if(uppI_[i]<lowI_[i]) {
      error("The Upper limit must be greater than the Lower limit!");
    }
    if(pfhold_[i]->isPeriodic()) {
      warning("INTERVAL is not used for periodic variables");
    } else {
      doInt_[i]=true;
    }
  }

  // parse selector stuff
  parse("SELECTOR", selector_);
  if(selector_.length()>0) {
    do_select_ = true;
    select_value_ = 0; // set defalt value or it might be not initialized if the user does not pass SELECTOR_ID
    parse("SELECTOR_ID", select_value_);
  }

  checkRead();

  log.printf("  Gaussian width ");
  if (adaptive_==FlexibleBin::diffusion) {
    log.printf(" (Note: The units of sigma are in timesteps) ");
  }
  if (adaptive_==FlexibleBin::geometry) {
    log.printf(" (Note: The units of sigma are in dist units) ");
  }
  for(unsigned i=0; i<sigma0_.size(); ++i) {
    log.printf(" %f",sigma0_[i]);
  }
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_);
  log.printf("  Gaussian files ");
  for(unsigned i=0; i<hillsfname_.size(); ++i) {
    log.printf("%s ",hillsfname_[i].c_str());
  }
  log.printf("\n");
  if(welltemp_) {
    log.printf("  Well-Tempered Bias Factor %f\n",biasf_);
    log.printf("  Hills relaxation time (tau) %f\n",tau);
    log.printf("  KbT %f\n",kbt_);
  }

  if(do_select_) {
    log.printf("  Add forces and update bias based on the value of SELECTOR %s\n",selector_.c_str());
    log.printf("  Id of the SELECTOR for this action %u\n", select_value_);
  }

  if(mw_n_>1) {
    if(walkers_mpi_) {
      error("MPI version of multiple walkers is not compatible with filesystem version of multiple walkers");
    }
    log.printf("  %d multiple walkers active\n",mw_n_);
    log.printf("  walker id %d\n",mw_id_);
    log.printf("  reading stride %d\n",mw_rstride_);
    if(mw_dir_!="") {
      log.printf("  directory with hills files %s\n",mw_dir_.c_str());
    }
  } else {
    if(walkers_mpi_) {
      log.printf("  Multiple walkers active using MPI communnication\n");
      if(mw_dir_!="") {
        log.printf("  directory with hills files %s\n",mw_dir_.c_str());
      }
      if(comm.Get_rank()==0) {
        // Only root of group can communicate with other walkers
        mpi_nw_ = multi_sim_comm.Get_size();
        mpi_id_ = multi_sim_comm.Get_rank();
      }
      // Communicate to the other members of the same group
      // info abount number of walkers and walker index
      comm.Bcast(mpi_nw_,0);
      comm.Bcast(mpi_id_,0);
    }
  }

  for(unsigned i=0; i<doInt_.size(); i++) {
    if(doInt_[i]) {
      log.printf("  Upper and Lower limits boundaries for the bias of CV %u are activated\n", i);
    }
  }
  if(grid_) {
    log.printf("  Grid min");
    for(unsigned i=0; i<gmin.size(); ++i) {
      log.printf(" %s",gmin[i].c_str() );
    }
    log.printf("\n");
    log.printf("  Grid max");
    for(unsigned i=0; i<gmax.size(); ++i) {
      log.printf(" %s",gmax[i].c_str() );
    }
    log.printf("\n");
    log.printf("  Grid bin");
    for(unsigned i=0; i<gbin.size(); ++i) {
      log.printf(" %u",gbin[i]);
    }
    log.printf("\n");
    if(spline) {
      log.printf("  Grid uses spline interpolation\n");
    }
    if(sparsegrid) {
      log.printf("  Grid uses sparse grid\n");
    }
    if(wgridstride_>0) {
      for(unsigned i=0; i<gridfilenames_.size(); ++i) {
        log.printf("  Grid is written on file %s with stride %d\n",gridfilenames_[i].c_str(),wgridstride_);
      }
    }
    if(gridreadfilenames_.size()>0) {
      for(unsigned i=0; i<gridreadfilenames_.size(); ++i) {
        log.printf("  Reading bias from grid in file %s \n",gridreadfilenames_[i].c_str());
      }
    }
  }

  // initializing vector of hills
  hills_.resize(pf_n_);

  // restart from external grid
  bool restartedFromGrid=false;

  // initializing and checking grid
  if(grid_) {
    // check for mesh and sigma size
    for(unsigned i=0; i<pf_n_; i++) {
      double a,b;
      int family = pfs_[i]; // point to families instead of arguments
      Tools::convert(gmin[family],a);
      Tools::convert(gmax[family],b);
      double mesh=(b-a)/((double)gbin[family]);
      if(adaptive_==FlexibleBin::none) {
        if(mesh>0.5*sigma0_[i]) {
          log<<"  WARNING: Using a PBMETAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
        }
      } else {
        if(mesh>0.5*sigma0min_[i]||sigma0min_[i]<0.) {
          log<<"  WARNING: to use a PBMETAD with a GRID and ADAPTIVE you need to set a Grid Spacing larger than half of the Gaussians \n";
        }
      }
    }
    std::string funcl=getLabel() + ".bias";
    for(unsigned i=0; i<pf_n_; ++i) {
      std::vector<Value*> args(1);
      args[0] = pfhold_[i];  //Use first argument in family for interactions.
      std::vector<std::string> gmin_t(1);
      std::vector<std::string> gmax_t(1);
      std::vector<unsigned>    gbin_t(1);
      gmin_t[0] = gmin[i];
      gmax_t[0] = gmax[i];
      gbin_t[0] = gbin[i];
      std::unique_ptr<GridBase> BiasGrid_;
      // Read grid from file
      if(gridreadfilenames_.size()>0) {
        IFile gridfile;
        gridfile.link(*this);
        if(gridfile.FileExist(gridreadfilenames_[i])) {
          gridfile.open(gridreadfilenames_[i]);
        } else {
          error("The GRID file you want to read: " + gridreadfilenames_[i] + ", cannot be found!");
        }
        std::string funcl = getLabel() + ".bias";
        BiasGrid_=GridBase::create(funcl, args, gridfile, gmin_t, gmax_t, gbin_t, sparsegrid, spline, true);
        if(BiasGrid_->getDimension() != args.size()) {
          error("mismatch between dimensionality of input grid and number of arguments");
        }
        if(pfhold_[i]->isPeriodic() != BiasGrid_->getIsPeriodic()[0]) {
          error("periodicity mismatch between arguments and input bias");
        }
        log.printf("  Restarting from %s:\n",gridreadfilenames_[i].c_str());
        if(getRestart()) {
          restartedFromGrid=true;
        }
      } else {
        if(!sparsegrid) {
          BiasGrid_=Tools::make_unique<Grid>(funcl,args,gmin_t,gmax_t,gbin_t,spline,true);
        } else           {
          BiasGrid_=Tools::make_unique<SparseGrid>(funcl,args,gmin_t,gmax_t,gbin_t,spline,true);
        }
        std::vector<std::string> actualmin=BiasGrid_->getMin();
        std::vector<std::string> actualmax=BiasGrid_->getMax();
        std::string is;
        Tools::convert(i,is);
        if(gmin_t[0]!=actualmin[0]) {
          error("GRID_MIN["+is+"] must be adjusted to "+actualmin[0]+" to fit periodicity");
        }
        if(gmax_t[0]!=actualmax[0]) {
          error("GRID_MAX["+is+"] must be adjusted to "+actualmax[0]+" to fit periodicity");
        }
      }
      BiasGrids_.emplace_back(std::move(BiasGrid_));
    }
  }



// creating vector of ifile* for hills reading
// open all files at the beginning and read Gaussians if restarting

  for(int j=0; j<mw_n_; ++j) {
    for(unsigned i=0; i<hillsfname_.size(); ++i) {
      unsigned k=j*hillsfname_.size()+i;
      std::string fname;
      if(mw_dir_!="") {
        if(mw_n_>1) {
          std::stringstream out;
          out << j;
          fname = mw_dir_+"/"+hillsfname_[i]+"."+out.str();
        } else if(walkers_mpi_) {
          fname = mw_dir_+"/"+hillsfname_[i];
        } else {
          fname = hillsfname_[i];
        }
      } else {
        if(mw_n_>1) {
          std::stringstream out;
          out << j;
          fname = hillsfname_[i]+"."+out.str();
        } else {
          fname = hillsfname_[i];
        }
      }
      ifiles_.emplace_back(Tools::make_unique<IFile>());
      // this is just a shortcut pointer to the last element:
      IFile *ifile = ifiles_.back().get();
      ifile->link(*this);
      ifilesnames_.push_back(fname);
      if(ifile->FileExist(fname)) {
        ifile->open(fname);
        if(getRestart()&&!restartedFromGrid) {
          log.printf("  Restarting from %s:",ifilesnames_[k].c_str());
          readGaussians(i,ifiles_[k].get());
        }
        ifiles_[k]->reset(false);
        // close only the walker own hills file for later writing
        if(j==mw_id_) {
          ifiles_[k]->close();
        }
      } else {
        // in case a file does not exist and we are restarting, complain that the file was not found
        if(getRestart()) {
          log<<"  WARNING: restart file "<<fname<<" not found\n";
        }
      }
    }
  }

  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi_) {
    multi_sim_comm.Barrier();
  }

  // open hills files for writing
  for(unsigned i=0; i<hillsfname_.size(); ++i) {
    auto ofile=Tools::make_unique<OFile>();
    ofile->link(*this);
    // if MPI multiple walkers, only rank 0 will write to file
    if(walkers_mpi_) {
      int r=0;
      if(comm.Get_rank()==0) {
        r=multi_sim_comm.Get_rank();
      }
      comm.Bcast(r,0);
      if(r>0) {
        ifilesnames_[mw_id_*hillsfname_.size()+i]="/dev/null";
      }
      ofile->enforceSuffix("");
    }
    if(mw_n_>1) {
      ofile->enforceSuffix("");
    }
    ofile->open(ifilesnames_[mw_id_*hillsfname_.size()+i]);
    if(fmt_.length()>0) {
      ofile->fmtField(fmt_);
    }
    ofile->addConstantField("multivariate");
    ofile->addConstantField("kerneltype");
    if(doInt_[i]) {
      ofile->addConstantField("lower_int").printField("lower_int",lowI_[i]);
      ofile->addConstantField("upper_int").printField("upper_int",uppI_[i]);
    }
    ofile->setHeavyFlush();
    // output periodicities of variables
    ofile->setupPrintValue( pfhold_[i] );  //assuming cvs in the same family have the same periodicity and boundaries.
    // push back
    hillsOfiles_.emplace_back(std::move(ofile));
  }

  // Dump grid to files
  if(wgridstride_ > 0) {
    for(unsigned i = 0; i < gridfilenames_.size(); ++i) {
      auto ofile=Tools::make_unique<OFile>();
      ofile->link(*this);
      std::string gridfname_tmp = gridfilenames_[i];
      if(walkers_mpi_) {
        int r = 0;
        if(comm.Get_rank() == 0) {
          r = multi_sim_comm.Get_rank();
        }
        comm.Bcast(r, 0);
        if(r>0) {
          gridfname_tmp = "/dev/null";
        }
        ofile->enforceSuffix("");
      }
      if(mw_n_>1) {
        ofile->enforceSuffix("");
      }
      ofile->open(gridfname_tmp);
      ofile->setHeavyFlush();
      gridfiles_.emplace_back(std::move(ofile));
    }
  }

  log<<"  Bibliography "<<plumed.cite("Pfaendtner and Bonomi. J. Chem. Theory Comput. 11, 5062 (2015)");
  if(doInt_[0])
    log<<plumed.cite(
         "Baftizadeh, Cossio, Pietrucci, and Laio, Curr. Phys. Chem. 2, 79 (2012)");
  if(mw_n_>1||walkers_mpi_)
    log<<plumed.cite(
         "Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");
  if(adaptive_!=FlexibleBin::none)
    log<<plumed.cite(
         "Branduardi, Bussi, and Parrinello, J. Chem. Theory Comput. 8, 2247 (2012)");
  if (do_pf_) {
    log<<plumed.cite("Prakash, Fu, Bonomi, and Pfaendtner, J. Chem. Theory Comput. 14, 4985 (2018)");
  }
  log<<"\n";


}

void PBMetaD::readGaussians(unsigned iarg, IFile *ifile) {
  std::vector<double> center(1);
  std::vector<double> sigma(1);
  double height;
  int nhills=0;
  bool multivariate=false;
  int family=pfs_[iarg];

  std::vector<Value> tmpvalues;
  tmpvalues.push_back( Value( this, pfhold_[family]->getName(), false ) );

  while(scanOneHill(iarg,ifile,tmpvalues,center,sigma,height,multivariate)) {
    ;
    nhills++;
    if(welltemp_) {
      height*=(biasf_-1.0)/biasf_;
    }
    addGaussian(family, Gaussian(center,sigma,height,multivariate));
  }
  log.printf("      %d Gaussians read\n",nhills);
}

void PBMetaD::writeGaussian(unsigned iarg, const Gaussian& hill, OFile *ofile) {
  int family=pfs_[iarg];
  ofile->printField("time",getTimeStep()*getStep());
  ofile->printField(pfhold_[family],hill.center[0]);

  ofile->printField("kerneltype","stretched-gaussian");
  if(hill.multivariate) {
    ofile->printField("multivariate","true");
    double lower = std::sqrt(1./hill.sigma[0]);
    ofile->printField("sigma_"+pfhold_[family]->getName()+"_"+
                      pfhold_[family]->getName(),lower);
  } else {
    ofile->printField("multivariate","false");
    ofile->printField("sigma_"+pfhold_[family]->getName(),hill.sigma[0]);
  }
  double height=hill.height;
  if(welltemp_) {
    height *= biasf_/(biasf_-1.0);
  }
  ofile->printField("height",height);
  ofile->printField("biasf",biasf_);
  if(mw_n_>1) {
    ofile->printField("clock",int(std::time(0)));
  }
  ofile->printField();
}

void PBMetaD::addGaussian(unsigned iarg, const Gaussian& hill) {
  if(!grid_) {
    hills_[iarg].push_back(hill);
  } else {
    std::vector<unsigned> nneighb=getGaussianSupport(iarg, hill);
    std::vector<Grid::index_t> neighbors=BiasGrids_[iarg]->getNeighbors(hill.center,nneighb);
    std::vector<double> der(1);
    std::vector<double> xx(1);
    if(comm.Get_size()==1) {
      for(unsigned i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        der[0]=0.0;
        BiasGrids_[iarg]->getPoint(ineigh,xx);
        double bias=evaluateGaussian(iarg,xx,hill,&der[0]);
        BiasGrids_[iarg]->addValueAndDerivatives(ineigh,bias,der);
      }
    } else {
      unsigned stride=comm.Get_size();
      unsigned rank=comm.Get_rank();
      std::vector<double> allder(neighbors.size(),0.0);
      std::vector<double> allbias(neighbors.size(),0.0);
      for(unsigned i=rank; i<neighbors.size(); i+=stride) {
        Grid::index_t ineigh=neighbors[i];
        BiasGrids_[iarg]->getPoint(ineigh,xx);
        allbias[i]=evaluateGaussian(iarg,xx,hill,&allder[i]);
      }
      comm.Sum(allbias);
      comm.Sum(allder);
      for(unsigned i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        der[0]=allder[i];
        BiasGrids_[iarg]->addValueAndDerivatives(ineigh,allbias[i],der);
      }
    }
  }
}

std::vector<unsigned> PBMetaD::getGaussianSupport(unsigned iarg, const Gaussian& hill) {
  std::vector<unsigned> nneigh;
  double cutoff;
  if(hill.multivariate) {
    double maxautoval=1./hill.sigma[0];
    cutoff=std::sqrt(2.0*dp2cutoff*maxautoval);
  } else {
    cutoff=std::sqrt(2.0*dp2cutoff)*hill.sigma[0];
  }

  if(doInt_[iarg]) {
    if(hill.center[0]+cutoff > uppI_[iarg] || hill.center[0]-cutoff < lowI_[iarg]) {
      // in this case, we updated the entire grid to avoid problems
      return BiasGrids_[iarg]->getNbin();
    } else {
      nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrids_[iarg]->getDx()[0])));
      return nneigh;
    }
  }

  nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrids_[iarg]->getDx()[0])) );

  return nneigh;
}

double PBMetaD::getBiasAndDerivatives(unsigned iarg, const std::vector<double>& cv, double* der) {
  double bias=0.0;
  int family = pfs_[iarg];
  if(!grid_) {
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    for(unsigned i=rank; i<hills_[family].size(); i+=stride) {
      bias += evaluateGaussian(iarg,cv,hills_[family][i],der);
    }
    comm.Sum(bias);
    if(der) {
      comm.Sum(der,1);
    }
  } else {
    if(der) {
      std::vector<double> vder(1);
      bias = BiasGrids_[family]->getValueAndDerivatives(cv,vder);
      der[0] = vder[0];
    } else {
      bias = BiasGrids_[family]->getValue(cv);
    }
  }

  return bias;
}

double PBMetaD::evaluateGaussian(unsigned iarg, const std::vector<double>& cv, const Gaussian& hill, double* der) {
  double bias=0.0;
// I use a pointer here because cv is const (and should be const)
// but when using doInt it is easier to locally replace cv[0] with
// the upper/lower limit in case it is out of range
  const double *pcv=NULL;
  double tmpcv[1]; // tmp array with cv (to be used with doInt_)
  tmpcv[0]=cv[0];
  bool isOutOfInt = false;
  if(doInt_[iarg]) {
    if(cv[0]<lowI_[iarg]) {
      tmpcv[0]=lowI_[iarg];
      isOutOfInt = true;
    } else if(cv[0]>uppI_[iarg]) {
      tmpcv[0]=uppI_[iarg];
      isOutOfInt = true;
    }
  }
  pcv=&(tmpcv[0]);

  if(hill.multivariate) {
    double dp  = difference(iarg, hill.center[0], pcv[0]);
    double dp2 = 0.5 * dp * dp * hill.sigma[0];
    if(dp2<dp2cutoff) {
      bias = hill.height*std::exp(-dp2);
      if(der && !isOutOfInt) {
        der[0] += -bias * dp * hill.sigma[0] * stretchA;
      }
      bias=stretchA*bias+hill.height*stretchB;
    }
  } else {
    double dp  = difference(iarg, hill.center[0], pcv[0]) * hill.invsigma[0];
    double dp2 = 0.5 * dp * dp;
    if(dp2<dp2cutoff) {
      bias = hill.height*std::exp(-dp2);
      if(der && !isOutOfInt) {
        der[0] += -bias * dp * hill.invsigma[0] * stretchA;
      }
      bias=stretchA*bias+hill.height*stretchB;
    }
  }

  return bias;
}

void PBMetaD::calculate() {
  // this is because presently there is no way to properly pass information
  // on adaptive hills (diff) after exchanges:
  if(adaptive_==FlexibleBin::diffusion && getExchangeStep()) {
    error("ADAPTIVE=DIFF is not compatible with replica exchange");
  }

  std::vector<double> cv(1);
  double der[1];
  std::vector<double> bias(getNumberOfArguments());
  std::vector<double> deriv(getNumberOfArguments());

  double ncv = (double) getNumberOfArguments();
  double bmin = 1.0e+19;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    cv[0]    = getArgument(i);
    der[0]   = 0.0;
    bias[i]  = getBiasAndDerivatives(i, cv, der);
    deriv[i] = der[0];
    if(bias[i] < bmin) {
      bmin = bias[i];
    }
  }
  double ene = 0.;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ene += std::exp((-bias[i]+bmin)/kbt_);
  }

  // set Forces - set them to zero if SELECTOR is active
  if(do_select_) {
    current_value_ = static_cast<unsigned>(plumed.passMap[selector_]);
  }

  if(!do_select_ || select_value_==current_value_) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      const double f = - std::exp((-bias[i]+bmin)/kbt_) / (ene) * deriv[i];
      setOutputForce(i, f);
    }
  }

  if(do_select_ && select_value_!=current_value_) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      setOutputForce(i, 0.0);
    }
  }

  // set bias
  ene = -kbt_ * (std::log(ene) - std::log(ncv)) + bmin;
  setBias(ene);
}

void PBMetaD::update() {
  bool multivariate;
  // adding hills criteria
  bool nowAddAHill;
  if(getStep()%stride_==0 && !isFirstStep_) {
    nowAddAHill=true;
  } else {
    nowAddAHill=false;
    isFirstStep_=false;
  }

  // if you use adaptive, call the FlexibleBin
  if(adaptive_!=FlexibleBin::none) {
    for(unsigned i=0; i<getNumberOfArguments(); i++) {
      flexbin_[i].update(nowAddAHill,i);
    }
    multivariate=true;
  } else {
    multivariate=false;
  }

  if(nowAddAHill && (!do_select_ || select_value_==current_value_)) {
    // get all biases and heights
    std::vector<double> cv(getNumberOfArguments());
    std::vector<double> bias(getNumberOfArguments());
    std::vector<double> thissigma(getNumberOfArguments());
    std::vector<double> height(getNumberOfArguments());
    std::vector<double> cv_tmp(1);
    std::vector<double> sigma_tmp(1);
    double norm = 0.0;
    double bmin = 1.0e+19;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      int family=pfs_[i];
      // get flex/sigmas for each family and assign them to this args sigma
      if(adaptive_!=FlexibleBin::none) {
        thissigma[i]=flexbin_[family].getInverseMatrix(i)[0];
      } else {
        thissigma[i]=sigma0_[family];
      }
      cv[i]     = getArgument(i);
      cv_tmp[0] = getArgument(i);
      bias[i] = getBiasAndDerivatives(i, cv_tmp);
      if(bias[i] < bmin) {
        bmin = bias[i];
      }
    }
    // calculate heights and norm
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      double h = std::exp((-bias[i]+bmin)/kbt_);
      norm += h;
      height[i] = h;
    }
    // normalize and apply welltemp correction
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      height[i] *=  height0_ / norm;
      if(welltemp_) {
        height[i] *= std::exp(-bias[i]/(kbt_*(biasf_-1.0)));
      }
    }

    // MPI Multiple walkers: share hills and add them all
    if(walkers_mpi_) {
      // Allocate arrays to store all walkers hills
      std::vector<double> all_cv(mpi_nw_*cv.size(), 0.0);
      std::vector<double> all_sigma(mpi_nw_*getNumberOfArguments(), 0.0);
      std::vector<double> all_height(mpi_nw_*height.size(), 0.0);
      if(comm.Get_rank()==0) {
        // fill in value
        for(unsigned i=0; i<getNumberOfArguments(); ++i) {
          unsigned j = mpi_id_ * getNumberOfArguments() + i;
          all_cv[j] = cv[i];
          all_sigma[j]  = thissigma[i];
          all_height[j] = height[i];
        }
        // Communicate (only root)
        multi_sim_comm.Sum(&all_cv[0], all_cv.size());
        multi_sim_comm.Sum(&all_sigma[0], all_sigma.size());
        multi_sim_comm.Sum(&all_height[0], all_height.size());
      }
      // Share info with group members
      comm.Sum(&all_cv[0], all_cv.size());
      comm.Sum(&all_sigma[0], all_sigma.size());
      comm.Sum(&all_height[0], all_height.size());
      // now add hills one by one
      for(unsigned j=0; j<mpi_nw_; ++j) {
        for(unsigned i=0; i<getNumberOfArguments(); ++i) {
          // Add CVs of same family together and write to same file
          int family = pfs_[i];
          cv_tmp[0]    = all_cv[j*cv.size()+i];
          double height_tmp = all_height[j*cv.size()+i];
          sigma_tmp[0] = all_sigma[j*cv.size()+i];
          Gaussian newhill = Gaussian(cv_tmp, sigma_tmp, height_tmp, multivariate);
          addGaussian(family, newhill);
          writeGaussian(i, newhill, hillsOfiles_[family].get());
        }
      }
      // just add your own hills
    } else {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        // Add CVs of same family together and write to same file
        int family = pfs_[i];
        cv_tmp[0] = cv[i];
        if(adaptive_!=FlexibleBin::none) {
          sigma_tmp[0]=thissigma[i];
        } else {
          sigma_tmp[0] = sigma0_[family];
        }
        Gaussian newhill = Gaussian(cv_tmp, sigma_tmp, height[i], multivariate);
        addGaussian(family, newhill);
        writeGaussian(i, newhill, hillsOfiles_[family].get());
      }
    }
  }

  // write grid files
  if(wgridstride_>0 && (getStep()%wgridstride_==0 || getCPT())) {
    int r = 0;
    if(walkers_mpi_) {
      if(comm.Get_rank()==0) {
        r=multi_sim_comm.Get_rank();
      }
      comm.Bcast(r,0);
    }
    if(r==0) {
      for(unsigned i=0; i<gridfiles_.size(); ++i) {
        gridfiles_[i]->rewind();
        BiasGrids_[i]->writeToFile(*gridfiles_[i]);
        gridfiles_[i]->flush();
      }
    }
  }

  // if multiple walkers and time to read Gaussians
  if(mw_n_>1 && getStep()%mw_rstride_==0) {
    for(int j=0; j<mw_n_; ++j) {
      for(unsigned i=0; i<hillsfname_.size(); ++i) {
        unsigned k=j*hillsfname_.size()+i;
        // don't read your own Gaussians
        if(j==mw_id_) {
          continue;
        }
        // if the file is not open yet
        if(!(ifiles_[k]->isOpen())) {
          // check if it exists now and open it!
          if(ifiles_[k]->FileExist(ifilesnames_[k])) {
            ifiles_[k]->open(ifilesnames_[k]);
            ifiles_[k]->reset(false);
          }
          // otherwise read the new Gaussians
        } else {
          log.printf("  Reading hills from %s:",ifilesnames_[k].c_str());
          readGaussians(i,ifiles_[k].get());
          ifiles_[k]->reset(false);
        }
      }
    }
  }

}

/// takes a pointer to the file and a template string with values v and gives back the next center, sigma and height
bool PBMetaD::scanOneHill(unsigned iarg, IFile *ifile, std::vector<Value> &tmpvalues, std::vector<double> &center, std::vector<double> &sigma, double &height, bool &multivariate) {
  double dummy;
  multivariate=false;
  Value* argPtr = pfhold_[pfs_[iarg]];
  if(ifile->scanField("time",dummy)) {
    ifile->scanField( &tmpvalues[0] );
    if( tmpvalues[0].isPeriodic() && ! argPtr->isPeriodic() ) {
      error("in hills file periodicity for variable " + tmpvalues[0].getName() + " does not match periodicity in input");
    } else if( tmpvalues[0].isPeriodic() ) {
      std::string imin, imax;
      tmpvalues[0].getDomain( imin, imax );
      std::string rmin, rmax;
      argPtr->getDomain( rmin, rmax );
      if( imin!=rmin || imax!=rmax ) {
        error("in hills file periodicity for variable " + tmpvalues[0].getName() + " does not match periodicity in input");
      }
    }
    center[0]=tmpvalues[0].get();
    std::string ktype="stretched-gaussian";
    if( ifile->FieldExist("kerneltype") ) {
      ifile->scanField("kerneltype",ktype);
    }

    if( ktype=="gaussian" ) {
      noStretchWarning();
    } else if( ktype!="stretched-gaussian") {
      error("non Gaussian kernels are not supported in MetaD");
    }

    std::string sss;
    ifile->scanField("multivariate",sss);
    if(sss=="true") {
      multivariate=true;
    } else if(sss=="false") {
      multivariate=false;
    } else {
      plumed_merror("cannot parse multivariate = "+ sss);
    }
    if(multivariate) {
      ifile->scanField("sigma_"+argPtr->getName()+"_"+
                       argPtr->getName(),sigma[0]);
      sigma[0] = 1./(sigma[0]*sigma[0]);
    } else {
      ifile->scanField("sigma_"+argPtr->getName(),sigma[0]);
    }
    ifile->scanField("height",height);
    ifile->scanField("biasf",dummy);
    if(ifile->FieldExist("clock")) {
      ifile->scanField("clock",dummy);
    }
    if(ifile->FieldExist("lower_int")) {
      ifile->scanField("lower_int",dummy);
    }
    if(ifile->FieldExist("upper_int")) {
      ifile->scanField("upper_int",dummy);
    }
    ifile->scanField();
    return true;
  } else {
    return false;
  }

}

bool PBMetaD::checkNeedsGradients()const {
  if(adaptive_==FlexibleBin::geometry) {
    if(getStep()%stride_==0 && !isFirstStep_) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

}
}

