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
#include "MetaDGrid.h"
#include "MetaDHills.h"
#include "MetaDKernel.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/FlexibleBin.h"
#include "tools/Exception.h"
#include "tools/Grid.h"
#include "tools/Matrix.h"
#include "tools/OpenMP.h"
#include "tools/Random.h"
#include "tools/File.h"
#include "tools/Communicator.h"
#include <algorithm>
#include <limits>
#include <sstream>

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS PBMETAD
/*
Used to performed Parallel Bias metadynamics.

This action provides an implementation of the Parallel Bias Metadynamics (PBMetaD) method that was introduced in the first paper cited below.
This variant of metadynamics (see [METAD](METAD.md)) uses multiple low-dimensional bias potentials that are applied in parallel.
In the current implementation, these have the form of mono-dimensional metadynamics bias
potentials:

$$
{V(s_1,t), ..., V(s_N,t)}
$$

where:

$$
V(s_i,t) = \sum_{ k \tau < t} W_i(k \tau)
\exp\left(
- \frac{(s_i-s_i^{(0)}(k \tau))^2}{2\sigma_i^2}
\right).
$$

To ensure the convergence of each mono-dimensional bias potential to the corresponding free energy,
at each deposition step the Gaussian heights are multiplied by the so-called conditional term:

$$
W_i(k \tau)=W_0 \frac{\exp\left(
- \frac{V(s_i,k \tau)}{k_B T}
\right)}{\sum_{i=1}^N
\exp\left(
- \frac{V(s_i,k \tau)}{k_B T}
\right)}
$$

where $W_0$ is the initial Gaussian height.

The PBMetaD bias potential is defined by:

$$
V_{PB}(\vec{s},t) = -k_B T \log{\sum_{i=1}^N
\exp\left(
- \frac{V(s_i,t)}{k_B T}
\right)}.
$$


Information on the Gaussian functions that build each bias potential are printed to
multiple HILLS files, which
are used both to restart the calculation and to reconstruct the mono-dimensional
free energies as a function of the corresponding CVs.
These can be reconstructed using the [sum_hills](sum_hills.md) utility because the final bias is given
by:

$$
V(s_i) = -F(s_i)
$$

Currently, only a subset of the [METAD](METAD.md) options are available in PBMetaD.

The bias potentials can be stored on a grid to increase performances of long PBMetaD simulations.
You should
provide either the number of bins for every collective variable (GRID_BIN) or
the desired grid spacing (GRID_SPACING). In case you provide both PLUMED will use
the most conservative choice (highest number of bins) for each dimension.
In case you do not provide any information about bin size (neither GRID_BIN nor GRID_SPACING)
and if Gaussian width is fixed PLUMED will use 1/5 of the Gaussian width as grid spacing.
This default choice should be reasonable for most applications.

As in [METAD](METAD.md), the grids are extended automatically, along the variables that are not periodic, whenever a
Gaussian does not fit in them: a grid always contains all the points where a Gaussian is not zero, so that the bias is
exactly zero outside it (with INTERVAL, the grid contains the interval and the bias beyond the grid is the constant value
beyond the limits).  The grid spacing does not change and the new points are on the same lattice as the initial ones.
Every extension is written on the log.  If you give GRID_SPACING without GRID_MIN and GRID_MAX each grid is placed
around its first Gaussian and then extended as needed.  GRID_SIZE_LIMIT sets the maximum number of points of each grid.
A grid read with GRID_RFILES can have any range, but it must be on the lattice of the grid given in input.

As in [METAD](METAD.md), no Gaussian is added at the first step of a simulation.  When a simulation is continued with
RESTART and the bias is read from the hills files, a Gaussian is added at the first step unless they already contain
it (i.e. the previous run had already performed this step).  With ADAPTIVE=DIFF and replica exchange, the running averages of the variables are
restarted when an exchange is accepted, and the variances are kept.

Another option that is available is the well-tempered metadynamics that is discussed in the second paper cited below. In this
variant of PBMetaD the heights of the Gaussian hills are scaled at each step by the
additional well-tempered metadynamics term.
This  ensures that each bias converges more smoothly. It should be noted that, in the case of well-tempered metadynamics, in
the output printed the Gaussian height is re-scaled using the bias factor.
Also notice that with well-tempered metadynamics the HILLS files do not contain the bias,
but the negative of the free-energy estimate. This choice has the advantage that
one can restart a simulation using a different value for the $\Delta T$. The applied bias will be scaled accordingly.

Note that you can also use here the flexible Gaussian approach that is discussed in the third paper cited below.
This methods allows you to adapt the Gaussian to the extent of Cartesian space covered by a variable or
to the space in collective variable covered in a given time. In this case the width of the deposited
Gaussian potential is denoted by one value only that is a Cartesian space (ADAPTIVE=GEOM) or a time
(ADAPTIVE=DIFF). Note that a specific integration technique for the deposited Gaussian kernels
should be used in this case. Check the documentation for utility sum_hills.

With the keyword INTERVAL one changes the metadynamics algorithm and sets the bias force equal to zero
outside boundary as discussed in the fourth paper cited below. If, for example, metadynamics is performed on a CV s and one is interested only
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
  be used together with a [UPPER_WALLS](UPPER_WALLS.md) or [LOWER_WALLS](LOWER_WALLS.md) at boundary.

For systems with multiple CVs that share identical properties, PBMetaD with partitioned families can be used
to group them under one bias potential that each contributes to as discussed in the penultimate paper in the references below. This is done with a list
of PF keywords, where each PF* argument contains the list of CVs from ARG to be placed in that family. Once
invoked, each CV in ARG must be placed in exactly one PF, even if it results in families containing only one CV.
Additionally, in cases where each of SIGMA or GRID entry would correspond to each ARG entry, they now correspond to
each PF and must be adjusted accordingly.

The multiple walkers that is discussed in the final paper cited below can also be used. See below the examples.

## Examples

The following input is for PBMetaD calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the PBMetaD bias potential are written to the COLVAR file every 100 steps.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 LABEL=pb FILE=HILLS_d1,HILLS_d2
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
```

(See also [DISTANCE](DISTANCE.md) and [PRINT](PRINT.md)).

If you use well-tempered metadynamics, you should specify a single bias factor and initial Gaussian height.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
... PBMETAD
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
```

Using partitioned families, each CV in ARG must be in exactly one family. Here,
the distance between atoms 1,2 is degenerate with 2,4, but not with the
distance between 3,5. Note that two SIGMA are provided to match the two PF.

```plumed
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
```

The following input enables the MPI version of multiple-walkers.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 LABEL=pb
FILE=HILLS_d1,HILLS_d2
WALKERS_MPI
... PBMETAD
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
```

The disk version of multiple-walkers can be
enabled by setting the number of walker used, the id of the
current walker which interprets the input file, the directory where the
hills containing files resides, and the frequency to read the other walkers.
Here is an example

```plumed
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
```

where  WALKERS_N is the total number of walkers, WALKERS_ID is the
id of the present walker (starting from 0 ) and the WALKERS_DIR is the directory
where all the walkers are located. WALKERS_RSTRIDE is the number of step between
one update and the other.

*/
//+ENDPLUMEDOC

class PBMetaD : public Bias {

private:
  using Gaussian = MetaDGaussian;
  /// A partitioned family: its arguments share a one dimensional bias, with its hills file and grid
  /// (without PF, there is one family per argument)
  struct Family {
    /// the first argument of the family, used for the hills file and the grid
    Value* arg=nullptr;
    /// evaluation of the hills (periodicity of the argument, INTERVAL)
    MetaDKernel kernel;
    /// the hills, without grid
    std::vector<Gaussian> hills;
    /// the hills, without grid, in contiguous arrays (for speed): center, inverse width (0 for multivariate hills),
    /// height, and whether they can be evaluated with these arrays (not multivariate, no INTERVAL)
    std::vector<double> center;
    std::vector<double> invsigma;
    std::vector<double> height;
    std::vector<char> diagonal;
    /// the value of hill h at cv (one dimensional), and its derivative (added to der) if der is given
    double evaluate(const std::vector<double>& cv, unsigned h, double* der) const {
      double dp[1];
      if(diagonal[h]) {
        return der ? kernel.evaluateDiagonalWithDerivatives(cv,&center[h],&invsigma[h],height[h],der,dp) :
               kernel.evaluateDiagonal(cv,&center[h],&invsigma[h],height[h]);
      }
      return der ? kernel.evaluateWithDerivatives(cv,hills[h],der,dp) : kernel.evaluate(cv,hills[h]);
    }
    /// the bias, with grid
    std::unique_ptr<MetaDBiasGrid> grid;
    /// INTERVAL
    bool interval=false;
    double lower=0.0;
    double upper=0.0;
    /// hills files that are read (own file when restarting, files of the other walkers), and writing of the hills
    std::unique_ptr<MetaDHillsFiles> hillsfiles;
    std::unique_ptr<OFile> hillsOfile;
    std::unique_ptr<OFile> gridOfile;
  };
  // general setup
  double kbt_=0.0;
  int stride_=0;
  // well-tempered MetaD
  bool welltemp_=false;
  double biasf_=1.0;
  // output files format
  std::string fmt_;
  /// the steps when hills are added (the first step of the run, the last hill read when restarting)
  MetaDHillSchedule schedule_;
  // Gaussian starting parameters
  double height0_=std::numeric_limits<double>::max();
  std::vector<double> sigma0_;
  std::vector<double> sigma0min_;
  std::vector<double> sigma0max_;
  int adaptive_=FlexibleBin::none;
  /// with ADAPTIVE, the widths of the hills of each argument
  std::vector<FlexibleBin> flexbin_;
  // Grids
  bool grid_=false;
  int wgridstride_=0;
  // Partitioned Families
  std::vector<Family> families_;
  /// the family of each argument
  std::vector<unsigned> pfs_;
  /// if partitioned families are enabled
  bool do_pf_=false;
  // multiple walkers
  int mw_n_=1;
  std::string mw_dir_;
  int mw_id_=0;
  int mw_rstride_=1;
  bool walkers_mpi_=false;
  unsigned mpi_nw_=0;
  unsigned mpi_id_=0;
  /// ADAPTIVE=DIFF and replica exchange: an exchange has been attempted at the last step, with the CVs of
  /// the partner (partner_cv) and the CVs at the previous update (last_cv)
  struct {
    bool attempted=false;
    std::vector<double> partner_cv;
    std::vector<double> last_cv;
  } exchange_;
  // variable for selector
  std::string selector_;
  bool  do_select_=false;
  unsigned select_value_=0;
  unsigned current_value_=0;

  /// read the hills in input file i of family f and add them to the bias; returns their number and
  /// (if lasttime is given and at least one hill is read) sets lasttime to the time of the last hill
  unsigned readGaussians(unsigned f, unsigned i, double* lasttime=nullptr);
  /// update the adaptive hill widths; with ADAPTIVE=DIFF also handle replica exchanges
  void   updateFlexibleBins(bool nowAddAHill);
  /// true if a hill has to be added at the current step
  bool   isHillStep() const;
  void   writeGaussian(unsigned f, const Gaussian&);
  void   addGaussian(unsigned f, const Gaussian&);
  /// the bias of each argument, and its derivative if der is given
  void   getBiases(const std::vector<double>& cv, std::vector<double>& bias, std::vector<double>* der);

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
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition, one for all biases");
  keys.add("optional","FILE","files in which the lists of added hills are stored, default names are assigned using arguments if FILE is not found");
  keys.add("optional","HEIGHT","the height of the Gaussian hills, one for all biases. Compulsory unless TAU, TEMP and BIASFACTOR are given");
  keys.add("optional","FMT","specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics with this bias factor, one for all biases.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","TAU","in well tempered metadynamics, sets height to (k_B Delta T*pace*timestep)/tau");
  keys.add("optional","GRID_MIN","the lower bounds for the grid (initial range: the grid is extended when needed)");
  keys.add("optional","GRID_MAX","the upper bounds for the grid (initial range: the grid is extended when needed)");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN). "
           "Without GRID_MIN and GRID_MAX it gives the grid spacing, and the grid is placed around the first Gaussian");
  keys.add("compulsory","GRID_SIZE_LIMIT","100000000","the maximum number of points of each grid, which is extended when a Gaussian does not fit in it");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","GRID_WSTRIDE", "frequency for dumping the grid");
  keys.add("optional","GRID_WFILES", "dump grid for the bias, default names are used if GRID_WSTRIDE is used without GRID_WFILES.");
  keys.add("optional","GRID_RFILES", "read grid for the bias");
  keys.add("optional","ADAPTIVE","use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or timestep dimensions");
  keys.add("optional","SIGMA_MAX","the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","SIGMA_MIN","the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.addInputKeyword("numbered","PF", "scalar", "specify which CVs belong in a partitioned family. Once a PF is specified, all CVs in ARG must be placed in a PF even if there is one CV per PF");
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
  keys.addDOI("10.1021/acs.jctc.5b00846");
  keys.addDOI("10.1103/PhysRevLett.100.020603");
  keys.addDOI("10.1021/ct3002464");
  keys.addDOI("10.2174/1877946811202010079");
  keys.addDOI("10.1021/acs.jctc.8b00448");
  keys.addDOI("10.1021/jp054359r");
}

PBMetaD::PBMetaD(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao) {

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
  std::vector<int> pfs(getNumberOfArguments(), -1);
  // the first argument of each family (for HILLS & GRID files)
  std::vector<Value*> pfhold;
  std::vector<Value*> familyargs;
  for(int i = 0;; i++) {
    parseArgumentList("PF", i, familyargs);
    if (familyargs.empty()) {
      break;
    }

    do_pf_ = true;
    pfhold.push_back(nullptr);
    log << "  Identified Partitioned Family " << i << ":";
    for (unsigned j = 0; j < familyargs.size(); j++) {
      log << " " << familyargs[j]->getName();
      // loop through the argument list to make sure it exists and assign it
      bool foundArg = false;
      for (unsigned argnum = 0; argnum < getNumberOfArguments(); argnum++) {
        if (familyargs[j]->getName() == getPntrToArgument(argnum)->getName()) {
          foundArg = true;
          if (pfs[argnum] != -1) {
            error(familyargs[j]->getName() + " already present in PF" + std::to_string(pfs[argnum]));
          }
          pfs[argnum] = i;  // store the pf# for each cv
          if (pfhold[i] == nullptr) {
            // if this is the first argument in the family, store a pointer for it (this is for HILLS & GRID files)
            pfhold[i] = getPntrToArgument(argnum);
          }
        }
      }
      if (!foundArg) {
        error(familyargs[j]->getName() + " in PF" + std::to_string(i) + " not found in ARG");
      }
    }
    log << "\n";
  }

  // if PF were specified, every argument gets treated as its own PF
  if (!do_pf_) {
    for(unsigned i=0; i < getNumberOfArguments(); i++) {
      pfhold.push_back(getPntrToArgument(i));
      pfs[i] = i;
    }
  } else {
    // If we are doing PF, make sure each argument got assigned to a family.
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      if (pfs[i] == -1) {
        error(getPntrToArgument(i)->getName() + " was not assigned a PF");
      }
    }
  }
  const unsigned pf_n=pfhold.size();
  pfs_.assign(pfs.begin(),pfs.end());
  families_.resize(pf_n);
  for(unsigned f=0; f<pf_n; f++) {
    families_[f].arg=pfhold[f];
    families_[f].kernel.setArguments(std::vector<Value*>(1,pfhold[f]));
  }

  // parse the sigma
  parseVector("SIGMA",sigma0_);
  if(adaptive_==FlexibleBin::none) {
    // if you use normal sigma you need one sigma per argument
    if( sigma0_.size()!=pf_n ) {
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
    if(sigma0min_.size()>0 && sigma0min_.size()!=pf_n) {
      error("the number of SIGMA_MIN values be the same of the number of the arguments/PF");
    } else if(sigma0min_.size()==0) {
      sigma0min_.assign(pf_n,-1.);
    }

    parseVector("SIGMA_MAX",sigma0max_);
    if(sigma0max_.size()>0 && sigma0max_.size()!=pf_n) {
      error("the number of SIGMA_MAX values be the same of the number of the arguments/PF");
    } else if(sigma0max_.size()==0) {
      sigma0max_.assign(pf_n,-1.);
    }

    // one for each argument, with the limits of its family
    for(unsigned i=0; i<getNumberOfArguments(); i++) {
      std::vector<double> tmp_smin(1,sigma0min_[pfs_[i]]);
      std::vector<double> tmp_smax(1,sigma0max_[pfs_[i]]);
      flexbin_.push_back(FlexibleBin(adaptive_,this,i,sigma0_[0],tmp_smin,tmp_smax));
    }
  }

  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  parse("PACE",stride_);
  if(stride_<=0) {
    error("frequency for hill addition is nonsensical");
  }

  std::vector<std::string> hillsfname;
  parseVector("FILE",hillsfname);
  if(hillsfname.size()==0) {
    for(unsigned i=0; i< pf_n; i++) {
      std::string hillName = do_pf_ ? "HILLS.PF"+std::to_string(i) : "HILLS."+getPntrToArgument(i)->getName();
      hillsfname.push_back(hillName);
    }
  }
  if( hillsfname.size()!=pf_n ) {
    error("number of FILE arguments does not match number of HILLS files");
  }

  parse("BIASFACTOR",biasf_);
  if( biasf_<1.0 ) {
    error("well tempered bias factor is nonsensical");
  }
  kbt_=getkBT();
  if(biasf_>1.0) {
    if(kbt_==0.0) {
      error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    }
    welltemp_=true;
  }
  // the biases of the families are combined, and the heights of the hills are computed, with kBT
  if(kbt_==0.0) {
    error("Unless the MD engine passes the temperature to plumed, you must specify it using TEMP");
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
  std::vector<std::string> gridfilenames;
  parseVector("GRID_WFILES",gridfilenames);
  if (wgridstride_ == 0 && gridfilenames.size() > 0) {
    error("frequency with which to output grid not specified use GRID_WSTRIDE");
  }
  if(gridfilenames.size()==0 && wgridstride_ > 0) {
    for(unsigned i=0; i<pf_n; i++) {
      std::string gridfname = do_pf_ ? "GRID.PF"+std::to_string(i) : "GRID."+getPntrToArgument(i)->getName();
      gridfilenames.push_back(gridfname);
    }
  }
  if(gridfilenames.size() > 0 && hillsfname.size() > 0 && gridfilenames.size() != hillsfname.size()) {
    error("number of GRID_WFILES arguments does not match number of HILLS files");
  }

  // Read grid
  std::vector<std::string> gridreadfilenames;
  parseVector("GRID_RFILES",gridreadfilenames);

  // Grid Stuff
  std::vector<std::string> gmin(pf_n);
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=pf_n && gmin.size()!=0) {
    error("not enough values for GRID_MIN");
  }
  std::vector<std::string> gmax(pf_n);
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=pf_n && gmax.size()!=0) {
    error("not enough values for GRID_MAX");
  }
  std::vector<unsigned> gbin(pf_n);
  std::vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=pf_n && gbin.size()!=0) {
    error("not enough values for GRID_BIN");
  }
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=pf_n && gspacing.size()!=0) {
    error("not enough values for GRID_SPACING");
  }
  const bool gridnobounds=MetaDBiasGrid::setupBins(*this,pfhold,gmin,gmax,gbin,gspacing,adaptive_!=FlexibleBin::none,sigma0_,sigma0min_);
  if(gbin.size()>0) {
    grid_=true;
  }
  double gridsizelimit=1.0e8;
  parse("GRID_SIZE_LIMIT",gridsizelimit);
  MetaDBiasGrid::checkSizeLimit(*this,gridsizelimit);

  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  if(!grid_&&gridfilenames.size() > 0) {
    error("To write a grid you need first to define it!");
  }
  if(!grid_&&gridreadfilenames.size() > 0) {
    error("To read a grid you need first to define it!");
  }
  if(gridreadfilenames.size()>0 && gridreadfilenames.size()!=pf_n) {
    error("number of GRID_RFILES arguments does not match number of HILLS files");
  }

  // Interval keyword
  std::vector<double> lowI;
  std::vector<double> uppI;
  parseVector("INTERVAL_MIN",lowI);
  parseVector("INTERVAL_MAX",uppI);
  // various checks
  if(lowI.size()!=uppI.size()) {
    error("both a lower and an upper limits must be provided with INTERVAL");
  }
  if(lowI.size()!=0 && lowI.size()!=pf_n) {
    error("check number of argument of INTERVAL");
  }
  for(unsigned i=0; i<lowI.size(); ++i) {
    if(uppI[i]<lowI[i]) {
      error("The Upper limit must be greater than the Lower limit!");
    }
    if(pfhold[i]->isPeriodic()) {
      warning("INTERVAL is not used for periodic variables");
    } else {
      families_[i].interval=true;
      families_[i].lower=lowI[i];
      families_[i].upper=uppI[i];
      families_[i].kernel.setInterval(lowI[i],uppI[i]);
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
  for(unsigned i=0; i<hillsfname.size(); ++i) {
    log.printf("%s ",hillsfname[i].c_str());
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
      mpi_nw_=getNumberOfReplicas();
      mpi_id_=getReplicaIndex();
    }
  }

  for(unsigned f=0; f<pf_n; f++) {
    if(families_[f].interval) {
      log.printf("  Upper and Lower limits boundaries for the bias of CV %u are activated\n", f);
    }
  }
  if(grid_ && gridnobounds) {
    log.printf("  Grid spacing");
    for(unsigned i=0; i<gspacing.size(); ++i) {
      log.printf(" %s",exactString(gspacing[i]).c_str());
    }
    log.printf(" (periodic variables: grid bin");
    for(unsigned i=0; i<gbin.size(); ++i) {
      if(pfhold[i]->isPeriodic()) {
        log.printf(" %u",gbin[i]);
      } else {
        log.printf(" -");
      }
    }
    log.printf(")\n  The grids are placed around the first Gaussian\n");
  } else if(grid_) {
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
  }
  if(grid_) {
    if(spline) {
      log.printf("  Grid uses spline interpolation\n");
    }
    if(sparsegrid) {
      log.printf("  Grid uses sparse grid\n");
    }
    if(wgridstride_>0) {
      for(unsigned i=0; i<gridfilenames.size(); ++i) {
        log.printf("  Grid is written on file %s with stride %d\n",gridfilenames[i].c_str(),wgridstride_);
      }
    }
    if(gridreadfilenames.size()>0) {
      for(unsigned i=0; i<gridreadfilenames.size(); ++i) {
        log.printf("  Reading bias from grid in file %s \n",gridreadfilenames[i].c_str());
      }
    }
    log.printf("  The grids are extended when a Gaussian does not fit in them, up to %s points\n",exactString(gridsizelimit).c_str());
  }

  // restart from external grid
  bool restartedFromGrid=false;

  // initializing and checking grid
  if(grid_) {
    std::string funcl=getLabel() + ".bias";
    for(unsigned f=0; f<pf_n; ++f) {
      // check for mesh and sigma size
      double a,b;
      Tools::convert(gmin[f],a);
      Tools::convert(gmax[f],b);
      double mesh=(b-a)/((double)gbin[f]);
      if(adaptive_==FlexibleBin::none) {
        if(mesh>0.5*sigma0_[f]) {
          log<<"  WARNING: Using a PBMETAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
        }
      } else {
        if(mesh>0.5*sigma0min_[f]||sigma0min_[f]<0.) {
          log<<"  WARNING: to use a PBMETAD with a GRID and ADAPTIVE you need to set a Grid Spacing lower than half of the Gaussians (SIGMA_MIN) \n";
        }
      }
      std::vector<Value*> args(1,pfhold[f]);
      std::vector<std::string> gmin_t(1,gmin[f]);
      std::vector<std::string> gmax_t(1,gmax[f]);
      std::vector<unsigned>    gbin_t(1,gbin[f]);
      families_[f].grid=Tools::make_unique<MetaDBiasGrid>(*this);
      MetaDBiasGrid& grid=*families_[f].grid;
      // Read grid from file
      if(gridreadfilenames.size()>0) {
        IFile gridfile;
        gridfile.link(*this);
        if(gridfile.FileExist(gridreadfilenames[f])) {
          gridfile.open(gridreadfilenames[f]);
        } else {
          error("The GRID file you want to read: " + gridreadfilenames[f] + ", cannot be found!");
        }
        grid.read(gridfile,gridreadfilenames[f],funcl,args,gmin_t,gmax_t,gbin_t,sparsegrid,spline,gridnobounds);
        log.printf("  Restarting from %s:\n",gridreadfilenames[f].c_str());
        if(getRestart()) {
          restartedFromGrid=true;
        }
      } else {
        grid.create(funcl,args,gmin_t,gmax_t,gbin_t,sparsegrid,spline,gridnobounds);
      }
      grid.setSizeLimit(gridsizelimit);
      // with INTERVAL the bias beyond the limits is constant: the grid contains the interval and at least one point
      // beyond each limit, which holds this constant value
      if(families_[f].interval) {
        grid.setInterval(families_[f].lower,families_[f].upper,getTime());
      }
    }
  }

  // open all files at the beginning and read Gaussians if restarting
  std::vector<std::string> ownpath(pf_n);  // actual path of the hills files of this walker, if they were read
  bool lastHillTimeKnown=false;  // time of the last hill in the hills files of this walker, if they were read
  double lastHillTime=0.0;
  for(unsigned f=0; f<pf_n; ++f) {
    families_[f].hillsfiles=Tools::make_unique<MetaDHillsFiles>(*this,comm,std::vector<Value*>(1,pfhold[f]));
  }
  for(int j=0; j<mw_n_; ++j) {
    for(unsigned f=0; f<pf_n; ++f) {
      std::string fname;
      if(mw_dir_!="") {
        if(mw_n_>1) {
          std::stringstream out;
          out << j;
          fname = mw_dir_+"/"+hillsfname[f]+"."+out.str();
        } else if(walkers_mpi_) {
          fname = mw_dir_+"/"+hillsfname[f];
        } else {
          fname = hillsfname[f];
        }
      } else {
        if(mw_n_>1) {
          std::stringstream out;
          out << j;
          fname = hillsfname[f]+"."+out.str();
        } else {
          fname = hillsfname[f];
        }
      }
      MetaDHillsFiles& hillsfiles=*families_[f].hillsfiles;
      const unsigned k=hillsfiles.addInput(fname);
      plumed_assert(k==static_cast<unsigned>(j));
      if(hillsfiles.isOpen(k)) {
        if(getRestart()&&!restartedFromGrid) {
          log.printf("  Restarting from %s:",hillsfiles.getInputName(k).c_str());
          if(j==mw_id_) {
            // the hills of all the families are written at the same time: the last one decides whether
            // a hill is added at the first step
            double lasttime=0.0;
            if(readGaussians(f,k,&lasttime)>0 && !(lastHillTimeKnown && lasttime<=lastHillTime)) {
              lastHillTimeKnown=true;
              lastHillTime=lasttime;
            }
            ownpath[f]=hillsfiles.getPath(k);
          } else {
            readGaussians(f,k);
          }
        }
        if(getRestart() && restartedFromGrid && j==mw_id_) {
          // new hills will be appended to this file also when the bias is read from a grid
          ownpath[f]=hillsfiles.getPath(k);
        }
        hillsfiles.reset(k);
        // close only the walker own hills file for later writing
        if(j==mw_id_) {
          hillsfiles.close(k);
        }
      } else {
        // in case a file does not exist and we are restarting, complain that the file was not found
        if(getRestart()) {
          log<<"  WARNING: restart file "<<fname<<" not found\n";
        }
      }
    }
  }

  if(lastHillTimeKnown) {
    schedule_.setLastHillTime(lastHillTime);
  }

  comm.Barrier();
  if(comm.Get_rank()==0 && walkers_mpi_) {
    multi_sim_comm.Barrier();
  }

  // all walkers have read the hills: the one that writes the files can now remove an incomplete
  // last line (left by a run that was interrupted while writing), before new hills are appended to them
  for(unsigned f=0; f<pf_n; ++f) {
    std::string message;
    if(ownpath[f].length()>0 && comm.Get_rank()==0 && !(walkers_mpi_ && getReplicaIndex()>0)) {
      message=removeIncompleteLastLine(ownpath[f]);
    }
    // the log must be written by all the processes
    if(comm.Get_size()>1) {
      bcastString(comm,message);
    }
    if(message.length()>0) {
      log<<message;
    }
  }

  // open hills files for writing
  for(unsigned f=0; f<pf_n; ++f) {
    auto ofile=Tools::make_unique<OFile>();
    ofile->link(*this);
    // if MPI multiple walkers, only rank 0 will write to file
    if(walkers_mpi_) {
      int r=getReplicaIndex();
      if(r>0) {
        families_[f].hillsfiles->setInputName(mw_id_,"/dev/null");
      }
      ofile->enforceSuffix("");
    }
    if(mw_n_>1) {
      ofile->enforceSuffix("");
    }
    ofile->open(families_[f].hillsfiles->getInputName(mw_id_));
    if(fmt_.length()>0) {
      ofile->fmtField(fmt_);
    }
    ofile->addConstantField("multivariate");
    ofile->addConstantField("kerneltype");
    if(families_[f].interval) {
      ofile->addConstantField("lower_int").printField("lower_int",families_[f].lower);
      ofile->addConstantField("upper_int").printField("upper_int",families_[f].upper);
    }
    ofile->setHeavyFlush();
    // output periodicities of variables
    ofile->setupPrintValue( pfhold[f] );  //assuming cvs in the same family have the same periodicity and boundaries.
    families_[f].hillsOfile=std::move(ofile);
  }

  // Dump grid to files
  if(wgridstride_ > 0) {
    for(unsigned f = 0; f < gridfilenames.size(); ++f) {
      auto ofile=Tools::make_unique<OFile>();
      ofile->link(*this);
      std::string gridfname_tmp = gridfilenames[f];
      if(walkers_mpi_) {
        int r=getReplicaIndex();
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
      families_[f].gridOfile=std::move(ofile);
    }
  }

  log<<"  Bibliography "<<plumed.cite("Pfaendtner and Bonomi. J. Chem. Theory Comput. 11, 5062 (2015)");
  bool interval=false;
  for(const auto & family : families_) {
    interval=interval || family.interval;
  }
  if(interval)
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

unsigned PBMetaD::readGaussians(unsigned f, unsigned i, double* lasttime) {
  Family& family=families_[f];
  // the callers have started a line of the log
  if(family.grid) {
    family.grid->setMidline(true);
  }
  const unsigned nhills=family.hillsfiles->read(i,[&](double time,bool multivariate,double height,const std::vector<double>& center,const std::vector<double>& sigma) {
    if(lasttime) {
      *lasttime=time;
    }
    if(welltemp_) {
      height*=(biasf_-1.0)/biasf_;
    }
    addGaussian(f,Gaussian(multivariate,height,center,sigma));
  });
  if(family.grid) {
    family.grid->setMidline(false);
  }
  log.printf("      %u Gaussians read\n",nhills);
  return nhills;
}

void PBMetaD::writeGaussian(unsigned f, const Gaussian& hill) {
  double height=hill.height;
  if(welltemp_) {
    height *= biasf_/(biasf_-1.0);
  }
  families_[f].hillsfiles->write(*families_[f].hillsOfile,getTimeStep()*getStep(),hill,height,biasf_,mw_n_>1);
}

void PBMetaD::addGaussian(unsigned f, const Gaussian& hill) {
  Family& family=families_[f];
  if(family.grid) {
    family.grid->addHill(hill,family.kernel,getTime());
  } else {
    family.hills.push_back(hill);
    family.center.push_back(hill.center[0]);
    family.invsigma.push_back(hill.multivariate ? 0.0 : hill.invsigma[0]);
    family.height.push_back(hill.height);
    family.diagonal.push_back(!hill.multivariate && !family.interval);
  }
}

void PBMetaD::getBiases(const std::vector<double>& cv, std::vector<double>& bias, std::vector<double>* der) {
  const unsigned nargs=getNumberOfArguments();
  bias.assign(nargs,0.0);
  if(der) {
    der->assign(nargs,0.0);
  }
  std::vector<double> cvi(1);
  std::vector<double> deri(1);
  if(grid_) {
    for(unsigned i=0; i<nargs; ++i) {
      const MetaDBiasGrid& grid=*families_[pfs_[i]].grid;
      cvi[0]=cv[i];
      if(der) {
        bias[i]=grid.getBiasAndDerivatives(cvi,deri);
        (*der)[i]=deri[0];
      } else {
        bias[i]=grid.getBias(cvi);
      }
    }
    return;
  }
  // the hills are split among the MPI processes and the OpenMP threads
  const unsigned stride=comm.Get_size();
  const unsigned rank=comm.Get_rank();
  const unsigned nt=OpenMP::getNumThreads();
  for(unsigned i=0; i<nargs; ++i) {
    const Family& family=families_[pfs_[i]];
    const unsigned nhills=family.hills.size();
    cvi[0]=cv[i];
    double b=0.0;
    double d=0.0;
    if(nt==1 || nhills<2*nt*stride) {
      for(unsigned h=rank; h<nhills; h+=stride) {
        b+=family.evaluate(cvi,h,der ? &d : nullptr);
      }
    } else {
      #pragma omp parallel num_threads(nt)
      {
        double dt=0.0;
        #pragma omp for reduction(+:b) nowait
        for(unsigned h=rank; h<nhills; h+=stride) {
          b+=family.evaluate(cvi,h,der ? &dt : nullptr);
        }
        #pragma omp atomic
        d+=dt;
      }
    }
    bias[i]=b;
    if(der) {
      (*der)[i]=d;
    }
  }
  if(stride>1) {
    // a single reduction for the biases and the derivatives
    std::vector<double> buffer(bias);
    if(der) {
      buffer.insert(buffer.end(),der->begin(),der->end());
    }
    comm.Sum(buffer);
    std::copy(buffer.begin(),buffer.begin()+nargs,bias.begin());
    if(der) {
      std::copy(buffer.begin()+nargs,buffer.end(),der->begin());
    }
  }
}

void PBMetaD::calculate() {
  schedule_.start(getStep());

  const unsigned nargs=getNumberOfArguments();
  std::vector<double> cv(nargs);
  for(unsigned i=0; i<nargs; ++i) {
    cv[i]=getArgument(i);
  }

  // during a replica exchange attempt the bias is computed on the coordinates of the partner.
  // The MD code does not tell whether the exchange is accepted: with ADAPTIVE=DIFF this is
  // detected at the next update, comparing the CVs with those of the partner (see updateFlexibleBins)
  if(adaptive_==FlexibleBin::diffusion && getExchangeStep()) {
    exchange_.attempted=true;
    exchange_.partner_cv=cv;
  }
  std::vector<double> bias;
  std::vector<double> deriv;
  getBiases(cv,bias,&deriv);

  double ncv = (double) nargs;
  double bmin = 1.0e+19;
  for(unsigned i=0; i<nargs; ++i) {
    if(bias[i] < bmin) {
      bmin = bias[i];
    }
  }
  double ene = 0.;
  for(unsigned i=0; i<nargs; ++i) {
    ene += std::exp((-bias[i]+bmin)/kbt_);
  }

  // set Forces - set them to zero if SELECTOR is active
  if(do_select_) {
    current_value_ = static_cast<unsigned>(plumed.passMap[selector_]);
  }

  if(!do_select_ || select_value_==current_value_) {
    for(unsigned i=0; i<nargs; ++i) {
      const double force = - std::exp((-bias[i]+bmin)/kbt_) / (ene) * deriv[i];
      setOutputForce(i, force);
    }
  }

  if(do_select_ && select_value_!=current_value_) {
    for(unsigned i=0; i<nargs; ++i) {
      setOutputForce(i, 0.0);
    }
  }

  // set bias
  ene = -kbt_ * (std::log(ene) - std::log(ncv)) + bmin;
  setBias(ene);
}

bool PBMetaD::isHillStep() const {
  return schedule_.isHillStep(getStep(),stride_,getTime(),0.5*stride_*getTimeStep());
}

void PBMetaD::updateFlexibleBins(bool nowAddAHill) {
  const unsigned nargs=getNumberOfArguments();
  if(adaptive_==FlexibleBin::diffusion) {
    std::vector<double> cv(nargs);
    for(unsigned i=0; i<nargs; ++i) {
      cv[i]=getArgument(i);
    }
    if(exchange_.attempted && exchange_.last_cv.size()==nargs) {
      // the exchange has been accepted if the CVs are now closer to those of the partner than
      // to their own value at the previous update (distances in units of the current fluctuations)
      double dpartner=0.0;
      double down=0.0;
      for(unsigned i=0; i<nargs; i++) {
        const std::vector<double> var=flexbin_[i].getMatrix();
        const double s2=(var.size()>0 && var[0]>0.0 ? var[0] : 1.0);
        const double dp=difference(i,exchange_.partner_cv[i],cv[i]);
        const double dl=difference(i,exchange_.last_cv[i],cv[i]);
        dpartner+=dp*dp/s2;
        down+=dl*dl/s2;
      }
      if(dpartner<down) {
        // the CVs jumped: restart the running averages, but keep the variances
        for(auto & fb : flexbin_) {
          fb.restartAverage();
        }
      }
    }
    exchange_.attempted=false;
    exchange_.last_cv=cv;
  }
  for(unsigned i=0; i<nargs; i++) {
    flexbin_[i].update(nowAddAHill,i);
  }
}

void PBMetaD::update() {
  // adding hills criteria
  bool nowAddAHill=isHillStep();

  const unsigned nargs=getNumberOfArguments();
  // if you use adaptive, call the FlexibleBin
  bool multivariate=false;
  if(adaptive_!=FlexibleBin::none) {
    updateFlexibleBins(nowAddAHill);
    multivariate=true;
  }

  // with ADAPTIVE=DIFF the width of a hill requires the fluctuations of at least two CV values
  // (this happens at the first update, e.g. after UPDATE_FROM or when restarting)
  if(nowAddAHill && adaptive_==FlexibleBin::diffusion) {
    for(const auto & fb : flexbin_) {
      if(fb.getNumberOfSamples()<2) {
        nowAddAHill=false;
      }
    }
  }

  if(nowAddAHill && (!do_select_ || select_value_==current_value_)) {
    // get all biases and heights
    std::vector<double> cv(nargs);
    std::vector<double> thissigma(nargs);
    std::vector<double> height(nargs);
    for(unsigned i=0; i<nargs; ++i) {
      // the width of the hills of each argument
      if(adaptive_!=FlexibleBin::none) {
        thissigma[i]=flexbin_[i].getInverseMatrix(i)[0];
      } else {
        thissigma[i]=sigma0_[pfs_[i]];
      }
      cv[i] = getArgument(i);
    }
    std::vector<double> bias;
    getBiases(cv,bias,nullptr);
    double norm = 0.0;
    double bmin = 1.0e+19;
    for(unsigned i=0; i<nargs; ++i) {
      if(bias[i] < bmin) {
        bmin = bias[i];
      }
    }
    // calculate heights and norm
    for(unsigned i=0; i<nargs; ++i) {
      double h = std::exp((-bias[i]+bmin)/kbt_);
      norm += h;
      height[i] = h;
    }
    // normalize and apply welltemp correction
    for(unsigned i=0; i<nargs; ++i) {
      height[i] *=  height0_ / norm;
      if(welltemp_) {
        height[i] *= std::exp(-bias[i]/(kbt_*(biasf_-1.0)));
      }
    }

    std::vector<double> cv_tmp(1);
    std::vector<double> sigma_tmp(1);
    // MPI Multiple walkers: share hills and add them all
    if(walkers_mpi_) {
      // the hills of all walkers (packed in a single array: center, sigma, height of each argument,
      // so that one collective is enough)
      const unsigned recordsize=3*nargs;
      std::vector<double> record(recordsize);
      for(unsigned i=0; i<nargs; ++i) {
        record[i]=cv[i];
        record[nargs+i]=thissigma[i];
        record[2*nargs+i]=height[i];
      }
      std::vector<double> all_records(mpi_nw_*recordsize, 0.0);
      if(comm.Get_rank()==0) {
        // Communicate (only root)
        multi_sim_comm.Allgather(record,all_records);
      }
      // Share info with group members
      comm.Bcast(all_records,0);
      // now add hills one by one
      for(unsigned j=0; j<mpi_nw_; ++j) {
        const double* rec=&all_records[j*recordsize];
        for(unsigned i=0; i<nargs; ++i) {
          // Add CVs of same family together and write to same file
          cv_tmp[0]    = rec[i];
          sigma_tmp[0] = rec[nargs+i];
          Gaussian newhill(multivariate, rec[2*nargs+i], cv_tmp, sigma_tmp);
          addGaussian(pfs_[i], newhill);
          writeGaussian(pfs_[i], newhill);
        }
      }
      // just add your own hills
    } else {
      for(unsigned i=0; i<nargs; ++i) {
        // Add CVs of same family together and write to same file
        cv_tmp[0] = cv[i];
        sigma_tmp[0] = thissigma[i];
        Gaussian newhill(multivariate, height[i], cv_tmp, sigma_tmp);
        addGaussian(pfs_[i], newhill);
        writeGaussian(pfs_[i], newhill);
      }
    }
  }

  // with multiple walkers the hills must reach the files before the other walkers read them
  // (this should be outside of the if block in case WALKERS_RSTRIDE is not a multiple of PACE)
  if(mw_n_>1 && getStep()%mw_rstride_==0) {
    for(auto & family : families_) {
      family.hillsOfile->flush();
    }
  }

  // write grid files
  if(wgridstride_>0 && (getStep()%wgridstride_==0 || getCPT())) {
    int r = 0;
    if(walkers_mpi_) {
      r=getReplicaIndex();
    }
    if(r==0) {
      for(auto & family : families_) {
        if(family.gridOfile) {
          family.gridOfile->rewind();
          family.grid->getGrid().writeToFile(*family.gridOfile);
          family.gridOfile->flush();
        }
      }
    }
  }

  // if multiple walkers and time to read Gaussians
  if(mw_n_>1 && getStep()%mw_rstride_==0) {
    for(int j=0; j<mw_n_; ++j) {
      for(unsigned f=0; f<families_.size(); ++f) {
        // don't read your own Gaussians
        if(j==mw_id_) {
          continue;
        }
        MetaDHillsFiles& hillsfiles=*families_[f].hillsfiles;
        // if the file is not open yet
        if(!hillsfiles.isOpen(j)) {
          // check if it exists now and open it!
          if(hillsfiles.openIfExists(j)) {
            hillsfiles.reset(j);
          }
          // otherwise read the new Gaussians
        } else {
          log.printf("  Reading hills from %s:",hillsfiles.getInputName(j).c_str());
          readGaussians(f,j);
          hillsfiles.reset(j);
        }
      }
    }
  }

}

bool PBMetaD::checkNeedsGradients()const {
  return adaptive_==FlexibleBin::geometry && isHillStep();
}

}
}
