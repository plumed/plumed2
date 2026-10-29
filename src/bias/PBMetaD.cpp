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
#include "core/PlumedMain.h"
#include "core/FlexibleBin.h"
#include "tools/Exception.h"
#include "tools/OpenMP.h"
#include "tools/File.h"
#include "tools/Communicator.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>
#include <limits>

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
V_{PB}(\vec{s},t) = -k_B T \log{\frac{1}{N}\sum_{i=1}^N
\exp\left(
- \frac{V(s_i,t)}{k_B T}
\right)}.
$$

PBMETAD needs the temperature: give it with TEMP, unless the MD engine passes it to PLUMED.
As in [METAD](METAD.md), the Gaussians are truncated and stretched so that they go to zero at the cutoff
(the hills files contain stretched Gaussians, which is what [sum_hills](sum_hills.md) assumes).


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
exactly zero outside it (with INTERVAL_MIN and INTERVAL_MAX, the grid contains the interval, and beyond the limits the bias
is the constant value at the limits).  The grid spacing does not change and the new points are on the same lattice as the initial ones.
Every extension is written on the log.  If you give GRID_SPACING without GRID_MIN and GRID_MAX each grid is placed
around its first Gaussian and then extended as needed.  GRID_SIZE_LIMIT sets the maximum number of points of each grid.
A grid read with GRID_RFILES can have any range, but it must be on the lattice of the grid given in input.

As in [METAD](METAD.md), no Gaussian is added at the first step of a simulation.  When a simulation is continued with
RESTART and the bias is read from the hills files, a Gaussian is added at the first step unless they already contain
it (i.e. the previous run had already performed this step).  When the bias is read from the grids (GRID_RFILES), the hills files are still read,
without adding their hills to the bias, to know the time of their last hill.  With WALKERS_N, the grids written by a walker
contain also the hills it has read from the files of the other walkers (their number is written in a comment line after
the grid), and only the hills of the other walkers that are not in the grids are added.

With ADAPTIVE=DIFF, as in [METAD](METAD.md), when restarting the statistics of the fluctuations of each variable
(running average and variance) are initialized from its last hills in the hills file of its family, so that the hills
deposited after the restart have the same width they would have had without the restart (use DIFF_FROM_SCRATCH to
accumulate the statistics from scratch).  With replica exchange, the running averages of the variables are
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
This method allows you to adapt the Gaussian to the extent of Cartesian space covered by a variable or
to the space in collective variable covered in a given time. In this case the width of the deposited
Gaussian potential is denoted by one value only that is a Cartesian space (ADAPTIVE=GEOM) or a time
(ADAPTIVE=DIFF), and SIGMA_MIN and SIGMA_MAX (one value for each family) limit the widths. Note that a specific
integration technique for the deposited Gaussian kernels should be used in this case. Check the documentation for utility
sum_hills.

With the keywords INTERVAL_MIN and INTERVAL_MAX (one value for each family; the values for periodic variables are
ignored) one changes the metadynamics algorithm and sets the bias force equal to zero
outside the boundaries as discussed in the fourth paper cited below. If, for example, metadynamics is performed on a CV s and one is interested only
to the free energy for s > boundary, the history dependent potential is still updated according to the above
equations but the metadynamics force is set to zero for s < boundary. Notice that Gaussian kernels are added also
if s < boundary, as the tails of these Gaussian kernels influence VG in the relevant region s > boundary. In this way, the
force on the system in the region s > boundary comes from both metadynamics and the force field, in the region
s < boundary only from the latter. This approach allows obtaining a history-dependent bias potential VG that
fluctuates around a stable estimator, equal to the negative of the free energy far enough from the
boundaries. Note that:

- It works only for one-dimensional biases;
- It works both with and without GRID;
- The interval limit boundary should be set in a region where the free energy derivative is not large;
- If in the region outside the limit boundary the system has a free energy minimum, INTERVAL_MIN and INTERVAL_MAX should
  be used together with a [UPPER_WALLS](UPPER_WALLS.md) or [LOWER_WALLS](LOWER_WALLS.md) at boundary.

For systems with multiple CVs that share identical properties, PBMetaD with partitioned families can be used
to group them under one bias potential, to which each of them contributes, as discussed in the penultimate paper in the
references below. This is done with a list of numbered PF keywords (PF0, PF1, ...), each of which contains the list of
CVs from ARG to be placed in that family. Once invoked, each CV in ARG must be placed in exactly one PF, even if it
results in families containing only one CV, and the CVs in a family must have the same periodicity.
Keywords that take one value for each argument without families (SIGMA, SIGMA_MIN, SIGMA_MAX, FILE, the GRID keywords,
INTERVAL_MIN, INTERVAL_MAX) take one value for each family instead.  The hills of all the CVs of a family are written in
the hills file of the family, as hills of its first CV.

With SELECTOR, the bias is applied and the hills are added only when the value of the selector (see
[SELECTOR](SELECTOR.md)) is equal to SELECTOR_ID.

The multiple walkers that is discussed in the final paper cited below can also be used. See below the examples.

## Examples

The following input is for PBMetaD calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the PBMetaD bias potential are written to the COLVAR file every 100 steps.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 TEMP=300 LABEL=pb FILE=HILLS_d1,HILLS_d2
PRINT ARG=d1,d2,pb.bias STRIDE=100 FILE=COLVAR
```

(See also [DISTANCE](DISTANCE.md) and [PRINT](PRINT.md)).

If you use well-tempered metadynamics, you should specify a single bias factor and initial Gaussian height.

```plumed
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
PBMETAD ...
ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3
PACE=500 BIASFACTOR=8 TEMP=300 LABEL=pb
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
PACE=500 BIASFACTOR=8 TEMP=300 LABEL=pb
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
PACE=500 BIASFACTOR=8 TEMP=300 LABEL=pb
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
PACE=500 BIASFACTOR=8 TEMP=300 LABEL=pb
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
    /// ADAPTIVE=DIFF: the last hills of the own hills file, for the statistics of the fluctuations at restart
    MetaDDiffusionRestart diffrestart;
    /// the arguments in the family
    std::vector<unsigned> args;
    /// WALKERS_N, restart from a grid: the number of hills of each of the other walkers in the grid (-1 for this walker);
    /// empty if not known (grid files written by older versions)
    std::vector<long long> gridWalkerHills;
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

  /// read the hills in input file i of family f and add them to the bias, except the first skipcount ones and those with
  /// time up to skipupto (they are already in the bias), or all of them if addtobias is false; returns their number and
  /// (if lasttime is given and at least one hill is read) sets lasttime to the time of the last hill
  unsigned readGaussians(unsigned f, unsigned i, double* lasttime=nullptr, bool addtobias=true,
                         double skipupto=-std::numeric_limits<double>::infinity(), unsigned long long skipcount=0);
  /// ADAPTIVE=DIFF: initialize the statistics of the fluctuations from the hills read when restarting
  void   initializeDiffusionStatistics();
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
  keys.add("optional","TEMP","the system temperature - it is needed unless the MD engine passes it to PLUMED");
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
  keys.add("optional","GRID_WFILES", "dump grid for the bias (one file for each family; METAD: GRID_WFILE), default names are used if GRID_WSTRIDE is used without GRID_WFILES.");
  keys.add("optional","GRID_RFILES", "read grid for the bias (one file for each family; METAD: GRID_RFILE)");
  keys.add("optional","ADAPTIVE","use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or timestep dimensions");
  keys.add("optional","SIGMA_MAX","the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","SIGMA_MIN","the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.addFlag("DIFF_FROM_SCRATCH",false,"with ADAPTIVE=DIFF, when restarting accumulate the statistics of the CV fluctuations from scratch, instead of initializing them from the last hills of the hills files");
  keys.addInputKeyword("numbered","PF", "scalar", "specify which CVs belong in a partitioned family. Once a PF is specified, all CVs in ARG must be placed in a PF even if there is one CV per PF");
  keys.add("optional","SELECTOR", "add forces and do update based on the value of SELECTOR");
  keys.add("optional","SELECTOR_ID", "value of SELECTOR");
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers - not compatible with WALKERS_* options other than WALKERS_DIR");
  keys.add("optional","INTERVAL_MIN","one dimensional lower limits (one for each family), outside the limits the system will not feel the biasing force (METAD: INTERVAL).");
  keys.add("optional","INTERVAL_MAX","one dimensional upper limits (one for each family), outside the limits the system will not feel the biasing force (METAD: INTERVAL).");
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
  // the arguments of a family share a one dimensional bias (and its hills file and grid)
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    const Value* arg=getPntrToArgument(i);
    const Value* first=pfhold[pfs_[i]];
    bool same=(arg->isPeriodic()==first->isPeriodic());
    if(same && arg->isPeriodic()) {
      std::string amin,amax,fmin,fmax;
      arg->getDomain(amin,amax);
      first->getDomain(fmin,fmax);
      same=(amin==fmin && amax==fmax);
    }
    if(!same) {
      error("the arguments in a partitioned family must have the same periodicity: "+arg->getName()+" and "+first->getName()+" do not");
    }
  }
  families_.resize(pf_n);
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    families_[pfs_[i]].args.push_back(i);
  }
  for(unsigned f=0; f<pf_n; f++) {
    families_[f].arg=pfhold[f];
    families_[f].kernel.setArguments(std::vector<Value*>(1,pfhold[f]));
    // stretched Gaussians (they go to zero at the cutoff), as in METAD and as written in the hills files
    if(!dp2cutoffNoStretch()) {
      families_[f].kernel.setStretch(dp2cutoffA,dp2cutoffB);
    }
  }

  // parse the sigma
  parseVector("SIGMA",sigma0_);
  for(const auto s : sigma0_) {
    if(!(s>0.0)) {
      error("SIGMA should be positive");
    }
  }
  if(adaptive_==FlexibleBin::none) {
    // if you use normal sigma you need one sigma per family
    if( sigma0_.size()!=pf_n ) {
      std::string fields = do_pf_ ? "PFs" : "arguments";
      error("number of " + fields + " does not match number of SIGMA parameters");
    }
    std::vector<double> smin,smax;
    parseVector("SIGMA_MIN",smin);
    parseVector("SIGMA_MAX",smax);
    if(smin.size()>0 || smax.size()>0) {
      error("SIGMA_MIN and SIGMA_MAX can only be used with ADAPTIVE");
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
    // the limits of the widths (one value for each family)
    parseVector("SIGMA_MIN",sigma0min_);
    if(sigma0min_.size()>0 && sigma0min_.size()!=pf_n) {
      error("the number of SIGMA_MIN values should be the same as the number of arguments (PFs)");
    } else if(sigma0min_.size()==0) {
      sigma0min_.assign(pf_n,-1.);
    }

    parseVector("SIGMA_MAX",sigma0max_);
    if(sigma0max_.size()>0 && sigma0max_.size()!=pf_n) {
      error("the number of SIGMA_MAX values should be the same as the number of arguments (PFs)");
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
  bool diff_from_scratch=false;
  parseFlag("DIFF_FROM_SCRATCH",diff_from_scratch);
  if(diff_from_scratch && adaptive_!=FlexibleBin::diffusion) {
    error("DIFF_FROM_SCRATCH can only be used with ADAPTIVE=DIFF");
  }
  // the hills are kept only if needed (restart with ADAPTIVE=DIFF), see readGaussians
  for(auto & family : families_) {
    family.diffrestart.pending=(adaptive_==FlexibleBin::diffusion && getRestart() && !diff_from_scratch);
  }
  if(adaptive_==FlexibleBin::diffusion && getRestart()) {
    if(diff_from_scratch) {
      log.printf("  The statistics for the hills width are accumulated from scratch after the restart (DIFF_FROM_SCRATCH)\n");
    } else {
      log.printf("  The statistics for the hills width are initialized from the last hills of the hills files, if available\n");
    }
  }

  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  if(height0_<0.0) {
    error("HEIGHT should not be negative");
  }
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
    error("the number of FILE names does not match the number of arguments (PFs)");
  }

  parse("BIASFACTOR",biasf_);
  if( biasf_<1.0 ) {
    error("well tempered bias factor is nonsensical");
  }
  kbt_=getkBT();
  if(biasf_>1.0) {
    welltemp_=true;
  }
  // the biases of the families are combined, and the heights of the hills are computed, with kBT
  if(kbt_==0.0) {
    error("Unless the MD engine passes the temperature to plumed, you must specify it using TEMP");
  }
  double tau=0.0;
  parse("TAU",tau);
  if(tau<0.0) {
    error("TAU should be positive");
  }
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
  if(mw_n_<1) {
    error("the number of walkers should be positive");
  }
  if(mw_id_<0 || mw_n_<=mw_id_) {
    error("walker ID should be a numerical value between 0 and the total number of walkers minus one");
  }
  parse("WALKERS_DIR",mw_dir_);
  parse("WALKERS_RSTRIDE",mw_rstride_);
  if(mw_rstride_<=0) {
    error("WALKERS_RSTRIDE should be positive");
  }

  // MPI version
  parseFlag("WALKERS_MPI",walkers_mpi_);

  // Grid file
  parse("GRID_WSTRIDE",wgridstride_);
  if(wgridstride_<0) {
    error("GRID_WSTRIDE should not be negative");
  }
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
  // (parseVector checks the number of values: the vectors are empty if the keywords are not given)
  std::vector<std::string> gmin(pf_n);
  parseVector("GRID_MIN",gmin);
  std::vector<std::string> gmax(pf_n);
  parseVector("GRID_MAX",gmax);
  std::vector<unsigned> gbin(pf_n);
  parseVector("GRID_BIN",gbin);
  std::vector<double> gspacing(pf_n);
  parseVector("GRID_SPACING",gspacing);
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
    error("both the lower and the upper limits must be provided with INTERVAL_MIN and INTERVAL_MAX");
  }
  if(lowI.size()!=0 && lowI.size()!=pf_n) {
    error("the number of INTERVAL_MIN and INTERVAL_MAX values should be the same as the number of arguments (PFs)");
  }
  for(unsigned i=0; i<lowI.size(); ++i) {
    if(uppI[i]<lowI[i]) {
      error("The Upper limit must be greater than the Lower limit!");
    }
    // (a value must be given for every family: the values for the periodic ones are not used)
    if(pfhold[i]->isPeriodic()) {
      warning("INTERVAL is not used for the periodic variable "+pfhold[i]->getName());
      continue;
    }
    families_[i].interval=true;
    families_[i].lower=lowI[i];
    families_[i].upper=uppI[i];
    families_[i].kernel.setInterval(lowI[i],uppI[i]);
  }

  // parse selector stuff
  parse("SELECTOR", selector_);
  if(selector_.length()>0) {
    do_select_ = true;
    select_value_ = 0; // default value, used if SELECTOR_ID is not given
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
  log.printf("\n  Gaussian height %f\n",height0_);
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
      log.printf("  Multiple walkers active using MPI communication\n");
      if(mw_dir_!="") {
        log.printf("  directory with hills files %s\n",mw_dir_.c_str());
      }
      mpi_nw_=getNumberOfReplicas();
      mpi_id_=getReplicaIndex();
    }
  }

  for(unsigned f=0; f<pf_n; f++) {
    if(families_[f].interval) {
      log.printf("  Upper and Lower limits boundaries for the bias of %s%s are activated at %f - %f\n",
                 do_pf_ ? ("PF"+std::to_string(f)+": ").c_str() : "",pfhold[f]->getName().c_str(),families_[f].lower,families_[f].upper);
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

  // with WALKERS_MPI, the replicas must take the same decision; the replicas can have different directories
  auto checkReplicas=[&](bool restarted,const std::string& what) {
    if(getRestart() && walkers_mpi_) {
      std::vector<int> all(mpi_nw_,0);
      if(comm.Get_rank()==0) {
        multi_sim_comm.Allgather(int(restarted),all);
      }
      comm.Bcast(all,0);
      const unsigned n=std::accumulate(all.begin(),all.end(),0);
      if(n!=0 && n!=mpi_nw_) {
        error("in this WALKERS_MPI run some replicas have restarted from "+what+" while others have not");
      }
    }
  };
  // the directory of the first replica, so that all the walkers read the same files
  auto firstReplicaDirectory=[&]() {
    std::string dir=std::filesystem::current_path().string();
    // the paths of the replicas can have different lengths
    if(comm.Get_rank()==0) {
      bcastString(multi_sim_comm,dir);
    }
    bcastString(comm,dir);
    return dir;
  };
  if(walkers_mpi_ && gridreadfilenames.size()>0) {
    const std::string dir=firstReplicaDirectory();
    for(auto & name : gridreadfilenames) {
      if(name.at(0)!='/') {
        name=dir+"/"+name;
      }
    }
  }

  // initializing and checking grid
  if(grid_) {
    std::string funcl=getLabel() + ".bias";
    for(unsigned f=0; f<pf_n; ++f) {
      // check for mesh and sigma size
      MetaDBiasGrid::checkSpacing(*this,std::vector<std::string>(1,gmin[f]),std::vector<std::string>(1,gmax[f]),std::vector<unsigned>(1,gbin[f]),
                                  adaptive_!=FlexibleBin::none,std::vector<double>(1,sigma0_[adaptive_!=FlexibleBin::none ? 0 : f]),
                                  std::vector<double>(1,adaptive_!=FlexibleBin::none ? sigma0min_[f] : -1.0));
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
        if(walkers_mpi_ || mw_n_>1) {
          // multiple walkers write their files without the suffix of the replica: with WALKERS_MPI all the walkers
          // read the same file
          gridfile.enforceSuffix("");
        }
        if(gridfile.FileExist(gridreadfilenames[f])) {
          gridfile.open(gridreadfilenames[f]);
        } else {
          error("The GRID file you want to read: " + gridreadfilenames[f] + ", cannot be found!");
        }
        grid.read(gridfile,gridreadfilenames[f],funcl,args,gmin_t,gmax_t,gbin_t,sparsegrid,spline,gridnobounds);
        if(mw_n_>1) {
          // the hills of the other walkers in the grid
          families_[f].gridWalkerHills=MetaDHillsFiles::readGridComment(comm,gridfile.getPath(),mw_n_);
        }
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
  checkReplicas(restartedFromGrid,"GRID");

  if(walkers_mpi_ && mw_dir_=="") {
    bool relative=true;
    for(const auto & name : hillsfname) {
      relative=relative && name.at(0)!='/';
    }
    if(relative) {
      mw_dir_=firstReplicaDirectory();
    }
  }

  // open all files at the beginning and read Gaussians if restarting
  bool restartedFromHills=false;
  std::vector<std::string> ownpath(pf_n);  // actual path of the hills files of this walker, if they were read
  bool lastHillTimeKnown=false;  // time of the last hill in the hills files of this walker, if they were read
  double lastHillTime=0.0;
  for(unsigned f=0; f<pf_n; ++f) {
    families_[f].hillsfiles=Tools::make_unique<MetaDHillsFiles>(*this,comm,std::vector<Value*>(1,pfhold[f]));
    // ADAPTIVE=DIFF: hills of the last ~6 averaging times (older ones have weight < exp(-6)), for each argument
    // of the family; with MPI walkers the file contains the hills of all walkers
    families_[f].diffrestart.maxhills=(static_cast<unsigned>(std::ceil(6.0*sigma0_[0]/stride_))+2)*
                                      std::max(1u,mpi_nw_)*families_[f].args.size();
  }
  for(int j=0; j<mw_n_; ++j) {
    for(unsigned f=0; f<pf_n; ++f) {
      const std::string fname=MetaDHillsFiles::getWalkerFileName(hillsfname[f],mw_dir_,mw_n_,j,walkers_mpi_);
      // with multiple walkers the files are written without the suffix of the replica
      const unsigned k=families_[f].hillsfiles->addInput(fname,mw_n_>1 || walkers_mpi_);
      plumed_assert(k==static_cast<unsigned>(j));
    }
  }
  // When the bias is read from the grids, the own hills files are read first: the time of their last hill tells
  // which hills of the other walkers (WALKERS_N) are already in the grids
  std::vector<int> order;
  if(getRestart() && restartedFromGrid) {
    order.push_back(mw_id_);
  }
  for(int j=0; j<mw_n_; ++j) {
    if(order.empty() || order[0]!=j) {
      order.push_back(j);
    }
  }
  for(const int j : order) {
    for(unsigned f=0; f<pf_n; ++f) {
      MetaDHillsFiles& hillsfiles=*families_[f].hillsfiles;
      const unsigned k=j;
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
          restartedFromHills=true;
        }
        if(getRestart() && restartedFromGrid && j==mw_id_) {
          // new hills will be appended to this file also when the bias is read from a grid
          ownpath[f]=hillsfiles.getPath(k);
          // the bias is read from the grid, but the hills are still needed for the time of the last hill (to decide
          // if a hill is added at the first step, and which hills of the other walkers are in the grid) and for the
          // statistics of ADAPTIVE=DIFF
          log.printf("  Reading hills from %s:",hillsfiles.getInputName(k).c_str());
          double lasttime=0.0;
          try {
            if(readGaussians(f,k,&lasttime,false)>0 && !(lastHillTimeKnown && lasttime<=lastHillTime)) {
              lastHillTimeKnown=true;
              lastHillTime=lasttime;
            }
          } catch(const std::exception& e) {
            // the bias comes from the grid: if the hills file cannot be used, its information is not used
            log.printf("\n  WARNING: the hills file cannot be used for the time of the last hill and the statistics of ADAPTIVE=DIFF"
                       " (new hills are appended to it):\n%s\n",e.what());
            families_[f].diffrestart.clear();
          }
        }
        if(getRestart() && restartedFromGrid && j!=mw_id_) {
          // the hills of the other walkers that are not in the grid are added
          const std::vector<long long>& ingridcount=families_[f].gridWalkerHills;
          const double ingrid=(lastHillTimeKnown ? lastHillTime : std::numeric_limits<double>::infinity());
          if(ingridcount.size()==static_cast<unsigned>(mw_n_) && ingridcount[j]>=0) {
            // the grid file contains the number of hills read from each walker
            log.printf("  Reading hills from %s (the first %lld are in the grid):",hillsfiles.getInputName(k).c_str(),ingridcount[j]);
            readGaussians(f,k,nullptr,true,-std::numeric_limits<double>::infinity(),ingridcount[j]);
          } else if(lastHillTimeKnown) {
            // grid written by an older version: those up to the last hill of this walker are assumed to be in the grid
            log.printf("  Reading hills from %s (those up to time %f are assumed to be in the grid):",hillsfiles.getInputName(k).c_str(),ingrid);
            readGaussians(f,k,nullptr,true,ingrid);
          } else {
            log.printf("  Reading hills from %s (they are assumed to be in the grid):",hillsfiles.getInputName(k).c_str());
            readGaussians(f,k,nullptr,true,ingrid);
          }
        }
        hillsfiles.reset(k);
        // close only the walker own hills file for later writing
        if(j==mw_id_) {
          hillsfiles.close(k);
        }
      } else {
        // in case a file does not exist and we are restarting, complain that the file was not found
        if(getRestart()&&!restartedFromGrid) {
          error("restart file "+hillsfiles.getInputName(k)+" not found");
        }
      }
    }
  }

  checkReplicas(restartedFromHills,"FILE");

  if(walkers_mpi_) {
    // the decision to add a hill at the first step must be the same for all the walkers (they share the hills)
    std::vector<double> last {lastHillTimeKnown ? 1.0 : 0.0, lastHillTime};
    if(comm.Get_rank()==0) {
      multi_sim_comm.Bcast(last,0);
    }
    comm.Bcast(last,0);
    lastHillTimeKnown=(last[0]!=0.0);
    lastHillTime=last[1];
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

unsigned PBMetaD::readGaussians(unsigned f, unsigned i, double* lasttime, bool addtobias, double skipupto, unsigned long long skipcount) {
  Family& family=families_[f];
  // the callers have started a line of the log
  if(family.grid) {
    family.grid->setMidline(true);
  }
  const std::string path=family.hillsfiles->getPath(i);
  unsigned long long nskipped=0;
  unsigned nhills=0;
  try {
    nhills=family.hillsfiles->read(i,[&](double time,bool multivariate,double height,const std::vector<double>& center,const std::vector<double>& sigma) {
      if(lasttime) {
        *lasttime=time;
      }
      if(lasttime && family.diffrestart.pending) {
        // the own hills file is being read: keep the last hills for the statistics of ADAPTIVE=DIFF
        family.diffrestart.add(path,time,multivariate,center,sigma);
      }
      if(!addtobias) {
        return;
      }
      // times are written with a finite precision (FMT): times closer than half a PACE are the same time
      if(nskipped<skipcount || time<=skipupto+0.5*stride_*getTimeStep()) {
        nskipped++;
        return;
      }
      if(welltemp_) {
        height*=(biasf_-1.0)/biasf_;
      }
      addGaussian(f,Gaussian(multivariate,height,center,sigma));
    });
  } catch(...) {
    if(family.grid) {
      family.grid->setMidline(false);
    }
    throw;
  }
  if(family.grid) {
    family.grid->setMidline(false);
  }
  if(!addtobias) {
    log.printf("      %u Gaussians scanned (the bias is read from the grid)\n",nhills);
  } else if(nskipped>0) {
    log.printf("      %u Gaussians read, %llu of them already in the grid\n",nhills,nskipped);
  } else {
    log.printf("      %u Gaussians read\n",nhills);
  }
  return nhills;
}

void PBMetaD::initializeDiffusionStatistics() {
  const double dt=getTimeStep();
  for(auto & family : families_) {
    MetaDDiffusionRestart& diffrestart=family.diffrestart;
    if(!diffrestart.pending) {
      continue;
    }
    diffrestart.pending=false;
    // In the hills file of the family, each deposition is a block with a hill for each argument of the family, in
    // their order (and with MPI walkers for each walker, in the order of the replicas)
    const unsigned nf=family.args.size();
    for(unsigned p=0; p<nf; p++) {
      const unsigned i=family.args[p];
      std::vector<double> average;
      std::vector<double> variance;
      double lasttime=0.0;
      // times are written with a finite precision (FMT): times closer than half a PACE are the same time
      const unsigned n=diffrestart.compute(0.5*stride_*dt,getTime(),(walkers_mpi_ ? mpi_nw_ : 1)*nf,(walkers_mpi_ ? mpi_id_ : 0)*nf+p,
                                           1.0/sigma0_[0],dt*getStride(),
      [&](unsigned,double a,double b) {
        return difference(i,a,b);
      },
      [&](unsigned,double x) {
        return bringBackInPbc(i,x);
      },
      average,variance,lasttime);
      if(n==0) {
        log.printf("  ADAPTIVE=DIFF: no hills available for %s, the statistics for the hills width start from scratch\n",
                   getPntrToArgument(i)->getName().c_str());
        continue;
      }
      flexbin_[i].setStatistics(average,variance,2);
      log.printf("  ADAPTIVE=DIFF: the statistics for the hills width of %s are initialized from %u hills of %s (last one at time %f)\n",
                 getPntrToArgument(i)->getName().c_str(),n,diffrestart.getFile().c_str(),lasttime);
    }
    diffrestart.clear();
  }
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

  // the biases are combined as -kT log(1/N sum exp(-V_i/kT)), shifting the exponents by the smallest bias
  const double ncv = static_cast<double>(nargs);
  const double bmin = *std::min_element(bias.begin(),bias.end());
  double ene = 0.;
  for(unsigned i=0; i<nargs; ++i) {
    ene += std::exp((-bias[i]+bmin)/kbt_);
  }

  // set Forces - set them to zero if SELECTOR is active
  if(do_select_) {
    const auto selector=plumed.passMap.find(selector_);
    if(selector==plumed.passMap.end()) {
      error("SELECTOR "+selector_+" has not been defined");
    }
    current_value_ = static_cast<unsigned>(selector->second);
  }

  const bool selected=(!do_select_ || select_value_==current_value_);
  for(unsigned i=0; i<nargs; ++i) {
    const double force = (selected ? - std::exp((-bias[i]+bmin)/kbt_) / (ene) * deriv[i] : 0.0);
    setOutputForce(i, force);
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
    initializeDiffusionStatistics();
  }
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

  // with SELECTOR, hills are added only when the selector has the value SELECTOR_ID; with WALKERS_MPI all the walkers
  // take part in the exchange of the hills (also those that do not add hills)
  const bool selected=(!do_select_ || select_value_==current_value_);
  if(nowAddAHill && (selected || walkers_mpi_)) {
    // get all biases and heights
    std::vector<double> cv(nargs);
    std::vector<double> thissigma(nargs);
    std::vector<double> height(nargs);
    for(unsigned i=0; i<nargs && selected; ++i) {
      // the width of the hills of each argument
      if(adaptive_!=FlexibleBin::none) {
        thissigma[i]=flexbin_[i].getInverseMatrix(i)[0];
      } else {
        thissigma[i]=sigma0_[pfs_[i]];
      }
      cv[i] = getArgument(i);
    }
    if(selected) {
      std::vector<double> bias;
      getBiases(cv,bias,nullptr);
      const double bmin=*std::min_element(bias.begin(),bias.end());
      // calculate heights and norm
      double norm = 0.0;
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
    }

    std::vector<double> cv_tmp(1);
    std::vector<double> sigma_tmp(1);
    // MPI Multiple walkers: share hills and add them all
    if(walkers_mpi_) {
      // the hills of all walkers (packed in a single array: center, sigma, height of each argument, and whether the
      // walker adds hills, so that one collective is enough)
      const unsigned recordsize=3*nargs+1;
      std::vector<double> record(recordsize,0.0);
      for(unsigned i=0; i<nargs; ++i) {
        record[i]=cv[i];
        record[nargs+i]=thissigma[i];
        record[2*nargs+i]=height[i];
      }
      record[3*nargs]=(selected ? 1.0 : 0.0);
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
        if(rec[3*nargs]==0.0) {
          continue;
        }
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
          if(mw_n_>1) {
            // the hills of the other walkers in the grid (see the constructor), in a comment after the grid
            family.gridOfile->printf("%s",family.hillsfiles->getGridComment(mw_id_).c_str());
          }
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
