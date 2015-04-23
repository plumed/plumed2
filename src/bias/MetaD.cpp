/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "tools/KernelFunctions.h"
#include "core/FlexibleBin.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include "time.h"
#include <iostream>
#include <limits>
#include <deque>

#define DP2CUTOFF 6.25

using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS METAD
/*
Used to performed MetaDynamics on one or more collective variables.

In a metadynamics simulations a history dependent bias composed of
intermittently added Gaussian functions is added to the potential \cite metad.

\f[
V(\vec{s},t) = \sum_{ k \tau < t} W(k \tau)
\exp\left(
-\sum_{i=1}^{d} \frac{(s_i-s_i^{(0)}(k \tau))^2}{2\sigma_i^2}
\right).
\f]

This potential forces the system away from the kinetic traps in the potential energy surface
and out into the unexplored parts of the energy landscape. Information on the Gaussian
functions from which this potential is composed is output to a file called HILLS, which
is used both the restart the calculation and to reconstruct the free energy as a function of the CVs.
The free energy can be reconstructed from a metadynamics calculation because the final bias is given
by:

\f[
V(\vec{s}) = -F(\vec(s))
\f]

During post processing the free energy can be calculated in this way using the \ref sum_hills
utility.

In the simplest possible implementation of a metadynamics calculation the expense of a metadynamics
calculation increases with the length of the simulation as one has to, at every step, evaluate
the values of a larger and larger number of Gaussians. To avoid this issue you can in plumed 2.0
store the bias on a grid.  This approach is similar to that proposed in \cite babi+08jcp but has the
advantage that the grid spacing is independent on the Gaussian width.
Notice that you should
provide either the number of bins for every collective variable (GRID_BIN) or
the desired grid spacing (GRID_SPACING). In case you provide both PLUMED will use
the most conservative choice (highest number of bins) for each dimension.
In case you do not provide any information about bin size (neither GRID_BIN nor GRID_SPACING)
and if Gaussian width is fixed PLUMED will use 1/5 of the Gaussian width as grid spacing.
This default choice should be reasonable for most applications.

Another option that is available in plumed 2.0 is well-tempered metadynamics \cite Barducci:2008. In this
varient of metadynamics the heights of the Gaussian hills are rescaled at each step so the bias is now
given by:

\f[
V({s},t)= \sum_{t'=0,\tau_G,2\tau_G,\dots}^{t'<t} W e^{-V({s}({q}(t'),t')/\Delta T} \exp\left(
-\sum_{i=1}^{d} \frac{(s_i({q})-s_i({q}(t'))^2}{2\sigma_i^2}
\right),
\f]

This method ensures that the bias converges more smoothly. It should be noted that, in the case of well-tempered metadynamics, in
the output printed the Gaussian height is re-scaled using the bias factor.
Also notice that with well-tempered metadynamics the HILLS file does not contain the bias,
but the negative of the free-energy estimate. This choice has the advantage that
one can restart a simulation using a different value for the \f$\Delta T\f$. The applied bias will be scaled accordingly.

Note that you can use here also the flexible gaussian approach  \cite Branduardi:2012dl
in which you can adapt the gaussian to the extent of Cartesian space covered by a variable or
to the space in collective variable covered in a given time. In this case the width of the deposited
gaussian potential is denoted by one value only that is a Cartesian space (ADAPTIVE=GEOM) or a time
(ADAPTIVE=DIFF). Note that a specific integration technique for the deposited gaussians
should be used in this case. Check the documentation for utility sum_hills.

With the keyword INTERVAL one changes the metadynamics algorithm setting the bias force equal to zero
outside boundary \cite baftizadeh2012protein. If, for example, metadynamics is performed on a CV s and one is interested only
to the free energy for s > sw, the history dependent potential is still updated according to the above
equations but the metadynamics force is set to zero for s < sw. Notice that Gaussians are added also
if s < sw, as the tails of these Gaussians influence VG in the relevant region s > sw. In this way, the
force on the system in the region s > sw comes from both metadynamics and the force field, in the region
s < sw only from the latter. This approach allows obtaining a history-dependent bias potential VG that
fluctuates around a stable estimator, equal to the negative of the free energy far enough from the
boundaries. Note that:
- It works only for one-dimensional biases;
- It works both with and without GRID;
- The interval limit sw in a region where the free energy derivative is not large;
- If in the region outside the limit sw the system has a free energy minimum, the INTERVAL keyword should
  be used together with a \ref UPPER_WALLS or \ref LOWER_WALLS at sw.

As a final note, since version 2.0.2 when the system is outside of the selected interval the force
is set to zero and the bias value to the value at the corresponding boundary. This allows acceptances
for replica exchange methods to be computed correctly.

Multiple walkers  \cite multiplewalkers can also be used. See below the examples.

Additional material and examples can be also found in the tutorials:

- \ref belfast-6
- \ref belfast-7
- \ref belfast-8

\par Examples
The following input is for a standard metadynamics calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the metadynamics bias potential are written to the COLVAR file every 100 steps.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 LABEL=restraint
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endverbatim
(See also \ref DISTANCE \ref PRINT).

\par
If you use adaptive Gaussians, with diffusion scheme where you use
a Gaussian that should cover the space of 20 timesteps in collective variables.
Note that in this case the histogram correction is needed when summing up hills.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=20 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=DIFF
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endverbatim

\par
If you use adaptive Gaussians, with geometrical scheme where you use
a Gaussian that should cover the space of 0.05 nm in Cartesian space.
Note that in this case the histogram correction is needed when summing up hills.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=GEOM
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endverbatim

\par
When using adaptive Gaussians you might want to limit how the hills width can change.
You can use SIGMA_MIN and SIGMA_MAX keywords.
The sigmas should specified in terms of CV so you should use the CV units.
Note that if you use a negative number, this means that the limit is not set.
Note also that in this case the histogram correction is needed when summing up hills.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ...
  ARG=d1,d2 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=GEOM
  SIGMA_MIN=0.2,0.1 SIGMA_MAX=0.5,1.0
... METAD
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endverbatim

\par
Multiple walkers can be also use as in  \cite multiplewalkers
These are enabled by setting the number of walker used, the id of the
current walker which interprets the input file, the directory where the
hills containing files resides, and the frequency to read the other walkers.
Here is an example
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
METAD ...
   ARG=d1 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint
   WALKERS_N=10
   WALKERS_ID=3
   WALKERS_DIR=../
   WALKERS_RSTRIDE=100
... METAD
\endverbatim
where  WALKERS_N is the total number of walkers, WALKERS_ID is the
id of the present walker (starting from 0 ) and the WALKERS_DIR is the directory
where all the walkers are located. WALKERS_RSTRIDE is the number of step between
one update and the other.

\par
The kinetics of the transitions between basins can also be analysed on the fly as
in \cite PRL230602. The flag ACCELERATION turn on accumulation of the acceleration
factor that can then be used to determine the rate. This method can be used together
with \ref COMMITTOR analysis to stop the simulation when the system get to the target basin.
It must be used together with a defined temperature from the simulation engine or the TEMP 
keyword.

*/
//+ENDPLUMEDOC

class MetaD : public Bias {

 private:
  struct Gaussian {
    /// Hill center in CV coordinates
    vector<double> center;
    /// The height of the hill
    double height;
    /// Whether the hill is diagonal (false) or generically symmetric (true)
    bool multivariate;
    /// Hill shape: if diagonal, the width along each coordinate, otherwise a symmetric matrix.
    vector<double> sigma;
    vector<double> invsigma;
    Gaussian(const vector<double> &center, const vector<double> &sigma, double height, bool multivariate):
      center(center), sigma(sigma), height(height), multivariate(multivariate), invsigma(sigma) {
      // to avoid troubles from zero element in flexible hills
      for (unsigned i = 0; i < invsigma.size(); ++i) {
        abs(invsigma[i]) > 1.e-20 ? invsigma[i] = 1.0 / invsigma[i] : 0.;
      }
    }
  };
  vector<double> sigma0_;
  vector<double> sigma0min_;
  vector<double> sigma0max_;
  vector<Gaussian> hills_;
  OFile hillsOfile_;
  OFile gridfile_;
  Grid* BiasGrid_;
  Grid* ExtGrid_;
  bool storeOldGrids_;
  std::string gridfilename_, gridreadfilename_;
  int wgridstride_;
  bool grid_, hasextgrid_;
  double height0_;
  double biasf_;
  double kbt_;
  int stride_;
  bool welltempered_;
  double wt_biasthreshold_;
  bool transitiontempered_;
  double tt_biasf_;
  double tt_biasthreshold_;
  vector<vector<double> > transitionwells_;
  double tt_alpha_;
  bool benthic_toleration_;
  double benthic_tol_number_;
  bool benthic_erosion_;
  double benthic_erosive_time_;
  double last_benthic_erosion_;
  int benthic_histo_stride_;
  vector<double> benthic_histo_bandwidth_;
  Grid *BenthicHistogram_;
  bool use_domains_;
  bool scale_new_hills_;
  bool use_whole_grid_domain_;
  string domainsreadfilename_;
  string domains_histo_readfilename_;
  unsigned n_domains_;
  vector<unsigned> domain_ids_;
  vector<vector<unsigned> > domain_boundary_pts_;
  vector<double> domain_boundary_grid_levels_;
  double domain_implicit_bias_level_;
  vector<Grid*> HillScalingGrids_;
  vector<OFile> regionfiles_;
  bool use_adaptive_domains_;
  bool delay_adaptive_domains_;
  enum AdaptiveDomainRefType {kTransitionRef, kMinRef};
  AdaptiveDomainRefType adaptive_domains_reftype_;
  double adaptive_domains_eoffset_;
  int adaptive_domains_stride_;
  double adaptive_domains_eincrement_;
  double adaptive_domains_last_elevel_;
  double adaptive_domains_downsampling_;
  int domains_histo_stride_;
  vector<double> domains_histo_bandwidth_;
  Grid *DomainsHistogram_;
  string whistofilename_;
  OFile whistofile_;
  bool print_domains_scaling_;
  vector<OFile*> DomainsScalingFilePs_;
  bool print_adaptive_domains_energies_;
  OFile HistEnergyFile_;
  OFile HistBiasEnergyFile_;
  double* dp_;
  int adaptive_;
  FlexibleBin *flexbin;
  int mw_n_;
  string mw_dir_;
  int mw_id_;
  int mw_rstride_;
  bool walkers_mpi;
  bool acceleration;
  bool calc_average_bias_coft_;
  double average_bias_coft_;
  double acc;
  vector<IFile*> ifiles;
  vector<string> ifilesnames;
  double uppI_;
  double lowI_;
  bool doInt_;
  bool isFirstStep;

  void   readGaussians(IFile*);
  bool   readChunkOfGaussians(IFile *ifile, unsigned n);
  void   writeGaussian(const Gaussian &, OFile &);
  void   addGaussian(const Gaussian &);
  void   addGaussianToGrid(const Gaussian &, Grid * const);
  double getHeight(const vector<double> &);
  double getTransitionBarrierBias();
  void   defineDomains(const Grid * const);
  void   filterDomains();
  void   spreadDomain(unsigned, unsigned, const Grid * const);
  void   applyDomainBoundaryBias(unsigned, double);

  void   createScalingGrids();
  bool   shouldAdaptDomainsNow();
  void   adaptDomains();
  double getBiasAndDerivatives(const vector<double> &, double* der = NULL);
  double evaluateGaussian(const vector<double> &, const Gaussian &, double* der = NULL);
  void   finiteDifferenceGaussian(const vector<double> &, const Gaussian &);
  vector<unsigned> getGaussianSupport(const Gaussian &);
  bool   scanOneHill(IFile *ifile,  vector<Value> &v, vector<double> &center, vector<double>  &sigma, double &height, bool &multivariate);
  std::string fmt;

  void   dumpBias();
  void   dumpGrid(Grid *, OFile &);

 public:
  MetaD(const ActionOptions &);
  ~MetaD();
  void calculate();
  void update();
  static void registerKeywords(Keywords &keys);
  bool checkNeedsGradients()const {
    if (adaptive_ == FlexibleBin::geometry) {
      return true;
    } else {
      return false;
    }
  }
};

PLUMED_REGISTER_ACTION(MetaD, "METAD")

void MetaD::registerKeywords(Keywords &keys) {
  Bias::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias", "default", "the instantaneous value of the bias potential");
  keys.addOutputComponent("acc", "ACCELERATION", "the metadynamics acceleration factor");
  keys.addOutputComponent("coft", "CALC_AVERAGE_BIAS", "the metadynamics average bias c(t)");
  keys.use("ARG");
  keys.add("compulsory", "SIGMA", "the widths of the Gaussian hills");
  keys.add("compulsory", "PACE", "the frequency for hill addition");
  keys.add("compulsory", "FILE", "HILLS", "a file in which the list of added hills is stored");
  keys.add("optional", "HEIGHT", "the heights of the Gaussian hills. Compulsory unless TAU, TEMP and BIASFACTOR are given");
  keys.add("optional", "FMT", "specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional", "BIASFACTOR", "use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional", "WTBIASTHRESHOLD", "use well tempered metadynamics with this bias threshold.  Please note you must also specify BIASFACTOR");
  keys.add("optional", "TTBIASFACTOR", "use transition tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional", "TTBIASTHRESHOLD", "use transition tempered metadynamics with this bias threshold.  Please note you must also specify TTBIASFACTOR");
  keys.add("numbered", "TRANSITIONWELL", "This keyword appears multiple times as TRANSITIONWELLx with x=0,1,2,...,n. Each specifies the coordinates for one well in transition-tempered metadynamics. At least one must be provided.");
  keys.add("optional", "TTALPHA", "use transition tempered metadynamics with this decay shape parameter value between 0.5 and 1.0 (default 0.5).  Please note you must also specify TTBIASFACTOR");
  keys.add("optional", "BENTHIC_TOLERATION", "use benthic metadynamics with this number of mistakes tolerated in transition states");
  keys.add("optional", "BENTHIC_FILTER_STRIDE", "use benthic metadynamics accumulating filter samples on this timescale in units of simulation time.  Please note you must also specify BENTHIC_TOLERATION");
  keys.add("optional", "BENTHIC_EROSION", "use benthic metadynamics with erosion on this boosted timescale in units of simulation time.  Please note you must also specify BENTHIC_TOLERATION");
  keys.addFlag("USE_DOMAINS", false, "use metabasin metadynamics");
  keys.add("optional", "REGION_RFILE", "use metabasin metadynamics with this file defining areas to flatten");
  keys.addFlag("WHOLE_GRID_DOMAIN", false, "use metabasin metadynamics with the entire grid as the domain to flatten");
  keys.addFlag("USE_ADAPTIVE_DOMAINS", false, "use metabasin metadynamics with adaptively set regions");
  keys.add("optional", "ADAPTIVE_DOMAINS_HISTOGRAM_RFILE", "read in a histogram file for restarting an adaptive domain calculation");
  keys.add("optional", "ADAPTIVE_DOMAINS_HISTOGRAM_WFILE", "print histogram files for restarting an adaptive domain calculation");
  keys.addFlag("DELAY_ADAPTIVE_DOMAINS", false, "use adaptive metabasin metadynamics only after the reference bias crosses the offset.");
  keys.add("optional", "ADAPTIVE_DOMAINS_REFERENCE", "set metabasin metadynamics adaptive regions with reference to either the 'transition' free energy or the 'minimum' free energy");
  keys.add("optional", "ADAPTIVE_DOMAINS_ENERGY_OFFSET", "use adaptive metabasin metadynamics with regions below a free energy that is the reference free energy plus this value");
  keys.add("optional", "ADAPTIVE_DOMAINS_STRIDE", "use adaptive metabasin metadynamics with regions adapted every this number of steps");  
  keys.add("optional", "ADAPTIVE_DOMAINS_ENERGY_INCREMENT", "use adaptive metabasin metadynamics with regions adapted every this increase of the bias level");
  keys.add("optional", "ADAPTIVE_DOMAINS_DOWNSAMPLING", "when creating scaling grids, add this ratio fewer hills to minimize performance costs");
  keys.addFlag("PRINT_DOMAINS_SCALING", false, "print out the scaling functions for each metabasin metadynamics domain");
  keys.addFlag("PRINT_ADAPTIVE_DOMAINS_ENERGIES", false, "print out the scaling functions for each metabasin metadynamics domain");
  keys.add("optional", "TEMP", "the system temperature - this is only needed if you are doing well-tempered metadynamics, transition-tempered metadynamics, acceleration, or adaptive metabasin metadynamics");
  keys.add("optional", "TAU", "in well tempered metadynamics, sets height to (kb*DeltaT*pace*timestep)/tau");
  keys.add("optional", "GRID_MIN", "the lower bounds for the grid");
  keys.add("optional", "GRID_MAX", "the upper bounds for the grid");
  keys.add("optional", "GRID_BIN", "the number of bins for the grid");
  keys.add("optional", "GRID_SPACING", "the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.addFlag("GRID_SPARSE", false, "use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE", false, "don't use spline interpolation with grids");
  keys.add("optional", "GRID_WSTRIDE", "write the grid to a file every N steps");
  keys.add("optional", "GRID_WFILE", "the file on which to write the grid");
  keys.addFlag("STORE_GRIDS", false, "store all the grid files the calculation generates. They will be deleted if this keyword is not present");
  keys.add("optional", "ADAPTIVE", "use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or timestep dimensions");
  keys.add("optional", "WALKERS_ID", "walker id");
  keys.add("optional", "WALKERS_N", "number of walkers");
  keys.add("optional", "WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional", "WALKERS_RSTRIDE", "stride for reading hills files");
  keys.add("optional", "INTERVAL", "monodimensional lower and upper limits, outside the limits the system will not feel the biasing force.");
  keys.add("optional", "GRID_RFILE", "a grid file from which the bias should be read at the initial step of the simulation");
  keys.add("optional", "SIGMA_MAX", "the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional", "SIGMA_MIN", "the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.addFlag("WALKERS_MPI", false, "Switch on MPI version of multiple walkers - not compatible with other WALKERS_* options");
  keys.addFlag("ACCELERATION", false, "Set to TRUE if you want to compute the metadynamics acceleration factor.");
  keys.addFlag("CALC_AVERAGE_BIAS", false, "Set to TRUE if you want to compute the metadynamics average bias, c(t).");
}

MetaD::~MetaD() {
  if (flexbin) {
    delete flexbin;
  }
  if (BiasGrid_) {
    delete BiasGrid_;
  }
  if (ExtGrid_) {
    delete ExtGrid_;
  }
  if (benthic_toleration_) {
    delete BenthicHistogram_;
  }

  if (use_domains_) {
    if (scale_new_hills_) {
      for (unsigned i = 0; i < n_domains_; i++) {
        if (HillScalingGrids_[i]) delete HillScalingGrids_[i];
      }
    }
    if (print_domains_scaling_) {
      for (unsigned i = 0; i < DomainsScalingFilePs_.size(); i++) {
        DomainsScalingFilePs_[i]->close();
        delete DomainsScalingFilePs_[i];
      }
    }
    if (use_adaptive_domains_) {
      delete DomainsHistogram_;
    }
    if (wgridstride_ > 0 && whistofilename_.size() > 0) {
      whistofile_.close();
    }
    if (print_adaptive_domains_energies_) {
      HistEnergyFile_.close();
      HistBiasEnergyFile_.close();
    }
  }
  hillsOfile_.close();
  if (wgridstride_ > 0) {
    gridfile_.close();
  }
  delete [] dp_;
  // Close open hills files for each other walker
  for (int i = 0; i < mw_n_; ++i) {
    if (ifiles[i]->isOpen()) {
      ifiles[i]->close();
    }
    delete ifiles[i];
  }
}

MetaD::MetaD(const ActionOptions &ao):
  PLUMED_BIAS_INIT(ao),
  // Grid stuff default initialization
  BiasGrid_(NULL), ExtGrid_(NULL), wgridstride_(0), grid_(false), hasextgrid_(false),
  // Metadynamics basic parameters
  height0_(std::numeric_limits<double>::max()), biasf_(1.0), kbt_(0.0),
  stride_(0), welltempered_(false),
  wt_biasthreshold_(0.0),
  transitiontempered_(false),
  tt_biasf_(1.0),
  tt_biasthreshold_(0.0),
  tt_alpha_(0.5),
  benthic_toleration_(false),
  benthic_tol_number_(0.0),
  benthic_erosion_(false),
  benthic_erosive_time_(0.0),
  last_benthic_erosion_(0.0),
  benthic_histo_stride_(0),
  BenthicHistogram_(NULL),
  use_domains_(false),
  scale_new_hills_(false),
  use_whole_grid_domain_(false),
  n_domains_(0),
  domain_implicit_bias_level_(0.0),
  use_adaptive_domains_(false),
  delay_adaptive_domains_(false),
  adaptive_domains_reftype_(kMinRef),
  adaptive_domains_eoffset_(0.0),
  adaptive_domains_stride_(0),
  adaptive_domains_eincrement_(0.0),
  adaptive_domains_last_elevel_(0.0),
  adaptive_domains_downsampling_(1.0),
  domains_histo_stride_(0),
  DomainsHistogram_(NULL),
  print_domains_scaling_(false),
  DomainsScalingFilePs_(vector<OFile*>()),
  print_adaptive_domains_energies_(false),
  // Other stuff
  dp_(NULL), adaptive_(FlexibleBin::none),
  flexbin(NULL),
  // Multiple walkers initialization
  mw_n_(1), mw_dir_("./"), mw_id_(0), mw_rstride_(1),
  walkers_mpi(false),
  acceleration(false), acc(0.0),
  calc_average_bias_coft_(false), average_bias_coft_(0.0),
  // Interval initialization
  uppI_(-1), lowI_(-1), doInt_(false),
  // Event clock initialization
  isFirstStep(true)
{
  // Set the hills file name.
  string hillsfname = "HILLS";
  parse("FILE", hillsfname);

  // Set the hill adaptivity type
  string adaptiveoption;
  adaptiveoption = "NONE";
  parse("ADAPTIVE", adaptiveoption);
  if (adaptiveoption == "NONE") {
    adaptive_ = FlexibleBin::none;
  } else if (adaptiveoption == "GEOM") {
    log.printf("  Uses Geometry-based hills width: sigma must be in distance units and only one sigma is needed\n");
    adaptive_ = FlexibleBin::geometry;
  } else if (adaptiveoption == "DIFF") {
    log.printf("  Uses Diffusion-based hills width: sigma must be in timesteps and only one sigma is needed\n");
    adaptive_ = FlexibleBin::diffusion;
  } else {
    error("I do not know this type of adaptive scheme");
  }

  // Set the base sigma, i.e. hill shape, parameters
  // Pay careful attention to the definitions in the adaptive hill case--
  // these are not necessarily CV-valued widths.
  parseVector("SIGMA", sigma0_);
  parse("FMT", fmt);
  if (adaptive_ == FlexibleBin::none) {
    // if you use normal sigma you need one sigma per argument
    if (sigma0_.size() != getNumberOfArguments()) {
      error("number of arguments does not match number of SIGMA parameters");
    }
  } else {
    // If you use flexible hills you need one sigma in an atomistic length unit.
    if (sigma0_.size() != 1) {
      error("If you choose ADAPTIVE you need only one sigma according to your choice of type (GEOM/DIFF)");
    }
    // If adaptive diffusion-based hills are used then the number 
    // must be an integer number of timesteps.
    if (adaptive_ == FlexibleBin::diffusion) {
      if (int(sigma0_[0]) - sigma0_[0] > 1.e-9 || int(sigma0_[0]) - sigma0_[0] < -1.e-9 || int(sigma0_[0]) < 1) {
        error("In case of adaptive hills with diffusion, the sigma must be an integer which is the number of timesteps\n");
      }
    }
    // here evtl parse the sigma min and max values
    parseVector("SIGMA_MIN", sigma0min_);
    if (sigma0min_.size() > 0 && sigma0min_.size() < getNumberOfArguments()) {
      error("the number of SIGMA_MIN values be at least the number of the arguments");
    } else if (sigma0min_.size() == 0) {
      sigma0min_.resize(getNumberOfArguments());
      for (unsigned i = 0; i < getNumberOfArguments(); i++) {
        sigma0min_[i] = -1.;
      }
    }
    parseVector("SIGMA_MAX", sigma0max_);
    if (sigma0max_.size() > 0 && sigma0max_.size() < getNumberOfArguments()) {
      error("the number of SIGMA_MAX values be at least the number of the arguments");
    } else if (sigma0max_.size() == 0) {
      sigma0max_.resize(getNumberOfArguments());
      for (unsigned i = 0; i < getNumberOfArguments(); i++) {
        sigma0max_[i] = -1.;
      }
    }
    flexbin = new FlexibleBin(adaptive_, this, sigma0_[0], sigma0min_, sigma0max_);
  }
  
  // Set the temperature.
  double temp = 0.0;
  parse("TEMP", temp);
  // If a temp is specified, use it.
  if (temp > 0.0) {
    kbt_ = plumed.getAtoms().getKBoltzmann() * temp;
  // Otherwise ask for the temperature from the MD engine.
  // If the MD engine does not have a set temperature, this
  // quantity will remain zero and no error will be recorded.
  } else {
    kbt_ = plumed.getAtoms().getKbT();
  }

  // Set well tempering parameters.
  parse("BIASFACTOR", biasf_);
  if (biasf_ < 1.0) {
    error("well tempered bias factor is nonsensical");
  }
  if (biasf_ > 1.0) {
    if (kbt_ == 0.0) {
      error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    }
    welltempered_ = true;
    parse("WTBIASTHRESHOLD", wt_biasthreshold_);
    if (wt_biasthreshold_ < 0.0) {
      error("well tempered bias threshold is nonsensical");
    }
  }

  // Set transition tempering parameters.
  parse("TTBIASFACTOR", tt_biasf_);
  if (tt_biasf_ < 1.0) {
    error("transition tempered bias factor is nonsensical");
  }
  if (tt_biasf_ > 1.0) {
    if (kbt_ == 0.0) {
      error("Unless the MD engine passes the temperature to plumed, with transition-tempered metad you must specify it using TEMP");
    }
    transitiontempered_ = true;
    parse("TTBIASTHRESHOLD", tt_biasthreshold_);
    if (tt_biasthreshold_ < 0.0) {
      error("transition tempered bias threshold is nonsensical");
    }
    parse("TTALPHA", tt_alpha_);
    if (tt_alpha_ < 0.5 || tt_alpha_ > 1.0) {
      error("transition tempered decay shape parameter alpha is nonsensical");
    }
    vector<double> tempcoords(getNumberOfArguments());    
    for (unsigned i = 0; ; i++) {
      if (!parseNumberedVector("TRANSITIONWELL", i, tempcoords) ) break;
      if (tempcoords.size() != getNumberOfArguments()) {
        error("incorrect number of coordinates for transition tempering well");
      }
      transitionwells_.push_back(tempcoords);
    }
  }

  // Check for a benthic metadynamics toleration threshold.
  // Wait to set the energy threshold until the hill height is parsed.
  parse("BENTHIC_TOLERATION", benthic_tol_number_);
  if (benthic_tol_number_ > 0.0) {
    benthic_toleration_ = true;
    // Check for a benthic metadynamics filter time.
    parse("BENTHIC_FILTER_STRIDE", benthic_histo_stride_);
    // Check for a benthic metadynamics erosion time.
    parse("BENTHIC_EROSION", benthic_erosive_time_);
    if (benthic_erosive_time_ > 0.0) {
      benthic_erosion_ = true;
    }
  }

  // Set initial bias deposition rate parameters.
  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT", height0_);
  parse("PACE", stride_);
  if (stride_ <= 0) {
    error("frequency for hill addition is nonsensical");
  }
  double tau = 0.0;
  parse("TAU", tau);
  if (tau == 0.0) {
    if (height0_ == std::numeric_limits<double>::max()) {
      error("At least one between HEIGHT and TAU should be specified");
    }
    // if tau is not set, we compute it here from the other input parameters
    if (welltempered_) {
      tau = (kbt_ * (biasf_ - 1.0)) / height0_ * getTimeStep() * stride_;
    }
  } else {
    if (!welltempered_) {
      error("TAU only makes sense in tempered metadynamics");
    }
    if (height0_ != std::numeric_limits<double>::max()) {
      error("At most one between HEIGHT and TAU should be specified");
    }
    height0_ = (kbt_ * (biasf_ - 1.0)) / tau * getTimeStep() * stride_;
  }

  // Set grid parameters.
  vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN", gmin);
  if (gmin.size() != getNumberOfArguments() && gmin.size() != 0) {
    error("not enough values for GRID_MIN");
  }
  vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX", gmax);
  if (gmax.size() != getNumberOfArguments() && gmax.size() != 0) {
    error("not enough values for GRID_MAX");
  }
  vector<unsigned> gbin(getNumberOfArguments());
  vector<double>   gspacing;
  parseVector("GRID_BIN", gbin);
  if (gbin.size() != getNumberOfArguments() && gbin.size() != 0) {
    error("not enough values for GRID_BIN");
  }
  parseVector("GRID_SPACING", gspacing);
  if (gspacing.size() != getNumberOfArguments() && gspacing.size() != 0) {
    error("not enough values for GRID_SPACING");
  }
  if (gmin.size() != gmax.size()) {
    error("GRID_MAX and GRID_MIN should be either present or absent");
  }
  if (gspacing.size() != 0 && gmin.size() == 0) {
    error("If GRID_SPACING is present also GRID_MIN should be present");
  }
  if (gbin.size() != 0     && gmin.size() == 0) {
    error("If GRID_SPACING is present also GRID_MIN should be present");
  }
  if (gmin.size() != 0) {
    if (gbin.size() == 0 && gspacing.size() == 0) {
      if (adaptive_ == FlexibleBin::none) {
        log << "  Binsize not spacified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0_.size() == getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for (unsigned i = 0; i < gspacing.size(); i++) {
          gspacing[i] = 0.2 * sigma0_[i];
        }
      } else {
        // with adaptive hills and grid a sigma min must be specified
        if (sigma0min_.size() == 0) {
          error("When using Adaptive Gaussians on a grid SIGMA_MIN must be specified");
        }
        log << "  Binsize not spacified, 1/5 of sigma_min will be be used\n";
        plumed_assert(sigma0_.size() == getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for (unsigned i = 0; i < gspacing.size(); i++) {
          gspacing[i] = 0.2 * sigma0min_[i];
        }
        //error("At least one among GRID_BIN and GRID_SPACING should be used");
      }
    } else if (gspacing.size() != 0 && gbin.size() == 0) {
      log << "  The number of bins will be estimated from GRID_SPACING\n";
    } else if (gspacing.size() != 0 && gbin.size() != 0) {
      log << "  You specified both GRID_BIN and GRID_SPACING\n";
      log << "  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if (gbin.size() == 0) {
      gbin.assign(getNumberOfArguments(), 1);
    }
    if (gspacing.size() != 0) {
      for (unsigned i = 0; i < getNumberOfArguments(); i++) {
        double a, b;
        Tools::convert(gmin[i], a);
        Tools::convert(gmax[i], b);
        unsigned n = ((b - a) / gspacing[i]);
        if (gbin[i] < n) {
          gbin[i] = n;
        }
      }
    }
  }
  bool sparsegrid = false;
  parseFlag("GRID_SPARSE", sparsegrid);
  bool nospline = false;
  parseFlag("GRID_NOSPLINE", nospline);
  bool spline = !nospline;
  if (gbin.size() > 0) {
    grid_ = true;
  }
  parse("GRID_WSTRIDE", wgridstride_);
  parse("GRID_WFILE", gridfilename_);
  parseFlag("STORE_GRIDS", storeOldGrids_);
  if (grid_ && gridfilename_.length() > 0) {
    if (wgridstride_ == 0) {
      error("frequency with which to output grid not specified use GRID_WSTRIDE");
    }
  }
  if (grid_ && wgridstride_ > 0) {
    if (gridfilename_.length() == 0) {
      error("grid filename not specified use GRID_WFILE");
    }
  }
  parse("GRID_RFILE", gridreadfilename_);
  
  // Set metabasin metadynamics parameters
  parseFlag("USE_DOMAINS", use_domains_);
  if (use_domains_) {
    parse("REGION_RFILE", domainsreadfilename_);
    parseFlag("WHOLE_GRID_DOMAIN", use_whole_grid_domain_);
    parseFlag("USE_ADAPTIVE_DOMAINS", use_adaptive_domains_);
    if (use_adaptive_domains_) {
      if (kbt_ == 0.0) {
        error("Unless the MD engine passes the temperature to plumed, with adaptive domains metabasin metad you must specify it using TEMP");
      }
      parse("ADAPTIVE_DOMAINS_HISTOGRAM_RFILE", domains_histo_readfilename_);
      parse("ADAPTIVE_DOMAINS_HISTOGRAM_WFILE", whistofilename_);
      parseFlag("DELAY_ADAPTIVE_DOMAINS", delay_adaptive_domains_);
      string refstring;
      parse("ADAPTIVE_DOMAINS_REFERENCE", refstring);
      if (refstring == "minimum") {
        adaptive_domains_reftype_ = kMinRef;
      } else if (refstring == "transition") {
        adaptive_domains_reftype_ = kTransitionRef;
      } else {
        error("unrecognized adaptive domains reference type for metabasin metadynamics");
      }
      parse("ADAPTIVE_DOMAINS_ENERGY_OFFSET", adaptive_domains_eoffset_);
      if (adaptive_domains_reftype_ == kMinRef && adaptive_domains_eoffset_ <= 0.0) {
        error("adaptive domains energy offset must be positive for minimum-referenced to make sense");
      }
      parse("ADAPTIVE_DOMAINS_STRIDE", adaptive_domains_stride_);
      if (adaptive_domains_stride_ <= 0) {
        error("adaptive domains requires a stride value greater than zero");
      }
      parse("ADAPTIVE_DOMAINS_ENERGY_INCREMENT", adaptive_domains_eincrement_);
      if (adaptive_domains_eincrement_ < 0.0) {
        error("adaptive domains requires an energy increment value greater than zero");
      }
      parse("ADAPTIVE_DOMAINS_DOWNSAMPLING", adaptive_domains_downsampling_);
      if (adaptive_domains_downsampling_ < 1.0) {
        error("adaptive domains downsampling must be at least one");
      }
      parseFlag("PRINT_ADAPTIVE_DOMAINS_ENERGIES", print_adaptive_domains_energies_);
      if (print_adaptive_domains_energies_) {
        std::string histogram_energy_filename = "histogram_energy.dat";
        HistEnergyFile_.link(*this);
        HistEnergyFile_.open(histogram_energy_filename);
        std::string histbias_energy_filename = "histogram-bias_energy.dat";
        HistBiasEnergyFile_.link(*this);
        HistBiasEnergyFile_.open(histbias_energy_filename);
      }
    }
    parseFlag("PRINT_DOMAINS_SCALING", print_domains_scaling_);
  }

  // Set multiple walkers parameters.
  parse("WALKERS_N", mw_n_);
  parse("WALKERS_ID", mw_id_);
  if (mw_n_ <= mw_id_) {
    error("walker ID should be a numerical value less than the total number of walkers");
  }
  parse("WALKERS_DIR", mw_dir_);
  parse("WALKERS_RSTRIDE", mw_rstride_);
  // MPI version
  parseFlag("WALKERS_MPI", walkers_mpi);
  
  // Set boundary conditions.
  vector<double> tmpI(2);
  parseVector("INTERVAL", tmpI);
  if (tmpI.size() != 2 && tmpI.size() != 0) {
    error("both a lower and an upper limits must be provided with INTERVAL");
  } else if (tmpI.size() == 2) {
    lowI_ = tmpI.at(0);
    uppI_ = tmpI.at(1);
    if (getNumberOfArguments() != 1) {
      error("INTERVAL limits correction works only for monodimensional metadynamics!");
    }
    if (uppI_ < lowI_) {
      error("The Upper limit must be greater than the Lower limit!");
    }
    doInt_ = true;
  }

  // Set special clock options.
  acceleration = false;
  parseFlag("ACCELERATION", acceleration);
  // Set to calculate the average bias.
  parseFlag("CALC_AVERAGE_BIAS", calc_average_bias_coft_);
  checkRead();
  
  // Log all of what has just been set.
  log.printf("  Gaussian width ");
  if (adaptive_ == FlexibleBin::diffusion) {
    log.printf(" (Note: The units of sigma are in timesteps) ");
  }
  if (adaptive_ == FlexibleBin::geometry) {
    log.printf(" (Note: The units of sigma are in dist units) ");
  }
  for (unsigned i = 0; i < sigma0_.size(); ++i) {
    log.printf(" %f", sigma0_[i]);
  }
  log.printf("  Gaussian height %f\n", height0_);
  log.printf("  Gaussian deposition pace %d\n", stride_);
  log.printf("  Gaussian file %s\n", hillsfname.c_str());
  if (welltempered_) {
    log.printf("  Well-Tempered bias factor %f\n", biasf_);
    log.printf("  Hills relaxation time (tau) %f\n", tau);
    log.printf("  KbT %f\n", kbt_);
    if (wt_biasthreshold_ != 0.0) {
      log.printf("  Well-Tempered bias threshold %f\n", wt_biasthreshold_);
    }
  }
  // Transition tempered metadynamics options
  if (transitiontempered_) {
    log.printf("  Transition-Tempered bias factor %f\n", tt_biasf_);
    log.printf("  Transition-Tempered bias threshold %f\n", tt_biasthreshold_);
    log.printf("  Transition-Tempered decay shape parameter alpha %f\n", tt_alpha_);
    log.printf("  Number of transition wells %d\n", transitionwells_.size());
    for (unsigned i = 0; i < transitionwells_.size(); i++) {
      log.printf("  Transition well %d at coordinate ", i);
      for (unsigned j = 0; j < getNumberOfArguments(); j++) {
        log.printf("%f ", transitionwells_[i][j]);
      }
      log.printf("\n", i);
    }
    log.printf("  KbT %f\n", kbt_);
    // Check that a grid is in use.
    if (!grid_) {
      error(" transition tempering requires a grid for the bias");
    }
    // For each dimension, check that the transition well coordinates are in the grid.
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      double max, min;
      Tools::convert(gmin[i], min);
      Tools::convert(gmax[i], max);
      for (unsigned j = 0; j < transitionwells_.size(); j++) {
        if (transitionwells_[j][i] < min || transitionwells_[j][i] > max) {
          error(" transition well is not in grid");
        }
      }
    }
  }
  // Benthic metadynamics options
  if (benthic_toleration_) {
    // Benthic erosion relies on using a grid and an acceleration factor.
    if (!grid_) {
      error(" benthic requires a grid for the bias");
    }
    if (!acceleration) {
      error(" benthic requires calculation of the acceleration factor");
    }
    log.printf("  Benthic number of samples to tolerate %d\n", benthic_tol_number_);
    if (benthic_erosion_) {
      // Log the timescale.
      log.printf("  Benthic erosion on timescale %f \n", benthic_erosive_time_);
    }
    // Set the histogramming parameters for the thresholding if not already
    // specified.
    if (benthic_histo_stride_ == 0) {
      benthic_histo_stride_ = max(1, stride_);
    }
    if (benthic_histo_bandwidth_.size() == 0) {
      benthic_histo_bandwidth_ = vector<double>(getNumberOfArguments());
      for (unsigned i = 0; i < getNumberOfArguments(); i++) {
        benthic_histo_bandwidth_[i] = sigma0_[i];
      }
    }
    log.printf("  Histogram update stride is %d and bandwiths are", domains_histo_stride_);
    for (unsigned i = 0; i < benthic_histo_bandwidth_.size(); i++) {
      log.printf(" %f", benthic_histo_bandwidth_[i]);
    }
    log.printf("\n");
  }

  if (doInt_) {
    log.printf("  Upper and Lower limits boundaries for the bias are activated at %f - %f\n", lowI_, uppI_);
  }
  if (grid_) {
    log.printf("  Grid min");
    for (unsigned i = 0; i < gmin.size(); ++i) {
      log.printf(" %s", gmin[i].c_str());
    }
    log.printf("\n");
    log.printf("  Grid max");
    for (unsigned i = 0; i < gmax.size(); ++i) {
      log.printf(" %s", gmax[i].c_str());
    }
    log.printf("\n");
    log.printf("  Grid bin");
    for (unsigned i = 0; i < gbin.size(); ++i) {
      log.printf(" %d", gbin[i]);
    }
    log.printf("\n");
    if (spline) {
      log.printf("  Grid uses spline interpolation\n");
    }
    if (sparsegrid) {
      log.printf("  Grid uses sparse grid\n");
    }
    if (wgridstride_ > 0) {
      log.printf("  Grid is written on file %s with stride %d\n", gridfilename_.c_str(), wgridstride_);
      if (whistofilename_.size() > 0) {
        log.printf("  Adaptive domains histogram is written on file %s with stride %d\n", whistofilename_.c_str(), wgridstride_);
      }
    }
  }
  if (gridreadfilename_.length() > 0) {
    log.printf("  Reading an additional bias from grid in file %s \n", gridreadfilename_.c_str());
  }

  if (use_domains_) {
    log.printf("  Using metabasin metadynamics \n");      
    if (!grid_ || sparsegrid) {
      error(" using domains requires a dense grid for the bias");
    }
    if (!nospline) {
      error(" using domains is currently numerically unstable when using splines for the bias");
    }
    if (adaptive_ != FlexibleBin::none) {
      error(" using domains only works for fixed, non-adaptive Gaussians");
    }
    if (!(use_whole_grid_domain_ || use_adaptive_domains_ || domainsreadfilename_.size())) {
      error(" must specify exactly one of WHOLE_GRID_DOMAIN, USE_ADAPTIVE_DOMAINS, or REGION_RFILE with USE_DOMAINS");
    }
    if ((use_whole_grid_domain_ && use_adaptive_domains_) || (use_whole_grid_domain_ && domainsreadfilename_.size()) || (domainsreadfilename_.size() && use_adaptive_domains_) ) {
      error(" must specify only one of WHOLE_GRID_DOMAIN, USE_ADAPTIVE_DOMAINS, or REGION_RFILE with USE_DOMAINS");
    }
    if (use_whole_grid_domain_) {
      log.printf("  Using the entire grid as the metabasin metadynamics domain \n");      
    }
    if (domainsreadfilename_.size() > 0) {
      log.printf("  Reading grid in file %s for initial metabasin metadynamics domains \n", domainsreadfilename_.c_str());
    }
    if (use_adaptive_domains_) {
      if (adaptive_domains_reftype_ == kTransitionRef) {
        if (!transitiontempered_) {
          error(" using transition-referenced adaptive domains requires transition tempering");
        } else {
          log.printf("  Using adaptive domain with a free energy limit of the transition barrier free energy");
        }
      } else if (adaptive_domains_reftype_ == kMinRef) {
        log.printf("  Using adaptive domains with a free energy limit of the minimum free energy");
      }
      log.printf(" plus %f \n", adaptive_domains_eoffset_);
      if (domains_histo_readfilename_.size() > 0) {
        if (gridreadfilename_.size() == 0) {
          error( " adaptive domains cannot be restart with only an input histogram, also provide an input bias");
        }
        log.printf("  Reading histogram in file %s for restarting metabasin metadynamics with adaptive domains \n", domains_histo_readfilename_.c_str());
      }
      log.printf("  Domains will be updated every %d steps\n", adaptive_domains_stride_);
      if (adaptive_domains_eincrement_ > 0) {
        log.printf("  After initialization domains will change only if the domains bias level has increased by %d energy units since last change\n", adaptive_domains_eincrement_);        
      }
      log.printf("  Scaling grids will be generated using 1 / %f of the grid points\n", adaptive_domains_downsampling_);
      // Set the histogramming parameters for the future free energy estimates
      // if they are not already specified.
      if (domains_histo_stride_ == 0) {
        domains_histo_stride_ = max(1, stride_ / 5);
      }
      if (domains_histo_bandwidth_.size() == 0) {
        domains_histo_bandwidth_ = vector<double>(getNumberOfArguments());
        for (unsigned i = 0; i < getNumberOfArguments(); i++) {
          domains_histo_bandwidth_[i] = sigma0_[i] / 3.0;
        }
      }
      log.printf("  Histogram update stride is %d and bandwiths are", domains_histo_stride_);
      for (unsigned i = 0; i < domains_histo_bandwidth_.size(); i++) {
        log.printf(" %f", domains_histo_bandwidth_[i]);
      }
      log.printf("\n");
    }
  }

  if (mw_n_ > 1) {
    if (walkers_mpi) {
      error("MPI version of multiple walkers is not compatible with filesystem version of multiple walkers");
    }
    log.printf("  %d multiple walkers active\n", mw_n_);
    log.printf("  walker id %d\n", mw_id_);
    log.printf("  reading stride %d\n", mw_rstride_);
    log.printf("  directory with hills files %s\n", mw_dir_.c_str());
  } else {
    if (walkers_mpi) {
      log.printf("  Multiple walkers active using MPI communnication\n");
    }
  }
  addComponent("bias");
  componentIsNotPeriodic("bias");
  if (acceleration) {
    if (kbt_ == 0.0) {
      error("The calculation of the acceleration works only if simulation temperature has been defined");
    }
    log.printf("  calculation on the fly of the acceleration factor\n");
    addComponent("acc");
    componentIsNotPeriodic("acc");
  }
  if (calc_average_bias_coft_) {
    if (kbt_ == 0.0) {
      error("The calculation of the average bias on the fly works only if simulation temperature has been defined");
    }
    if (!grid_) {
      error("Calculating the average bias on the fly works only with a grid");
    }
    log.printf("  calculation on the fly of the average bias c(t)\n");
    addComponent("coft");
    componentIsNotPeriodic("coft");
  }
  
  // Perform initializations based on the options just set and logged.
  // For performance, allocate memory for the hill distance parameter.
  dp_ = new double[getNumberOfArguments()];
  
  // Initialize and check grid.
  if (grid_) {
    // check for adaptive and sigma_min
    if (sigma0min_.size() == 0 && adaptive_ != FlexibleBin::none) {
      error("When using Adaptive Gaussians on a grid SIGMA_MIN must be specified");
    }
    // check for mesh and sigma size
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      double a, b;
      Tools::convert(gmin[i], a);
      Tools::convert(gmax[i], b);
      double mesh = (b - a) / ((double)gbin[i]);
      if (mesh > 0.5 * sigma0_[i]) {
        log << "  WARNING: Using a METAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
      }
    }
    std::string funcl = getLabel() + ".bias";
    if (!sparsegrid) {
      BiasGrid_ = new Grid(funcl, getArguments(), gmin, gmax, gbin, spline, true);
    } else {
      BiasGrid_ = new SparseGrid(funcl, getArguments(), gmin, gmax, gbin, spline, true);
    }
    std::vector<std::string> actualmin = BiasGrid_->getMin();
    std::vector<std::string> actualmax = BiasGrid_->getMax();
    for (unsigned i = 0; i < getNumberOfArguments(); i++) {
      if (gmin[i] != actualmin[i]) {
        log << "  WARNING: GRID_MIN[" << i << "] has been adjusted to " << actualmin[i] << " to fit periodicity\n";
      }
      if (gmax[i] != actualmax[i]) {
        log << "  WARNING: GRID_MAX[" << i << "] has been adjusted to " << actualmax[i] << " to fit periodicity\n";
      }
    }
  }
  if (wgridstride_ > 0) {
    gridfile_.link(*this);
    gridfile_.open(gridfilename_);
  }

  // Initialize and read an external grid if requested.
  // If the bias uses a grid itself, add the external grid bias to the
  // internal bias; this is necessary for methods that calculate 
  // bias-related quantities on-the-fly in calculation, such as TTMetaD
  // and metabasin metadynamics.
  if (gridreadfilename_.length() > 0) {
    hasextgrid_ = true;
    // Read the grid in input, find the keys.
    IFile gridfile;
    gridfile.open(gridreadfilename_);
    std::string funcl = getLabel() + ".bias";
    ExtGrid_ = Grid::create(funcl, getArguments(), gridfile, false, false, true);
    gridfile.close();
    // Check for consistency between the external grid and the MetaD input specs.
    if (ExtGrid_->getDimension() != getNumberOfArguments()) {
      error("mismatch between dimensionality of input grid and number of arguments");
    }
    for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
      if (getPntrToArgument(i)->isPeriodic() != ExtGrid_->getIsPeriodic()[i]) {
        error("periodicity mismatch between arguments and input bias");
      }
    }
    // Add this to the internal grid if one exists and erase it if so.
    if (grid_) {
      double temp_val = 0.0;
      vector<double> temp_der(getNumberOfArguments());
      vector<double> temp_point(getNumberOfArguments());
      for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
        BiasGrid_->getPoint(i, temp_point);
        temp_val = ExtGrid_->getValueAndDerivatives(temp_point, temp_der);
        if(temp_val != 0.0) {
          BiasGrid_->addValueAndDerivatives(i, temp_val, temp_der);
        }
      }
      delete ExtGrid_;
      ExtGrid_ = NULL;
      hasextgrid_ = false;
    }
  }

  // Initialize an auxiliary histogram for the benthic metadynamics control.
  if (benthic_toleration_) {
    // Use the same parameters as for the main grid, but don't use splines or derivatives.
    std::string funcl = getLabel() + ".number";
    BenthicHistogram_ = new Grid(funcl, getArguments(), gmin, gmax, gbin, false, false);
  }

  log.printf("  Entering USE_DOMAINS initialization.\n"); log.flush();
  if (use_domains_) {
    // Initialize domains based on an initial region if requested.
    if (use_whole_grid_domain_ || domainsreadfilename_.length() > 0) {
      // Initialize and set up a whole grid region if requested.
      if (use_whole_grid_domain_) {
        // Create a copy of the grid.
        Grid *whole_grid_region;
        std::string funcl = getLabel() + ".indicator";
        vector<unsigned> input_nbins = BiasGrid_->getNbin();
        for (unsigned j = 0; j < getNumberOfArguments(); ++j) {
          if (!getPntrToArgument(j)->isPeriodic()) {
            input_nbins[j] -= 1;
          }
        }
        whole_grid_region = new Grid(funcl, getArguments(), BiasGrid_->getMin(), BiasGrid_->getMax(), input_nbins, false, false);
        // Set all values to 1.0.
        for (unsigned i = 0; i < whole_grid_region->getMaxSize(); i++) {
          whole_grid_region->setValue(i, 1.0);
        }
        // Use the region grid to set up domains in the bias grid.
        defineDomains(whole_grid_region);
        // Get rid of the temporary.
        delete whole_grid_region;
      }
      // Initialize and read domains from a region grid if requested.
      if (domainsreadfilename_.length() > 0) {
        // Read the grid in input, find the keys.
        IFile gridfile;
        gridfile.open(domainsreadfilename_);
        std::string funcl = getLabel() + ".indicator";
        Grid* input_region = Grid::create(funcl, getArguments(), gridfile, false, false, false);
        gridfile.close();
        // Check for consistency between the region grid and the MetaD input specs.
        if (input_region->getDimension() != getNumberOfArguments()) {
          error("mismatch between dimensionality of region input grid and number of arguments");
        }
        for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
          if (getPntrToArgument(i)->isPeriodic() != input_region->getIsPeriodic()[i]) {
            error("periodicity mismatch between arguments and region input function");
          }
        }
        // Use the region grid to set up domains in the bias grid.
        defineDomains(input_region);
        // Get rid of the temporary.
        delete input_region;
        if (n_domains_ == 0) {
          error("no domains found in region file " + domainsreadfilename_);
        }
      }
      // Use the domains to set up scaling grids.
      createScalingGrids();
      scale_new_hills_ = true;
      if (print_domains_scaling_) {
        if (DomainsScalingFilePs_.size() < HillScalingGrids_.size()) {
          for (unsigned i = DomainsScalingFilePs_.size(); i < HillScalingGrids_.size(); i++) {
            std::ostringstream filename_stream;
            filename_stream << "domain_" << i << "_scaling.dat";
            std::string scaling_filename = filename_stream.str();
            OFile* scaling_filep = new OFile();
            scaling_filep->link(*this);
            scaling_filep->open(scaling_filename);
            DomainsScalingFilePs_.push_back(scaling_filep);
          }
        }
        for (unsigned i = 0; i < DomainsScalingFilePs_.size(); i++) DomainsScalingFilePs_[i]->rewind();
        for (unsigned i = 0; i < HillScalingGrids_.size(); i++) HillScalingGrids_[i]->writeToFile(*(DomainsScalingFilePs_[i]));
        for (unsigned i = 0; i < DomainsScalingFilePs_.size(); i++) DomainsScalingFilePs_[i]->flush();
      }
    }
  
    // Initialize an auxiliary histogram for adaptive domain metabasin metadynamics.
    // Use the same parameters as for the main grid, but don't use splines or derivatives.
    if (use_adaptive_domains_) {
      // Either initialize by reading in a past histogram or from scratch.
      if (domains_histo_readfilename_.size() > 0) {
        IFile histofile;
        histofile.open(domains_histo_readfilename_);
        std::string funcl = getLabel() + ".number";
        DomainsHistogram_ = Grid::create(funcl, getArguments(), histofile, false, false, false);
        histofile.close();
        // Check for consistency between the external grid and the MetaD input specs.
        if (DomainsHistogram_->getDimension() != getNumberOfArguments()) {
          error("mismatch between dimensionality of input histogram and number of arguments");
        }
        for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
          if (getPntrToArgument(i)->isPeriodic() != DomainsHistogram_->getIsPeriodic()[i]) {
            error("periodicity mismatch between arguments and input histogram");
          }
        }
        adaptDomains();
      } else {
        std::string funcl = getLabel() + ".number";
        DomainsHistogram_ = new Grid(funcl, getArguments(), gmin, gmax, gbin, false, false);
        for (unsigned i = 0; i < DomainsHistogram_->getMaxSize(); i++) {
          DomainsHistogram_->setValue(i, .3);
        }
      }
      if (wgridstride_ > 0 && whistofilename_.size() > 0) {
        whistofile_.link(*this);
        whistofile_.open(whistofilename_);
      }
    }
  }
  log.printf("  Exiting USE_DOMAINS initialization.\n"); log.flush();

  // Create vector of ifile* for hills reading.
  // Open all files at the beginning and read Gaussians if restarting.
  // For multiple walkers, keep all other walkers' files open but close
  // this walker's own file.
  for (int i = 0; i < mw_n_; ++i) {
    string fname;
    // Find the filename for the ith walker from the base filename.
    if (mw_n_ > 1) {
      stringstream out;
      out << i;
      fname = mw_dir_ + "/" + hillsfname + "." + out.str();
    } else {
      fname = hillsfname;
    }
    // Add to the list of filenames and the list of open files.
    IFile *ifile = new IFile();
    ifile->link(*this);
    ifiles.push_back(ifile);
    ifilesnames.push_back(fname);
    if (ifile->FileExist(fname)) {
      ifile->open(fname);
      // Read the Gaussians in the file if and only if this is a restart.
      if (plumed.getRestart()) {
        log.printf("  Restarting from %s:", ifilesnames[i].c_str());
        readGaussians(ifiles[i]);
      }
      ifiles[i]->reset(false);
      // Close only this walker's own hills file for later writing.
      if (i == mw_id_) {
        ifiles[i]->close();
      }
    }
  }
  // Reopen this walker's hills file for writing.
  hillsOfile_.link(*this);
  hillsOfile_.open(ifilesnames[mw_id_]);
  if (fmt.length() > 0) {
    hillsOfile_.fmtField(fmt);
  }
  hillsOfile_.addConstantField("multivariate");
  if (doInt_) {
    hillsOfile_.addConstantField("lower_int").printField("lower_int", lowI_);
    hillsOfile_.addConstantField("upper_int").printField("upper_int", uppI_);
  }
  hillsOfile_.setHeavyFlush();

  // Output periodicities of variables
  for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
    hillsOfile_.setupPrintValue(getPntrToArgument(i));
  }

  // Print citatation suggestions for the specific methods used to the log.
  log << "  Bibliography " << plumed.cite("Laio and Parrinello, PNAS 99, 12562 (2002)");
  if (welltempered_) log << plumed.cite(
                         "Barducci, Bussi, and Parrinello, Phys. Rev. Lett. 100, 020603 (2008)");
  if (transitiontempered_) log << plumed.cite(
                         "Dama, Rotskoff, Parrinello, and Voth, J. Chem. Theory Comput. 10, 3626 (2014)");
  if (mw_n_ > 1 || walkers_mpi) log << plumed.cite(
                                        "Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");
  if (adaptive_ != FlexibleBin::none) log << plumed.cite(
        "Branduardi, Bussi, and Parrinello, J. Chem. Theory Comput. 8, 2247 (2012)");
  if (doInt_) log << plumed.cite(
                      "Baftizadeh, Cossio, Pietrucci, and Laio, Curr. Phys. Chem. 2, 79 (2012)");
  if (acceleration) log << plumed.cite(
                            "Pratyush and Parrinello, Phys. Rev. Lett. 111, 230602 (2013)");
  if (calc_average_bias_coft_) log << plumed.cite(
                                       "Pratyush and Parrinello, J. Phys. Chem. B 119, 736–742 (2015)");
  log << "\n";
  log.flush();
}

void MetaD::readGaussians(IFile *ifile) {
  unsigned ncv = getNumberOfArguments();
  vector<double> center(ncv);
  vector<double> sigma(ncv);
  double height;
  int nhills = 0;
  bool multivariate = false;
  std::vector<Value> tmpvalues;
  for (unsigned j = 0; j < getNumberOfArguments(); ++j) {
    tmpvalues.push_back(Value(this, getPntrToArgument(j)->getName(), false));
  }
  while (scanOneHill(ifile, tmpvalues, center, sigma, height, multivariate)) {
    nhills++;
    if (welltempered_) {
      height *= (biasf_ - 1.0) / biasf_;
    }
    addGaussian(Gaussian(center, sigma, height, multivariate));
  }
  log.printf("      %d Gaussians read\n", nhills);
}

bool MetaD::readChunkOfGaussians(IFile *ifile, unsigned n) {
  unsigned ncv = getNumberOfArguments();
  vector<double> center(ncv);
  vector<double> sigma(ncv);
  double height;
  unsigned nhills = 0;
  bool multivariate = false;
  std::vector<Value> tmpvalues;
  for (unsigned j = 0; j < getNumberOfArguments(); ++j) {
    tmpvalues.push_back(Value(this, getPntrToArgument(j)->getName(), false));
  }
  while (scanOneHill(ifile, tmpvalues, center, sigma, height, multivariate)) {
    if (welltempered_) {
      height *= (biasf_ - 1.0) / biasf_;
    }
    addGaussian(Gaussian(center, sigma, height, multivariate));
    if (nhills == n) {
      log.printf("      %u Gaussians read\n", nhills);
      return true;
    }
    nhills++;
  }
  log.printf("      %u Gaussians read\n", nhills);
  return false;
}

// Write the specification of one Gaussian hill to the hills output file.

void MetaD::writeGaussian(const Gaussian &hill, OFile &file) {
  unsigned ncv = getNumberOfArguments();
  file.printField("time", getTimeStep()*getStep());
  for (unsigned i = 0; i < ncv; ++i) {
    file.printField(getPntrToArgument(i), hill.center[i]);
  }
  if (hill.multivariate) {
    file.printField("multivariate", "true");
    Matrix<double> mymatrix(ncv, ncv);
    unsigned k = 0;
    for (unsigned i = 0; i < ncv; i++) {
      for (unsigned j = i; j < ncv; j++) {
        mymatrix(i, j) = mymatrix(j, i) = hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    // Invert it.
    Matrix<double> invmatrix(ncv, ncv);
    Invert(mymatrix, invmatrix);
    // Enforce symmetry.
    for (unsigned i = 0; i < ncv; i++) {
      for (unsigned j = i; j < ncv; j++) {
        invmatrix(i, j) = invmatrix(j, i);
      }
    }
    // Do a Cholesky decomposition so as to have a "sigma like" number.
    // Sigma like means width-like, a square root of the multivariate Gaussian matrix.
    Matrix<double> lower(ncv, ncv);
    cholesky(invmatrix, lower); // Now this, in band form, is similar to the sigmas.
    // loop in band form
    for (unsigned i = 0; i < ncv; i++) {
      for (unsigned j = 0; j < ncv - i; j++) {
        file.printField("sigma_" + getPntrToArgument(j + i)->getName() + "_" + getPntrToArgument(j)->getName(), lower(j + i, j));
      }
    }
  } else {
    hillsOfile_.printField("multivariate", "false");
    for (unsigned i = 0; i < ncv; ++i) {
      file.printField("sigma_" + getPntrToArgument(i)->getName(), hill.sigma[i]);
    }
  }
  double height = hill.height;
  if (welltempered_) {
    height *= biasf_ / (biasf_ - 1.0);
  }
  file.printField("height", height).printField("biasf", biasf_);
  if (mw_n_ > 1) {
    file.printField("clock", int(time(0)));
  }
  file.printField();
}

double tame_scaling(double scale) {
  if (scale > 1.0) {
    return scale;
  } else {
    return scale + .5 * (1 - scale) * (1 - scale);
  }
}

void MetaD::addGaussian(const Gaussian &hill) {
  // If the height is zero, add no hill.
  if (hill.height == 0.0) {
    return;
  }
  // Add the hill to the list of Gaussians if there is no grid.
  if (!grid_) {
    hills_.push_back(hill);
  // Otherwise, update the bias grid.
  } else {
    unsigned ncv = getNumberOfArguments();
    vector<double> der(ncv);
    unsigned hilldomain = 0;
    if (scale_new_hills_) {
      hilldomain = domain_ids_[BiasGrid_->getIndex(hill.center)];
    }
    vector<double> xx(ncv);
    // Determine which grid points might be updated.
    vector<unsigned> neighbors = BiasGrid_->getNeighbors(hill.center, getGaussianSupport(hill));
    // Evaluate the hill over the grid points serially if single-core.
    if (comm.Get_size() == 1) {
      // On each grid point that may change,
      for (unsigned i = 0; i < neighbors.size(); ++i) {
        // Get the CV coordinate of the grid point.
        unsigned ineigh = neighbors[i];
        for (unsigned j = 0; j < ncv; ++j) {
          der[j] = 0.0;
        }
        BiasGrid_->getPoint(ineigh, xx);
        // Evaluate the Gaussian.
        double bias = evaluateGaussian(xx, hill, &der[0]);
        if (scale_new_hills_ && hilldomain != 0) {
          bias /= tame_scaling(HillScalingGrids_[hilldomain - 1]->getValue(ineigh));
        }
        // Transfer the result to the grid.
        BiasGrid_->addValueAndDerivatives(ineigh, bias, der);
      }
    // Otherwise, evaluate the hill over the grid points in parallel.
    } else {
      unsigned stride = comm.Get_size();
      unsigned rank = comm.Get_rank();
      vector<double> allder(ncv * neighbors.size(), 0.0);
      vector<double> allbias(neighbors.size(), 0.0);
      // For a specific partition set of the grid points,
      for (unsigned i = rank; i < neighbors.size(); i += stride) {
        // Get the CV coordinate of the grid point.
        unsigned ineigh = neighbors[i];
        BiasGrid_->getPoint(ineigh, xx);
        // Evaluate the Gaussian at that CV coordinate.
        allbias[i] = evaluateGaussian(xx, hill, &allder[ncv * i]);
        if (scale_new_hills_ && hilldomain != 0) {
          allbias[i] /= tame_scaling(HillScalingGrids_[hilldomain - 1]->getValue(ineigh));
        }
      }
      // Combine all the Gaussian evaluations together.
      comm.Sum(allbias);
      comm.Sum(allder);
      // Transfer all the evaluations to the grid.
      for (unsigned i = 0; i < neighbors.size(); ++i) {
        unsigned ineigh = neighbors[i];
        for (unsigned j = 0; j < ncv; ++j) {
          der[j] = allder[ncv * i + j];
        }
        BiasGrid_->addValueAndDerivatives(ineigh, allbias[i], der);
      }
    }
    // After adding a hill, correct the derivatives and the 
    // boundary bias.
    if (scale_new_hills_ && hilldomain != 0) {
      // Correct the boundary bias
      if (!use_whole_grid_domain_) applyDomainBoundaryBias(hilldomain - 1, 0.125 * height0_);
      // Set all derivatives.
      for (unsigned i = 0; i < neighbors.size(); ++i) {
        BiasGrid_->setDerivFromValues(neighbors[i]);
      }
    }
  }
}

void MetaD::addGaussianToGrid(const Gaussian &hill, Grid * const grid) {
  unsigned ncv = getNumberOfArguments();
  vector<double> der(ncv);
  vector<double> xx(ncv);
  // Determine which grid points might be updated.
  vector<unsigned> neighbors = grid->getNeighbors(hill.center, getGaussianSupport(hill));
  // Evaluate the hill over the grid points serially if single-core.
  if (comm.Get_size() == 1) {
    // On each grid point that may change,
    for (unsigned i = 0; i < neighbors.size(); ++i) {
      // Get the CV coordinate of the grid point.
      unsigned ineigh = neighbors[i];
      for (unsigned j = 0; j < ncv; ++j) {
        der[j] = 0.0;
      }
      grid->getPoint(ineigh, xx);
      // Evaluate the Gaussian.
      double bias = evaluateGaussian(xx, hill, &der[0]);
      // Transfer the result to the grid.
      grid->addValueAndDerivatives(ineigh, bias, der);
    }
  // Otherwise, evaluate the hill over the grid points in parallel.
  } else {
    unsigned stride = comm.Get_size();
    unsigned rank = comm.Get_rank();
    vector<double> allder(ncv * neighbors.size(), 0.0);
    vector<double> allbias(neighbors.size(), 0.0);
    // For a specific partition set of the grid points,
    for (unsigned i = rank; i < neighbors.size(); i += stride) {
      // Get the CV coordinate of the grid point.
      unsigned ineigh = neighbors[i];
      grid->getPoint(ineigh, xx);
      // Evaluate the Gaussian at that CV coordinate.
      allbias[i] = evaluateGaussian(xx, hill, &allder[ncv * i]);
    }
    // Combine all the Gaussian evaluations together.
    comm.Sum(allbias);
    comm.Sum(allder);
    // Transfer all the evaluations to the grid.
    for (unsigned i = 0; i < neighbors.size(); ++i) {
      unsigned ineigh = neighbors[i];
      for (unsigned j = 0; j < ncv; ++j) {
        der[j] = allder[ncv * i + j];
      }
      grid->addValueAndDerivatives(ineigh, allbias[i], der);
    }
  }
}

vector<unsigned> MetaD::getGaussianSupport(const Gaussian &hill) {
  vector<unsigned> nneigh;
  if (doInt_) {
    double cutoff = sqrt(2.0 * DP2CUTOFF) * hill.sigma[0];
    if (hill.center[0] + cutoff > uppI_ || hill.center[0] - cutoff < lowI_) {
      // in this case, we updated the entire grid to avoid problems
      return BiasGrid_->getNbin();
    } else {
      nneigh.push_back(static_cast<unsigned>(ceil(cutoff / BiasGrid_->getDx()[0])));
      return nneigh;
    }
  }
  // If using a diagonal hill, a good bounding box for the cutoff surface is
  // the number of grid binwidths required to reach the cutoff radius along 
  // each axis.
  if (!hill.multivariate) {
    for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
      double cutoff = sqrt(2.0 * DP2CUTOFF) * hill.sigma[i];
      nneigh.push_back(static_cast<unsigned>(ceil(cutoff / BiasGrid_->getDx()[i])));
    }
  // For general hills the hill axes may not line up with the coordinate axes, 
  // so a more involved calculation is necessary to find the bounding box for
  // the cutoff surface.
  } else {
    unsigned ncv = getNumberOfArguments();
    unsigned k = 0;
    //log<<"------- GET GAUSSIAN SUPPORT --------\n";
    // Unpack the hill's sigma vector into a symmetric matrix.
    Matrix<double> mymatrix(ncv, ncv);
    for (unsigned i = 0; i < ncv; i++) {
      for (unsigned j = i; j < ncv; j++) {
        mymatrix(i, j) = mymatrix(j, i) = hill.sigma[k];
        k++;
      }
    }
    // Reinvert so to have the ellipses
    Matrix<double> myinv(ncv, ncv);
    Invert(mymatrix, myinv);
    Matrix<double> myautovec(ncv, ncv);
    vector<double> myautoval(ncv); //should I take this or their square root?
    diagMat(myinv, myautoval, myautovec);
    double maxautoval;
    maxautoval = 0.;
    unsigned ind_maxautoval;
    ind_maxautoval = ncv;
    for (unsigned i = 0; i < ncv; i++) {
      if (myautoval[i] > maxautoval) {
        maxautoval = myautoval[i];
        ind_maxautoval = i;
      }
    }
    for (unsigned i = 0; i < ncv; i++) {
      double cutoff = sqrt(2.0 * DP2CUTOFF) * abs(sqrt(maxautoval) * myautovec(i, ind_maxautoval));
      //log<<"AUTOVAL "<<myautoval[0]<<" COMP "<<abs(myautoval[0]*myautovec(i,0)) <<" CUTOFF "<<cutoff<<"\n";
      nneigh.push_back(static_cast<unsigned>(ceil(cutoff / BiasGrid_->getDx()[i])));
    }
  }
  //log<<"------- END GET GAUSSIAN SUPPORT --------\n";
  return nneigh;
}

double MetaD::getBiasAndDerivatives(const vector<double> &cv, double* der) {
  double bias = 0.0;
  if (!grid_) {
    unsigned stride = comm.Get_size();
    unsigned rank = comm.Get_rank();
    for (unsigned i = rank; i < hills_.size(); i += stride) {
      bias += evaluateGaussian(cv, hills_[i], der);
      //finite difference test
      //finiteDifferenceGaussian(cv,hills_[i]);
    }
    comm.Sum(bias);
    if (der) {
      comm.Sum(der, getNumberOfArguments());
    }
  } else {
    if (der) {
      vector<double> vder(getNumberOfArguments());
      bias = BiasGrid_->getValueAndDerivatives(cv, vder);
      for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
        der[i] = vder[i];
      }
    } else {
      bias = BiasGrid_->getValue(cv);
    }
  }
  if (hasextgrid_) {
    if (der) {
      vector<double> vder(getNumberOfArguments());
      bias += ExtGrid_->getValueAndDerivatives(cv, vder);
      for (unsigned i = 0; i < getNumberOfArguments(); ++i) {
        der[i] += vder[i];
      }
    } else {
      bias += ExtGrid_->getValue(cv);
    }
  }

  // Add in the overall bias level if domains are being used.
  if (use_domains_) {
    bias += domain_implicit_bias_level_;
  }

  return bias;
}

double MetaD::evaluateGaussian(const vector<double> &cv, const Gaussian &hill, double* der) {
  double dp2 = 0.0;
  double bias = 0.0;
  // I use a pointer here because cv is const (and should be const)
  // but when using doInt it is easier to locally replace cv[0] with
  // the upper/lower limit in case it is out of range
  const double *pcv = NULL; // pointer to cv
  double tmpcv[1]; // tmp array with cv (to be used with doInt_)
  if (cv.size() > 0) {
    pcv = &cv[0];
  }
  if (doInt_) {
    plumed_dbg_assert(cv.size() == 1);
    pcv = &(tmpcv[0]);
    tmpcv[0] = cv[0];
    if (cv[0] < lowI_) {
      tmpcv[0] = lowI_;
    }
    if (cv[0] > uppI_) {
      tmpcv[0] = uppI_;
    }
  }
  if (hill.multivariate) {
    unsigned k = 0;
    unsigned ncv = cv.size();
    // recompose the full sigma from the upper diag cholesky
    Matrix<double> mymatrix(ncv, ncv);
    for (unsigned i = 0; i < ncv; i++) {
      for (unsigned j = i; j < ncv; j++) {
        mymatrix(i, j) = mymatrix(j, i) = hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    for (unsigned i = 0; i < cv.size(); ++i) {
      double dp_i = difference(i, hill.center[i], pcv[i]);
      dp_[i] = dp_i;
      for (unsigned j = i; j < cv.size(); ++j) {
        if (i == j) {
          dp2 += dp_i * dp_i * mymatrix(i, j) * 0.5;
        } else {
          double dp_j = difference(j, hill.center[j], pcv[j]);
          dp2 += dp_i * dp_j * mymatrix(i, j) ;
        }
      }
    }
    if (dp2 < DP2CUTOFF) {
      bias = hill.height * exp(-dp2);
      if (der) {
        for (unsigned i = 0; i < cv.size(); ++i) {
          double tmp = 0.0;
          k = i;
          for (unsigned j = 0; j < cv.size(); ++j) {
            tmp +=   dp_[j] * mymatrix(i, j) * bias;
          }
          der[i] -= tmp;
        }
      }
    }
  } else {
    for (unsigned i = 0; i < cv.size(); ++i) {
      double dp = difference(i, hill.center[i], pcv[i]) * hill.invsigma[i];
      dp2 += dp * dp;
      dp_[i] = dp;
    }
    dp2 *= 0.5;
    if (dp2 < DP2CUTOFF) {
      bias = hill.height * exp(-dp2);
      if (der) {
        for (unsigned i = 0; i < cv.size(); ++i) {
          der[i] += -bias * dp_[i] * hill.invsigma[i];
        }
      }
    }
  }
  if (doInt_) {
    if ((cv[0] < lowI_ || cv[0] > uppI_) && der) for (unsigned i = 0; i < cv.size(); ++i) {
      der[i] = 0;
    }
  }
  return bias;
}

double MetaD::getHeight(const vector<double> &cv) {
  double height = height0_;
  if (welltempered_) {
    double vbias = getBiasAndDerivatives(cv);
    height *= exp(-max(0.0, vbias - wt_biasthreshold_) / (kbt_ * (biasf_ - 1.0)));
  }
  if (transitiontempered_) {
    double vbarrier = getTransitionBarrierBias();
    if (tt_alpha_ == 1.0) {
      height *= exp(-max(0.0, vbarrier - tt_biasthreshold_) / (kbt_ * (tt_biasf_ - 1.0)));
    } else {
      height *= pow(1 + max(0.0, vbarrier - tt_biasthreshold_) / (kbt_ * (tt_biasf_ - 1.0)), - tt_alpha_ / (1 - tt_alpha_));
    }
  }
  if (benthic_toleration_) {
    if (BenthicHistogram_->getValue(BiasGrid_->getIndex(cv)) < benthic_tol_number_) {
      height = 0.0;
    }
  }
  if (use_domains_ && scale_new_hills_) {
    if (domain_ids_[BiasGrid_->getIndex(cv)] == 0) {
      height = 0.0;
    }
  }
  return height;
}

double MetaD::getTransitionBarrierBias() {
  
  // If there is only one well of interest, return the bias at that well point.
  if (transitionwells_.size() == 1) {
    double tb_bias = getBiasAndDerivatives(transitionwells_[0], NULL);
    return tb_bias;
  
  // Otherwise, check for the least barrier bias between all pairs of wells.
  // Note that because the paths can be considered edges between the wells' nodes
  // to make a graph and the path barriers satisfy certain cycle inequalities, it
  // is sufficient to look at paths corresponding to a minimal spanning tree of the
  // overall graph rather than examining every edge in the graph.
  // For simplicity, I chose the star graph with center well 0 as the spanning tree.
  // It is most efficient to start the path searches from the wells that are
  // expected to be sampled last, so transitionwell_[0] should correspond to the
  // starting well. With this choice the searches will terminate in one step until
  // transitionwell_[1] is sampled.
  } else {
    double least_transition_bias, curr_transition_bias;
    vector<double> sink = transitionwells_[0];
    vector<double> source = transitionwells_[1];
    least_transition_bias = BiasGrid_->findMaximalPathMinimum(source, sink);
    for (unsigned i = 2; i < transitionwells_.size(); i++) {
      if (least_transition_bias == 0.0) {
          break;
      }
      source = transitionwells_[i];
      curr_transition_bias = BiasGrid_->findMaximalPathMinimum(source, sink);
      least_transition_bias = fmin(curr_transition_bias, least_transition_bias);
    }
    if (use_domains_) {
      least_transition_bias += domain_implicit_bias_level_;
    }
    return least_transition_bias;
  }
}

void MetaD::spreadDomain(unsigned i, unsigned domain, const Grid * const input_region) {
  deque<unsigned> point_stack = deque<unsigned>();
  point_stack.push_back(i);
  unsigned curr_point;
  vector<double> temp_point(getNumberOfArguments());
  while (!point_stack.empty()) {
    curr_point = point_stack.front();
    vector<unsigned> neighs = BiasGrid_->getNearestNeighbors(curr_point);
    for (unsigned j = 0; j < neighs.size(); j++) {
      BiasGrid_->getPoint(neighs[j], temp_point);
      bool should_be_biased = (input_region->getValue(temp_point) == 1.0);
      bool domain_is_not_set = (domain_ids_[neighs[j]] == 0);
      if (should_be_biased && domain_is_not_set) {
        domain_ids_[neighs[j]] = domain;
        point_stack.push_back(neighs[j]);
      }
    }
    point_stack.pop_front();
  }
}

void MetaD::applyDomainBoundaryBias(unsigned idomain, double threshold) {
  // Set the overall bias level as needed.
  double boundary_ave = 0.0;
  // The bias level is set so that the current domain's 
  // average boundary bias matches the exterior bias.
  // First calculate the average boundary bias.
  for (unsigned i = 0; i < domain_boundary_pts_[idomain].size(); i++) {
     boundary_ave += BiasGrid_->getValue(domain_boundary_pts_[idomain][i]);
  }
  boundary_ave /= (double) domain_boundary_pts_[idomain].size();
  // If the boundary bias differs from zero appreciably, then
  // correct the bias level by 1) adding to the level and 2)
  // subtracting off a bias with the same boundary level that 
  // could be added through flat sampling of the domain. 
  if ((boundary_ave - domain_boundary_grid_levels_[idomain]) >= threshold) {
    // Add to the level.
    domain_implicit_bias_level_ += boundary_ave;
    // Subtract bias from the domain and its immediate surroundings.
    for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
      if (HillScalingGrids_[idomain]->getValue(i) > 0.0) {
        if (domain_ids_[i] == idomain) {
          BiasGrid_->addValue(i, -boundary_ave); 
        } else {
          double scaling = HillScalingGrids_[idomain]->getValue(i);
          double bias_corr = -boundary_ave * scaling / tame_scaling(scaling);
          BiasGrid_->addValue(i, bias_corr);
        }
      }
    }
  }
}

void MetaD::defineDomains(const Grid * const input_region) {
  
  // Create a new gridded integer function with the same exact
  // points as the bias grid. Clear it if it already exists.
  n_domains_ = 0;
  domain_ids_ = vector<unsigned>(BiasGrid_->getMaxSize(), 0);

  // Define the connected domains given the input region.
  // Search over all points to find ones that are supposed to be biased but
  // are not yet assigned to any domain.
  vector<double> temp_point(getNumberOfArguments());
  for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
    BiasGrid_->getPoint(i, temp_point);
    bool should_be_biased = (input_region->getValue(temp_point) == 1.0);
    bool domain_is_not_set = (domain_ids_[i] == 0);
    // If a point satisfies both conditions in this outer loop, it is the first point
    // discovered in a new connected domain. 
    if (should_be_biased && domain_is_not_set) {
      n_domains_++;
      // Set its domain value.
      domain_ids_[i] = n_domains_;
      // Find all other points to be biased that are connected to this one by other
      // points that should also be biased and set their domain values. Recursive.
      spreadDomain(i, n_domains_, input_region);
    }
  }

  // Find and store the boundary points for each domain.
  domain_boundary_pts_ = vector< vector<unsigned> >(n_domains_);
  for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
    unsigned idomain = domain_ids_[i];
    if (idomain != 0) {
      vector<unsigned> neighs = BiasGrid_->getNearestNeighbors(i);
      bool is_boundary_pt = false;
      for (unsigned j = 0; j < neighs.size(); j++) {
        if (idomain != domain_ids_[neighs[j]]) {
          is_boundary_pt = true;
          break;
        }
      }
      domain_boundary_pts_[idomain - 1].push_back(i);
    }
  }
  
  // Filter out domains that don't meet minimum size criteria if
  // they weren't specifically entered by hand.
  if (use_adaptive_domains_) {
    filterDomains();
  }

  // Calculate the level targets for each domain boundary
  // to keep the level unchanged after future hill additions.
  domain_boundary_grid_levels_ = vector<double>(n_domains_, 0.0);
  for (unsigned i = 0; i < n_domains_; i++) {
    // Find the average of the grid bias over the domain boundary.
    double boundary_ave = 0.0;
    for (unsigned j = 0; j < domain_boundary_pts_[i].size(); j++) {
       boundary_ave += BiasGrid_->getValue(domain_boundary_pts_[i][j]);
    }
    boundary_ave /= (double) domain_boundary_pts_[i].size();  
    domain_boundary_grid_levels_[i] = boundary_ave;
  }
}

void MetaD::filterDomains() {
  // Only do anything if there are actually domains to consider.
  if (n_domains_ == 0) {
    return;
  }
  // Set up. Allocate tags recording if each domain should be
  // filtered out or not and calculate the number of points
  // required to be in a domain for it to be retained after
  // filtering.
  vector<bool> eliminate_domain(n_domains_, false);
  int min_points = 1;
  
  // Start with the volume of a unit sphere.
  double vol = 1.0;
  double half_args = ((double) getNumberOfArguments()) / 2.0;
  vol *= pow(M_PI, half_args) / tgamma(1.0 + half_args);
  // Rescale by hill size (in units of grid points) along each dimension.
  for (unsigned i = 0; i < getNumberOfArguments(); i++) {
    vol *= sigma0_[i] / BiasGrid_->getDx()[i];
  }
  // Set the minimum number of points based on the volume.
  // If it is 1 or 0, there is no filtering to be done at all.
  if (floor(vol) > 1) {
    min_points = floor(vol);
  } else {
    return;
  }

  // For every domain currently defined, check if it is large
  // enough to survive the filter.
  for (unsigned i = 0; i < n_domains_; i++) {
    // Find the number of points in this domain
    int npoints = 0;
    for (unsigned j = 0; j < BiasGrid_->getMaxSize(); j++) {
      if (domain_ids_[j] == (i + 1)) {
        npoints++;
      }
    }
    // Record if it is above or below threshold.
    if (npoints < min_points) {
      eliminate_domain[i] = true;
    }
  }

  // Eliminate the domains that fell below the cutoff size.
  // Loop over all domains, assuming none will make the cut.
  int new_n_domains = 0;
  vector<vector<unsigned> > new_domain_boundary_pts;
  // Reassign all domain IDs for each domain.
  for (unsigned i = 0; i < n_domains_; i++) {
    if (eliminate_domain[i]) {
      // Eliminate any domain that should not exist by setting
      // the IDs to zero.
      for (unsigned j = 0; j < BiasGrid_->getMaxSize(); j++) {
        if (domain_ids_[j] == (i + 1)) {
          domain_ids_[j] = 0;
        }
      }
    } else if (!eliminate_domain[i]) {
      // Record that another domain made the cut.
      new_n_domains++;
      // Reset the domain's points' IDs
      for (unsigned j = 0; j < BiasGrid_->getMaxSize(); j++) {
        if (domain_ids_[j] == (i + 1)) {
          domain_ids_[j] = new_n_domains;
        }
      }
      // The new boundary points are the same as before, just
      // in a different place in the lists.
      new_domain_boundary_pts.push_back(domain_boundary_pts_[i]);
    }
  }
  domain_boundary_pts_ = new_domain_boundary_pts;
  n_domains_ = new_n_domains;
}

void MetaD::createScalingGrids() {
  plumed_dbg_assert(scale_new_hills_ == false);
  // Erase existing scaling grids if any are present.
  for (unsigned i = 0; i < HillScalingGrids_.size(); i++) {
    if (HillScalingGrids_[i] != NULL) {
      delete HillScalingGrids_[i];
    }
  }
  HillScalingGrids_ = vector<Grid *>();
  // Tell the program to scale future hills if a domain exists.
  if (n_domains_ > 0) {
    scale_new_hills_ = true;
  } else {
    scale_new_hills_ = false;
  }
  // For every domain, create a new scaling grid and evaluate a scaling function
  // for hills added in that domain. Only executes if a domain exists.
  for (unsigned i = 0; i < n_domains_; i++) {
    // Prepare to create a new scaling
    Grid * newScalingGrid;
    std::ostringstream label_stream;
    label_stream << getLabel() << ".scaling" << i;
    string funcl = label_stream.str();
    vector<unsigned> input_nbins = BiasGrid_->getNbin();
    for (unsigned j = 0; j < getNumberOfArguments(); ++j) {
      if (!getPntrToArgument(j)->isPeriodic()) {
        input_nbins[j] -= 1;
      }
    }
    // If grid size equals max grid size, create dense scaling grids.
    if (BiasGrid_->getMaxSize() == BiasGrid_->getSize()) {
      newScalingGrid = new Grid(funcl, getArguments(), BiasGrid_->getMin(), BiasGrid_->getMax(), input_nbins, false, true);
    // Otherwise create sparse scaling grids.
    } else {
      newScalingGrid = new SparseGrid(funcl, getArguments(), BiasGrid_->getMin(), BiasGrid_->getMax(), input_nbins, false, true);
    }
    // Create the base scaling function by adding a Gaussian to the grid for
    // each point in the domain.
    int pts_seen = 0;
    int hills_added = 0;
    vector<double> temp_point(getNumberOfArguments());
    for (unsigned j = 0; j < BiasGrid_->getMaxSize(); j++) {
      if (domain_ids_[j] == (i + 1)) {
        pts_seen++;
        if (pts_seen >= adaptive_domains_downsampling_ * hills_added) {
          hills_added++;
          BiasGrid_->getPoint(j, temp_point);
          Gaussian newhill = Gaussian(temp_point, sigma0_, height0_, false);
          addGaussianToGrid(newhill, newScalingGrid);
        }
      }
    }

    // Rescale the base scaling function using the min of the boundary values.
    // This sets the scaling function to have a minimum of one on the boundary.
    double boundary_min = newScalingGrid->getValue(domain_boundary_pts_[i][0]);
    for (unsigned j = 1; j < domain_boundary_pts_[i].size(); j++) {
      boundary_min = fmin(boundary_min, newScalingGrid->getValue(domain_boundary_pts_[i][j]));
    }
    for (unsigned j = 0; j < newScalingGrid->getMaxSize(); j++) {
      newScalingGrid->setValue(j, newScalingGrid->getValue(j) / boundary_min);
    }
    for (unsigned j = 0; j < newScalingGrid->getMaxSize(); j++) {
      newScalingGrid->setDerivFromValues(j);
    }
    // Save the completed scaling grid.
    HillScalingGrids_.push_back(newScalingGrid);
  }
}

bool MetaD::shouldAdaptDomainsNow() {
  
  // Check that this step is a candidate.
  if (getStep() % adaptive_domains_stride_ != 0 || isFirstStep) {
    return false;
  }
  if (adaptive_domains_eincrement_ > 0.0) {
    if (scale_new_hills_) { 
      if (domain_implicit_bias_level_ < adaptive_domains_eincrement_ + adaptive_domains_last_elevel_) {
        return false;
      }
    }
  }

  // Check the delay condition if necessary.
  bool adapt_domains_now = true;
  if (delay_adaptive_domains_) {
    // For the min-referenced adaptive domains, delay until the 
    // bias (scaled into free energy by the biasfactor) is above
    // the offset somewhere.
    double trigger_bias = 0.0;
    if (adaptive_domains_reftype_ == kMinRef) {
      double max_bias = 0.0;
      for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
        max_bias = max(max_bias, BiasGrid_->getValue(i));
      }
      trigger_bias = max_bias + domain_implicit_bias_level_;
    // For the transition-referenced adaptive domains, delay until
    // the transition barrier bias is above the offset.
    } else if (adaptive_domains_reftype_ == kTransitionRef) {
      trigger_bias = getTransitionBarrierBias();
    }
    // Rescale the bias to match a free energy estimate as needed.
    if (welltempered_) {
      trigger_bias *= biasf_ / (biasf_ - 1.0);
    }
    // If the trigger is not passed, continue the delay.
    if (trigger_bias <= adaptive_domains_eoffset_) {
      adapt_domains_now = false;
    // Otherwise allow the adaptation and end the delay.
    } else {
      delay_adaptive_domains_ = false;
    }
  }
  return adapt_domains_now;
}

void MetaD::adaptDomains() {

  // Ensure that the boundary biases are perfectly consistent
  // before redefinition to avoid drift.
  if (scale_new_hills_) {
    for (unsigned i = 0; i < n_domains_; i++) {
      applyDomainBoundaryBias(i, 0.0);
    }
  }

  // Calculate a new free energy estimate.
  Grid *new_region;
  std::string funcl = getLabel() + ".indicator";
  vector<unsigned> input_nbins = BiasGrid_->getNbin();
  for (unsigned j = 0; j < getNumberOfArguments(); ++j) {
    if (!getPntrToArgument(j)->isPeriodic()) {
      input_nbins[j] -= 1;
    }
  }
  new_region = new Grid(funcl, getArguments(), BiasGrid_->getMin(), BiasGrid_->getMax(), input_nbins, false, false);
  // Copy the histogram.
  for (unsigned i = 0; i < new_region->getMaxSize(); i++) {
    new_region->setValue(i, DomainsHistogram_->getValue(i));
  }
  // Take a log scaled by kbT to convert the histogram into an energy.
  new_region->logAllValuesAndDerivatives(-kbt_);
  if (print_adaptive_domains_energies_) {
    HistEnergyFile_.rewind();
    new_region->writeToFile(HistEnergyFile_);
    HistEnergyFile_.flush();
  }

  // Subtract the current bias.
  for (unsigned i = 0; i < new_region->getMaxSize(); i++) {
    new_region->addValue(i, -BiasGrid_->getValue(i));
  }
  if (print_adaptive_domains_energies_) {
    HistBiasEnergyFile_.rewind();
    new_region->writeToFile(HistBiasEnergyFile_);
    HistBiasEnergyFile_.flush();
  }

  // Find the threshold to use. 
  // First find the adaptive reference value.
  double reference_energy = 0.0;
  if (adaptive_domains_reftype_ == kMinRef) {
    reference_energy = new_region->getValue(0);
    for (unsigned i = 1; i < new_region->getMaxSize(); i++) {
      reference_energy = min(reference_energy, new_region->getValue(i));
    }
  } else if (adaptive_domains_reftype_ == kTransitionRef) {
    // Trick the transition barrier routine for the bias into providing the
    // transition barrier on the new free energy surface.
    // Store bias data.
    Grid *temp_grid = BiasGrid_;
    double temp_level = domain_implicit_bias_level_;
    // Replace the bias data with the new free energy estimate data.
    BiasGrid_ = new_region;
    domain_implicit_bias_level_ = 0.0;
    // Flip the energy estimate.
    for (unsigned i = 0; i < new_region->getMaxSize(); i++) {
      new_region->setValue(i, -new_region->getValue(i));
    }
    // Perform the search and take the negative of the output value.
    reference_energy = -getTransitionBarrierBias();
    // Flip the energy estimate back.
    for (unsigned i = 0; i < new_region->getMaxSize(); i++) {
      new_region->setValue(i, -new_region->getValue(i));
    }
    // Reset the bias data to its true values.
    domain_implicit_bias_level_ = temp_level;
    BiasGrid_ = temp_grid;
  }
  // The threshold is now a simple offset from that value.
  double threshold = reference_energy + adaptive_domains_eoffset_;
  
  // The new region is the region where the free energy is below
  // that threshold.
  for (unsigned i = 0; i < new_region->getMaxSize(); i++) {
    if (new_region->getValue(i) < threshold) {
      new_region->setValue(i, 1.0);
    } else {
      new_region->setValue(i, 0.0);
    }
  }

  // Update the domains and the scaling grid.
  scale_new_hills_ = false;
  defineDomains(new_region);
  delete new_region;
  createScalingGrids();
  if (n_domains_ > 0) {
    // Print for examination if desired.
    if (print_domains_scaling_) {
      if (DomainsScalingFilePs_.size() < HillScalingGrids_.size()) {
        for (unsigned i = DomainsScalingFilePs_.size(); i < HillScalingGrids_.size(); i++) {
          std::ostringstream filename_stream;
          filename_stream << "domain_" << i << "_scaling.dat";
          std::string scaling_filename = filename_stream.str();
          OFile* scaling_filep = new OFile();
          scaling_filep->link(*this);
          scaling_filep->open(scaling_filename);
          DomainsScalingFilePs_.push_back(scaling_filep);
        }
      }
      for (unsigned i = 0; i < DomainsScalingFilePs_.size(); i++) DomainsScalingFilePs_[i]->rewind();
      for (unsigned i = 0; i < HillScalingGrids_.size(); i++) HillScalingGrids_[i]->writeToFile(*(DomainsScalingFilePs_[i]));
      for (unsigned i = 0; i < DomainsScalingFilePs_.size(); i++) DomainsScalingFilePs_[i]->flush();
    }
  }
  // Record the bias level for the purpose of 
  // timing the next update.
  adaptive_domains_last_elevel_ = domain_implicit_bias_level_;
}

/// Calculate the bias and bias force at the current step, update the
/// acceleration factor, and communicate the bias force to the simulation.

void MetaD::calculate() {
// this is because presently there is no way to properly pass information
// on adaptive hills (diff) after exchanges:
  if (adaptive_ == FlexibleBin::diffusion && getExchangeStep()) {
    error("ADAPTIVE=DIFF is not compatible with replica exchange");
  }
  unsigned ncv = getNumberOfArguments();
  vector<double> cv(ncv);
  for (unsigned i = 0; i < ncv; ++i) {
    cv[i] = getArgument(i);
  }
  double* der = new double[ncv];
  for (unsigned i = 0; i < ncv; ++i) {
    der[i] = 0.0;
  }
  double ene = getBiasAndDerivatives(cv, der);
  getPntrToComponent("bias")->set(ene);
  // calculate the acceleration factor
  if (acceleration && !isFirstStep) {
    if (!use_domains_) {
      acc += exp(ene / (kbt_));
    } else {
      acc += exp((ene - domain_implicit_bias_level_) / (kbt_));
    }
    double mean_acc = acc / ((double) getStep());
    getPntrToComponent("acc")->set(mean_acc);
  } else if (acceleration && isFirstStep) {
    getPntrToComponent("acc")->set(1.0);
  }
  // Set the average bias
  if (calc_average_bias_coft_) {
    getPntrToComponent("coft")->set(average_bias_coft_);  
  } 
  // set Forces
  for (unsigned i = 0; i < ncv; ++i) {
    const double f = -der[i];
    setOutputForce(i, f);
  }
  delete [] der;
}



void MetaD::update() {
  vector<double> cv(getNumberOfArguments());
  vector<double> thissigma;
  bool multivariate;
  // adding hills criteria (could be more complex though)
  bool nowAddAHill;
  if (getStep() % stride_ == 0 && !isFirstStep) {
    nowAddAHill = true;
  } else {
    nowAddAHill = false;
    isFirstStep = false;
  }
  for (unsigned i = 0; i < cv.size(); ++i) {
    cv[i] = getArgument(i);
  }
  // if you use adaptive, call the FlexibleBin
  if (adaptive_ != FlexibleBin::none) {
    flexbin->update(nowAddAHill);
    multivariate = true;
  } else {
    multivariate = false;
  };
  // When using adaptive domains metabasin metadynamics, record a histogram
  // in addition to the usual hills. These are not normed because this will
  // only be used through the log and normalization thus corresponds to an
  // irrelevant constant. 
  if (use_adaptive_domains_) {
    if (getStep() % domains_histo_stride_ == 0 && !isFirstStep) {
      KernelFunctions kernel(cv, domains_histo_bandwidth_, "gaussian", false, 1.0, false);
      DomainsHistogram_->addKernel(kernel);
    }
    // Check if the delay conditions are satisfied.
    if (shouldAdaptDomainsNow()) {
      // If so calculate a new region and new domains.
      adaptDomains();
    }
  }
  if (benthic_toleration_) {
    if (getStep() % benthic_histo_stride_ == 0 && !isFirstStep) {
      KernelFunctions kernel(cv, benthic_histo_bandwidth_, "gaussian", false, 1.0, false);
      BenthicHistogram_->addKernel(kernel);
    }
    // When using benthic erosion, erosion can't affect anything until a
    // new hill is added--so only perform erosion just before adding hills.
    // Hills are added directly or from reading other walkers' hills.
    if (benthic_erosion_ && (nowAddAHill || (mw_n_ > 1 && getStep() % mw_rstride_ == 0))) {
      // The erosion time is a real time scale, so it is compared to the
      // boosted time rather than the plain simulation time.
      if ( acc * getTimeStep() > (last_benthic_erosion_ + benthic_erosive_time_)) {
        for (unsigned i = 0; i < BenthicHistogram_->getMaxSize(); i++) {
          if (BenthicHistogram_->getValue(i) <= benthic_tol_number_) {
            BenthicHistogram_->setValue(i, 0.0);
          }
        }
        last_benthic_erosion_ = acc * getTimeStep();
      }
    }
  }
  if (nowAddAHill) { // probably this can be substituted with a signal
    // add a Gaussian
    double height = getHeight(cv);
    // use normal sigma or matrix form?
    if (adaptive_ != FlexibleBin::none) {
      thissigma = flexbin->getInverseMatrix(); // returns upper diagonal inverse
      //cerr<<"ADDING HILLS "<<endl;
    } else {
      thissigma = sigma0_;  // returns normal sigma
    }
    // In case we use walkers_mpi, it is now necessary to communicate with other replicas.
    if (walkers_mpi) {
      int nw = 0;
      int mw = 0;
      if (comm.Get_rank() == 0) {
        // Only root of group can communicate with other walkers
        nw = multi_sim_comm.Get_size();
        mw = multi_sim_comm.Get_rank();
      }
      // Communicate to the other members of the same group
      // info abount number of walkers and walker index
      comm.Bcast(nw, 0);
      comm.Bcast(mw, 0);
      // Allocate arrays to store all walkers hills
      std::vector<double> all_cv(nw * cv.size(), 0.0);
      std::vector<double> all_sigma(nw * thissigma.size(), 0.0);
      std::vector<double> all_height(nw, 0.0);
      std::vector<int>    all_multivariate(nw, 0);
      if (comm.Get_rank() == 0) {
        // Communicate (only root)
        multi_sim_comm.Allgather(cv, all_cv);
        multi_sim_comm.Allgather(thissigma, all_sigma);
        multi_sim_comm.Allgather(height, all_height);
        multi_sim_comm.Allgather(int(multivariate), all_multivariate);
      }
      // Share info with group members
      comm.Bcast(all_cv, 0);
      comm.Bcast(all_sigma, 0);
      comm.Bcast(all_height, 0);
      comm.Bcast(all_multivariate, 0);
      for (int i = 0; i < nw; i++) {
        // actually add hills one by one
        std::vector<double> cv_now(cv.size());
        std::vector<double> sigma_now(thissigma.size());
        for (unsigned j = 0; j < cv.size(); j++) {
          cv_now[j] = all_cv[i * cv.size() + j];
        }
        for (unsigned j = 0; j < thissigma.size(); j++) {
          sigma_now[j] = all_sigma[i * thissigma.size() + j];
        }
        Gaussian newhill = Gaussian(cv_now, sigma_now, all_height[i], all_multivariate[i]);
        addGaussian(newhill);
        writeGaussian(newhill, hillsOfile_);
      }
    } else {
      Gaussian newhill = Gaussian(cv, thissigma, height, multivariate);
      addGaussian(newhill);
      // print on HILLS file
      writeGaussian(newhill, hillsOfile_);
    }
  }
  // dump grid on file
  if (wgridstride_ > 0 && getStep() % wgridstride_ == 0) {
    dumpBias();
    if (use_adaptive_domains_ && whistofilename_.size() > 0) {
      dumpGrid(DomainsHistogram_, whistofile_);
    }
  }
  // if multiple walkers and time to read Gaussians
  if (mw_n_ > 1 && getStep() % mw_rstride_ == 0) {
    for (int i = 0; i < mw_n_; ++i) {
      // don't read your own Gaussians
      if (i == mw_id_) {
        continue;
      }
      // if the file is not open yet
      if (!(ifiles[i]->isOpen())) {
        // check if it exists now and open it!
        if (ifiles[i]->FileExist(ifilesnames[i])) {
          ifiles[i]->open(ifilesnames[i]);
          ifiles[i]->reset(false);
        }
        // otherwise read the new Gaussians
      } else {
        log.printf("  Reading hills from %s:", ifilesnames[i].c_str());
        readGaussians(ifiles[i]);
        ifiles[i]->reset(false);
      }
    }
  }
  // Calculate the new average bias after any bias update.
  // This follows the Tiwary and Parrinello JPCB paper.
  if (calc_average_bias_coft_ && (nowAddAHill || (mw_n_ > 1 && getStep() % mw_rstride_ == 0))) {
    // Calc sums rather than integrals because the normalization
    // is irrelevant.
    double exp_free_energy_sum = 0.0;
    double exp_biased_free_energy_sum = 0.0;
    // The formula depends on how the final free energy 
    // should be inferred from the bias.
    // For reasons I don't understand, the first branch of the if
    // statement can't be followed unless I flush the log first.
    log.flush();
    if (biasf_ == 1.0) {
      for (unsigned i; i < BiasGrid_->getMaxSize(); i++) {
        double pt_bias = BiasGrid_->getValue(i);
        exp_free_energy_sum += exp(pt_bias / kbt_);
        exp_biased_free_energy_sum += 1.0;
      }
    } else if (biasf_ > 1.0) {
      for (unsigned i; i < BiasGrid_->getMaxSize(); i++) {
        double pt_bias = BiasGrid_->getValue(i);
        exp_free_energy_sum += exp(biasf_ * pt_bias / (kbt_  * (biasf_ - 1)));
        exp_biased_free_energy_sum += exp(pt_bias / (kbt_ * (biasf_ - 1)));
      }
    }
    average_bias_coft_ = kbt_ * ( std::log(exp_free_energy_sum) - std::log(exp_biased_free_energy_sum));
    getPntrToComponent("coft")->set(average_bias_coft_);
  }
}

void MetaD::finiteDifferenceGaussian(const vector<double> &cv, const Gaussian &hill) {
  log << "--------- finiteDifferenceGaussian: size " << cv.size() << "------------\n";
// for each cv
// first get the bias and the derivative
  vector<double> oldder(cv.size());
  vector<double> der(cv.size());
  vector<double> mycv(cv.size());
  mycv = cv;
  double step = 1.e-6;
  Random random;
// just displace a tiny bit
  for (unsigned i = 0; i < cv.size(); i++) {
    log << "CV " << i << " V " << mycv[i] << "\n";
  }
  for (unsigned i = 0; i < cv.size(); i++) {
    mycv[i] += 1.e-2 * 2 * (random.RandU01() - 0.5);
  }
  for (unsigned i = 0; i < cv.size(); i++) {
    log << "NENEWWCV " << i << " V " << mycv[i] << "\n";
  }
  double oldbias = evaluateGaussian(mycv, hill, &oldder[0]);
  for (unsigned i = 0; i < mycv.size(); i++) {
    double delta = step * 2 * (random.RandU01() - 0.5);
    mycv[i] += delta;
    double newbias = evaluateGaussian(mycv, hill, &der[0]);
    log << "CV " << i;
    log << " ANAL " << oldder[i] << " NUM " << (newbias - oldbias) / delta << " DIFF " << (oldder[i] - (newbias - oldbias) / delta) << "\n";
    mycv[i] -= delta;
  }
  log << "--------- END finiteDifferenceGaussian ------------\n";
}

/// Read a single hill specification from file. In order to parse the line this function requires the number
/// of CVs, passed through the size of the tmpvalues vector. If there is no hill on this line, the function returns false.
/// The center coordinate of the hill is written into the vector center, the vector of widths for a diagonal hill or the
/// flattened matrix for a nondiagonal hill is written into sigma, the height is written to height, and whether it is 
/// non-diagonal (true) or diagonal (false) is written to multivariate.

bool MetaD::scanOneHill(IFile *ifile,  vector<Value> &tmpvalues, vector<double> &center, vector<double>  &sigma, double &height , bool &multivariate) {
  double dummy;
  multivariate = false;
  if (ifile->scanField("time", dummy)) {
    unsigned ncv;
    ncv = tmpvalues.size();
    for (unsigned i = 0; i < ncv; ++i) {
      ifile->scanField(&tmpvalues[i]);
      if (tmpvalues[i].isPeriodic() && ! getPntrToArgument(i)->isPeriodic()) {
        error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
      } else if (tmpvalues[i].isPeriodic()) {
        std::string imin, imax;
        tmpvalues[i].getDomain(imin, imax);
        std::string rmin, rmax;
        getPntrToArgument(i)->getDomain(rmin, rmax);
        if (imin != rmin || imax != rmax) {
          error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
        }
      }
      center[i] = tmpvalues[i].get();
    }
    // scan for multivariate label: record the actual file position so to eventually rewind
    std::string sss;
    ifile->scanField("multivariate", sss);
    if (sss == "true") {
      multivariate = true;
    } else if (sss == "false") {
      multivariate = false;
    } else {
      plumed_merror("cannot parse multivariate = " + sss);
    }
    if (multivariate) {
      sigma.resize(ncv * (ncv + 1) / 2);
      Matrix<double> upper(ncv, ncv);
      Matrix<double> lower(ncv, ncv);
      for (unsigned i = 0; i < ncv; i++) {
        for (unsigned j = 0; j < ncv - i; j++) {
          ifile->scanField("sigma_" + getPntrToArgument(j + i)->getName() + "_" + getPntrToArgument(j)->getName(), lower(j + i, j));
          upper(j, j + i) = lower(j + i, j);
        }
      }
      Matrix<double> mymult(ncv, ncv);
      Matrix<double> invmatrix(ncv, ncv);
      //log<<"Lower \n";
      //matrixOut(log,lower);
      //log<<"Upper \n";
      //matrixOut(log,upper);
      mult(lower, upper, mymult);
      //log<<"Mult \n";
      //matrixOut(log,mymult);
      // now invert and get the sigmas
      Invert(mymult, invmatrix);
      //log<<"Invert \n";
      //matrixOut(log,invmatrix);
      // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form)
      unsigned k = 0;
      for (unsigned i = 0; i < ncv; i++) {
        for (unsigned j = i; j < ncv; j++) {
          sigma[k] = invmatrix(i, j);
          k++;
        }
      }
    } else {
      for (unsigned i = 0; i < ncv; ++i) {
        ifile->scanField("sigma_" + getPntrToArgument(i)->getName(), sigma[i]);
      }
    }
    ifile->scanField("height", height);
    ifile->scanField("biasf", dummy);
    if (ifile->FieldExist("clock")) {
      ifile->scanField("clock", dummy);
    }
    if (ifile->FieldExist("lower_int")) {
      ifile->scanField("lower_int", dummy);
    }
    if (ifile->FieldExist("upper_int")) {
      ifile->scanField("upper_int", dummy);
    }
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}

void MetaD::dumpBias() {
  // Add the internally used domain bias level.
  if (use_domains_ && domain_implicit_bias_level_ > 0.0) {
    for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
      BiasGrid_->addValue(i, domain_implicit_bias_level_);
    }
  }
  dumpGrid(BiasGrid_, gridfile_);
  // Subtract it again.
  if (use_domains_ && domain_implicit_bias_level_ > 0.0) {
    for (unsigned i = 0; i < BiasGrid_->getMaxSize(); i++) {
      BiasGrid_->addValue(i, -domain_implicit_bias_level_);
    }
  }
}

void MetaD::dumpGrid(Grid *grid, OFile &gridfile) {
  // in case old grids are stored, a sequence of grids should appear
  // this call results in a repetition of the header:
  if (storeOldGrids_) {
    gridfile.clearFields();
  // in case only latest grid is stored, file should be rewound
  // this will overwrite previously written grids
  } else {
    gridfile.rewind();
  }
  grid->writeToFile(gridfile);
  // if a single grid is stored, it is necessary to flush it, otherwise
  // the file might stay empty forever (when a single grid is not large enough to
  // trigger flushing from the operating system).
  // on the other hand, if grids are stored one after the other this is
  // no necessary, and we leave the flushing control to the user as usual
  // (with FLUSH keyword)
  if (!storeOldGrids_) {
    gridfile.flush();
  }
}

}
}