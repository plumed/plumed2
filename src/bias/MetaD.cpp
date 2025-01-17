/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

//+PLUMEDOC BIAS METAD
/*
Used to performed metadynamics on one or more collective variables.

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
V(\vec{s}) = -F(\vec{s})
\f]

During post processing the free energy can be calculated in this way using the \ref sum_hills
utility.

In the simplest possible implementation of a metadynamics calculation the expense of a metadynamics
calculation increases with the length of the simulation as one has to, at every step, evaluate
the values of a larger and larger number of Gaussian kernels. To avoid this issue you can
store the bias on a grid.  This approach is similar to that proposed in \cite babi08jcp but has the
advantage that the grid spacing is independent on the Gaussian width.
Notice that you should provide the grid boundaries (GRID_MIN and GRID_MAX) and either the number of bins
for every collective variable (GRID_BIN) or the desired grid spacing (GRID_SPACING).
In case you provide both PLUMED will use the most conservative choice (highest number of bins) for each dimension.
In case you do not provide any information about bin size (neither GRID_BIN nor GRID_SPACING)
PLUMED will use 1/5 of the Gaussian width (SIGMA) as grid spacing if the width is fixed or 1/5 of the minimum
Gaussian width (SIGMA_MIN) if the width is variable. This default choice should be reasonable for most applications.

Alternatively to the use of grids, it is possible to use a neighbor list to decrease the cost of evaluating the bias,
this can be enabled using NLIST. NLIST can be beneficial with more than 2 collective variables, where GRID becomes
expensive and memory consuming. The neighbor list will be updated everytime the CVs go farther than a cut-off value
from the position they were at last neighbor list update. Gaussians are added to the neigbhor list if their center
is within 6.*DP2CUTOFF*sigma*sigma. While the list is updated if the CVs are farther from the center than 0.5 of the
standard deviation of the Gaussian center distribution of the list. These parameters (6 and 0.5) can be modified using
NLIST_PARAMETERS. Note that the use of neighbor list does not provide the exact bias.

Metadynamics can be restarted either from a HILLS file as well as from a GRID, in this second
case one can first save a GRID using GRID_WFILE (and GRID_WSTRIDE) and at a later stage read
it using GRID_RFILE.

The work performed by the METAD bias can be calculated using CALC_WORK, note that this is expensive when not using grids.

Another option that is available in plumed is well-tempered metadynamics \cite Barducci:2008. In this
variant of metadynamics the heights of the Gaussian hills are scaled at each step so the bias is now
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

As a final note, since version 2.0.2 when the system is outside of the selected interval the force
is set to zero and the bias value to the value at the corresponding boundary. This allows acceptances
for replica exchange methods to be computed correctly.

Multiple walkers  \cite multiplewalkers can also be used. See below the examples.


The \f$c(t)\f$ reweighting factor can also be calculated on the fly using the equations
presented in \cite Tiwary_jp504920s.
The expression used to calculate \f$c(t)\f$ follows directly from Eq. 3 in \cite Tiwary_jp504920s,
where \f$F(\vec{s})=-\gamma/(\gamma-1) V(\vec{s})\f$.
This gives smoother results than equivalent Eqs. 13 and Eqs. 14 in that paper.
The \f$c(t)\f$ is given by the rct component while the bias
normalized by \f$c(t)\f$ is given by the rbias component (rbias=bias-rct) which can be used
to obtain a reweighted histogram.
The calculation of \f$c(t)\f$ is enabled by using the keyword CALC_RCT.
By default \f$c(t)\f$ is updated every time the bias changes, but if this slows down the simulation
the keyword RCT_USTRIDE can be set to a value higher than 1.
This option requires that a grid is used.

Additional material and examples can be also found in the tutorials:

- \ref lugano-3

Concurrent metadynamics
as done e.g. in Ref. \cite gil2015enhanced . This indeed can be obtained by using the METAD
action multiple times in the same input file.

\par Examples

The following input is for a standard metadynamics calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the metadynamics bias potential are written to the COLVAR file every 100 steps.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 PACE=500 LABEL=restraint
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endplumedfile
(See also \ref DISTANCE \ref PRINT).

\par
If you use adaptive Gaussian kernels, with diffusion scheme where you use
a Gaussian that should cover the space of 20 time steps in collective variables.
Note that in this case the histogram correction is needed when summing up hills.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=20 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=DIFF
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endplumedfile

\par
If you use adaptive Gaussian kernels, with geometrical scheme where you use
a Gaussian that should cover the space of 0.05 nm in Cartesian space.
Note that in this case the histogram correction is needed when summing up hills.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=GEOM
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endplumedfile

\par
When using adaptive Gaussian kernels you might want to limit how the hills width can change.
You can use SIGMA_MIN and SIGMA_MAX keywords.
The sigmas should specified in terms of CV so you should use the CV units.
Note that if you use a negative number, this means that the limit is not set.
Note also that in this case the histogram correction is needed when summing up hills.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ...
  ARG=d1,d2 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint ADAPTIVE=GEOM
  SIGMA_MIN=0.2,0.1 SIGMA_MAX=0.5,1.0
... METAD
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endplumedfile

\par
Multiple walkers can be also use as in  \cite multiplewalkers
These are enabled by setting the number of walker used, the id of the
current walker which interprets the input file, the directory where the
hills containing files resides, and the frequency to read the other walkers.
Here is an example
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
METAD ...
   ARG=d1 SIGMA=0.05 HEIGHT=0.3 PACE=500 LABEL=restraint
   WALKERS_N=10
   WALKERS_ID=3
   WALKERS_DIR=../
   WALKERS_RSTRIDE=100
... METAD
\endplumedfile
where  WALKERS_N is the total number of walkers, WALKERS_ID is the
id of the present walker (starting from 0 ) and the WALKERS_DIR is the directory
where all the walkers are located. WALKERS_RSTRIDE is the number of step between
one update and the other. Since version 2.2.5, hills files are automatically
flushed every WALKERS_RSTRIDE steps.

\par
The \f$c(t)\f$ reweighting factor can be calculated on the fly using the equations
presented in \cite Tiwary_jp504920s as described above.
This is enabled by using the keyword CALC_RCT,
and can be done only if the bias is defined on a grid.
\plumedfile
phi: TORSION ATOMS=1,2,3,4
psi: TORSION ATOMS=5,6,7,8

METAD ...
 LABEL=metad
 ARG=phi,psi SIGMA=0.20,0.20 HEIGHT=1.20 BIASFACTOR=5 TEMP=300.0 PACE=500
 GRID_MIN=-pi,-pi GRID_MAX=pi,pi GRID_BIN=150,150
 CALC_RCT
 RCT_USTRIDE=10
... METAD
\endplumedfile
Here we have asked that the calculation is performed every 10 hills deposition by using
RCT_USTRIDE keyword. If this keyword is not given, the calculation will
by default be performed every time the bias changes. The \f$c(t)\f$ reweighting factor will be given
in the rct component while the instantaneous value of the bias potential
normalized using the \f$c(t)\f$ reweighting factor is given in the rbias component
[rbias=bias-rct] which can be used to obtain a reweighted histogram or
free energy surface using the \ref HISTOGRAM analysis.

\par
The kinetics of the transitions between basins can also be analyzed on the fly as
in \cite PRL230602. The flag ACCELERATION turn on accumulation of the acceleration
factor that can then be used to determine the rate. This method can be used together
with \ref COMMITTOR analysis to stop the simulation when the system get to the target basin.
It must be used together with Well-Tempered Metadynamics. If restarting from a previous
metadynamics you need to use the ACCELERATION_RFILE keyword to give the name of the
data file from which the previous value of the acceleration factor should be read, otherwise the
calculation of the acceleration factor will be wrong.

\par
By using the flag FREQUENCY_ADAPTIVE the frequency adaptive scheme introduced in \cite Wang-JCP-2018
is turned on. The frequency for hill addition then changes dynamically based on the acceleration factor
according to the following equation
\f[
\tau_{\mathrm{dep}}(t) =
\min\left[
\tau_0 \cdot
\max\left[\frac{\alpha(t)}{\theta},1\right]
,\tau_{c}
\right]
\f]
where \f$\tau_0\f$ is the initial hill addition frequency given by the PACE keyword,
\f$\tau_{c}\f$ is the maximum allowed frequency given by the FA_MAX_PACE keyword,
\f$\alpha(t)\f$ is the instantaneous acceleration factor at time \f$t\f$,
and \f$\theta\f$ is a threshold value that acceleration factor has to reach before
triggering a change in the hill addition frequency given by the FA_MIN_ACCELERATION keyword.
The frequency for updating the hill addition frequency according to this equation is
given by the FA_UPDATE_FREQUENCY keyword, by default it is the same as the value given
in PACE. The hill hill addition frequency increase monotonously such that if the
instantaneous acceleration factor is lower than in the previous updating step the
previous \f$\tau_{\mathrm{dep}}\f$ is kept rather than updating it to a lower value.
The instantaneous hill addition frequency \f$\tau_{\mathrm{dep}}(t)\f$ is outputted
to pace component. Note that if restarting from a previous metadynamics run you need to
use the ACCELERATION_RFILE keyword to read in the acceleration factors from the
previous run, otherwise the hill addition frequency will start from the initial
frequency.


\par
You can also provide a target distribution using the keyword TARGET
\cite white2015designing
\cite marinelli2015ensemble
\cite gil2016empirical
The TARGET should be a grid containing a free-energy (i.e. the -\f$k_B\f$T*log of the desired target distribution).
Gaussian kernels will then be scaled by a factor
\f[
e^{\beta(\tilde{F}(s)-\tilde{F}_{max})}
\f]
Here \f$\tilde{F}(s)\f$ is the free energy defined on the grid and \f$\tilde{F}_{max}\f$ its maximum value.
Notice that we here used the maximum value as in ref \cite gil2016empirical
This choice allows to avoid exceedingly large Gaussian kernels to be added. However,
it could make the Gaussian too small. You should always choose carefully the HEIGHT parameter
in this case.
The grid file should be similar to other PLUMED grid files in that it should contain
both the target free-energy and its derivatives.

Notice that if you wish your simulation to converge to the target free energy you should use
the DAMPFACTOR command to provide a global tempering \cite dama2014well
Alternatively, if you use a BIASFACTOR your simulation will converge to a free
energy that is a linear combination of the target free energy and of the intrinsic free energy
determined by the original force field.

\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
METAD ...
 LABEL=t1
 ARG=d1 SIGMA=0.05 TAU=200 DAMPFACTOR=100 PACE=250
 GRID_MIN=1.14 GRID_MAX=1.32 GRID_BIN=6
 TARGET=dist.grid
... METAD

PRINT ARG=d1,t1.bias STRIDE=100 FILE=COLVAR
\endplumedfile

The file dist.dat for this calculation would read:

\auxfile{dist.grid}
#! FIELDS d1 t1.target der_d1
#! SET min_d1 1.14
#! SET max_d1 1.32
#! SET nbins_d1  6
#! SET periodic_d1 false
   1.1400   0.0031   0.1101
   1.1700   0.0086   0.2842
   1.2000   0.0222   0.6648
   1.2300   0.0521   1.4068
   1.2600   0.1120   2.6873
   1.2900   0.2199   4.6183
   1.3200   0.3948   7.1055
\endauxfile

Notice that BIASFACTOR can also be chosen as equal to 1. In this case one will perform
unbiased sampling. Instead of using HEIGHT, one should provide the TAU parameter.
\plumedfile
d: DISTANCE ATOMS=3,5
METAD ARG=d SIGMA=0.1 TAU=4.0 TEMP=300 PACE=100 BIASFACTOR=1.0
\endplumedfile
The HILLS file obtained will still work with `plumed sum_hills` so as to plot a free-energy.
The case where this makes sense is probably that of RECT simulations.

Regarding RECT simulations, you can also use the RECT keyword so as to avoid using multiple input files.
For instance, a single input file will be
\plumedfile
d: DISTANCE ATOMS=3,5
METAD ARG=d SIGMA=0.1 TAU=4.0 TEMP=300 PACE=100 RECT=1.0,1.5,2.0,3.0
\endplumedfile
The number of elements in the RECT array should be equal to the number of replicas.

*/
//+ENDPLUMEDOC

class MetaD : public Bias {

private:
  struct Gaussian {
    bool   multivariate; // this is required to discriminate the one dimensional case
    double height;
    std::vector<double> center;
    std::vector<double> sigma;
    std::vector<double> invsigma;
    Gaussian(const bool m, const double h, const std::vector<double>& c, const std::vector<double>& s):
      multivariate(m),height(h),center(c),sigma(s),invsigma(s) {
      // to avoid troubles from zero element in flexible hills
      for(unsigned i=0; i<invsigma.size(); ++i) {
        if(std::abs(invsigma[i])>1.e-20) invsigma[i]=1.0/invsigma[i] ;
        else invsigma[i]=0.0;
      }
    }
  };
  struct TemperingSpecs {
    bool is_active;
    std::string name_stem;
    std::string name;
    double biasf;
    double threshold;
    double alpha;
    inline TemperingSpecs(bool is_active, const std::string &name_stem, const std::string &name, double biasf, double threshold, double alpha) :
      is_active(is_active), name_stem(name_stem), name(name), biasf(biasf), threshold(threshold), alpha(alpha)
    {}
  };
  // general setup
  double kbt_;
  int stride_;
  bool calc_work_;
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
  std::vector<Gaussian> hills_;
  std::unique_ptr<FlexibleBin> flexbin_;
  int adaptive_;
  OFile hillsOfile_;
  std::vector<std::unique_ptr<IFile>> ifiles_;
  std::vector<std::string> ifilesnames_;
  // Grids
  bool grid_;
  std::unique_ptr<GridBase> BiasGrid_;
  OFile gridfile_;
  bool storeOldGrids_;
  int wgridstride_;
  // multiple walkers
  int mw_n_;
  std::string mw_dir_;
  int mw_id_;
  int mw_rstride_;
  bool walkers_mpi_;
  unsigned mpi_nw_;
  // flying gaussians
  bool flying_;
  // kinetics from metadynamics
  bool acceleration_;
  double acc_;
  double acc_restart_mean_;
  // transition-tempering metadynamics
  bool calc_max_bias_;
  double max_bias_;
  bool calc_transition_bias_;
  double transition_bias_;
  std::vector<std::vector<double> > transitionwells_;
  static const size_t n_tempering_options_ = 1;
  static const std::string tempering_names_[1][2];
  double dampfactor_;
  struct TemperingSpecs tt_specs_;
  std::string targetfilename_;
  std::unique_ptr<GridBase> TargetGrid_;
  // frequency adaptive metadynamics
  int current_stride_;
  bool freq_adaptive_;
  int fa_update_frequency_;
  int fa_max_stride_;
  double fa_min_acceleration_;
  // intervals
  double uppI_;
  double lowI_;
  bool doInt_;
  // reweighting
  bool calc_rct_;
  double reweight_factor_;
  unsigned rct_ustride_;
  // work
  double work_;
  // neighbour list stuff
  bool nlist_;
  bool nlist_update_;
  unsigned nlist_steps_;
  std::array<double,2> nlist_param_;
  std::vector<Gaussian> nlist_hills_;
  std::vector<double> nlist_center_;
  std::vector<double> nlist_dev2_;

  double stretchA=1.0;
  double stretchB=0.0;

  bool noStretchWarningDone=false;

  void noStretchWarning() {
    if(!noStretchWarningDone) {
      log<<"\nWARNING: you are using a HILLS file with Gaussian kernels, PLUMED 2.8 uses stretched Gaussians by default\n";
    }
    noStretchWarningDone=true;
  }

  static void registerTemperingKeywords(const std::string &name_stem, const std::string &name, Keywords &keys);
  void   readTemperingSpecs(TemperingSpecs &t_specs);
  void   logTemperingSpecs(const TemperingSpecs &t_specs);
  void   readGaussians(IFile*);
  void   writeGaussian(const Gaussian&,OFile&);
  void   addGaussian(const Gaussian&);
  double getHeight(const std::vector<double>&);
  void   temperHeight(double &height, const TemperingSpecs &t_specs, const double tempering_bias);
  double getBias(const std::vector<double>&);
  double getBiasAndDerivatives(const std::vector<double>&, std::vector<double>&);
  double evaluateGaussian(const std::vector<double>&, const Gaussian&);
  double evaluateGaussianAndDerivatives(const std::vector<double>&, const Gaussian&,std::vector<double>&,std::vector<double>&);
  double getGaussianNormalization(const Gaussian&);
  std::vector<unsigned> getGaussianSupport(const Gaussian&);
  bool   scanOneHill(IFile* ifile, std::vector<Value>& v, std::vector<double>& center, std::vector<double>& sigma, double& height, bool& multivariate);
  void   computeReweightingFactor();
  double getTransitionBarrierBias();
  void   updateFrequencyAdaptiveStride();
  void   updateNlist();

public:
  explicit MetaD(const ActionOptions&);
  void calculate() override;
  void update() override;
  static void registerKeywords(Keywords& keys);
  bool checkNeedsGradients()const override;
};

PLUMED_REGISTER_ACTION(MetaD,"METAD")

void MetaD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.addOutputComponent("rbias","CALC_RCT","the instantaneous value of the bias normalized using the c(t) reweighting factor [rbias=bias-rct]."
                          "This component can be used to obtain a reweighted histogram.");
  keys.addOutputComponent("rct","CALC_RCT","the reweighting factor \\f$c(t)\\f$.");
  keys.addOutputComponent("work","CALC_WORK","accumulator for work");
  keys.addOutputComponent("acc","ACCELERATION","the metadynamics acceleration factor");
  keys.addOutputComponent("maxbias", "CALC_MAX_BIAS", "the maximum of the metadynamics V(s, t)");
  keys.addOutputComponent("transbias", "CALC_TRANSITION_BIAS", "the metadynamics transition bias V*(t)");
  keys.addOutputComponent("pace","FREQUENCY_ADAPTIVE","the hill addition frequency when employing frequency adaptive metadynamics");
  keys.addOutputComponent("nlker","NLIST","number of hills in the neighbor list");
  keys.addOutputComponent("nlsteps","NLIST","number of steps from last neighbor list update");
  keys.use("ARG");
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","FILE","HILLS","a file in which the list of added hills is stored");
  keys.add("optional","HEIGHT","the heights of the Gaussian hills. Compulsory unless TAU and either BIASFACTOR or DAMPFACTOR are given");
  keys.add("optional","FMT","specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this bias factor.  Please note you must also specify temp");
  keys.addFlag("CALC_WORK",false,"calculate the total accumulated work done by the bias since last restart");
  keys.add("optional","RECT","list of bias factors for all the replicas");
  keys.add("optional","DAMPFACTOR","damp hills with exp(-max(V)/(kT*DAMPFACTOR)");
  for (size_t i = 0; i < n_tempering_options_; i++) {
    registerTemperingKeywords(tempering_names_[i][0], tempering_names_[i][1], keys);
  }
  keys.add("optional","TARGET","target to a predefined distribution");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","TAU","in well tempered metadynamics, sets height to (k_B Delta T*pace*timestep)/tau");
  keys.addFlag("CALC_RCT",false,"calculate the c(t) reweighting factor and use that to obtain the normalized bias [rbias=bias-rct]."
               "This method is not compatible with metadynamics not on a grid.");
  keys.add("optional","RCT_USTRIDE","the update stride for calculating the \\f$c(t)\\f$ reweighting factor."
           "The default 1, so \\f$c(t)\\f$ is updated every time the bias is updated.");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","GRID_WSTRIDE","write the grid to a file every N steps");
  keys.add("optional","GRID_WFILE","the file on which to write the grid");
  keys.add("optional","GRID_RFILE","a grid file from which the bias should be read at the initial step of the simulation");
  keys.addFlag("STORE_GRIDS",false,"store all the grid files the calculation generates. They will be deleted if this keyword is not present");
  keys.addFlag("NLIST",false,"Use neighbor list for kernels summation, faster but experimental");
  keys.add("optional", "NLIST_PARAMETERS","(default=6.,0.5) the two cutoff parameters for the Gaussians neighbor list");
  keys.add("optional","ADAPTIVE","use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or time step dimensions");
  keys.add("optional","SIGMA_MAX","the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","SIGMA_MIN","the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers - not compatible with WALKERS_* options other than WALKERS_DIR");
  keys.add("optional","INTERVAL","one dimensional lower and upper limits, outside the limits the system will not feel the biasing force.");
  keys.addFlag("FLYING_GAUSSIAN",false,"Switch on flying Gaussian method, must be used with WALKERS_MPI");
  keys.addFlag("ACCELERATION",false,"Set to TRUE if you want to compute the metadynamics acceleration factor.");
  keys.add("optional","ACCELERATION_RFILE","a data file from which the acceleration should be read at the initial step of the simulation");
  keys.addFlag("CALC_MAX_BIAS", false, "Set to TRUE if you want to compute the maximum of the metadynamics V(s, t)");
  keys.addFlag("CALC_TRANSITION_BIAS", false, "Set to TRUE if you want to compute a metadynamics transition bias V*(t)");
  keys.add("numbered", "TRANSITIONWELL", "This keyword appears multiple times as TRANSITIONWELL followed by an integer. Each specifies the coordinates for one well as in transition-tempered metadynamics. At least one must be provided.");
  keys.addFlag("FREQUENCY_ADAPTIVE",false,"Set to TRUE if you want to enable frequency adaptive metadynamics such that the frequency for hill addition to change dynamically based on the acceleration factor.");
  keys.add("optional","FA_UPDATE_FREQUENCY","the frequency for updating the hill addition pace in frequency adaptive metadynamics, by default this is equal to the value given in PACE");
  keys.add("optional","FA_MAX_PACE","the maximum hill addition frequency allowed in frequency adaptive metadynamics. By default there is no maximum value.");
  keys.add("optional","FA_MIN_ACCELERATION","only update the hill addition pace in frequency adaptive metadynamics after reaching the minimum acceleration factor given here. By default it is 1.0.");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

const std::string MetaD::tempering_names_[1][2] = {{"TT", "transition tempered"}};

void MetaD::registerTemperingKeywords(const std::string &name_stem, const std::string &name, Keywords &keys) {
  keys.add("optional", name_stem + "BIASFACTOR", "use " + name + " metadynamics with this bias factor.  Please note you must also specify temp");
  keys.add("optional", name_stem + "BIASTHRESHOLD", "use " + name + " metadynamics with this bias threshold.  Please note you must also specify " + name_stem + "BIASFACTOR");
  keys.add("optional", name_stem + "ALPHA", "use " + name + " metadynamics with this hill size decay exponent parameter.  Please note you must also specify " + name_stem + "BIASFACTOR");
}

MetaD::MetaD(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao),
  kbt_(0.0),
  stride_(0),
  calc_work_(false),
  welltemp_(false),
  biasf_(-1.0),
  isFirstStep_(true),
  height0_(std::numeric_limits<double>::max()),
  adaptive_(FlexibleBin::none),
  grid_(false),
  wgridstride_(0),
  mw_n_(1), mw_dir_(""), mw_id_(0), mw_rstride_(1),
  walkers_mpi_(false), mpi_nw_(0),
  flying_(false),
  acceleration_(false), acc_(0.0), acc_restart_mean_(0.0),
  calc_max_bias_(false), max_bias_(0.0),
  calc_transition_bias_(false), transition_bias_(0.0),
  dampfactor_(0.0),
  tt_specs_(false, "TT", "Transition Tempered", -1.0, 0.0, 1.0),
  current_stride_(0),
  freq_adaptive_(false),
  fa_update_frequency_(0),
  fa_max_stride_(0),
  fa_min_acceleration_(1.0),
  uppI_(-1), lowI_(-1), doInt_(false),
  calc_rct_(false),
  reweight_factor_(0.0),
  rct_ustride_(1),
  work_(0),
  nlist_(false),
  nlist_update_(false),
  nlist_steps_(0)
{
  if(!dp2cutoffNoStretch()) {
    stretchA=dp2cutoffA;
    stretchB=dp2cutoffB;
  }
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

  // parse the sigma
  parseVector("SIGMA",sigma0_);
  if(adaptive_==FlexibleBin::none) {
    // if you use normal sigma you need one sigma per argument
    if( sigma0_.size()!=getNumberOfArguments() ) error("number of arguments does not match number of SIGMA parameters");
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
    if(sigma0min_.size()>0 && sigma0min_.size()!=getNumberOfArguments()) {
      error("the number of SIGMA_MIN values be the same of the number of the arguments");
    } else if(sigma0min_.size()==0) {
      sigma0min_.resize(getNumberOfArguments());
      for(unsigned i=0; i<getNumberOfArguments(); i++) {sigma0min_[i]=-1.;}
    }

    parseVector("SIGMA_MAX",sigma0max_);
    if(sigma0max_.size()>0 && sigma0max_.size()!=getNumberOfArguments()) {
      error("the number of SIGMA_MAX values be the same of the number of the arguments");
    } else if(sigma0max_.size()==0) {
      sigma0max_.resize(getNumberOfArguments());
      for(unsigned i=0; i<getNumberOfArguments(); i++) {sigma0max_[i]=-1.;}
    }

    flexbin_=Tools::make_unique<FlexibleBin>(adaptive_,this,sigma0_[0],sigma0min_,sigma0max_);
  }

  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  parse("PACE",stride_);
  if(stride_<=0) error("frequency for hill addition is nonsensical");
  current_stride_ = stride_;
  std::string hillsfname="HILLS";
  parse("FILE",hillsfname);

  // Manually set to calculate special bias quantities
  // throughout the course of simulation. (These are chosen due to
  // relevance for tempering and event-driven logic as well.)
  parseFlag("CALC_MAX_BIAS", calc_max_bias_);
  parseFlag("CALC_TRANSITION_BIAS", calc_transition_bias_);

  std::vector<double> rect_biasf_;
  parseVector("RECT",rect_biasf_);
  if(rect_biasf_.size()>0) {
    int r=0;
    if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
    comm.Bcast(r,0);
    biasf_=rect_biasf_[r];
    log<<"  You are using RECT\n";
  } else {
    parse("BIASFACTOR",biasf_);
  }
  if( biasf_<1.0  && biasf_!=-1.0) error("well tempered bias factor is nonsensical");
  parse("DAMPFACTOR",dampfactor_);
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
  if(biasf_>=1.0) {
    if(kbt_==0.0) error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    welltemp_=true;
  }
  if(dampfactor_>0.0) {
    if(kbt_==0.0) error("Unless the MD engine passes the temperature to plumed, with damped metad you must specify it using TEMP");
  }

  parseFlag("CALC_WORK",calc_work_);

  // Set transition tempering parameters.
  // Transition wells are read later via calc_transition_bias_.
  readTemperingSpecs(tt_specs_);
  if (tt_specs_.is_active) calc_transition_bias_ = true;

  // If any previous option specified to calculate a transition bias,
  // now read the transition wells for that quantity.
  if (calc_transition_bias_) {
    std::vector<double> tempcoords(getNumberOfArguments());
    for (unsigned i = 0; ; i++) {
      if (!parseNumberedVector("TRANSITIONWELL", i, tempcoords) ) break;
      if (tempcoords.size() != getNumberOfArguments()) {
        error("incorrect number of coordinates for transition tempering well");
      }
      transitionwells_.push_back(tempcoords);
    }
  }

  parse("TARGET",targetfilename_);
  if(targetfilename_.length()>0 && kbt_==0.0)  error("with TARGET temperature must be specified");
  double tau=0.0;
  parse("TAU",tau);
  if(tau==0.0) {
    if(height0_==std::numeric_limits<double>::max()) error("At least one between HEIGHT and TAU should be specified");
    // if tau is not set, we compute it here from the other input parameters
    if(welltemp_) tau=(kbt_*(biasf_-1.0))/height0_*getTimeStep()*stride_;
    else if(dampfactor_>0.0) tau=(kbt_*dampfactor_)/height0_*getTimeStep()*stride_;
  } else {
    if(height0_!=std::numeric_limits<double>::max()) error("At most one between HEIGHT and TAU should be specified");
    if(welltemp_) {
      if(biasf_!=1.0) height0_=(kbt_*(biasf_-1.0))/tau*getTimeStep()*stride_;
      else           height0_=kbt_/tau*getTimeStep()*stride_; // special case for gamma=1
    }
    else if(dampfactor_>0.0) height0_=(kbt_*dampfactor_)/tau*getTimeStep()*stride_;
    else error("TAU only makes sense in well-tempered or damped metadynamics");
  }

  // Grid Stuff
  std::vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  std::vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  std::vector<unsigned> gbin(getNumberOfArguments());
  std::vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) error("not enough values for GRID_SPACING");
  if(gmin.size()!=gmax.size()) error("GRID_MAX and GRID_MIN should be either present or absent");
  if(gspacing.size()!=0 && gmin.size()==0) error("If GRID_SPACING is present also GRID_MIN and GRID_MAX should be present");
  if(gbin.size()!=0     && gmin.size()==0) error("If GRID_BIN is present also GRID_MIN and GRID_MAX should be present");
  if(gmin.size()!=0) {
    if(gbin.size()==0 && gspacing.size()==0) {
      if(adaptive_==FlexibleBin::none) {
        log<<"  Binsize not specified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0_.size()==getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for(unsigned i=0; i<gspacing.size(); i++) gspacing[i]=0.2*sigma0_[i];
      } else {
        // with adaptive hills and grid a sigma min must be specified
        for(unsigned i=0; i<sigma0min_.size(); i++) if(sigma0min_[i]<=0) error("When using ADAPTIVE Gaussians on a grid SIGMA_MIN must be specified");
        log<<"  Binsize not specified, 1/5 of sigma_min will be be used\n";
        gspacing.resize(getNumberOfArguments());
        for(unsigned i=0; i<gspacing.size(); i++) gspacing[i]=0.2*sigma0min_[i];
      }
    } else if(gspacing.size()!=0 && gbin.size()==0) {
      log<<"  The number of bins will be estimated from GRID_SPACING\n";
    } else if(gspacing.size()!=0 && gbin.size()!=0) {
      log<<"  You specified both GRID_BIN and GRID_SPACING\n";
      log<<"  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if(gbin.size()==0) gbin.assign(getNumberOfArguments(),1);
    if(gspacing.size()!=0) for(unsigned i=0; i<getNumberOfArguments(); i++) {
        double a,b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        unsigned n=std::ceil(((b-a)/gspacing[i]));
        if(gbin[i]<n) gbin[i]=n;
      }
  }
  if(gbin.size()>0) grid_=true;

  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  parse("GRID_WSTRIDE",wgridstride_);
  std::string gridfilename_;
  parse("GRID_WFILE",gridfilename_);
  parseFlag("STORE_GRIDS",storeOldGrids_);
  if(grid_ && gridfilename_.length()>0) {
    if(wgridstride_==0 ) error("frequency with which to output grid not specified use GRID_WSTRIDE");
  }
  if(grid_ && wgridstride_>0) {
    if(gridfilename_.length()==0) error("grid filename not specified use GRID_WFILE");
  }

  std::string gridreadfilename_;
  parse("GRID_RFILE",gridreadfilename_);

  if(!grid_&&gridfilename_.length()> 0) error("To write a grid you need first to define it!");
  if(!grid_&&gridreadfilename_.length()>0) error("To read a grid you need first to define it!");

  /*setup neighbor list stuff*/
  parseFlag("NLIST", nlist_);
  nlist_center_.resize(getNumberOfArguments());
  nlist_dev2_.resize(getNumberOfArguments());
  if(nlist_&&grid_) error("NLIST and GRID cannot be combined!");
  std::vector<double> nlist_param;
  parseVector("NLIST_PARAMETERS",nlist_param);
  if(nlist_param.size()==0)
  {
    nlist_param_[0]=6.0;//*DP2CUTOFF -> max distance of neighbors
    nlist_param_[1]=0.5;//*nlist_dev2_[i] -> condition for rebuilding
  }
  else
  {
    plumed_massert(nlist_param.size()==2,"two cutoff parameters are needed for the neighbor list");
    plumed_massert(nlist_param[0]>1.0,"the first of NLIST_PARAMETERS must be greater than 1. The smaller the first, the smaller should be the second as well");
    const double min_PARAM_1=(1.-1./std::sqrt(nlist_param[0]/2))+0.16;
    plumed_massert(nlist_param[1]>0,"the second of NLIST_PARAMETERS must be greater than 0");
    plumed_massert(nlist_param[1]<=min_PARAM_1,"the second of NLIST_PARAMETERS must be smaller to avoid systematic errors. Largest suggested value is: 1.16-1/sqrt(PARAM_0/2) = "+std::to_string(min_PARAM_1));
    nlist_param_[0]=nlist_param[0];
    nlist_param_[1]=nlist_param[1];
  }

  // Reweighting factor rct
  parseFlag("CALC_RCT",calc_rct_);
  if (calc_rct_) plumed_massert(grid_,"CALC_RCT is supported only if bias is on a grid");
  parse("RCT_USTRIDE",rct_ustride_);

  if(dampfactor_>0.0) {
    if(!grid_) error("With DAMPFACTOR you should use grids");
  }

  // Multiple walkers
  parse("WALKERS_N",mw_n_);
  parse("WALKERS_ID",mw_id_);
  if(mw_n_<=mw_id_) error("walker ID should be a numerical value less than the total number of walkers");
  parse("WALKERS_DIR",mw_dir_);
  parse("WALKERS_RSTRIDE",mw_rstride_);

  // MPI version
  parseFlag("WALKERS_MPI",walkers_mpi_);

  //If this Action is not compiled with MPI the user is informed and we exit gracefully
  if(walkers_mpi_) {
    plumed_assert(Communicator::plumedHasMPI()) << "Invalid walkers configuration: WALKERS_MPI flag requires MPI compilation";
    plumed_assert(Communicator::initialized()) << "Invalid walkers configuration: WALKERS_MPI needs the communicator correctly initialized.";
  }

  // Flying Gaussian
  parseFlag("FLYING_GAUSSIAN", flying_);

  // Inteval keyword
  std::vector<double> tmpI(2);
  parseVector("INTERVAL",tmpI);
  if(tmpI.size()!=2&&tmpI.size()!=0) error("both a lower and an upper limits must be provided with INTERVAL");
  else if(tmpI.size()==2) {
    lowI_=tmpI.at(0);
    uppI_=tmpI.at(1);
    if(getNumberOfArguments()!=1) error("INTERVAL limits correction works only for monodimensional metadynamics!");
    if(uppI_<lowI_) error("The Upper limit must be greater than the Lower limit!");
    if(getPntrToArgument(0)->isPeriodic()) error("INTERVAL cannot be used with periodic variables!");
    doInt_=true;
  }

  parseFlag("ACCELERATION",acceleration_);
  // Check for a restart acceleration if acceleration is active.
  std::string acc_rfilename;
  if(acceleration_) {
    parse("ACCELERATION_RFILE", acc_rfilename);
  }

  freq_adaptive_=false;
  parseFlag("FREQUENCY_ADAPTIVE",freq_adaptive_);
  //
  fa_update_frequency_=0;
  parse("FA_UPDATE_FREQUENCY",fa_update_frequency_);
  if(fa_update_frequency_!=0 && !freq_adaptive_) {
    plumed_merror("It doesn't make sense to use the FA_MAX_PACE keyword if frequency adaptive METAD hasn't been activated by using the FREQUENCY_ADAPTIVE flag");
  }
  if(fa_update_frequency_==0 && freq_adaptive_) {
    fa_update_frequency_=stride_;
  }
  //
  fa_max_stride_=0;
  parse("FA_MAX_PACE",fa_max_stride_);
  if(fa_max_stride_!=0 && !freq_adaptive_) {
    plumed_merror("It doesn't make sense to use the FA_MAX_PACE keyword if frequency adaptive METAD hasn't been activated by using the FREQUENCY_ADAPTIVE flag");
  }
  //
  fa_min_acceleration_=1.0;
  parse("FA_MIN_ACCELERATION",fa_min_acceleration_);
  if(fa_min_acceleration_!=1.0 && !freq_adaptive_) {
    plumed_merror("It doesn't make sense to use the FA_MIN_ACCELERATION keyword if frequency adaptive METAD hasn't been activated by using the FREQUENCY_ADAPTIVE flag");
  }

  checkRead();

  log.printf("  Gaussian width ");
  if (adaptive_==FlexibleBin::diffusion)log.printf(" (Note: The units of sigma are in timesteps) ");
  if (adaptive_==FlexibleBin::geometry)log.printf(" (Note: The units of sigma are in dist units) ");
  for(unsigned i=0; i<sigma0_.size(); ++i) log.printf(" %f",sigma0_[i]);
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_);
  log.printf("  Gaussian file %s\n",hillsfname.c_str());
  if(welltemp_) {
    log.printf("  Well-Tempered Bias Factor %f\n",biasf_);
    log.printf("  Hills relaxation time (tau) %f\n",tau);
    log.printf("  KbT %f\n",kbt_);
  }

  // Transition tempered metadynamics options
  if (tt_specs_.is_active) {
    logTemperingSpecs(tt_specs_);
    // Check that the appropriate transition bias quantity is calculated.
    // (Should never trip, given that the flag is automatically set.)
    if (!calc_transition_bias_) {
      error(" transition tempering requires calculation of a transition bias");
    }
  }

  // Overall tempering sanity check (this gets tricky when multiple are active).
  // When multiple temperings are active, it's fine to have one tempering attempt
  // to increase hill size with increasing bias, so long as the others can shrink
  // the hills faster than it increases their size in the long-time limit.
  // This set of checks ensures that the hill sizes eventually decay to zero as c(t)
  // diverges to infinity.
  // The alpha parameter allows hills to decay as 1/t^alpha instead of 1/t,
  // a slower decay, so as t -> infinity, only the temperings with the largest
  // alphas govern the final asymptotic decay. (Alpha helps prevent false convergence.)
  if (welltemp_ || dampfactor_ > 0.0 || tt_specs_.is_active) {
    // Determine the number of active temperings.
    int n_active = 0;
    if (welltemp_) n_active++;
    if (dampfactor_ > 0.0) n_active++;
    if (tt_specs_.is_active) n_active++;
    // Find the greatest alpha.
    double greatest_alpha = 0.0;
    if (welltemp_) greatest_alpha = std::max(greatest_alpha, 1.0);
    if (dampfactor_ > 0.0) greatest_alpha = std::max(greatest_alpha, 1.0);
    if (tt_specs_.is_active) greatest_alpha = std::max(greatest_alpha, tt_specs_.alpha);
    // Find the least alpha.
    double least_alpha = 1.0;
    if (welltemp_) least_alpha = std::min(least_alpha, 1.0);
    if (dampfactor_ > 0.0) least_alpha = std::min(least_alpha, 1.0);
    if (tt_specs_.is_active) least_alpha = std::min(least_alpha, tt_specs_.alpha);
    // Find the inverse harmonic average of the delta T parameters for all
    // of the temperings with the greatest alpha values.
    double total_governing_deltaT_inv = 0.0;
    if (welltemp_ && 1.0 == greatest_alpha && biasf_ != 1.0) total_governing_deltaT_inv += 1.0 / (biasf_ - 1.0);
    if (dampfactor_ > 0.0 && 1.0 == greatest_alpha) total_governing_deltaT_inv += 1.0 / (dampfactor_);
    if (tt_specs_.is_active && tt_specs_.alpha == greatest_alpha) total_governing_deltaT_inv += 1.0 / (tt_specs_.biasf - 1.0);
    // Give a newbie-friendly error message for people using one tempering if
    // only one is active.
    if (n_active == 1 && total_governing_deltaT_inv < 0.0) {
      error("for stable tempering, the bias factor must be greater than one");
      // Give a slightly more complex error message to users stacking multiple
      // tempering options at a time, but all with uniform alpha values.
    } else if (total_governing_deltaT_inv < 0.0 && greatest_alpha == least_alpha) {
      error("for stable tempering, the sum of the inverse Delta T parameters must be greater than zero!");
      // Give the most technical error message to users stacking multiple tempering
      // options with different alpha parameters.
    } else if (total_governing_deltaT_inv < 0.0 && greatest_alpha != least_alpha) {
      error("for stable tempering, the sum of the inverse Delta T parameters for the greatest asymptotic hill decay exponents must be greater than zero!");
    }
  }

  if(doInt_) log.printf("  Upper and Lower limits boundaries for the bias are activated at %f - %f\n", lowI_, uppI_);

  if(grid_) {
    log.printf("  Grid min");
    for(unsigned i=0; i<gmin.size(); ++i) log.printf(" %s",gmin[i].c_str() );
    log.printf("\n");
    log.printf("  Grid max");
    for(unsigned i=0; i<gmax.size(); ++i) log.printf(" %s",gmax[i].c_str() );
    log.printf("\n");
    log.printf("  Grid bin");
    for(unsigned i=0; i<gbin.size(); ++i) log.printf(" %u",gbin[i]);
    log.printf("\n");
    if(spline) {log.printf("  Grid uses spline interpolation\n");}
    if(sparsegrid) {log.printf("  Grid uses sparse grid\n");}
    if(wgridstride_>0) {log.printf("  Grid is written on file %s with stride %d\n",gridfilename_.c_str(),wgridstride_);}
  }

  if(mw_n_>1) {
    if(walkers_mpi_) error("MPI version of multiple walkers is not compatible with filesystem version of multiple walkers");
    log.printf("  %d multiple walkers active\n",mw_n_);
    log.printf("  walker id %d\n",mw_id_);
    log.printf("  reading stride %d\n",mw_rstride_);
    if(mw_dir_!="")log.printf("  directory with hills files %s\n",mw_dir_.c_str());
  } else {
    if(walkers_mpi_) {
      log.printf("  Multiple walkers active using MPI communnication\n");
      if(mw_dir_!="")log.printf("  directory with hills files %s\n",mw_dir_.c_str());
      if(comm.Get_rank()==0) {
        // Only root of group can communicate with other walkers
        mpi_nw_=multi_sim_comm.Get_size();
      }
      // Communicate to the other members of the same group
      // info abount number of walkers and walker index
      comm.Bcast(mpi_nw_,0);
    }
  }

  if(flying_) {
    if(!walkers_mpi_) error("Flying Gaussian method must be used with MPI version of multiple walkers");
    log.printf("  Flying Gaussian method with %d walkers active\n",mpi_nw_);
  }

  if(nlist_) {
    addComponent("nlker");
    componentIsNotPeriodic("nlker");
    addComponent("nlsteps");
    componentIsNotPeriodic("nlsteps");
  }

  if(calc_rct_) {
    addComponent("rbias"); componentIsNotPeriodic("rbias");
    addComponent("rct"); componentIsNotPeriodic("rct");
    log.printf("  The c(t) reweighting factor will be calculated every %u hills\n",rct_ustride_);
    getPntrToComponent("rct")->set(reweight_factor_);
  }

  if(calc_work_) { addComponent("work"); componentIsNotPeriodic("work"); }

  if(acceleration_) {
    if (kbt_ == 0.0) {
      error("The calculation of the acceleration works only if simulation temperature has been defined");
    }
    log.printf("  calculation on the fly of the acceleration factor\n");
    addComponent("acc"); componentIsNotPeriodic("acc");
    // Set the initial value of the the acceleration.
    // If this is not a restart, set to 1.0.
    if (acc_rfilename.length() == 0) {
      getPntrToComponent("acc")->set(1.0);
      if(getRestart()) {
        log.printf("  WARNING: calculating the acceleration factor in a restarted run without reading in the previous value will most likely lead to incorrect results.\n");
        log.printf("           You should use the ACCELERATION_RFILE keyword.\n");
      }
      // Otherwise, read and set the restart value.
    } else {
      // Restart of acceleration does not make sense if the restart timestep is zero.
      //if (getStep() == 0) {
      //  error("Restarting calculation of acceleration factors works only if simulation timestep is restarted correctly");
      //}
      // Open the ACCELERATION_RFILE.
      IFile acc_rfile;
      acc_rfile.link(*this);
      if(acc_rfile.FileExist(acc_rfilename)) {
        acc_rfile.open(acc_rfilename);
      } else {
        error("The ACCELERATION_RFILE file you want to read: " + acc_rfilename + ", cannot be found!");
      }
      // Read the file to find the restart acceleration.
      double acc_rmean=0.0;
      double acc_rtime=0.0;
      bool found=false;
      std::string acclabel = getLabel() + ".acc";
      acc_rfile.allowIgnoredFields();
      while(acc_rfile.scanField("time", acc_rtime)) {
        acc_rfile.scanField(acclabel, acc_rmean);
        acc_rfile.scanField();
        found=true;
      }
      if(!found) error("The ACCELERATION_RFILE file you want to read: " + acc_rfilename + ", does not contain a time field!");
      acc_restart_mean_ = acc_rmean;
      // Set component based on the read values.
      getPntrToComponent("acc")->set(acc_rmean);
      log.printf("  initial acceleration factor read from file %s: value of %f at time %f\n",acc_rfilename.c_str(),acc_rmean,acc_rtime);
    }
  }

  if (calc_max_bias_) {
    if (!grid_) error("Calculating the maximum bias on the fly works only with a grid");
    log.printf("  calculation on the fly of the maximum bias max(V(s,t)) \n");
    addComponent("maxbias");
    componentIsNotPeriodic("maxbias");
  }

  if (calc_transition_bias_) {
    if (!grid_) error("Calculating the transition bias on the fly works only with a grid");
    log.printf("  calculation on the fly of the transition bias V*(t)\n");
    addComponent("transbias");
    componentIsNotPeriodic("transbias");
    log<<"  Number of transition wells "<<transitionwells_.size()<<"\n";
    if (transitionwells_.size() == 0) error("Calculating the transition bias on the fly requires definition of at least one transition well");
    // Check that a grid is in use.
    if (!grid_) error(" transition barrier finding requires a grid for the bias");
    // Log the wells and check that they are in the grid.
    for (unsigned i = 0; i < transitionwells_.size(); i++) {
      // Log the coordinate.
      log.printf("  Transition well %d at coordinate ", i);
      for (unsigned j = 0; j < getNumberOfArguments(); j++) log.printf("%f ", transitionwells_[i][j]);
      log.printf("\n");
      // Check that the coordinate is in the grid.
      for (unsigned j = 0; j < getNumberOfArguments(); j++) {
        double max, min;
        Tools::convert(gmin[j], min);
        Tools::convert(gmax[j], max);
        if (transitionwells_[i][j] < min || transitionwells_[i][j] > max) error(" transition well is not in grid");
      }
    }
  }

  if(freq_adaptive_) {
    if(!acceleration_) {
      plumed_merror("Frequency adaptive metadynamics only works if the calculation of the acceleration factor is enabled with the ACCELERATION keyword\n");
    }
    if(walkers_mpi_) {
      plumed_merror("Combining frequency adaptive metadynamics with MPI multiple walkers is not allowed");
    }

    log.printf("  Frequency adaptive metadynamics enabled\n");
    if(getRestart() && acc_rfilename.length() == 0) {
      log.printf("  WARNING: using the frequency adaptive scheme in a restarted run without reading in the previous value of the acceleration factor will most likely lead to incorrect results.\n");
      log.printf("           You should use the ACCELERATION_RFILE keyword.\n");
    }
    log.printf("  The frequency for hill addition will change dynamically based on the metadynamics acceleration factor\n");
    log.printf("  The hill addition frequency will be updated every %d steps\n",fa_update_frequency_);
    if(fa_min_acceleration_>1.0) {
      log.printf("  The hill addition frequency will only be updated once the metadynamics acceleration factor becomes larger than %.1f \n",fa_min_acceleration_);
    }
    if(fa_max_stride_!=0) {
      log.printf("  The hill addition frequency will not become larger than %d steps\n",fa_max_stride_);
    }
    addComponent("pace"); componentIsNotPeriodic("pace");
    updateFrequencyAdaptiveStride();
  }

  // initializing and checking grid
  bool restartedFromGrid=false;  // restart from grid file
  if(grid_) {
    if(!(gridreadfilename_.length()>0)) {
      // check for mesh and sigma size
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        double a,b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        double mesh=(b-a)/((double)gbin[i]);
        if(adaptive_==FlexibleBin::none) {
          if(mesh>0.5*sigma0_[i]) log<<"  WARNING: Using a METAD with a Grid Spacing larger than half of the Gaussians width (SIGMA) can produce artifacts\n";
        } else {
          if(sigma0min_[i]<0.) error("When using ADAPTIVE Gaussians on a grid SIGMA_MIN must be specified");
          if(mesh>0.5*sigma0min_[i]) log<<"  WARNING: to use a METAD with a GRID and ADAPTIVE you need to set a Grid Spacing lower than half of the Gaussians (SIGMA_MIN) \n";
        }
      }
      std::string funcl=getLabel() + ".bias";
      if(!sparsegrid) {BiasGrid_=Tools::make_unique<Grid>(funcl,getArguments(),gmin,gmax,gbin,spline,true);}
      else {BiasGrid_=Tools::make_unique<SparseGrid>(funcl,getArguments(),gmin,gmax,gbin,spline,true);}
      std::vector<std::string> actualmin=BiasGrid_->getMin();
      std::vector<std::string> actualmax=BiasGrid_->getMax();
      for(unsigned i=0; i<getNumberOfArguments(); i++) {
        std::string is;
        Tools::convert(i,is);
        if(gmin[i]!=actualmin[i]) error("GRID_MIN["+is+"] must be adjusted to "+actualmin[i]+" to fit periodicity");
        if(gmax[i]!=actualmax[i]) error("GRID_MAX["+is+"] must be adjusted to "+actualmax[i]+" to fit periodicity");
      }
    } else {
      // read the grid in input, find the keys
#ifdef __PLUMED_HAS_GETCWD
      if(walkers_mpi_&&gridreadfilename_.at(0)!='/') {
        //if possible the root replica will share its current folder so that all walkers will read the same file
        char cwd[4096]= {0};
        const char* ret=getcwd(cwd,4096);
        plumed_assert(ret)<<"Name of current directory too long, increase buffer size";
        gridreadfilename_ = "/" + gridreadfilename_;
        gridreadfilename_ = ret + gridreadfilename_;
        if(comm.Get_rank()==0) multi_sim_comm.Bcast(gridreadfilename_,0);
        comm.Bcast(gridreadfilename_,0);
      }
#endif
      IFile gridfile;
      gridfile.link(*this);
      if(gridfile.FileExist(gridreadfilename_)) {
        gridfile.open(gridreadfilename_);
      } else {
        error("The GRID file you want to read: " + gridreadfilename_ + ", cannot be found!");
      }
      std::string funcl=getLabel() + ".bias";
      BiasGrid_=GridBase::create(funcl, getArguments(), gridfile, gmin, gmax, gbin, sparsegrid, spline, true);
      if(BiasGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
        if( getPntrToArgument(i)->isPeriodic()!=BiasGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias");
        double a, b;
        Tools::convert(gmin[i],a);
        Tools::convert(gmax[i],b);
        double mesh=(b-a)/((double)gbin[i]);
        if(mesh>0.5*sigma0_[i]) log<<"  WARNING: Using a METAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
      }
      log.printf("  Restarting from %s\n",gridreadfilename_.c_str());
      if(getRestart()) restartedFromGrid=true;
    }
  }

  // if we are restarting from GRID and using WALKERS_MPI we can check that all walkers have actually read the grid
  if(getRestart()&&walkers_mpi_) {
    std::vector<int> restarted(mpi_nw_,0);
    if(comm.Get_rank()==0) multi_sim_comm.Allgather(int(restartedFromGrid), restarted);
    comm.Bcast(restarted,0);
    int result = std::accumulate(restarted.begin(),restarted.end(),0);
    if(result!=0&&result!=mpi_nw_) error("in this WALKERS_MPI run some replica have restarted from GRID while other do not!");
  }

#ifdef __PLUMED_HAS_GETCWD
  if(walkers_mpi_&&mw_dir_==""&&hillsfname.at(0)!='/') {
    //if possible the root replica will share its current folder so that all walkers will read the same file
    char cwd[4096]= {0};
    const char* ret=getcwd(cwd,4096);
    plumed_assert(ret)<<"Name of current directory too long, increase buffer size";
    mw_dir_ = ret;
    mw_dir_ = mw_dir_ + "/";
    if(comm.Get_rank()==0) multi_sim_comm.Bcast(mw_dir_,0);
    comm.Bcast(mw_dir_,0);
  }
#endif

  // creating std::vector of ifile* for hills reading
  // open all files at the beginning and read Gaussians if restarting
  bool restartedFromHills=false;  // restart from hills files
  for(int i=0; i<mw_n_; ++i) {
    std::string fname;
    if(mw_dir_!="") {
      if(mw_n_>1) {
        std::stringstream out; out << i;
        fname = mw_dir_+"/"+hillsfname+"."+out.str();
      } else if(walkers_mpi_) {
        fname = mw_dir_+"/"+hillsfname;
      } else {
        fname = hillsfname;
      }
    } else {
      if(mw_n_>1) {
        std::stringstream out; out << i;
        fname = hillsfname+"."+out.str();
      } else {
        fname = hillsfname;
      }
    }
    ifiles_.emplace_back(Tools::make_unique<IFile>());
    // this is just a shortcut pointer to the last element:
    IFile *ifile = ifiles_.back().get();
    ifilesnames_.push_back(fname);
    ifile->link(*this);
    if(ifile->FileExist(fname)) {
      ifile->open(fname);
      if(getRestart()&&!restartedFromGrid) {
        log.printf("  Restarting from %s:",ifilesnames_[i].c_str());
        readGaussians(ifiles_[i].get());
        restartedFromHills=true;
      }
      ifiles_[i]->reset(false);
      // close only the walker own hills file for later writing
      if(i==mw_id_) ifiles_[i]->close();
    } else {
      // in case a file does not exist and we are restarting, complain that the file was not found
      if(getRestart()&&!restartedFromGrid) error("restart file "+fname+" not found");
    }
  }

  // if we are restarting from FILE and using WALKERS_MPI we can check that all walkers have actually read the FILE
  if(getRestart()&&walkers_mpi_) {
    std::vector<int> restarted(mpi_nw_,0);
    if(comm.Get_rank()==0) multi_sim_comm.Allgather(int(restartedFromHills), restarted);
    comm.Bcast(restarted,0);
    int result = std::accumulate(restarted.begin(),restarted.end(),0);
    if(result!=0&&result!=mpi_nw_) error("in this WALKERS_MPI run some replica have restarted from FILE while other do not!");
  }

  comm.Barrier();
  // this barrier is needed when using walkers_mpi
  // to be sure that all files have been read before
  // backing them up
  // it should not be used when walkers_mpi is false otherwise
  // it would introduce troubles when using replicas without METAD
  // (e.g. in bias exchange with a neutral replica)
  // see issue #168 on github
  if(comm.Get_rank()==0 && walkers_mpi_) multi_sim_comm.Barrier();

  if(targetfilename_.length()>0) {
    IFile gridfile; gridfile.open(targetfilename_);
    std::string funcl=getLabel() + ".target";
    TargetGrid_=GridBase::create(funcl,getArguments(),gridfile,false,false,true);
    if(TargetGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->isPeriodic()!=TargetGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias");
    }
  }

  if(getRestart()) {
    // if this is a restart the neighbor list should be immediately updated
    if(nlist_) nlist_update_=true;
    // Calculate the Tiwary-Parrinello reweighting factor if we are restarting from previous hills
    if(calc_rct_) computeReweightingFactor();
    // Calculate all special bias quantities desired if restarting with nonzero bias.
    if(calc_max_bias_) {
      max_bias_ = BiasGrid_->getMaxValue();
      getPntrToComponent("maxbias")->set(max_bias_);
    }
    if(calc_transition_bias_) {
      transition_bias_ = getTransitionBarrierBias();
      getPntrToComponent("transbias")->set(transition_bias_);
    }
  }

  // open grid file for writing
  if(wgridstride_>0) {
    gridfile_.link(*this);
    if(walkers_mpi_) {
      int r=0;
      if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
      comm.Bcast(r,0);
      if(r>0) gridfilename_="/dev/null";
      gridfile_.enforceSuffix("");
    }
    if(mw_n_>1) gridfile_.enforceSuffix("");
    gridfile_.open(gridfilename_);
  }

  // open hills file for writing
  hillsOfile_.link(*this);
  if(walkers_mpi_) {
    int r=0;
    if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
    comm.Bcast(r,0);
    if(r>0) ifilesnames_[mw_id_]="/dev/null";
    hillsOfile_.enforceSuffix("");
  }
  if(mw_n_>1) hillsOfile_.enforceSuffix("");
  hillsOfile_.open(ifilesnames_[mw_id_]);
  if(fmt_.length()>0) hillsOfile_.fmtField(fmt_);
  hillsOfile_.addConstantField("multivariate");
  hillsOfile_.addConstantField("kerneltype");
  if(doInt_) {
    hillsOfile_.addConstantField("lower_int").printField("lower_int",lowI_);
    hillsOfile_.addConstantField("upper_int").printField("upper_int",uppI_);
  }
  hillsOfile_.setHeavyFlush();
  // output periodicities of variables
  for(unsigned i=0; i<getNumberOfArguments(); ++i) hillsOfile_.setupPrintValue( getPntrToArgument(i) );

  bool concurrent=false;
  const ActionSet&actionSet(plumed.getActionSet());
  for(const auto & p : actionSet) if(dynamic_cast<MetaD*>(p.get())) { concurrent=true; break; }
  if(concurrent) log<<"  You are using concurrent metadynamics\n";
  if(rect_biasf_.size()>0) {
    if(walkers_mpi_) {
      log<<"  You are using RECT in its 'altruistic' implementation\n";
    }{
      log<<"  You are using RECT\n";
    }
  }

  log<<"  Bibliography "<<plumed.cite("Laio and Parrinello, PNAS 99, 12562 (2002)");
  if(welltemp_) log<<plumed.cite("Barducci, Bussi, and Parrinello, Phys. Rev. Lett. 100, 020603 (2008)");
  if(tt_specs_.is_active) {
    log << plumed.cite("Dama, Rotskoff, Parrinello, and Voth, J. Chem. Theory Comput. 10, 3626 (2014)");
    log << plumed.cite("Dama, Parrinello, and Voth, Phys. Rev. Lett. 112, 240602 (2014)");
  }
  if(mw_n_>1||walkers_mpi_) log<<plumed.cite("Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");
  if(adaptive_!=FlexibleBin::none) log<<plumed.cite("Branduardi, Bussi, and Parrinello, J. Chem. Theory Comput. 8, 2247 (2012)");
  if(doInt_) log<<plumed.cite("Baftizadeh, Cossio, Pietrucci, and Laio, Curr. Phys. Chem. 2, 79 (2012)");
  if(acceleration_) log<<plumed.cite("Pratyush and Parrinello, Phys. Rev. Lett. 111, 230602 (2013)");
  if(calc_rct_) log<<plumed.cite("Pratyush and Parrinello, J. Phys. Chem. B, 119, 736 (2015)");
  if(concurrent || rect_biasf_.size()>0) log<<plumed.cite("Gil-Ley and Bussi, J. Chem. Theory Comput. 11, 1077 (2015)");
  if(rect_biasf_.size()>0 && walkers_mpi_) log<<plumed.cite("Hosek, Toulcova, Bortolato, and Spiwok, J. Phys. Chem. B 120, 2209 (2016)");
  if(targetfilename_.length()>0) {
    log<<plumed.cite("White, Dama, and Voth, J. Chem. Theory Comput. 11, 2451 (2015)");
    log<<plumed.cite("Marinelli and Faraldo-Gómez,  Biophys. J. 108, 2779 (2015)");
    log<<plumed.cite("Gil-Ley, Bottaro, and Bussi, J. Chem. Theory Comput. 12, 2790 (2016)");
  }
  if(freq_adaptive_) log<<plumed.cite("Wang, Valsson, Tiwary, Parrinello, and Lindorff-Larsen, J. Chem. Phys. 149, 072309 (2018)");
  log<<"\n";
}

void MetaD::readTemperingSpecs(TemperingSpecs &t_specs)
{
  // Set global tempering parameters.
  parse(t_specs.name_stem + "BIASFACTOR", t_specs.biasf);
  if (t_specs.biasf != -1.0) {
    if (kbt_ == 0.0) {
      error("Unless the MD engine passes the temperature to plumed, with tempered metad you must specify it using TEMP");
    }
    if (t_specs.biasf == 1.0) {
      error("A bias factor of 1 corresponds to zero delta T and zero hill size, so it is not allowed.");
    }
    t_specs.is_active = true;
    parse(t_specs.name_stem + "BIASTHRESHOLD", t_specs.threshold);
    if (t_specs.threshold < 0.0) {
      error(t_specs.name + " bias threshold is nonsensical");
    }
    parse(t_specs.name_stem + "ALPHA", t_specs.alpha);
    if (t_specs.alpha <= 0.0 || t_specs.alpha > 1.0) {
      error(t_specs.name + " decay shape parameter alpha is nonsensical");
    }
  }
}

void MetaD::logTemperingSpecs(const TemperingSpecs &t_specs)
{
  log.printf("  %s bias factor %f\n", t_specs.name.c_str(), t_specs.biasf);
  log.printf("  KbT %f\n", kbt_);
  if (t_specs.threshold != 0.0) log.printf("  %s bias threshold %f\n", t_specs.name.c_str(), t_specs.threshold);
  if (t_specs.alpha != 1.0) log.printf("  %s decay shape parameter alpha %f\n", t_specs.name.c_str(), t_specs.alpha);
}

void MetaD::readGaussians(IFile *ifile)
{
  unsigned ncv=getNumberOfArguments();
  std::vector<double> center(ncv);
  std::vector<double> sigma(ncv);
  double height;
  int nhills=0;
  bool multivariate=false;

  std::vector<Value> tmpvalues;
  for(unsigned j=0; j<getNumberOfArguments(); ++j) tmpvalues.push_back( Value( this, getPntrToArgument(j)->getName(), false ) );

  while(scanOneHill(ifile,tmpvalues,center,sigma,height,multivariate))
  {
    nhills++;
    // note that for gamma=1 we store directly -F
    if(welltemp_ && biasf_>1.0) height*=(biasf_-1.0)/biasf_;
    addGaussian(Gaussian(multivariate,height,center,sigma));
  }
  log.printf("      %d Gaussians read\n",nhills);
}

void MetaD::writeGaussian(const Gaussian& hill, OFile&file)
{
  unsigned ncv=getNumberOfArguments();
  file.printField("time",getTimeStep()*getStep());
  for(unsigned i=0; i<ncv; ++i) {
    file.printField(getPntrToArgument(i),hill.center[i]);
  }
  hillsOfile_.printField("kerneltype","stretched-gaussian");
  if(hill.multivariate) {
    hillsOfile_.printField("multivariate","true");
    Matrix<double> mymatrix(ncv,ncv);
    unsigned k=0;
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        // recompose the full inverse matrix
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k];
        k++;
      }
    }
    // invert it
    Matrix<double> invmatrix(ncv,ncv);
    Invert(mymatrix,invmatrix);
    // enforce symmetry
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        invmatrix(i,j)=invmatrix(j,i);
      }
    }

    // do cholesky so to have a "sigma like" number
    Matrix<double> lower(ncv,ncv);
    cholesky(invmatrix,lower);
    // loop in band form
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=0; j<ncv-i; j++) {
        file.printField("sigma_"+getPntrToArgument(j+i)->getName()+"_"+getPntrToArgument(j)->getName(),lower(j+i,j));
      }
    }
  } else {
    hillsOfile_.printField("multivariate","false");
    for(unsigned i=0; i<ncv; ++i)
      file.printField("sigma_"+getPntrToArgument(i)->getName(),hill.sigma[i]);
  }
  double height=hill.height;
  // note that for gamma=1 we store directly -F
  if(welltemp_ && biasf_>1.0) height*=biasf_/(biasf_-1.0);
  file.printField("height",height).printField("biasf",biasf_);
  if(mw_n_>1) file.printField("clock",int(std::time(0)));
  file.printField();
}

void MetaD::addGaussian(const Gaussian& hill)
{
  if(grid_) {
    size_t ncv=getNumberOfArguments();
    std::vector<unsigned> nneighb=getGaussianSupport(hill);
    std::vector<Grid::index_t> neighbors=BiasGrid_->getNeighbors(hill.center,nneighb);
    std::vector<double> der(ncv);
    std::vector<double> xx(ncv);
    if(comm.Get_size()==1) {
      // for performance reasons and thread safety
      std::vector<double> dp(ncv);
      for(size_t i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        for(unsigned j=0; j<ncv; ++j) der[j]=0.0;
        BiasGrid_->getPoint(ineigh,xx);
        double bias=evaluateGaussianAndDerivatives(xx,hill,der,dp);
        BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
      }
    } else {
      unsigned stride=comm.Get_size();
      unsigned rank=comm.Get_rank();
      std::vector<double> allder(ncv*neighbors.size(),0.0);
      std::vector<double> n_der(ncv,0.0);
      std::vector<double> allbias(neighbors.size(),0.0);
      // for performance reasons and thread safety
      std::vector<double> dp(ncv);
      for(unsigned i=rank; i<neighbors.size(); i+=stride) {
        Grid::index_t ineigh=neighbors[i];
        for(unsigned j=0; j<ncv; ++j) n_der[j]=0.0;
        BiasGrid_->getPoint(ineigh,xx);
        allbias[i]=evaluateGaussianAndDerivatives(xx,hill,n_der,dp);
        for(unsigned j=0; j<ncv; j++) allder[ncv*i+j]=n_der[j];
      }
      comm.Sum(allbias);
      comm.Sum(allder);
      for(unsigned i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        for(unsigned j=0; j<ncv; ++j) der[j]=allder[ncv*i+j];
        BiasGrid_->addValueAndDerivatives(ineigh,allbias[i],der);
      }
    }
  } else hills_.push_back(hill);
}

std::vector<unsigned> MetaD::getGaussianSupport(const Gaussian& hill)
{
  std::vector<unsigned> nneigh;
  std::vector<double> cutoff;
  unsigned ncv=getNumberOfArguments();

  // traditional or flexible hill?
  if(hill.multivariate) {
    unsigned k=0;
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        // recompose the full inverse matrix
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k];
        k++;
      }
    }
    // Reinvert so to have the ellipses
    Matrix<double> myinv(ncv,ncv);
    Invert(mymatrix,myinv);
    Matrix<double> myautovec(ncv,ncv);
    std::vector<double> myautoval(ncv); //should I take this or their square root?
    diagMat(myinv,myautoval,myautovec);
    double maxautoval=0.;
    unsigned ind_maxautoval; ind_maxautoval=ncv;
    for(unsigned i=0; i<ncv; i++) {
      if(myautoval[i]>maxautoval) {maxautoval=myautoval[i]; ind_maxautoval=i;}
    }
    for(unsigned i=0; i<ncv; i++) {
      cutoff.push_back(std::sqrt(2.0*dp2cutoff)*std::abs(std::sqrt(maxautoval)*myautovec(i,ind_maxautoval)));
    }
  } else {
    for(unsigned i=0; i<ncv; ++i) {
      cutoff.push_back(std::sqrt(2.0*dp2cutoff)*hill.sigma[i]);
    }
  }

  if(doInt_) {
    if(hill.center[0]+cutoff[0] > uppI_ || hill.center[0]-cutoff[0] < lowI_) {
      // in this case, we updated the entire grid to avoid problems
      return BiasGrid_->getNbin();
    } else {
      nneigh.push_back( static_cast<unsigned>(ceil(cutoff[0]/BiasGrid_->getDx()[0])) );
      return nneigh;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      nneigh.push_back( static_cast<unsigned>(ceil(cutoff[i]/BiasGrid_->getDx()[i])) );
    }
  }

  return nneigh;
}

double MetaD::getBias(const std::vector<double>& cv)
{
  double bias=0.0;
  if(grid_) bias = BiasGrid_->getValue(cv);
  else {
    unsigned nt=OpenMP::getNumThreads();
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();

    if(!nlist_) {
      #pragma omp parallel num_threads(nt)
      {
        #pragma omp for reduction(+:bias) nowait
        for(unsigned i=rank; i<hills_.size(); i+=stride) bias+=evaluateGaussian(cv,hills_[i]);
      }
    } else {
      #pragma omp parallel num_threads(nt)
      {
        #pragma omp for reduction(+:bias) nowait
        for(unsigned i=rank; i<nlist_hills_.size(); i+=stride) bias+=evaluateGaussian(cv,nlist_hills_[i]);
      }
    }
    comm.Sum(bias);
  }

  return bias;
}

double MetaD::getBiasAndDerivatives(const std::vector<double>& cv, std::vector<double>& der)
{
  unsigned ncv=getNumberOfArguments();
  double bias=0.0;
  if(grid_) {
    std::vector<double> vder(ncv);
    bias=BiasGrid_->getValueAndDerivatives(cv,vder);
    for(unsigned i=0; i<ncv; i++) der[i]=vder[i];
  } else {
    unsigned nt=OpenMP::getNumThreads();
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();

    if(!nlist_) {
      if(hills_.size()<2*nt*stride||nt==1) {
        // for performance reasons and thread safety
        std::vector<double> dp(ncv);
        for(unsigned i=rank; i<hills_.size(); i+=stride) {
          bias+=evaluateGaussianAndDerivatives(cv,hills_[i],der,dp);
        }
      } else {
        #pragma omp parallel num_threads(nt)
        {
          std::vector<double> omp_deriv(ncv,0.);
          // for performance reasons and thread safety
          std::vector<double> dp(ncv);
          #pragma omp for reduction(+:bias) nowait
          for(unsigned i=rank; i<hills_.size(); i+=stride) {
            bias+=evaluateGaussianAndDerivatives(cv,hills_[i],omp_deriv,dp);
          }
          #pragma omp critical
          for(unsigned i=0; i<ncv; i++) der[i]+=omp_deriv[i];
        }
      }
    } else {
      if(hills_.size()<2*nt*stride||nt==1) {
        // for performance reasons and thread safety
        std::vector<double> dp(ncv);
        for(unsigned i=rank; i<nlist_hills_.size(); i+=stride) {
          bias+=evaluateGaussianAndDerivatives(cv,nlist_hills_[i],der,dp);
        }
      } else {
        #pragma omp parallel num_threads(nt)
        {
          std::vector<double> omp_deriv(ncv,0.);
          // for performance reasons and thread safety
          std::vector<double> dp(ncv);
          #pragma omp for reduction(+:bias) nowait
          for(unsigned i=rank; i<nlist_hills_.size(); i+=stride) {
            bias+=evaluateGaussianAndDerivatives(cv,nlist_hills_[i],omp_deriv,dp);
          }
          #pragma omp critical
          for(unsigned i=0; i<ncv; i++) der[i]+=omp_deriv[i];
        }
      }
    }
    comm.Sum(bias);
    comm.Sum(der);
  }

  return bias;
}

double MetaD::getGaussianNormalization(const Gaussian& hill)
{
  double norm=1;
  unsigned ncv=hill.center.size();

  if(hill.multivariate) {
    // recompose the full sigma from the upper diag cholesky
    unsigned k=0;
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
      double ldet; logdet( mymatrix, ldet );
      norm = std::exp( ldet );  // Not sure here if mymatrix is sigma or inverse
    }
  } else {
    for(unsigned i=0; i<hill.sigma.size(); i++) norm*=hill.sigma[i];
  }

  return norm*std::pow(2*pi,static_cast<double>(ncv)/2.0);
}

double MetaD::evaluateGaussian(const std::vector<double>& cv, const Gaussian& hill)
{
  unsigned ncv=cv.size();

  // I use a pointer here because cv is const (and should be const)
  // but when using doInt it is easier to locally replace cv[0] with
  // the upper/lower limit in case it is out of range
  double tmpcv[1];
  const double *pcv=NULL; // pointer to cv
  if(ncv>0) pcv=&cv[0];
  if(doInt_) {
    plumed_assert(ncv==1);
    tmpcv[0]=cv[0];
    if(cv[0]<lowI_) tmpcv[0]=lowI_;
    if(cv[0]>uppI_) tmpcv[0]=uppI_;
    pcv=&(tmpcv[0]);
  }

  double dp2=0.0;
  if(hill.multivariate) {
    unsigned k=0;
    // recompose the full sigma from the upper diag cholesky
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    for(unsigned i=0; i<ncv; i++) {
      double dp_i=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_i*dp_i*mymatrix(i,j)*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_i*dp_j*mymatrix(i,j);
        }
      }
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      double dp=difference(i,hill.center[i],pcv[i])*hill.invsigma[i];
      dp2+=dp*dp;
    }
    dp2*=0.5;
  }

  double bias=0.0;
  if(dp2<dp2cutoff) bias=hill.height*(stretchA*std::exp(-dp2)+stretchB);

  return bias;
}

double MetaD::evaluateGaussianAndDerivatives(const std::vector<double>& cv, const Gaussian& hill, std::vector<double>& der, std::vector<double>& dp_)
{
  unsigned ncv=cv.size();

  // I use a pointer here because cv is const (and should be const)
  // but when using doInt it is easier to locally replace cv[0] with
  // the upper/lower limit in case it is out of range
  const double *pcv=NULL; // pointer to cv
  double tmpcv[1]; // tmp array with cv (to be used with doInt_)
  if(ncv>0) pcv=&cv[0];
  if(doInt_) {
    plumed_assert(ncv==1);
    tmpcv[0]=cv[0];
    if(cv[0]<lowI_) tmpcv[0]=lowI_;
    if(cv[0]>uppI_) tmpcv[0]=uppI_;
    pcv=&(tmpcv[0]);
  }

  bool int_der=false;
  if(doInt_) {
    if(cv[0]<lowI_ || cv[0]>uppI_) int_der=true;
  }

  double dp2=0.0;
  double bias=0.0;
  if(hill.multivariate) {
    unsigned k=0;
    // recompose the full sigma from the upper diag cholesky
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    for(unsigned i=0; i<ncv; i++) {
      dp_[i]=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_[i]*dp_[i]*mymatrix(i,j)*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_[i]*dp_j*mymatrix(i,j);
        }
      }
    }
    if(dp2<dp2cutoff) {
      bias=hill.height*std::exp(-dp2);
      if(!int_der) {
        for(unsigned i=0; i<ncv; i++) {
          double tmp=0.0;
          for(unsigned j=0; j<ncv; j++) tmp += dp_[j]*mymatrix(i,j)*bias;
          der[i]-=tmp*stretchA;
        }
      } else {
        for(unsigned i=0; i<ncv; i++) der[i]=0.;
      }
      bias=stretchA*bias+hill.height*stretchB;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      dp_[i]=difference(i,hill.center[i],pcv[i])*hill.invsigma[i];
      dp2+=dp_[i]*dp_[i];
    }
    dp2*=0.5;
    if(dp2<dp2cutoff) {
      bias=hill.height*std::exp(-dp2);
      if(!int_der) {
        for(unsigned i=0; i<ncv; i++) der[i]-=bias*dp_[i]*hill.invsigma[i]*stretchA;
      } else {
        for(unsigned i=0; i<ncv; i++) der[i]=0.;
      }
      bias=stretchA*bias+hill.height*stretchB;
    }
  }

  return bias;
}

double MetaD::getHeight(const std::vector<double>& cv)
{
  double height=height0_;
  if(welltemp_) {
    double vbias = getBias(cv);
    if(biasf_>1.0) {
      height = height0_*std::exp(-vbias/(kbt_*(biasf_-1.0)));
    } else {
      // notice that if gamma=1 we store directly -F
      height = height0_*std::exp(-vbias/kbt_);
    }
  }
  if(dampfactor_>0.0) {
    plumed_assert(BiasGrid_);
    double m=BiasGrid_->getMaxValue();
    height*=std::exp(-m/(kbt_*(dampfactor_)));
  }
  if (tt_specs_.is_active) {
    double vbarrier = transition_bias_;
    temperHeight(height, tt_specs_, vbarrier);
  }
  if(TargetGrid_) {
    double f=TargetGrid_->getValue(cv)-TargetGrid_->getMaxValue();
    height*=std::exp(f/kbt_);
  }
  return height;
}

void MetaD::temperHeight(double& height, const TemperingSpecs& t_specs, const double tempering_bias)
{
  if (t_specs.alpha == 1.0) {
    height *= std::exp(-std::max(0.0, tempering_bias - t_specs.threshold) / (kbt_ * (t_specs.biasf - 1.0)));
  } else {
    height *= std::pow(1 + (1 - t_specs.alpha) / t_specs.alpha * std::max(0.0, tempering_bias - t_specs.threshold) / (kbt_ * (t_specs.biasf - 1.0)), - t_specs.alpha / (1 - t_specs.alpha));
  }
}

void MetaD::calculate()
{
  // this is because presently there is no way to properly pass information
  // on adaptive hills (diff) after exchanges:
  if(adaptive_==FlexibleBin::diffusion && getExchangeStep()) error("ADAPTIVE=DIFF is not compatible with replica exchange");

  const unsigned ncv=getNumberOfArguments();
  std::vector<double> cv(ncv);
  for(unsigned i=0; i<ncv; ++i) cv[i]=getArgument(i);

  if(nlist_) {
    nlist_steps_++;
    if(getExchangeStep()) nlist_update_=true;
    else {
      for(unsigned i=0; i<ncv; ++i) {
        double d = difference(i, cv[i], nlist_center_[i]);
        double nk_dist2 = d*d/nlist_dev2_[i];
        if(nk_dist2>nlist_param_[1]) {nlist_update_=true; break;}
      }
    }
    if(nlist_update_) updateNlist();
  }

  double ene = 0.;
  std::vector<double> der(ncv,0.);
  if(biasf_!=1.0) ene = getBiasAndDerivatives(cv,der);
  setBias(ene);
  for(unsigned i=0; i<ncv; i++) setOutputForce(i,-der[i]);

  if(calc_work_) getPntrToComponent("work")->set(work_);
  if(calc_rct_) getPntrToComponent("rbias")->set(ene - reweight_factor_);
  // calculate the acceleration factor
  if(acceleration_&&!isFirstStep_) {
    acc_ += static_cast<double>(getStride()) * std::exp(ene/(kbt_));
    const double mean_acc = acc_/((double) getStep());
    getPntrToComponent("acc")->set(mean_acc);
  } else if (acceleration_ && isFirstStep_ && acc_restart_mean_ > 0.0) {
    acc_ = acc_restart_mean_ * static_cast<double>(getStep());
    if(freq_adaptive_) {
      // has to be done here if restarting, as the acc is not defined before
      updateFrequencyAdaptiveStride();
    }
  }
}

void MetaD::update()
{
  // adding hills criteria (could be more complex though)
  bool nowAddAHill;
  if(getStep()%current_stride_==0 && !isFirstStep_) nowAddAHill=true;
  else {
    nowAddAHill=false;
    isFirstStep_=false;
  }

  unsigned ncv=getNumberOfArguments();
  std::vector<double> cv(ncv);
  for(unsigned i=0; i<ncv; ++i) cv[i] = getArgument(i);

  double vbias=0.;
  if(calc_work_) vbias=getBias(cv);

  // if you use adaptive, call the FlexibleBin
  bool multivariate=false;
  if(adaptive_!=FlexibleBin::none) {
    flexbin_->update(nowAddAHill);
    multivariate=true;
  }

  std::vector<double> thissigma;
  if(nowAddAHill) {
    // add a Gaussian
    double height=getHeight(cv);
    // returns upper diagonal inverse
    if(adaptive_!=FlexibleBin::none) thissigma=flexbin_->getInverseMatrix();
    // returns normal sigma
    else thissigma=sigma0_;

    // In case we use walkers_mpi, it is now necessary to communicate with other replicas.
    if(walkers_mpi_) {
      // Allocate arrays to store all walkers hills
      std::vector<double> all_cv(mpi_nw_*ncv,0.0);
      std::vector<double> all_sigma(mpi_nw_*thissigma.size(),0.0);
      std::vector<double> all_height(mpi_nw_,0.0);
      std::vector<int>    all_multivariate(mpi_nw_,0);
      if(comm.Get_rank()==0) {
        // Communicate (only root)
        multi_sim_comm.Allgather(cv,all_cv);
        multi_sim_comm.Allgather(thissigma,all_sigma);
        // notice that if gamma=1 we store directly -F so this scaling is not necessary:
        multi_sim_comm.Allgather(height*(biasf_>1.0?biasf_/(biasf_-1.0):1.0),all_height);
        multi_sim_comm.Allgather(int(multivariate),all_multivariate);
      }
      // Share info with group members
      comm.Bcast(all_cv,0);
      comm.Bcast(all_sigma,0);
      comm.Bcast(all_height,0);
      comm.Bcast(all_multivariate,0);

      // Flying Gaussian
      if (flying_) {
        hills_.clear();
        comm.Barrier();
      }

      for(unsigned i=0; i<mpi_nw_; i++) {
        // actually add hills one by one
        std::vector<double> cv_now(ncv);
        std::vector<double> sigma_now(thissigma.size());
        for(unsigned j=0; j<ncv; j++) cv_now[j]=all_cv[i*ncv+j];
        for(unsigned j=0; j<thissigma.size(); j++) sigma_now[j]=all_sigma[i*thissigma.size()+j];
        // notice that if gamma=1 we store directly -F so this scaling is not necessary:
        double fact=(biasf_>1.0?(biasf_-1.0)/biasf_:1.0);
        Gaussian newhill=Gaussian(all_multivariate[i],all_height[i]*fact,cv_now,sigma_now);
        addGaussian(newhill);
        if(!flying_) writeGaussian(newhill,hillsOfile_);
      }
    } else {
      Gaussian newhill=Gaussian(multivariate,height,cv,thissigma);
      addGaussian(newhill);
      writeGaussian(newhill,hillsOfile_);
    }

    // this is to update the hills neighbor list
    if(nlist_) nlist_update_=true;
  }

  // this should be outside of the if block in case
  // mw_rstride_ is not a multiple of stride_
  if(mw_n_>1 && getStep()%mw_rstride_==0) hillsOfile_.flush();

  if(calc_work_) {
    if(nlist_) updateNlist();
    double vbias1=getBias(cv);
    work_+=vbias1-vbias;
  }

  // dump grid on file
  if(wgridstride_>0&&(getStep()%wgridstride_==0||getCPT())) {
    // in case old grids are stored, a sequence of grids should appear
    // this call results in a repetition of the header:
    if(storeOldGrids_) gridfile_.clearFields();
    // in case only latest grid is stored, file should be rewound
    // this will overwrite previously written grids
    else {
      int r = 0;
      if(walkers_mpi_) {
        if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
        comm.Bcast(r,0);
      }
      if(r==0) gridfile_.rewind();
    }
    BiasGrid_->writeToFile(gridfile_);
    // if a single grid is stored, it is necessary to flush it, otherwise
    // the file might stay empty forever (when a single grid is not large enough to
    // trigger flushing from the operating system).
    // on the other hand, if grids are stored one after the other this is
    // no necessary, and we leave the flushing control to the user as usual
    // (with FLUSH keyword)
    if(!storeOldGrids_) gridfile_.flush();
  }

  // if multiple walkers and time to read Gaussians
  if(mw_n_>1 && getStep()%mw_rstride_==0) {
    for(int i=0; i<mw_n_; ++i) {
      // don't read your own Gaussians
      if(i==mw_id_) continue;
      // if the file is not open yet
      if(!(ifiles_[i]->isOpen())) {
        // check if it exists now and open it!
        if(ifiles_[i]->FileExist(ifilesnames_[i])) {
          ifiles_[i]->open(ifilesnames_[i]);
          ifiles_[i]->reset(false);
        }
        // otherwise read the new Gaussians
      } else {
        log.printf("  Reading hills from %s:",ifilesnames_[i].c_str());
        readGaussians(ifiles_[i].get());
        ifiles_[i]->reset(false);
      }
    }
    // this is to update the hills neighbor list
    if(nlist_) nlist_update_=true;
  }

  // Recalculate special bias quantities whenever the bias has been changed by the update.
  bool bias_has_changed = (nowAddAHill || (mw_n_ > 1 && getStep() % mw_rstride_ == 0));
  if (calc_rct_ && bias_has_changed && getStep()%(stride_*rct_ustride_)==0) computeReweightingFactor();
  if (calc_max_bias_ && bias_has_changed) {
    max_bias_ = BiasGrid_->getMaxValue();
    getPntrToComponent("maxbias")->set(max_bias_);
  }
  if (calc_transition_bias_ && bias_has_changed) {
    transition_bias_ = getTransitionBarrierBias();
    getPntrToComponent("transbias")->set(transition_bias_);
  }

  // Frequency adaptive metadynamics - update hill addition frequency
  if(freq_adaptive_ && getStep()%fa_update_frequency_==0) {
    updateFrequencyAdaptiveStride();
  }
}

/// takes a pointer to the file and a template std::string with values v and gives back the next center, sigma and height
bool MetaD::scanOneHill(IFile* ifile, std::vector<Value>& tmpvalues, std::vector<double>& center, std::vector<double>& sigma, double& height, bool& multivariate)
{
  double dummy;
  multivariate=false;
  if(ifile->scanField("time",dummy)) {
    unsigned ncv=tmpvalues.size();
    for(unsigned i=0; i<ncv; ++i) {
      ifile->scanField( &tmpvalues[i] );
      if( tmpvalues[i].isPeriodic() && ! getPntrToArgument(i)->isPeriodic() ) {
        error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
      } else if( tmpvalues[i].isPeriodic() ) {
        std::string imin, imax; tmpvalues[i].getDomain( imin, imax );
        std::string rmin, rmax; getPntrToArgument(i)->getDomain( rmin, rmax );
        if( imin!=rmin || imax!=rmax ) {
          error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
        }
      }
      center[i]=tmpvalues[i].get();
    }
    // scan for kerneltype
    std::string ktype="stretched-gaussian";
    if( ifile->FieldExist("kerneltype") ) ifile->scanField("kerneltype",ktype);
    if( ktype=="gaussian" ) {
      noStretchWarning();
    } else if( ktype!="stretched-gaussian") {
      error("non Gaussian kernels are not supported in MetaD");
    }
    // scan for multivariate label: record the actual file position so to eventually rewind
    std::string sss;
    ifile->scanField("multivariate",sss);
    if(sss=="true") multivariate=true;
    else if(sss=="false") multivariate=false;
    else plumed_merror("cannot parse multivariate = "+ sss);
    if(multivariate) {
      sigma.resize(ncv*(ncv+1)/2);
      Matrix<double> upper(ncv,ncv);
      Matrix<double> lower(ncv,ncv);
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=0; j<ncv-i; j++) {
          ifile->scanField("sigma_"+getPntrToArgument(j+i)->getName()+"_"+getPntrToArgument(j)->getName(),lower(j+i,j));
          upper(j,j+i)=lower(j+i,j);
        }
      }
      Matrix<double> mymult(ncv,ncv);
      Matrix<double> invmatrix(ncv,ncv);
      mult(lower,upper,mymult);
      // now invert and get the sigmas
      Invert(mymult,invmatrix);
      // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form)
      unsigned k=0;
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) {
          sigma[k]=invmatrix(i,j);
          k++;
        }
      }
    } else {
      for(unsigned i=0; i<ncv; ++i) {
        ifile->scanField("sigma_"+getPntrToArgument(i)->getName(),sigma[i]);
      }
    }

    ifile->scanField("height",height);
    ifile->scanField("biasf",dummy);
    if(ifile->FieldExist("clock")) ifile->scanField("clock",dummy);
    if(ifile->FieldExist("lower_int")) ifile->scanField("lower_int",dummy);
    if(ifile->FieldExist("upper_int")) ifile->scanField("upper_int",dummy);
    ifile->scanField();
    return true;
  } else {
    return false;
  }
}

void MetaD::computeReweightingFactor()
{
  if(biasf_==1.0) { // in this case we have no bias, so reweight factor is 0.0
    getPntrToComponent("rct")->set(0.0);
    return;
  }

  double Z_0=0; //proportional to the integral of exp(-beta*F)
  double Z_V=0; //proportional to the integral of exp(-beta*(F+V))
  double minusBetaF=biasf_/(biasf_-1.)/kbt_;
  double minusBetaFplusV=1./(biasf_-1.)/kbt_;
  if (biasf_==-1.0) { //non well-tempered case
    minusBetaF=1./kbt_;
    minusBetaFplusV=0;
  }
  max_bias_=BiasGrid_->getMaxValue(); //to avoid exp overflow

  const unsigned rank=comm.Get_rank();
  const unsigned stride=comm.Get_size();
  for (Grid::index_t t=rank; t<BiasGrid_->getSize(); t+=stride) {
    const double val=BiasGrid_->getValue(t);
    Z_0+=std::exp(minusBetaF*(val-max_bias_));
    Z_V+=std::exp(minusBetaFplusV*(val-max_bias_));
  }
  comm.Sum(Z_0);
  comm.Sum(Z_V);

  reweight_factor_=kbt_*std::log(Z_0/Z_V)+max_bias_;
  getPntrToComponent("rct")->set(reweight_factor_);
}

double MetaD::getTransitionBarrierBias()
{
  // If there is only one well of interest, return the bias at that well point.
  if (transitionwells_.size() == 1) {
    double tb_bias = getBias(transitionwells_[0]);
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
    double least_transition_bias;
    std::vector<double> sink = transitionwells_[0];
    std::vector<double> source = transitionwells_[1];
    least_transition_bias = BiasGrid_->findMaximalPathMinimum(source, sink);
    for (unsigned i = 2; i < transitionwells_.size(); i++) {
      if (least_transition_bias == 0.0) {
        break;
      }
      source = transitionwells_[i];
      double curr_transition_bias = BiasGrid_->findMaximalPathMinimum(source, sink);
      least_transition_bias = fmin(curr_transition_bias, least_transition_bias);
    }
    return least_transition_bias;
  }
}

void MetaD::updateFrequencyAdaptiveStride()
{
  plumed_massert(freq_adaptive_,"should only be used if frequency adaptive metadynamics is enabled");
  plumed_massert(acceleration_,"frequency adaptive metadynamics can only be used if the acceleration factor is calculated");
  const double mean_acc = acc_/((double) getStep());
  int tmp_stride= stride_*floor((mean_acc/fa_min_acceleration_)+0.5);
  if(mean_acc >= fa_min_acceleration_) {
    if(tmp_stride > current_stride_) {current_stride_ = tmp_stride;}
  }
  if(fa_max_stride_!=0 && current_stride_>fa_max_stride_) {
    current_stride_=fa_max_stride_;
  }
  getPntrToComponent("pace")->set(current_stride_);
}

bool MetaD::checkNeedsGradients()const
{
  if(adaptive_==FlexibleBin::geometry) {
    if(getStep()%stride_==0 && !isFirstStep_) return true;
    else return false;
  } else return false;
}

void MetaD::updateNlist()
{
  // no need to check for neighbors
  if(hills_.size()==0) return;

  // here we generate the neighbor list
  nlist_hills_.clear();
  std::vector<Gaussian> local_flat_nl;
  unsigned nt=OpenMP::getNumThreads();
  if(hills_.size()<2*nt) nt=1;
  #pragma omp parallel num_threads(nt)
  {
    std::vector<Gaussian> private_flat_nl;
    #pragma omp for nowait
    for(unsigned k=0; k<hills_.size(); k++)
    {
      double dist2=0;
      for(unsigned i=0; i<getNumberOfArguments(); i++)
      {
        const double d=difference(i,getArgument(i),hills_[k].center[i])/hills_[k].sigma[i];
        dist2+=d*d;
      }
      if(dist2<=nlist_param_[0]*dp2cutoff) private_flat_nl.push_back(hills_[k]);
    }
    #pragma omp critical
    local_flat_nl.insert(local_flat_nl.end(), private_flat_nl.begin(), private_flat_nl.end());
  }
  nlist_hills_ = local_flat_nl;

  // here we set some properties that are used to decide when to update it again
  for(unsigned i=0; i<getNumberOfArguments(); i++) nlist_center_[i]=getArgument(i);
  std::vector<double> dev2;
  dev2.resize(getNumberOfArguments(),0);
  for(unsigned k=0; k<nlist_hills_.size(); k++)
  {
    for(unsigned i=0; i<getNumberOfArguments(); i++)
    {
      const double d=difference(i,getArgument(i),nlist_hills_[k].center[i]);
      dev2[i]+=d*d;
    }
  }
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    if(dev2[i]>0.) nlist_dev2_[i]=dev2[i]/static_cast<double>(nlist_hills_.size());
    else nlist_dev2_[i]=hills_.back().sigma[i]*hills_.back().sigma[i];
  }

  // we are done
  getPntrToComponent("nlker")->set(nlist_hills_.size());
  getPntrToComponent("nlsteps")->set(nlist_steps_);
  nlist_steps_=0;
  nlist_update_=false;
}

}
}
