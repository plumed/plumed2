/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "tools/Grid.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Exception.h"
#include "core/FlexibleBin.h"
#include "tools/Matrix.h"
#include "tools/Random.h"
#include <string>
#include <cstring>
#include "tools/File.h"
#include <iostream>
#include <limits>
#include <ctime>

#define DP2CUTOFF 6.25

using namespace std;


namespace PLMD{
namespace bias{

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
the values of a larger and larger number of Gaussians. To avoid this issue you can
store the bias on a grid.  This approach is similar to that proposed in \cite babi08jcp but has the 
advantage that the grid spacing is independent on the Gaussian width.
Notice that you should
provide either the number of bins for every collective variable (GRID_BIN) or
the desired grid spacing (GRID_SPACING). In case you provide both PLUMED will use
the most conservative choice (highest number of bins) for each dimension.
In case you do not provide any information about bin size (neither GRID_BIN nor GRID_SPACING)
and if Gaussian width is fixed PLUMED will use 1/5 of the Gaussian width as grid spacing.
This default choice should be reasonable for most applications.

Another option that is available in plumed is well-tempered metadynamics \cite Barducci:2008. In this
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


The c(t) reweighting factor can also be calculated on the fly using the equations 
presented in \cite Tiwary_jp504920s. 
The expression used to calculate c(t) follows directly from using Eq. 12 in 
Eq. 3 in \cite Tiwary_jp504920s and gives smoother results than equivalent Eqs. 13 
and Eqs. 14 in that paper. The c(t) is given by the rct component while the bias 
normalized by c(t) is given by the rbias component (rbias=bias-ct) which can be used 
to obtain a reweighted histogram.
The calculation of c(t) is enabled by using the keyword REWEIGHTING_NGRID where the grid used for the 
calculation is specified.   This grid should have a size that is equal or larger than the grid given in GRID_BIN./
By default c(t) is updated every 50 Gaussian hills but this 
can be changed by using the REWEIGHTING_NHILLS keyword. 
This option can only be employed together with Well-Tempered Metadynamics and requires that 
a grid is used.

Additional material and examples can be also found in the tutorials: 

- \ref belfast-6
- \ref belfast-7
- \ref belfast-8

Notice that at variance with PLUMED 1.3 it is now straightforward to apply concurrent metadynamics
as done e.g. in Ref. \cite gil2015enhanced . This indeed can be obtained by using the METAD
action multiple times in the same input file.

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
The c(t) reweighting factor can be calculated on the fly using the equations
presented in \cite Tiwary_jp504920s as described above.
This is enabled by using the keyword REWEIGHTING_NGRID where the grid used for
the calculation is set. The number of grid points given in REWEIGHTING_NGRID
should be equal or larger than the number of grid points given in GRID_BIN.
\verbatim
METAD ...
 LABEL=metad
 ARG=phi,psi SIGMA=0.20,0.20 HEIGHT=1.20 BIASFACTOR=5 TEMP=300.0 PACE=500
 GRID_MIN=-pi,-pi GRID_MAX=pi,pi GRID_BIN=150,150
 REWEIGHTING_NGRID=150,150
 REWEIGHTING_NHILLS=20
... METAD
\endverbatim
Here we have asked that the calculation is performed every 20 hills by using
REWEIGHTING_NHILLS keyword. If this keyword is not given the calculation will
by default be performed every 50 hills. The c(t) reweighting factor will be given
in the rct component while the instantaneous value of the bias potential
normalized using the c(t) reweighting factor is given in the rbias component
[rbias=bias-c(t)] which can be used to obtain a reweighted histogram or
free energy surface using the \ref HISTOGRAM analysis.

\par
The kinetics of the transitions between basins can also be analysed on the fly as
in \cite PRL230602. The flag ACCELERATION turn on accumulation of the acceleration
factor that can then be used to determine the rate. This method can be used together
with \ref COMMITTOR analysis to stop the simulation when the system get to the target basin.
It must be used together with Well-Tempered Metadynamics.


*/
//+ENDPLUMEDOC

class MetaD : public Bias{

private:
  struct Gaussian {
   vector<double> center;
   vector<double> sigma;
   double height;
   bool   multivariate; // this is required to discriminate the one dimensional case 
   vector<double> invsigma;
   Gaussian(const vector<double> & center,const vector<double> & sigma,double height, bool multivariate ):
     center(center),sigma(sigma),height(height),multivariate(multivariate),invsigma(sigma){
       // to avoid troubles from zero element in flexible hills
       for(unsigned i=0;i<invsigma.size();++i)abs(invsigma[i])>1.e-20?invsigma[i]=1.0/invsigma[i]:0.; 
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
  std::string gridfilename_,gridreadfilename_;
  int wgridstride_; 
  bool grid_,hasextgrid_;
  double height0_;
  double biasf_;
  double kbt_;
  int stride_;
  bool welltemp_;
  double* dp_;
  int adaptive_;
  FlexibleBin *flexbin;
  int mw_n_;
  string mw_dir_;
  int mw_id_;
  int mw_rstride_;
  bool walkers_mpi;
  bool acceleration;
  double acc;
  vector<IFile*> ifiles;
  vector<string> ifilesnames;
  double uppI_;
  double lowI_;
  bool doInt_;
  bool isFirstStep;
  double reweight_factor;
  std::vector<unsigned> rewf_grid_; 
  unsigned int rewf_ustride_;
/// accumulator for work
  double work_;
  long int last_step_warn_grid;
  
 
  void   readGaussians(IFile*);
  bool   readChunkOfGaussians(IFile *ifile, unsigned n);
  void   writeGaussian(const Gaussian&,OFile&);
  void   addGaussian(const Gaussian&);
  double getHeight(const vector<double>&);
  double getBiasAndDerivatives(const vector<double>&,double* der=NULL);
  double evaluateGaussian(const vector<double>&, const Gaussian&,double* der=NULL);
  void   finiteDifferenceGaussian(const vector<double>&, const Gaussian&);
  double getGaussianNormalization( const Gaussian& );
  vector<unsigned> getGaussianSupport(const Gaussian&);
  bool   scanOneHill(IFile *ifile,  vector<Value> &v, vector<double> &center, vector<double>  &sigma, double &height, bool &multivariate  );
  void   computeReweightingFactor();
  std::string fmt;

public:
  explicit MetaD(const ActionOptions&);
  ~MetaD();
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
  bool checkNeedsGradients()const{if(adaptive_==FlexibleBin::geometry){return true;}else{return false;}}
};

PLUMED_REGISTER_ACTION(MetaD,"METAD")

void MetaD::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("rbias","REWEIGHTING_NGRID","the instantaneous value of the bias normalized using the \\f$c(t)\\f$ reweighting factor [rbias=bias-c(t)]. This component can be used to obtain a reweighted histogram.");
  keys.addOutputComponent("rct","REWEIGHTING_NGRID","the reweighting factor \\f$c(t)\\f$.");
  keys.addOutputComponent("work","default","accumulator for work");
  keys.addOutputComponent("acc","ACCELERATION","the metadynamics acceleration factor");
  keys.use("ARG");
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","FILE","HILLS","a file in which the list of added hills is stored");
  keys.add("optional","HEIGHT","the heights of the Gaussian hills. Compulsory unless TAU, TEMP and BIASFACTOR are given");
  keys.add("optional","FMT","specify format for HILLS files (useful for decrease the number of digits in regtests)");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","TAU","in well tempered metadynamics, sets height to (kb*DeltaT*pace*timestep)/tau");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("optional","REWEIGHTING_NGRID","calculate the c(t) reweighting factor and use that to obtain the normalized bias [rbias=bias-c(t)]. Here you should specify the number of grid points required in each dimension. The number of grid points should be equal or larger to the number of grid points given in GRID_BIN." "This method is not compatible with metadynamics not on a grid.");
  keys.add("optional","REWEIGHTING_NHILLS","how many Gaussian hills should be deposited between calculating the c(t) reweighting factor. The default is to do this every 50 hills.");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","GRID_WSTRIDE","write the grid to a file every N steps");
  keys.add("optional","GRID_WFILE","the file on which to write the grid");
  keys.addFlag("STORE_GRIDS",false,"store all the grid files the calculation generates. They will be deleted if this keyword is not present");
  keys.add("optional","ADAPTIVE","use a geometric (=GEOM) or diffusion (=DIFF) based hills width scheme. Sigma is one number that has distance units or timestep dimensions");
  keys.add("optional","WALKERS_ID", "walker id");
  keys.add("optional","WALKERS_N", "number of walkers");
  keys.add("optional","WALKERS_DIR", "shared directory with the hills files from all the walkers");
  keys.add("optional","WALKERS_RSTRIDE","stride for reading hills files");
  keys.add("optional","INTERVAL","monodimensional lower and upper limits, outside the limits the system will not feel the biasing force.");
  keys.add("optional","GRID_RFILE","a grid file from which the bias should be read at the initial step of the simulation");
  keys.add("optional","SIGMA_MAX","the upper bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.add("optional","SIGMA_MIN","the lower bounds for the sigmas (in CV units) when using adaptive hills. Negative number means no bounds ");
  keys.addFlag("WALKERS_MPI",false,"Switch on MPI version of multiple walkers - not compatible with other WALKERS_* options");
  keys.addFlag("ACCELERATION",false,"Set to TRUE if you want to compute the metadynamics acceleration factor.");  
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

MetaD::~MetaD(){
  if(flexbin) delete flexbin;
  if(BiasGrid_) delete BiasGrid_;
  hillsOfile_.close();
  if(wgridstride_>0) gridfile_.close();
  delete [] dp_;
  // close files
  for(int i=0;i<mw_n_;++i){
   if(ifiles[i]->isOpen()) ifiles[i]->close();
   delete ifiles[i];
  }
}

MetaD::MetaD(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
// Grid stuff initialization
BiasGrid_(NULL),ExtGrid_(NULL), wgridstride_(0), grid_(false), hasextgrid_(false),
// Metadynamics basic parameters
height0_(std::numeric_limits<double>::max()), biasf_(1.0), kbt_(0.0),
stride_(0), welltemp_(false),
// Other stuff
dp_(NULL), adaptive_(FlexibleBin::none),
flexbin(NULL),
// Multiple walkers initialization
mw_n_(1), mw_dir_("./"), mw_id_(0), mw_rstride_(1),
walkers_mpi(false),
acceleration(false), acc(0.0),
// Interval initialization
uppI_(-1), lowI_(-1), doInt_(false),
isFirstStep(true),
reweight_factor(0.0),
rewf_ustride_(1),
work_(0),
last_step_warn_grid(0)
{
  // parse the flexible hills
  string adaptiveoption;
  adaptiveoption="NONE";
  parse("ADAPTIVE",adaptiveoption);
  if (adaptiveoption=="GEOM"){
		  log.printf("  Uses Geometry-based hills width: sigma must be in distance units and only one sigma is needed\n");
		  adaptive_=FlexibleBin::geometry;	
  }else if (adaptiveoption=="DIFF"){
		  log.printf("  Uses Diffusion-based hills width: sigma must be in timesteps and only one sigma is needed\n");
		  adaptive_=FlexibleBin::diffusion;	
  }else if (adaptiveoption=="NONE"){
		  adaptive_=FlexibleBin::none;	
  }else{
		  error("I do not know this type of adaptive scheme");	
  }
  // parse the sigma
  parseVector("SIGMA",sigma0_);

  parse("FMT",fmt);

  if (adaptive_==FlexibleBin::none){
         // if you use normal sigma you need one sigma per argument 
         if( sigma0_.size()!=getNumberOfArguments() ) error("number of arguments does not match number of SIGMA parameters");
  }else{
         // if you use flexible hills you need one sigma  
         if(sigma0_.size()!=1){
        	 error("If you choose ADAPTIVE you need only one sigma according to your choice of type (GEOM/DIFF)");
         } 
	 // if adaptive then the number must be an integer
 	 if(adaptive_==FlexibleBin::diffusion){
		if(int(sigma0_[0])-sigma0_[0]>1.e-9 || int(sigma0_[0])-sigma0_[0] <-1.e-9 || int(sigma0_[0])<1 ){
		 	error("In case of adaptive hills with diffusion, the sigma must be an integer which is the number of timesteps\n");	
		} 
	 } 
	 // here evtl parse the sigma min and max values
	 
	 parseVector("SIGMA_MIN",sigma0min_);
	 if(sigma0min_.size()>0 && sigma0min_.size()<getNumberOfArguments()) {
           error("the number of SIGMA_MIN values be at least the number of the arguments");
	 } else if(sigma0min_.size()==0) { 
           sigma0min_.resize(getNumberOfArguments());
           for(unsigned i=0;i<getNumberOfArguments();i++){sigma0min_[i]=-1.;}	
         } 

	 parseVector("SIGMA_MAX",sigma0max_);
	 if(sigma0max_.size()>0 && sigma0max_.size()<getNumberOfArguments()) {
           error("the number of SIGMA_MAX values be at least the number of the arguments"); 
         } else if(sigma0max_.size()==0) { 
           sigma0max_.resize(getNumberOfArguments());
           for(unsigned i=0;i<getNumberOfArguments();i++){sigma0max_[i]=-1.;}	
         } 

         flexbin=new FlexibleBin(adaptive_,this,sigma0_[0],sigma0min_,sigma0max_);
  }
  // note: HEIGHT is not compulsory, since one could use the TAU keyword, see below
  parse("HEIGHT",height0_);
  parse("PACE",stride_);
  if(stride_<=0 ) error("frequency for hill addition is nonsensical");
  string hillsfname="HILLS";
  parse("FILE",hillsfname);
  parse("BIASFACTOR",biasf_);
  if( biasf_<1.0 ) error("well tempered bias factor is nonsensical");
  double temp=0.0;
  parse("TEMP",temp);
  if(temp>0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
  else kbt_=plumed.getAtoms().getKbT();
  if(biasf_>1.0){
    if(kbt_==0.0) error("Unless the MD engine passes the temperature to plumed, with well-tempered metad you must specify it using TEMP");
    welltemp_=true;
  }
  double tau=0.0;
  parse("TAU",tau);
  if(tau==0.0){
    if(height0_==std::numeric_limits<double>::max()) error("At least one between HEIGHT and TAU should be specified");
    // if tau is not set, we compute it here from the other input parameters
    if(welltemp_) tau=(kbt_*(biasf_-1.0))/height0_*getTimeStep()*stride_;
  } else {
    if(!welltemp_)error("TAU only makes sense in well-tempered metadynamics");
    if(height0_!=std::numeric_limits<double>::max()) error("At most one between HEIGHT and TAU should be specified");
    height0_=(kbt_*(biasf_-1.0))/tau*getTimeStep()*stride_;
  }
  
  // Grid Stuff
  vector<std::string> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  if(gmin.size()!=getNumberOfArguments() && gmin.size()!=0) error("not enough values for GRID_MIN");
  vector<std::string> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  if(gmax.size()!=getNumberOfArguments() && gmax.size()!=0) error("not enough values for GRID_MAX");
  vector<unsigned> gbin(getNumberOfArguments());
  vector<double>   gspacing;
  parseVector("GRID_BIN",gbin);
  if(gbin.size()!=getNumberOfArguments() && gbin.size()!=0) error("not enough values for GRID_BIN");
  parseVector("GRID_SPACING",gspacing);
  if(gspacing.size()!=getNumberOfArguments() && gspacing.size()!=0) error("not enough values for GRID_SPACING");
  if(gmin.size()!=gmax.size()) error("GRID_MAX and GRID_MIN should be either present or absent");
  if(gspacing.size()!=0 && gmin.size()==0) error("If GRID_SPACING is present also GRID_MIN should be present");
  if(gbin.size()!=0     && gmin.size()==0) error("If GRID_SPACING is present also GRID_MIN should be present");
  if(gmin.size()!=0){
    if(gbin.size()==0 && gspacing.size()==0){
      if(adaptive_==FlexibleBin::none){
        log<<"  Binsize not specified, 1/5 of sigma will be be used\n";
        plumed_assert(sigma0_.size()==getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for(unsigned i=0;i<gspacing.size();i++) gspacing[i]=0.2*sigma0_[i];
      } else {
        // with adaptive hills and grid a sigma min must be specified
        if(sigma0min_.size()==0) error("When using Adaptive Gaussians on a grid SIGMA_MIN must be specified");
        log<<"  Binsize not specified, 1/5 of sigma_min will be be used\n";
        plumed_assert(sigma0_.size()==getNumberOfArguments());
        gspacing.resize(getNumberOfArguments());
        for(unsigned i=0;i<gspacing.size();i++) gspacing[i]=0.2*sigma0min_[i];
        //error("At least one among GRID_BIN and GRID_SPACING should be used");
      }
    } else if(gspacing.size()!=0 && gbin.size()==0){
      log<<"  The number of bins will be estimated from GRID_SPACING\n";
    } else if(gspacing.size()!=0 && gbin.size()!=0){
      log<<"  You specified both GRID_BIN and GRID_SPACING\n";
      log<<"  The more conservative (highest) number of bins will be used for each variable\n";
    }
    if(gbin.size()==0) gbin.assign(getNumberOfArguments(),1);
    if(gspacing.size()!=0) for(unsigned i=0;i<getNumberOfArguments();i++){
      double a,b;
      Tools::convert(gmin[i],a);
      Tools::convert(gmax[i],b);
      unsigned n=((b-a)/gspacing[i])+1;
      if(gbin[i]<n) gbin[i]=n;
    }
  }
  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  if(gbin.size()>0){grid_=true;}
  parse("GRID_WSTRIDE",wgridstride_);
  parse("GRID_WFILE",gridfilename_); 
  parseFlag("STORE_GRIDS",storeOldGrids_);
  if(grid_ && gridfilename_.length()>0){
    if(wgridstride_==0 ) error("frequency with which to output grid not specified use GRID_WSTRIDE");
  }

  if(grid_ && wgridstride_>0){
    if(gridfilename_.length()==0) error("grid filename not specified use GRID_WFILE"); 
  }

  parse("GRID_RFILE",gridreadfilename_);

  if(grid_){ 
     parseVector("REWEIGHTING_NGRID",rewf_grid_); 
     if( rewf_grid_.size()>0 && rewf_grid_.size()!=getNumberOfArguments() ){
         error("size mismatch for REWEIGHTING_NGRID keyword");
     } else if( rewf_grid_.size()==getNumberOfArguments() ){
         for(unsigned j=0;j<getNumberOfArguments();++j){
            if( !getPntrToArgument(j)->isPeriodic() ) rewf_grid_[j] += 1; 
         }
     }
     if( adaptive_==FlexibleBin::diffusion || adaptive_==FlexibleBin::geometry) warning("reweighting has not been proven to work with adaptive Gaussians");
     rewf_ustride_=50; parse("REWEIGHTING_NHILLS",rewf_ustride_);
  }

  // Multiple walkers
  parse("WALKERS_N",mw_n_);
  parse("WALKERS_ID",mw_id_);
  if(mw_n_<=mw_id_) error("walker ID should be a numerical value less than the total number of walkers");
  parse("WALKERS_DIR",mw_dir_);
  parse("WALKERS_RSTRIDE",mw_rstride_);

  // MPI version
  parseFlag("WALKERS_MPI",walkers_mpi);

  // Inteval keyword
  vector<double> tmpI(2);
  parseVector("INTERVAL",tmpI);
  if(tmpI.size()!=2&&tmpI.size()!=0) error("both a lower and an upper limits must be provided with INTERVAL");
  else if(tmpI.size()==2) {
    lowI_=tmpI.at(0);
    uppI_=tmpI.at(1);
    if(getNumberOfArguments()!=1) error("INTERVAL limits correction works only for monodimensional metadynamics!");
    if(uppI_<lowI_) error("The Upper limit must be greater than the Lower limit!");
    doInt_=true;
  }

  acceleration=false;
  parseFlag("ACCELERATION",acceleration);

  checkRead();

  log.printf("  Gaussian width ");
  if (adaptive_==FlexibleBin::diffusion)log.printf(" (Note: The units of sigma are in timesteps) ");
  if (adaptive_==FlexibleBin::geometry)log.printf(" (Note: The units of sigma are in dist units) ");
  for(unsigned i=0;i<sigma0_.size();++i) log.printf(" %f",sigma0_[i]);
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_); 
  log.printf("  Gaussian file %s\n",hillsfname.c_str());
  if(welltemp_){
    log.printf("  Well-Tempered Bias Factor %f\n",biasf_);
    log.printf("  Hills relaxation time (tau) %f\n",tau);
    log.printf("  KbT %f\n",kbt_);
  }
  if(doInt_) log.printf("  Upper and Lower limits boundaries for the bias are activated at %f - %f\n", lowI_, uppI_);
  if(grid_){
   log.printf("  Grid min");
   for(unsigned i=0;i<gmin.size();++i) log.printf(" %s",gmin[i].c_str() );
   log.printf("\n");
   log.printf("  Grid max");
   for(unsigned i=0;i<gmax.size();++i) log.printf(" %s",gmax[i].c_str() );
   log.printf("\n");
   log.printf("  Grid bin");
   for(unsigned i=0;i<gbin.size();++i) log.printf(" %u",gbin[i]);
   log.printf("\n");
   if(spline){log.printf("  Grid uses spline interpolation\n");}
   if(sparsegrid){log.printf("  Grid uses sparse grid\n");}
   if(wgridstride_>0){log.printf("  Grid is written on file %s with stride %d\n",gridfilename_.c_str(),wgridstride_);} 
  }
  if(gridreadfilename_.length()>0){
	   log.printf("  Reading an additional bias from grid in file %s \n",gridreadfilename_.c_str());
  }


  if(mw_n_>1){
   if(walkers_mpi) error("MPI version of multiple walkers is not compatible with filesystem version of multiple walkers");
   log.printf("  %d multiple walkers active\n",mw_n_);
   log.printf("  walker id %d\n",mw_id_);
   log.printf("  reading stride %d\n",mw_rstride_);
   log.printf("  directory with hills files %s\n",mw_dir_.c_str());
  } else {
   if(walkers_mpi) log.printf("  Multiple walkers active using MPI communnication\n"); 
  }

  addComponent("bias"); componentIsNotPeriodic("bias");
  if( rewf_grid_.size()>0 ){ 
   addComponent("rbias"); componentIsNotPeriodic("rbias");
   addComponent("rct"); componentIsNotPeriodic("rct"); 
  }
  addComponent("work"); componentIsNotPeriodic("work");

  if(acceleration) {
    if(!welltemp_) error("The calculation of the acceleration works only if Well-Tempered Metadynamics is on"); 
    log.printf("  calculation on the fly of the acceleration factor");
    addComponent("acc"); componentIsNotPeriodic("acc");
  }

// for performance
  dp_ = new double[getNumberOfArguments()];

// initializing and checking grid
  if(grid_){
   // check for adaptive and sigma_min
   if(sigma0min_.size()==0&&adaptive_!=FlexibleBin::none) error("When using Adaptive Gaussians on a grid SIGMA_MIN must be specified");
   // check for mesh and sigma size
   for(unsigned i=0;i<getNumberOfArguments();i++) {
     double a,b;
     Tools::convert(gmin[i],a);
     Tools::convert(gmax[i],b);
     double mesh=(b-a)/((double)gbin[i]);
     if(mesh>0.5*sigma0_[i]) log<<"  WARNING: Using a METAD with a Grid Spacing larger than half of the Gaussians width can produce artifacts\n";
   }
   std::string funcl=getLabel() + ".bias";
   if(!sparsegrid){BiasGrid_=new Grid(funcl,getArguments(),gmin,gmax,gbin,spline,true);}
   else{BiasGrid_=new SparseGrid(funcl,getArguments(),gmin,gmax,gbin,spline,true);}
   std::vector<std::string> actualmin=BiasGrid_->getMin();
   std::vector<std::string> actualmax=BiasGrid_->getMax();
   for(unsigned i=0;i<getNumberOfArguments();i++){
     if(gmin[i]!=actualmin[i]) log<<"  WARNING: GRID_MIN["<<i<<"] has been adjusted to "<<actualmin[i]<<" to fit periodicity\n";
     if(gmax[i]!=actualmax[i]) log<<"  WARNING: GRID_MAX["<<i<<"] has been adjusted to "<<actualmax[i]<<" to fit periodicity\n";
   }
  }

  if(wgridstride_>0){
    gridfile_.link(*this);
    gridfile_.open(gridfilename_);
  }

// initializing external grid
  if(gridreadfilename_.length()>0){
   hasextgrid_=true;
   // read the grid in input, find the keys
   IFile gridfile; gridfile.open(gridreadfilename_);
   std::string funcl=getLabel() + ".bias";
   ExtGrid_=Grid::create(funcl,getArguments(),gridfile,false,false,true);
   gridfile.close();
   if(ExtGrid_->getDimension()!=getNumberOfArguments()) error("mismatch between dimensionality of input grid and number of arguments");
   for(unsigned i=0;i<getNumberOfArguments();++i){
     if( getPntrToArgument(i)->isPeriodic()!=ExtGrid_->getIsPeriodic()[i] ) error("periodicity mismatch between arguments and input bias");
   }
  }

// Tiwary-Parrinello reweighting factor
  if(rewf_grid_.size()>0){
   log.printf("  the c(t) reweighting factor will be calculated every %u hills\n",rewf_ustride_);
   getPntrToComponent("rct")->set(reweight_factor);
  }

// creating vector of ifile* for hills reading 
// open all files at the beginning and read Gaussians if restarting
  for(int i=0;i<mw_n_;++i){
   string fname;
   if(mw_n_>1) {
    stringstream out; out << i;
    fname = mw_dir_+"/"+hillsfname+"."+out.str();
   } else {
    fname = hillsfname;
   }
   IFile *ifile = new IFile();
   ifile->link(*this);
   ifiles.push_back(ifile);                                                             
   ifilesnames.push_back(fname);
   if(ifile->FileExist(fname)){
    ifile->open(fname);
    if(getRestart()){
     log.printf("  Restarting from %s:",ifilesnames[i].c_str());                  
     readGaussians(ifiles[i]);                                                    
    }
    ifiles[i]->reset(false);
    // close only the walker own hills file for later writing
    if(i==mw_id_) ifiles[i]->close();
   }
  }

  comm.Barrier();

// this barrier is needed when using walkers_mpi
// to be sure that all files have been read before
// backing them up
// it should not be used when walkers_mpi is false otherwise
// it would introduce troubles when using replicas without METAD
// (e.g. in bias exchange with a neutral replica)
// see issue #168 on github
  if(comm.Get_rank()==0 && walkers_mpi) multi_sim_comm.Barrier();

// Calculate the Tiwary-Parrinello reweighting factor if we are restarting from previous hills
  if(plumed.getRestart() && rewf_grid_.size()>0 ){computeReweightingFactor();}

// open hills file for writing
  hillsOfile_.link(*this);
  if(walkers_mpi){
    int r=0;
    if(comm.Get_rank()==0) r=multi_sim_comm.Get_rank();
    comm.Bcast(r,0);
    if(r>0) ifilesnames[mw_id_]="/dev/null";
    hillsOfile_.enforceSuffix("");
  }
  hillsOfile_.open(ifilesnames[mw_id_]);
  if(fmt.length()>0) hillsOfile_.fmtField(fmt);
  hillsOfile_.addConstantField("multivariate");
  if(doInt_) {
    hillsOfile_.addConstantField("lower_int").printField("lower_int",lowI_);
    hillsOfile_.addConstantField("upper_int").printField("upper_int",uppI_);
  }
  hillsOfile_.setHeavyFlush();
// output periodicities of variables
  for(unsigned i=0;i<getNumberOfArguments();++i) hillsOfile_.setupPrintValue( getPntrToArgument(i) );

  bool concurrent=false;

  const ActionSet&actionSet(plumed.getActionSet());
  for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p) if(dynamic_cast<MetaD*>(*p)){ concurrent=true; break; }
  if(concurrent) log<<"  You are using concurrent metadynamics\n";

  log<<"  Bibliography "<<plumed.cite("Laio and Parrinello, PNAS 99, 12562 (2002)");
  if(welltemp_) log<<plumed.cite(
    "Barducci, Bussi, and Parrinello, Phys. Rev. Lett. 100, 020603 (2008)");
  if(mw_n_>1||walkers_mpi) log<<plumed.cite(
    "Raiteri, Laio, Gervasio, Micheletti, and Parrinello, J. Phys. Chem. B 110, 3533 (2006)");
  if(adaptive_!=FlexibleBin::none) log<<plumed.cite(
    "Branduardi, Bussi, and Parrinello, J. Chem. Theory Comput. 8, 2247 (2012)");
  if(doInt_) log<<plumed.cite(
     "Baftizadeh, Cossio, Pietrucci, and Laio, Curr. Phys. Chem. 2, 79 (2012)");
  if(acceleration) log<<plumed.cite(
     "Pratyush and Parrinello, Phys. Rev. Lett. 111, 230602 (2013)");
  if(rewf_grid_.size()>0) log<<plumed.cite(
     "Pratyush and Parrinello, J. Phys. Chem. B, 119, 736 (2015)");
  if(concurrent) log<<plumed.cite(
     "Gil-Ley and Bussi, J. Chem. Theory Comput. 11, 1077 (2015)");
  log<<"\n";

}

void MetaD::readGaussians(IFile *ifile)
{
 unsigned ncv=getNumberOfArguments();
 vector<double> center(ncv);
 vector<double> sigma(ncv);
 double height;
 int nhills=0; 
 bool multivariate=false;

 std::vector<Value> tmpvalues;
 for(unsigned j=0;j<getNumberOfArguments();++j) tmpvalues.push_back( Value( this, getPntrToArgument(j)->getName(), false ) ); 

 while(scanOneHill(ifile,tmpvalues,center,sigma,height,multivariate)){;
  nhills++;
  if(welltemp_){height*=(biasf_-1.0)/biasf_;}
  addGaussian(Gaussian(center,sigma,height,multivariate));
 }     
 log.printf("      %d Gaussians read\n",nhills);
}

bool MetaD::readChunkOfGaussians(IFile *ifile, unsigned n)
{
 unsigned ncv=getNumberOfArguments();
 vector<double> center(ncv);
 vector<double> sigma(ncv);
 double height;
 unsigned nhills=0;
 bool multivariate=false;
 std::vector<Value> tmpvalues;
 for(unsigned j=0;j<getNumberOfArguments();++j) tmpvalues.push_back( Value( this, getPntrToArgument(j)->getName(), false ) ); 

 while(scanOneHill(ifile,tmpvalues,center,sigma,height,multivariate)){;
  if(welltemp_){height*=(biasf_-1.0)/biasf_;}
  addGaussian(Gaussian(center,sigma,height,multivariate));
  if(nhills==n){
      log.printf("      %u Gaussians read\n",nhills);
      return true;
  }
  nhills++;
 }     
 log.printf("      %u Gaussians read\n",nhills);
 return false;
}

void MetaD::writeGaussian(const Gaussian& hill, OFile&file){
  unsigned ncv=getNumberOfArguments();
  file.printField("time",getTimeStep()*getStep());
  for(unsigned i=0;i<ncv;++i){
    file.printField(getPntrToArgument(i),hill.center[i]);
  }
  if(hill.multivariate){
    hillsOfile_.printField("multivariate","true");
         Matrix<double> mymatrix(ncv,ncv);
         unsigned k=0;
         for(unsigned i=0;i<ncv;i++){
                for(unsigned j=i;j<ncv;j++){
                        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
                        k++;
                }
         }
         // invert it 
         Matrix<double> invmatrix(ncv,ncv);
         Invert(mymatrix,invmatrix);
         // enforce symmetry
         for(unsigned i=0;i<ncv;i++){
                for(unsigned j=i;j<ncv;j++){
                        invmatrix(i,j)=invmatrix(j,i);
                }
         }

         // do cholesky so to have a "sigma like" number
         Matrix<double> lower(ncv,ncv);
         cholesky(invmatrix,lower); // now this , in band form , is similar to the sigmas
         // loop in band form 
         for (unsigned i=0;i<ncv;i++){
              for (unsigned j=0;j<ncv-i;j++){
                      file.printField("sigma_"+getPntrToArgument(j+i)->getName()+"_"+getPntrToArgument(j)->getName(),lower(j+i,j));
              }
         }
  } else {
    hillsOfile_.printField("multivariate","false");
    for(unsigned i=0;i<ncv;++i)
      file.printField("sigma_"+getPntrToArgument(i)->getName(),hill.sigma[i]);
  }
  double height=hill.height;
  if(welltemp_){height*=biasf_/(biasf_-1.0);}
  file.printField("height",height).printField("biasf",biasf_);
  if(mw_n_>1) file.printField("clock",int(std::time(0)));
  file.printField();
}

void MetaD::addGaussian(const Gaussian& hill){

 if(!grid_){hills_.push_back(hill);} 
 else{
  unsigned ncv=getNumberOfArguments();
  vector<unsigned> nneighb=getGaussianSupport(hill);
  vector<Grid::index_t> neighbors=BiasGrid_->getNeighbors(hill.center,nneighb);
  vector<double> der(ncv);
  vector<double> xx(ncv);
  if(comm.Get_size()==1){
    for(unsigned i=0;i<neighbors.size();++i){
     Grid::index_t ineigh=neighbors[i];
     for(unsigned j=0;j<ncv;++j){der[j]=0.0;}
     BiasGrid_->getPoint(ineigh,xx);
     double bias=evaluateGaussian(xx,hill,&der[0]);
     BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
    } 
  } else {
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    vector<double> allder(ncv*neighbors.size(),0.0);
    vector<double> allbias(neighbors.size(),0.0);
    for(unsigned i=rank;i<neighbors.size();i+=stride){
     Grid::index_t ineigh=neighbors[i];
     BiasGrid_->getPoint(ineigh,xx);
     allbias[i]=evaluateGaussian(xx,hill,&allder[ncv*i]);
    }
    comm.Sum(allbias);
    comm.Sum(allder);
    for(unsigned i=0;i<neighbors.size();++i){
     Grid::index_t ineigh=neighbors[i];
     for(unsigned j=0;j<ncv;++j){der[j]=allder[ncv*i+j];}
     BiasGrid_->addValueAndDerivatives(ineigh,allbias[i],der);
    }
  }
 }
}

vector<unsigned> MetaD::getGaussianSupport(const Gaussian& hill)
{
 vector<unsigned> nneigh;
 if(doInt_){
   double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[0];
   if(hill.center[0]+cutoff > uppI_ || hill.center[0]-cutoff < lowI_) { 
     // in this case, we updated the entire grid to avoid problems
     return BiasGrid_->getNbin();
   } else {
     nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrid_->getDx()[0])) );
     return nneigh;
   }
 }

 // traditional or flexible hill? 
 if(hill.multivariate){
	unsigned ncv=getNumberOfArguments();
	unsigned k=0;
	//log<<"------- GET GAUSSIAN SUPPORT --------\n"; 
	Matrix<double> mymatrix(ncv,ncv);
	for(unsigned i=0;i<ncv;i++){
		for(unsigned j=i;j<ncv;j++){
			mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
			k++;
		}
	}
        // Reinvert so to have the ellipses 
	Matrix<double> myinv(ncv,ncv);
	Invert(mymatrix,myinv);
	Matrix<double> myautovec(ncv,ncv);
	vector<double> myautoval(ncv); //should I take this or their square root? 
	diagMat(myinv,myautoval,myautovec);
	double maxautoval;maxautoval=0.;
        unsigned ind_maxautoval;ind_maxautoval=ncv; 
	for (unsigned i=0;i<ncv;i++){
		if(myautoval[i]>maxautoval){maxautoval=myautoval[i];ind_maxautoval=i;}
        }  
	for (unsigned i=0;i<ncv;i++){
		double cutoff=sqrt(2.0*DP2CUTOFF)*abs(sqrt(maxautoval)*myautovec(i,ind_maxautoval));
		//log<<"AUTOVAL "<<myautoval[0]<<" COMP "<<abs(myautoval[0]*myautovec(i,0)) <<" CUTOFF "<<cutoff<<"\n";
	  	nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrid_->getDx()[i])) );
        }
 }else{
	 for(unsigned i=0;i<getNumberOfArguments();++i){
	  double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[i];
	  nneigh.push_back( static_cast<unsigned>(ceil(cutoff/BiasGrid_->getDx()[i])) );
 	}
 }
	//log<<"------- END GET GAUSSIAN SUPPORT --------\n"; 
 return nneigh;
}

double MetaD::getBiasAndDerivatives(const vector<double>& cv, double* der)
{
 double bias=0.0;
 if(!grid_){
  if(hills_.size()>10000 && (getStep()-last_step_warn_grid)>10000){
    std::string msg;
    Tools::convert(hills_.size(),msg);
    msg="You have accumulated "+msg+" hills, you should enable GRIDs to avoid serious performance hits";
    warning(msg);
    last_step_warn_grid=getStep();
  }
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  for(unsigned i=rank;i<hills_.size();i+=stride){
   bias+=evaluateGaussian(cv,hills_[i],der);
   //finite difference test 
   //finiteDifferenceGaussian(cv,hills_[i]);
  }
  comm.Sum(bias);
  if(der) comm.Sum(der,getNumberOfArguments());
 }else{
  if(der){
   vector<double> vder(getNumberOfArguments());
   bias=BiasGrid_->getValueAndDerivatives(cv,vder);
   for(unsigned i=0;i<getNumberOfArguments();++i) {der[i]=vder[i];}
  }else{
   bias=BiasGrid_->getValue(cv);
  }
 }
 if(hasextgrid_){
	  if(der){
	   vector<double> vder(getNumberOfArguments());
	   bias+=ExtGrid_->getValueAndDerivatives(cv,vder);
	   for(unsigned i=0;i<getNumberOfArguments();++i) {der[i]+=vder[i];}
	  }else{
	   bias+=ExtGrid_->getValue(cv);
	  }
 }
 return bias;
}

double MetaD::getGaussianNormalization( const Gaussian& hill ){
  double norm=1; unsigned ncv=hill.center.size();

  if(hill.multivariate){
  // recompose the full sigma from the upper diag cholesky 
    unsigned k=0; 
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0;i<ncv;i++){
        for(unsigned j=i;j<ncv;j++){
            mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
            k++;
        }
        double ldet; logdet( mymatrix, ldet );
        norm = exp( ldet );  // Not sure here if mymatrix is sigma or inverse
    }
  } else {
    for(unsigned i=0;i<hill.sigma.size();++i) norm*=hill.sigma[i]; 
  }

  return norm*pow(2*pi,static_cast<double>(ncv)/2.0);
}

double MetaD::evaluateGaussian
 (const vector<double>& cv, const Gaussian& hill, double* der)
{
 double dp2=0.0;
 double bias=0.0;
// I use a pointer here because cv is const (and should be const)
// but when using doInt it is easier to locally replace cv[0] with
// the upper/lower limit in case it is out of range
 const double *pcv=NULL; // pointer to cv
 double tmpcv[1]; // tmp array with cv (to be used with doInt_)
 if(cv.size()>0) pcv=&cv[0];
 if(doInt_){
   plumed_assert(cv.size()==1);
   pcv=&(tmpcv[0]);
   tmpcv[0]=cv[0];
   if(cv[0]<lowI_) tmpcv[0]=lowI_;
   if(cv[0]>uppI_) tmpcv[0]=uppI_;
 }
 if(hill.multivariate){ 
    unsigned k=0;
    unsigned ncv=cv.size(); 
    // recompose the full sigma from the upper diag cholesky 
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0;i<ncv;i++){
		for(unsigned j=i;j<ncv;j++){
			mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
			k++;
		}
    }

    for(unsigned i=0;i<cv.size();++i){
	double dp_i=difference(i,hill.center[i],pcv[i]);
        dp_[i]=dp_i;
    	for(unsigned j=i;j<cv.size();++j){
                  if(i==j){ 
              	  	 dp2+=dp_i*dp_i*mymatrix(i,j)*0.5; 
                  }else{ 
   		 	 double dp_j=difference(j,hill.center[j],pcv[j]);
              	  	 dp2+=dp_i*dp_j*mymatrix(i,j) ;
                  }
        }
    } 
    if(dp2<DP2CUTOFF){
       bias=hill.height*exp(-dp2);
       if(der){
        for(unsigned i=0;i<cv.size();++i){
                double tmp=0.0;
                k=i;
                for(unsigned j=0;j<cv.size();++j){
                                tmp+=   dp_[j]*mymatrix(i,j)*bias;
                        }
                        der[i]-=tmp;
                }   
       }
    }
 }else{
    for(unsigned i=0;i<cv.size();++i){
     double dp=difference(i,hill.center[i],pcv[i])*hill.invsigma[i];
     dp2+=dp*dp;
     dp_[i]=dp;
    }
    dp2*=0.5;
    if(dp2<DP2CUTOFF){
       bias=hill.height*exp(-dp2);
       if(der){
        for(unsigned i=0;i<cv.size();++i){der[i]+=-bias*dp_[i]*hill.invsigma[i];}
       }
    }
 }
 if(doInt_){
   if((cv[0]<lowI_ || cv[0]>uppI_) && der ) for(unsigned i=0;i<cv.size();++i)der[i]=0;
  }
 return bias;
}

double MetaD::getHeight(const vector<double>& cv)
{
 double height=height0_;
 if(welltemp_){
    double vbias=getBiasAndDerivatives(cv);
    height=height0_*exp(-vbias/(kbt_*(biasf_-1.0)));
 } 
 return height;
}

void MetaD::calculate()
{

// this is because presently there is no way to properly pass information
// on adaptive hills (diff) after exchanges:
  if(adaptive_==FlexibleBin::diffusion && getExchangeStep()) error("ADAPTIVE=DIFF is not compatible with replica exchange");

  unsigned ncv=getNumberOfArguments();
  vector<double> cv(ncv);
  for(unsigned i=0;i<ncv;++i){cv[i]=getArgument(i);}

  double* der=new double[ncv];
  for(unsigned i=0;i<ncv;++i){der[i]=0.0;}
  double ene=getBiasAndDerivatives(cv,der);
  getPntrToComponent("bias")->set(ene);
  if( rewf_grid_.size()>0 ) getPntrToComponent("rbias")->set(ene - reweight_factor);
// calculate the acceleration factor
  if(acceleration&&!isFirstStep) {
    acc += exp(ene/(kbt_));
    double mean_acc = acc/((double) getStep());
    getPntrToComponent("acc")->set(mean_acc);
  }
  getPntrToComponent("work")->set(work_);
// set Forces 
  for(unsigned i=0;i<ncv;++i){
   const double f=-der[i];
   setOutputForce(i,f);
  }

  delete [] der;
}

void MetaD::update(){
  vector<double> cv(getNumberOfArguments());
  vector<double> thissigma;
  bool multivariate;

  // adding hills criteria (could be more complex though)
  bool nowAddAHill; 
  if(getStep()%stride_==0 && !isFirstStep ){nowAddAHill=true;}else{nowAddAHill=false;isFirstStep=false;}

  for(unsigned i=0;i<cv.size();++i){cv[i]=getArgument(i);}

  double vbias=getBiasAndDerivatives(cv);

  // if you use adaptive, call the FlexibleBin 
  if (adaptive_!=FlexibleBin::none){
       flexbin->update(nowAddAHill);
       multivariate=true;
  }else{
       multivariate=false;
  };

  if(nowAddAHill){ // probably this can be substituted with a signal
   // add a Gaussian
   double height=getHeight(cv);
   // use normal sigma or matrix form? 
   if (adaptive_!=FlexibleBin::none){	
	thissigma=flexbin->getInverseMatrix(); // returns upper diagonal inverse 
	//cerr<<"ADDING HILLS "<<endl;	
   }else{
	thissigma=sigma0_;    // returns normal sigma
   }
// In case we use walkers_mpi, it is now necessary to communicate with other replicas.
   if(walkers_mpi){
     int nw=0;
     int mw=0;
     if(comm.Get_rank()==0){
// Only root of group can communicate with other walkers
       nw=multi_sim_comm.Get_size();
       mw=multi_sim_comm.Get_rank();
     }
// Communicate to the other members of the same group
// info abount number of walkers and walker index
     comm.Bcast(nw,0);
     comm.Bcast(mw,0);
// Allocate arrays to store all walkers hills
     std::vector<double> all_cv(nw*cv.size(),0.0);
     std::vector<double> all_sigma(nw*thissigma.size(),0.0);
     std::vector<double> all_height(nw,0.0);
     std::vector<int>    all_multivariate(nw,0);
     if(comm.Get_rank()==0){
// Communicate (only root)
       multi_sim_comm.Allgather(cv,all_cv);
       multi_sim_comm.Allgather(thissigma,all_sigma);
       multi_sim_comm.Allgather(height,all_height);
       multi_sim_comm.Allgather(int(multivariate),all_multivariate);
     }
// Share info with group members
     comm.Bcast(all_cv,0);
     comm.Bcast(all_sigma,0);
     comm.Bcast(all_height,0);
     comm.Bcast(all_multivariate,0);
     for(int i=0;i<nw;i++){
// actually add hills one by one
       std::vector<double> cv_now(cv.size());
       std::vector<double> sigma_now(thissigma.size());
       for(unsigned j=0;j<cv.size();j++) cv_now[j]=all_cv[i*cv.size()+j];
       for(unsigned j=0;j<thissigma.size();j++) sigma_now[j]=all_sigma[i*thissigma.size()+j];
       Gaussian newhill=Gaussian(cv_now,sigma_now,all_height[i],all_multivariate[i]);
       addGaussian(newhill);
       writeGaussian(newhill,hillsOfile_);
     }
   } else {
     Gaussian newhill=Gaussian(cv,thissigma,height,multivariate);
     addGaussian(newhill);
// print on HILLS file
     writeGaussian(newhill,hillsOfile_);
   }
  }

  double vbias1=getBiasAndDerivatives(cv);
  work_+=vbias1-vbias;

  // dump grid on file
  if(wgridstride_>0&&getStep()%wgridstride_==0){
    // in case old grids are stored, a sequence of grids should appear
    // this call results in a repetition of the header:
    if(storeOldGrids_) gridfile_.clearFields();
    // in case only latest grid is stored, file should be rewound
    // this will overwrite previously written grids
    else gridfile_.rewind();
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
 if(mw_n_>1 && getStep()%mw_rstride_==0){
   for(int i=0;i<mw_n_;++i){
    // don't read your own Gaussians
    if(i==mw_id_) continue;
    // if the file is not open yet 
    if(!(ifiles[i]->isOpen())){
     // check if it exists now and open it!
     if(ifiles[i]->FileExist(ifilesnames[i])) {
       ifiles[i]->open(ifilesnames[i]);
       ifiles[i]->reset(false);
     }
    // otherwise read the new Gaussians 
    } else {
     log.printf("  Reading hills from %s:",ifilesnames[i].c_str());
     readGaussians(ifiles[i]);
     ifiles[i]->reset(false);
    }
   }
 } 

 if(getStep()%(stride_*rewf_ustride_)==0 && nowAddAHill && rewf_grid_.size()>0 ){computeReweightingFactor();}
}

void MetaD::finiteDifferenceGaussian
 (const vector<double>& cv, const Gaussian& hill)
{
 log<<"--------- finiteDifferenceGaussian: size "<<cv.size() <<"------------\n";
 // for each cv
 // first get the bias and the derivative
 vector<double> oldder(cv.size()); 
 vector<double> der(cv.size()); 
 vector<double> mycv(cv.size()); 
 mycv=cv; 
 double step=1.e-6;
 Random random; 
 // just displace a tiny bit
 for(unsigned i=0;i<cv.size();i++)log<<"CV "<<i<<" V "<<mycv[i]<<"\n";
 for(unsigned i=0;i<cv.size();i++)mycv[i]+=1.e-2*2*(random.RandU01()-0.5);
 for(unsigned i=0;i<cv.size();i++)log<<"NENEWWCV "<<i<<" V "<<mycv[i]<<"\n";
 double oldbias=evaluateGaussian(mycv,hill,&oldder[0]);
 for (unsigned i=0;i<mycv.size();i++){
               double delta=step*2*(random.RandU01()-0.5);
               mycv[i]+=delta;
               double newbias=evaluateGaussian(mycv,hill,&der[0]);		
               log<<"CV "<<i;
               log<<" ANAL "<<oldder[i]<<" NUM "<<(newbias-oldbias)/delta<<" DIFF "<<(oldder[i]-(newbias-oldbias)/delta)<<"\n";
               mycv[i]-=delta;
 }
 log<<"--------- END finiteDifferenceGaussian ------------\n";
}

/// takes a pointer to the file and a template string with values v and gives back the next center, sigma and height 
bool MetaD::scanOneHill(IFile *ifile,  vector<Value> &tmpvalues, vector<double> &center, vector<double>  &sigma, double &height , bool &multivariate  ){
  double dummy;
  multivariate=false;
  if(ifile->scanField("time",dummy)){
     unsigned ncv; ncv=tmpvalues.size();
     for(unsigned i=0;i<ncv;++i){
       ifile->scanField( &tmpvalues[i] );
       if( tmpvalues[i].isPeriodic() && ! getPntrToArgument(i)->isPeriodic() ){
          error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
       } else if( tmpvalues[i].isPeriodic() ){
          std::string imin, imax; tmpvalues[i].getDomain( imin, imax );
          std::string rmin, rmax; getPntrToArgument(i)->getDomain( rmin, rmax );
          if( imin!=rmin || imax!=rmax ){
            error("in hills file periodicity for variable " + tmpvalues[i].getName() + " does not match periodicity in input");
          }
       }
       center[i]=tmpvalues[i].get();
     }
     // scan for multivariate label: record the actual file position so to eventually rewind 
     std::string sss;
     ifile->scanField("multivariate",sss);
     if(sss=="true") multivariate=true;
     else if(sss=="false") multivariate=false;
     else plumed_merror("cannot parse multivariate = "+ sss);
     if(multivariate){
        sigma.resize(ncv*(ncv+1)/2);
        Matrix<double> upper(ncv,ncv);
        Matrix<double> lower(ncv,ncv);
   	for (unsigned i=0;i<ncv;i++){
                 for (unsigned j=0;j<ncv-i;j++){
                         ifile->scanField("sigma_"+getPntrToArgument(j+i)->getName()+"_"+getPntrToArgument(j)->getName(),lower(j+i,j));
                         upper(j,j+i)=lower(j+i,j);
                 }
        }
        Matrix<double> mymult(ncv,ncv);       
        Matrix<double> invmatrix(ncv,ncv);       
	//log<<"Lower \n";
        //matrixOut(log,lower); 
	//log<<"Upper \n";
        //matrixOut(log,upper); 
        mult(lower,upper,mymult);          
	//log<<"Mult \n";
        //matrixOut(log,mymult); 
        // now invert and get the sigmas
        Invert(mymult,invmatrix);
	//log<<"Invert \n";
        //matrixOut(log,invmatrix); 
        // put the sigmas in the usual order: upper diagonal (this time in normal form and not in band form) 
        unsigned k=0;
   	for (unsigned i=0;i<ncv;i++){
   	       for (unsigned j=i;j<ncv;j++){
   	       	sigma[k]=invmatrix(i,j);
   	       	k++;
   	       }
   	}
     }else{
     	for(unsigned i=0;i<ncv;++i){
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
  }else{ 
    return false; 
  }; 
}

void MetaD::computeReweightingFactor(){
 if( !welltemp_ ) error("cannot compute the c(t) reweighting factors for non well-tempered metadynamics");

 // Recover the minimum values for the grid
 unsigned ncv=getNumberOfArguments();
 double grid_size=1.0; unsigned ntotgrid=1;
 std::vector<double> dmin( ncv ),dmax( ncv ), grid_spacing( ncv ), vals( ncv ); 
 for(unsigned j=0;j<ncv;++j){
    Tools::convert( BiasGrid_->getMin()[j], dmin[j] );
    Tools::convert( BiasGrid_->getMax()[j], dmax[j] );
    grid_spacing[j] = ( dmax[j] - dmin[j] ) / static_cast<double>( rewf_grid_[j] );
    grid_size *= grid_spacing[j]; 
    if( !getPntrToArgument(j)->isPeriodic() ) dmax[j] += grid_spacing[j]; 
    ntotgrid *= rewf_grid_[j];
 }
       
 // Now sum over whole grid
 reweight_factor=0.0; double* der=new double[ncv]; std::vector<unsigned> t_index( ncv );
 double sum1=0.0; double sum2=0.0;
 double afactor = biasf_ / (kbt_*(biasf_-1.0)); double afactor2 = 1.0 / (kbt_*(biasf_-1.0));
 unsigned rank=comm.Get_rank(), stride=comm.Get_size(); 
 for(unsigned i=rank;i<ntotgrid;i+=stride){
     t_index[0]=(i%rewf_grid_[0]);
     unsigned kk=i;
     for(unsigned j=1;j<ncv-1;++j){ kk=(kk-t_index[j-1])/rewf_grid_[i-1]; t_index[j]=(kk%rewf_grid_[i]); }
     if( ncv>=2 ) t_index[ncv-1]=((kk-t_index[ncv-1])/rewf_grid_[ncv-2]);
      
     for(unsigned j=0;j<ncv;++j) vals[j]=dmin[j] + t_index[j]*grid_spacing[j]; 

     double currentb=getBiasAndDerivatives(vals,der);
     sum1 += exp( afactor*currentb );
     sum2 += exp( afactor2*currentb );
 }
 delete [] der;
 comm.Sum( sum1 ); comm.Sum( sum2 );
 reweight_factor = kbt_ * std::log( sum1/sum2 );
 getPntrToComponent("rct")->set(reweight_factor);
}

}
}
