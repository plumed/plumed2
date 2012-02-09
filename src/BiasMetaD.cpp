#include "Bias.h"
#include "ActionRegister.h"
#include "Grid.h"
#include "PlumedMain.h"
#include "Atoms.h"
#include "PlumedException.h"

#define DP2CUTOFF 6.25

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS METAD 
/**
MetaDynamics on one or more variables

In a metadynamics simulations a history dependent bias composed of 
intermittently added Gaussian functions is added to the potential \cite metad.
This potential forces the system away from the kinetic traps in the potential energy surface
and out into the unexplored parts of the energy landscape.

\par Syntax
\verbatim
METAD ...
  ARG=x1,x2,... 
  SIGMA=s1,s2,... 
  HEIGHT=w 
  PACE=s
  [FILE=filename]
  [RESTART]
  [BIASFACTOR=biasf]
  [TEMP=temp]
  [GRID_MIN=min1,min2]
  [GRID_MAX=max1,max2]
  [GRID_BIN=bin1,bin2]
  [GRID_NOSPLINE]
  [GRID_SPARSE]
  [GRID_WFILE=filename]
  [GRID_WSTRIDE=ws]
... METAD
\endverbatim
SIGMA specifies an array of Gaussian widths, one for each variable,
HEIGHT the Gaussian height, PACE the Gaussian deposition stride in steps,
FILE the name of the file where the Gaussians are written to (or read from), 
RESTART to restart the run, BIASFACTOR the bias factor of well-tempered metad, 
TEMP the temperature. To store the bias potential on a grid,
you need to specify the grid boundaries with GRID_MIN and GRID_MAX and the number
of bins with GRID_BIN. An experimental sparse grid can be activated with GRID_SPARSE.
The use of spline can be disabled with GRID_NOSPLINE. You can dump the grid on file
with GRID_WFILE every GRID_WSTRIDE steps.
\par Examples
The following input is for a standard metadynamics calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the metadynamics bias potential are written to the COLVAR file every 100 steps.
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
METAD ARG=d1,d2 SIGMA=0.2,0.2 HEIGHT=0.3 STRIDE=500 LABEL=restraint
PRINT ARG=d1,d2,restraint.bias STRIDE=100  FILE=COLVAR
\endverbatim
(See also \ref DISTANCE \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasMetaD : public Bias{

private:
  struct Gaussian {
   vector<double> center;
   vector<double> sigma;
   double height;
   vector<double> invsigma;
   Gaussian(const vector<double> & center,const vector<double> & sigma,double height):
     center(center),sigma(sigma),height(height),invsigma(sigma){
       for(unsigned i=0;i<invsigma.size();++i)invsigma[i]=1.0/invsigma[i];
     }
  };
  vector<double> sigma0_;
  vector<Gaussian> hills_;
  FILE* hillsfile_;
  Grid* BiasGrid_;
  FILE* gridfile_;
  double height0_;
  double biasf_;
  double temp_;
  double* dp_;
  int stride_;
  int wgridstride_; 
  bool welltemp_;
  bool restart_;
  bool grid_;
  
  void   readGaussians(FILE*);
  void   writeGaussian(const Gaussian&,FILE*);
  void   addGaussian(const Gaussian&);
  double getHeight(const vector<double>&);
  double getBiasAndDerivatives(const vector<double>&,double* der=NULL);
  double evaluateGaussian(const vector<double>&, const Gaussian&,double* der=NULL);
  vector<unsigned> getGaussianSupport(const Gaussian&);


public:
  BiasMetaD(const ActionOptions&);
  ~BiasMetaD();
  void calculate();
  void update();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasMetaD,"METAD")

void BiasMetaD::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.add("compulsory","SIGMA","the widths of the Gaussian hills");
  keys.add("compulsory","HEIGHT","the heights of the Gaussian hills");
  keys.add("compulsory","PACE","the frequency for hill addition");
  keys.add("compulsory","FILE","HILLS","a file in which the list of added hills is stored");
  keys.addFlag("RESTART",false,"restart the calculation from a previous metadynamics calculation.");
  keys.add("optional","BIASFACTOR","use well tempered metadynamics and use this biasfactor.  Please note you must also specify temp");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
  keys.add("optional","GRID_MIN","the lower bounds for the grid");
  keys.add("optional","GRID_MAX","the upper bounds for the grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.addFlag("GRID_SPARSE",false,"use a sparse grid to store hills");
  keys.addFlag("GRID_NOSPLINE",false,"don't use spline interpolation with grids");
  keys.add("optional","GRID_WSTRIDE","write the grid to a file every N steps");
  keys.add("optional","GRID_WFILE","the file on which to write the grid");
}

BiasMetaD::~BiasMetaD(){
  if(BiasGrid_) delete BiasGrid_;
  if(hillsfile_) fclose(hillsfile_);
  if(gridfile_) fclose(gridfile_);
  delete [] dp_;
}

BiasMetaD::BiasMetaD(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
sigma0_(getNumberOfArguments(),0.0),
hillsfile_(NULL),
BiasGrid_(NULL),
gridfile_(NULL),
height0_(0.0),
biasf_(1.0),
temp_(0.0),
dp_(NULL),
stride_(0),
wgridstride_(0),
welltemp_(false),
restart_(false),
grid_(false)
{
  parseVector("SIGMA",sigma0_);
  plumed_assert(sigma0_.size()==getNumberOfArguments());
  parse("HEIGHT",height0_);
  plumed_assert(height0_>0.0);
  parse("PACE",stride_);
  plumed_assert(stride_>0);
  string hillsfname="HILLS";
  parse("FILE",hillsfname);
  parseFlag("RESTART",restart_);
  parse("BIASFACTOR",biasf_);
  plumed_assert(biasf_>=1.0);
  parse("TEMP",temp_);
  if(biasf_>1.0){
   plumed_assert(temp_>0.0);
   welltemp_=true;
  }

  vector<double> gmin(getNumberOfArguments());
  parseVector("GRID_MIN",gmin);
  plumed_assert(gmin.size()==getNumberOfArguments() || gmin.size()==0);
  vector<double> gmax(getNumberOfArguments());
  parseVector("GRID_MAX",gmax);
  plumed_assert(gmax.size()==getNumberOfArguments() || gmax.size()==0);
  vector<unsigned> gbin(getNumberOfArguments());
  parseVector("GRID_BIN",gbin);
  plumed_assert(gbin.size()==getNumberOfArguments() || gbin.size()==0);
  plumed_assert(gmin.size()==gmax.size() && gmin.size()==gbin.size());
  bool sparsegrid=false;
  parseFlag("GRID_SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("GRID_NOSPLINE",nospline);
  bool spline=!nospline;
  if(gbin.size()>0){grid_=true;}
  parse("GRID_WSTRIDE",wgridstride_);
  string gridfname;
  parse("GRID_WFILE",gridfname); 

  if(grid_&&gridfname.length()>0){plumed_assert(wgridstride_>0);}
  if(grid_&&wgridstride_>0){plumed_assert(gridfname.length()>0);}

  checkRead();

  log.printf("  Gaussian width");
  for(unsigned i=0;i<sigma0_.size();++i) log.printf(" %f",sigma0_[i]);
  log.printf("\n");
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_); 
  log.printf("  Gaussian file %s\n",hillsfname.c_str());
  if(welltemp_){log.printf("  Well-Tempered Bias Factor %f\n",biasf_);}
  if(grid_){
   log.printf("  Grid min");
   for(unsigned i=0;i<gmin.size();++i) log.printf(" %f",gmin[i]);
   log.printf("\n");
   log.printf("  Grid max");
   for(unsigned i=0;i<gmax.size();++i) log.printf(" %f",gmax[i]);
   log.printf("\n");
   log.printf("  Grid bin");
   for(unsigned i=0;i<gbin.size();++i) log.printf(" %d",gbin[i]);
   log.printf("\n");
   if(spline){log.printf("  Grid uses spline interpolation\n");}
   if(sparsegrid){log.printf("  Grid uses sparse grid\n");}
   if(wgridstride_>0){log.printf("  Grid is written on file %s with stride %d\n",gridfname.c_str(),wgridstride_);} 
 }
  
  addValue("bias");

// for performance
   dp_ = new double[getNumberOfArguments()];

// initializing grid
  if(grid_){
   vector<bool> pbc;
   for(unsigned i=0;i<getNumberOfArguments();++i){
    pbc.push_back(getArguments()[i]->isPeriodic());
// if periodic, use CV domain for grid boundaries
    if(pbc[i]){
     double dmin,dmax;
     getArguments()[i]->getDomain(dmin,dmax);
     gmin[i]=dmin;
     gmax[i]=dmax;
    }
   }
   if(!sparsegrid){BiasGrid_=new Grid(gmin,gmax,gbin,pbc,spline,true);}
   else{BiasGrid_=new SparseGrid(gmin,gmax,gbin,pbc,spline,true);}
// open file for grid writing 
   if(wgridstride_>0){gridfile_=fopen(gridfname.c_str(),"w");}
  }

// restarting from HILLS file
  if(restart_){
   hillsfile_=fopen(hillsfname.c_str(),"a+");
   log.printf("  Restarting from %s:",hillsfname.c_str());
   readGaussians(hillsfile_);
  }else{
   hillsfile_=fopen(hillsfname.c_str(),"w");
  } 
}

void BiasMetaD::readGaussians(FILE* file)
{
 unsigned ncv=getNumberOfArguments();
 double dummy;
 vector<double> center(ncv);
 vector<double> sigma(ncv);
 double height;
 int nhills=0;
 rewind(file);
 while(1){
  if(fscanf(file, "%lf", &dummy)!=1){break;}
  for(unsigned i=0;i<ncv;++i){fscanf(file, "%lf", &(center[i]));}
  for(unsigned i=0;i<ncv;++i){fscanf(file, "%lf", &(sigma[i]));}
  fscanf(file, "%lf", &height);
  fscanf(file, "%lf", &dummy);
  nhills++;
  if(welltemp_){height*=(biasf_-1.0)/biasf_;}
  addGaussian(Gaussian(center,sigma,height));
 }     
 log.printf("  %d Gaussians read\n",nhills);
}

void BiasMetaD::writeGaussian(const Gaussian& hill, FILE* file)
{
 unsigned ncv=getNumberOfArguments();
 fprintf(hillsfile_, "%10.3f   ", getTimeStep()*getStep());
 for(unsigned i=0;i<ncv;++i){fprintf(file, "%14.9f   ", hill.center[i]);}
 for(unsigned i=0;i<ncv;++i){fprintf(file, "%14.9f   ", hill.sigma[i]);}
 double height=hill.height;
 if(welltemp_){height*=biasf_/(biasf_-1.0);}
 fprintf(file, "%14.9f   %4.3f \n",height,biasf_);
}

void BiasMetaD::addGaussian(const Gaussian& hill)
{
 if(!grid_){hills_.push_back(hill);} 
 else{
  unsigned ncv=getNumberOfArguments();
  vector<unsigned> nneighb=getGaussianSupport(hill);
  vector<unsigned> neighbors=BiasGrid_->getNeighbors(hill.center,nneighb);
  vector<double> der(ncv);
  vector<double> xx(ncv);
  for(unsigned i=0;i<neighbors.size();++i){
   unsigned ineigh=neighbors[i];
   for(unsigned j=0;j<ncv;++j){der[j]=0.0;}
   BiasGrid_->getPoint(ineigh,xx);   
   double bias=evaluateGaussian(xx,hill,&der[0]);
   BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
  }
 }
}

vector<unsigned> BiasMetaD::getGaussianSupport(const Gaussian& hill)
{
 vector<unsigned> nneigh;
 for(unsigned i=0;i<getNumberOfArguments();++i){
  double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[i];
  nneigh.push_back((unsigned)ceil(cutoff/BiasGrid_->getDx()[i]));
 }
 return nneigh;
}

double BiasMetaD::getBiasAndDerivatives(const vector<double>& cv, double* der)
{
 double bias=0.0;
 if(!grid_){
  for(unsigned i=0;i<hills_.size();++i){
   bias+=evaluateGaussian(cv,hills_[i],der);
  }
 }else{
  if(der){
   vector<double> vder(getNumberOfArguments());
   bias=BiasGrid_->getValueAndDerivatives(cv,vder);
   for(unsigned i=0;i<getNumberOfArguments();++i){der[i]=vder[i];}
  }else{
   bias=BiasGrid_->getValue(cv);
  }
 }
 return bias;
}

double BiasMetaD::evaluateGaussian
 (const vector<double>& cv, const Gaussian& hill, double* der)
{
 double dp2=0.0;
 double bias=0.0;
 for(unsigned i=0;i<cv.size();++i){
  double dp=difference(i,hill.center[i],cv[i])*hill.invsigma[i];
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
 return bias;
}

double BiasMetaD::getHeight(const vector<double>& cv)
{
 double height=height0_;
 if(welltemp_){
    double vbias=getBiasAndDerivatives(cv);
    height=height0_*exp(-vbias/(plumed.getAtoms().getKBoltzmann()*temp_*(biasf_-1.0)));
 } 
 return height;
}

void BiasMetaD::calculate()
{
  unsigned ncv=getNumberOfArguments();
  vector<double> cv(ncv);
  for(unsigned i=0;i<ncv;++i){cv[i]=getArgument(i);}

  double* der=new double[ncv];
  for(unsigned i=0;i<ncv;++i){der[i]=0.0;}
  double ene=getBiasAndDerivatives(cv,der);

  Value* value=getValue("bias");
  setValue(value,ene);

// set Forces 
  for(unsigned i=0;i<ncv;++i){
   const double f=-der[i];
   setOutputForces(i,f);
  }

  delete [] der;
}

void BiasMetaD::update(){
  vector<double> cv(getNumberOfArguments());
  for(unsigned i=0;i<cv.size();++i){cv[i]=getArgument(i);}
  if(getStep()%stride_==0){
// add a Gaussian
   double height=getHeight(cv);
   Gaussian newhill=Gaussian(cv,sigma0_,height);
   addGaussian(newhill);
// print on HILLS file
   writeGaussian(newhill,hillsfile_);
  }
// dump grid on file
  if(wgridstride_>0&&getStep()%wgridstride_==0){
   BiasGrid_->writeToFile(gridfile_); 
  }
}

}
