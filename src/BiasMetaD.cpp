#include "Bias.h"
#include "ActionRegister.h"
#include "Grid.h"

#include <cassert>
#include <iostream>

#define DP2CUTOFF 6.25

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS METAD 
/**
MetaDynamics on one or more variables

\par Syntax
\verbatim
METAD ...
  ARG=x1,x2,... 
  SIGMA=s1,s2,... 
  HEIGHT=w 
  PACE=s
  RESTART
  [BIASFACTOR=biasf]
  [TEMP=temp]
  [GRIDMIN=min1,min2]
  [GRIDMAX=max1,max2]
  [GRIDBIN=bin1,bin2]
  [NOSPLINE]
  [SPARSEGRID]
... METAD
\endverbatim
SIGMA specifies an array of Gaussian widths, one for each variable,
HEIGHT the Gaussian height, PACE the Gaussian deposition stride in steps,
RESTART to restart the run, BIASFACTOR the bias factor of well-tempered metad, 
TEMP the temperature. To activate the use of Grid to store the bias potential,
you need to specify the grid boundaries with GRIDMIN and GRIDMAX and the number
of bins with GRIDBIN. An experimental sparse grid can be activated with SPARSEGRID.
The use of spline can be disabled with NOSPLINE.
\par Example
The following input is for a standard metadynamics calculation using as
collective variables the distance between atoms 3 and 5
and the distance between atoms 2 and 4. The value of the CVs and
the metadynamics bias potential is written to the COLVAR file every 100 steps.
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
  };
  vector<double> sigma0_;
  vector<Gaussian> hills_;
  FILE* hillsfile_;
  Grid* BiasGrid_;
  double height0_;
  double biasf_;
  double temp_;
  int stride_;
  bool welltemp_;
  bool restart_;
  bool grid_;
  
  void   readGaussians(FILE*);
  void   writeGaussian(Gaussian,FILE*);
  void   addGaussian(Gaussian);
  double getHeight(vector<double>);
  double getBiasAndDerivatives(vector<double>,double*);
  double evaluateGaussian(vector<double>,Gaussian,double*);
  vector<unsigned> getGaussianSupport(Gaussian);


public:
  BiasMetaD(const ActionOptions&);
  void calculate();
};

PLUMED_REGISTER_ACTION(BiasMetaD,"METAD")

BiasMetaD::BiasMetaD(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
sigma0_(getNumberOfArguments(),0.0),
height0_(0.0),
biasf_(1.0),
temp_(0.0),
stride_(0),
welltemp_(false),
restart_(false),
grid_(false)
{
  parseVector("SIGMA",sigma0_);
  assert(sigma0_.size()==getNumberOfArguments());
  parse("HEIGHT",height0_);
  assert(height0_>0.0);
  parse("PACE",stride_);
  assert(stride_>0); 
  parseFlag("RESTART",restart_);
  parse("BIASFACTOR",biasf_);
  assert(biasf_>=1.0);
  parse("TEMP",temp_);
  if(biasf_>1.0){
   assert(temp_>0.0);
   welltemp_=true;
  }
  vector<double> gmin;
  parseVector("GRIDMIN",gmin);
  assert(gmin.size()==getNumberOfArguments() || gmin.size()==0);
  vector<double> gmax;
  parseVector("GRIDMAX",gmax);
  assert(gmax.size()==getNumberOfArguments() || gmax.size()==0);
  vector<unsigned> gbin;
  parseVector("GRIDBIN",gbin);
  assert(gbin.size()==getNumberOfArguments() || gbin.size()==0);
  assert(gmin.size()==gmax.size() && gmin.size()==gbin.size());
  bool sparsegrid=false;
  parseFlag("SPARSEGRID",sparsegrid);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  bool spline=!nospline;
  if(gbin.size()>0){grid_=true;}

  checkRead();

  log.printf("  Gaussian width");
  for(unsigned i=0;i<sigma0_.size();++i) log.printf(" %f",sigma0_[i]);
  log.printf("\n");
  log.printf("  Gaussian height %f\n",height0_);
  log.printf("  Gaussian deposition pace %d\n",stride_); 
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
  }
  
  addValue("bias");

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
  }

// open HILLS file
  hillsfile_=fopen("HILLS","a+");
// restarting from HILLS file
  if(restart_){
   log.printf("  Restarting");
   readGaussians(hillsfile_);
  } 
}

void BiasMetaD::readGaussians(FILE* file)
{
 int ncv=getNumberOfArguments();
 double dummy;
 vector<double> center;
 vector<double> sigma;
 center.resize(ncv);
 sigma.resize(ncv);
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
  Gaussian newhill={center,sigma,height}; 
  addGaussian(newhill);
 }     
 log.printf("  %d Gaussians read\n",nhills);
}

void BiasMetaD::writeGaussian(Gaussian hill, FILE* file)
{
 int ncv=getNumberOfArguments();
 fprintf(hillsfile_, "%10.3f   ", getTimeStep()*getStep());
 for(unsigned i=0;i<ncv;++i){fprintf(file, "%14.9f   ", hill.center[i]);}
 for(unsigned i=0;i<ncv;++i){fprintf(file, "%14.9f   ", hill.sigma[i]);}
 double height=hill.height;
 if(welltemp_){height*=biasf_/(biasf_-1.0);}
 fprintf(file, "%14.9f   %4.3f \n",height,biasf_);
}

void BiasMetaD::addGaussian(Gaussian hill)
{
 if(!grid_){hills_.push_back(hill);} 
 else{
  int ncv=getNumberOfArguments();
  vector<unsigned> nneighb=getGaussianSupport(hill);
  vector<unsigned> neighbors=BiasGrid_->getNeighbors(hill.center,nneighb);
  double* der=new double[ncv];
  for(unsigned i=0;i<neighbors.size();++i){
   unsigned ineigh=neighbors[i];
   for(unsigned j=0;j<ncv;++j){der[j]=0.0;}
   double bias=evaluateGaussian(BiasGrid_->getPoint(ineigh),hill,der);
   vector<double> vder(der, der+ncv);
   BiasGrid_->addValueAndDerivatives(ineigh,bias,vder);
  }
  delete(der);
 }
}

vector<unsigned> BiasMetaD::getGaussianSupport(Gaussian hill)
{
 vector<unsigned> nneigh;
 for(unsigned i=0;i<getNumberOfArguments();++i){
  double cutoff=sqrt(2.0*DP2CUTOFF)*hill.sigma[i];
  nneigh.push_back((unsigned)ceil(cutoff/BiasGrid_->getDx()[i]));
 }
 return nneigh;
}

double BiasMetaD::getBiasAndDerivatives(vector<double> cv, double* der)
{
 double bias=0.0;
 if(!grid_){
  vector<Gaussian>::iterator it;
  for(it=hills_.begin();it!=hills_.end();it++){
   bias+=evaluateGaussian(cv,*it,der);
  }
 }else{
  if(der!=NULL){
   vector<double> vder;
   vder.resize(getNumberOfArguments());  
   bias=BiasGrid_->getValueAndDerivatives(cv,vder);
   for(unsigned i=0;i<getNumberOfArguments();++i){der[i]=vder[i];}
  }else{
   bias=BiasGrid_->getValue(cv);
  }
 }
 return bias;
}

double BiasMetaD::evaluateGaussian
 (vector<double> cv, Gaussian hill, double* der)
{
 vector<double> dp;
 double dp2=0.0;
 double bias=0.0;
 for(unsigned i=0;i<cv.size();++i){
  dp.push_back(difference(i,hill.center[i],cv[i])/hill.sigma[i]);
  dp2+=dp[i]*dp[i];
 }
 dp2*=0.5;
 if(dp2<DP2CUTOFF){
  bias=hill.height*exp(-dp2);
  if(der!=NULL){
   for(unsigned i=0;i<cv.size();++i){der[i]+=-bias*dp[i]/hill.sigma[i];}
  }
 }
 return bias;
}

double BiasMetaD::getHeight(vector<double> cv)
{
 double height=height0_;
 if(welltemp_){
    double vbias=getBiasAndDerivatives(cv,NULL);
    height=height0_*exp(-vbias/(temp_*(biasf_-1.0)));
 } 
 return height;
}

void BiasMetaD::calculate()
{
  vector<double> cv;
  int ncv=getNumberOfArguments();
  for(unsigned i=0;i<ncv;++i){cv.push_back(getArgument(i));}

  if(getStep()%stride_==0){
// add a Gaussian
   double height=getHeight(cv);
   Gaussian newhill={cv,sigma0_,height}; 
   addGaussian(newhill);
// print on HILLS file
   writeGaussian(newhill,hillsfile_);
  }

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
  delete(der);
}

}
