#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

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
  STRIDE=s
  [BIASFACTOR=biasf]
  [TEMP=temp]
... METAD
\endverbatim
SIGMA specifies an array of Gaussian widths, one for each variable,
HEIGHT the Gaussian height, STRIDE the Gaussian deposition pace in steps,
BIASFACTOR the bias factor for well-tempered metadynamics, 
TEMP the temperature.
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

struct Gaussian {
  vector<double> center;
  vector<double> sigma;
  double height; 
};

class BiasMetaD : public Bias{

private:
  vector<double> sigma_;
  vector<Gaussian> hills_;
  double w0_;
  double biasf_;
  double temp_;
  int stride_;
  bool welltemp_;
  bool grid_;
  double getHeight(vector<double>);
  void addGaussian(vector<double>,vector<double>,double);
  void addGaussianToList(vector<double>,vector<double>,double);
  double getBiasAndForces(vector<double>,double*);
  double getBiasAndForcesFromList(vector<double>,double*);

public:
  BiasMetaD(const ActionOptions&);
  void calculate();
};

PLUMED_REGISTER_ACTION(BiasMetaD,"METAD")

BiasMetaD::BiasMetaD(const ActionOptions& ao):
PLUMED_BIAS_INIT(ao),
sigma_(getNumberOfArguments(),0.0),
w0_(0.0),
biasf_(1.0),
temp_(0.0),
stride_(0),
welltemp_(false),
grid_(false)
{
  parseVector("SIGMA",sigma_);
  assert(sigma_.size()==getNumberOfArguments());
  parse("HEIGHT",w0_);
  assert(w0_>0.0);
  parse("PACE",stride_);
  assert(stride_>0); 
  parse("BIASFACTOR",biasf_);
  assert(biasf_>=1.0);
  parse("TEMP",temp_);
  if(biasf_>1.0){
   assert(temp_>0.0);
   welltemp_=true;
  }  
  checkRead();

  log.printf("  Gaussian width");
  for(unsigned i=0;i<sigma_.size();++i) log.printf(" %f",sigma_[i]);
  log.printf("\n");
  log.printf("  Gaussian height %f\n",w0_);
  log.printf("  Gaussian deposition pace %d\n",stride_); 
  if(welltemp_){log.printf("  Well-Tempered Bias Factor %f\n",biasf_);}

  addValue("bias");
}

void BiasMetaD::addGaussian(vector<double> center,
 vector<double> sigma, double height)
{
 if(!grid_){addGaussianToList(center,sigma,height);} 
}

void BiasMetaD::addGaussianToList(vector<double> center,
 vector<double> sigma, double height)
{
 Gaussian newhill;
 newhill.center=center;
 newhill.sigma=sigma;
 newhill.height=height;
 hills_.push_back(newhill);
}

double BiasMetaD::getBiasAndForces(vector<double> cv, double* ff)
{
 if(!grid_){return getBiasAndForcesFromList(cv,ff);}
}

double BiasMetaD::getBiasAndForcesFromList(vector<double> cv, double* ff)
{
 double bias=0.0;
 vector<Gaussian>::iterator it;
// cycle on the Gaussians deposited
 for(it=hills_.begin();it!=hills_.end();it++)
 {
  vector<double> dp;
  double dp2=0.0;
  for(unsigned i=0;i<cv.size();++i){
   dp.push_back(difference(i,cv[i],(*it).center[i])/(*it).sigma[i]);
   dp2+=dp[i]*dp[i];
  }
// put DP2CUTOFF here
  double newbias=(*it).height*exp(-0.5*dp2);
  bias+=newbias;
  if(ff!=NULL){
   for(unsigned i=0;i<cv.size();++i){
    ff[i]+=newbias*dp[i]/(*it).sigma[i];
   }
  }
 }
 return bias;
}

double BiasMetaD::getHeight(vector<double> cv)
{
 double ww=w0_;
 if(welltemp_){
    double vbias=getBiasAndForces(cv,NULL);
    ww=w0_*exp(-vbias/(temp_*(biasf_-1.0)));
 } 
 return ww;
}

void BiasMetaD::calculate(){
  vector<double> cv;
  int ncv=getNumberOfArguments();
  for(unsigned i=0;i<ncv;++i){cv.push_back(getArgument(i));}

  if(getStep()%stride_==0){
// add a Gaussian
   double ww=getHeight(cv);
   addGaussian(cv,sigma_,ww);
// print on HILLS file
   //printGaussian(cv,sigma_,ww);
  }

  double* ff=new double[ncv];
  for(unsigned i=0;i<ncv;++i){ff[i]=0.0;}
  double ene=getBiasAndForces(cv,ff);

  Value* value=getValue("bias");
  setValue(value,ene);

// set Forces 
  for(unsigned i=0;i<ncv;++i){
   const double f=ff[i];
   setOutputForces(i,f);
  }
  delete(ff);
}

}
