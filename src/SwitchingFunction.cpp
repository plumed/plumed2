#include "SwitchingFunction.h"
#include <vector>
#include <limits>
#include "Tools.h"

using namespace std;
using namespace PLMD;

void SwitchingFunction::set(const std::string & definition){
  vector<string> data=Tools::getWords(definition);
  plumed_assert(data.size()>=1);
  string name=data[0];
  data.erase(data.begin());
  invr0=0.0;
  d0=0.0;
  dmax=std::numeric_limits<double>::max();
  init=true;

  double r0;
  bool found_r0=Tools::parse(data,"R_0",r0);
  plumed_massert(found_r0,"R_0 is needed");
  invr0=1.0/r0;
  Tools::parse(data,"D_0",d0);
  Tools::parse(data,"D_MAX",dmax);

  if(name=="SPLINE"){
    type=spline;
    nn=6;
    mm=12;
    Tools::parse(data,"NN",nn);
    Tools::parse(data,"MM",mm);
  } else if(name=="EXP") type=exponential;
  else if(name=="GAUSSIAN") type=gaussian;
  else plumed_merror("Cannot understand switching function type '"+name+"'");
  plumed_massert(data.size()==0,"Error reading");
}

double SwitchingFunction::calculate(double distance,double&dfunc)const{
  plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  const double rdist = (distance-d0)*invr0;
  double result;
  if(rdist<=0.){
     result=1.;
     dfunc=0.0;
  }else if(rdist>dmax){
     result=0.;
     dfunc=0.0;
  }else{
    if(type==spline){
      if(rdist>(1.-100.0*epsilon) && rdist<(1+100.0*epsilon)){
         result=nn/mm;
         dfunc=0.5*nn*(nn-mm)/mm;
      }else{
         double rNdist=rdist;
         double rMdist=rdist;
    // this is a naive optimization
    // we probably have to implement some generic, fast pow(double,int)
        if(nn>2) for(int i=0;i<nn-2;i++) rNdist*=rdist;
        else rNdist = pow(rdist, nn-1);
        if(mm>2) for(int i=0;i<mm-2;i++) rMdist*=rdist;
        else rMdist = pow(rdist, mm-1);
         double num = 1.-rNdist*rdist;
         double iden = 1./(1.-rMdist*rdist);
         double func = num*iden;
         result = func;
         dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist));
      }
    }else if(type==exponential){
      result=exp(-rdist);
      dfunc=-result;
    }else if(type==gaussian){
      result=exp(-0.5*rdist*rdist);
      dfunc=-rdist*result;
    }else plumed_merror("Unknown switching function type");
// this is for the chain rule:
    dfunc*=invr0;
// this is because calculate() sets dfunc to the derivative divided times the distance.
// (I think this is misleading and I would like to modify it - GB)
    dfunc/=distance;
  }
  return result;
}

SwitchingFunction::SwitchingFunction():
  init(false){
}

void SwitchingFunction::set(int nn,int mm,double r0,double d0){
  init=true;
  type=spline;
  this->nn=nn;
  this->mm=mm;
  this->invr0=1.0/r0;
  this->d0=d0;
  this->dmax=pow(0.00001,1./(nn-mm));
}




