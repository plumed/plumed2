#ifndef __PLUMED_SwitchingFunction_h
#define __PLUMED_SwitchingFunction_h

#include <cmath>
#include "PlumedException.h"

namespace PLMD {

//+DEVELDOC TOOLBOX SwitchingFunction
/**
A class for calculating the switching function : \f$\frac{1 - \left( \frac{r-d_0}{r_0} \right)^n}{1 - \left( \frac{r-d_0}{r_0} \right)^m}\f$
*/
//+ENDDEVELDOC

/// Small class to compure switching functions.
/// In the future we might extend it so as to be set using
/// a string:
/// void set(std::string);
/// which can then be parsed for more complex stuff, e.g. exponentials
/// tabulated functions from file, matheval, etc...
class SwitchingFunction{
  bool init;
  int nn,mm;
  double r_0,d_0;
  double epsilon;
  double threshold;
public:
  SwitchingFunction();
  void set(int nn,int mm,double r_0,double d_0);
  double calculate(double x,double&df)const;
};

inline
SwitchingFunction::SwitchingFunction():
  init(false){
}

inline
void SwitchingFunction::set(int nn,int mm,double r_0,double d_0){
  init=true;
  this->nn=nn;
  this->mm=mm;
  this->r_0=r_0;
  this->d_0=d_0;
  epsilon=1e-10;
  threshold=pow(0.00001,1./(nn-mm));
}

inline
double SwitchingFunction::calculate(double distance,double&dfunc)const{
  plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  const double rdist = (distance-d_0)/r_0;
  double result=0.;
  if(rdist<=0.){
     result+=1.;
     dfunc=0.;
  }else if(rdist>(1.-epsilon) && rdist<(1+epsilon)){
     result+=nn/mm;
     dfunc=0.5*nn*(nn-mm)/mm;
  }else if(rdist>threshold){
     dfunc=0.;
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
     result += func;
     dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist))/(distance*r_0);
  }
  return result;
}

}

#endif

