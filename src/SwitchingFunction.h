#ifndef __PLUMED_SwitchingFunction_h
#define __PLUMED_SwitchingFunction_h

#include <cmath>
#include <string>
#include "Keywords.h"
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
  enum {spline,exponential,gaussian} type;
  int nn,mm;
  double invr0,d0,dmax;
public:
  static std::string documentation();
  SwitchingFunction();
  void set(int nn,int mm,double r_0,double d_0);
  void set(const std::string& definition, std::string& errormsg);
  std::string description() const ;
  double calculate(double x,double&df)const;
  double get_r0() const;
  void printKeywords( Log& log ) const ;
};

}

#endif

