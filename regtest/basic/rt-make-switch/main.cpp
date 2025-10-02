
#include "plumed/tools/SwitchingFunction.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <utility>


using PLMD::SwitchingFunction;
//using PLMD::SwitchingFunctionAccelerable;


struct callSQR {
  template<typename T>
  static std::pair <double,double> call(const T& sw, double point) {
    double deriv;
    point *= point;
    double value = sw.calculateSqr(point,deriv);
    return std::make_pair(value,deriv);
  }
};

struct callplain {
  template<typename T>
  static std::pair <double,double> call(const T& sw, const double point) {
    double deriv;
    double value = sw.calculate(point,deriv);
    return std::make_pair(value,deriv);
  }
};

template<typename caller,typename SW>
std::ostream& calculateAndSubtract(const SW& sw,
                                   const double point,
                                   const double value,
                                   const double derivative,
                                   std::ostream& out) {
  double mydev;
  double myval;
  std::tie(myval,mydev) = caller::call(sw,point);
  out << (myval -value)<< " " << (mydev - derivative)<<" " ;
  return out;
}


void test(const std::string& name, std::string definition,std::string custom,bool stretch=true) {
  if (!stretch) {
    definition += " NOSTRETCH ";
    custom += " NOSTRETCH ";
  }

  std::ofstream os("out_"+name+(stretch ? "" : "_nostretch"));
  os <<std::fixed<<std::setprecision(6);
  SwitchingFunction sw;
  std::string error;
  sw.set(definition,error);
  if (!error.empty()) {
    std::cerr<<error<<"\n";
  }
  error.clear();

  SwitchingFunction swc;
  swc.set(custom,error);
  if (!error.empty()) {
    std::cerr<<error<<"\n";
  }
  error.clear();
/*
  SwitchingFunctionAccelerable swa;
  swa.set(definition,error);
  if (!error.empty()) {
    std::cerr<<error<<"\n";
  }
*/
  os << "point : value deriv vsq_delta dsq_delta "
     "vAcc_delta dAcc_delta "
     //"vAccsq_delta dAccsq_delta "
     "vcus_delta dcus_delta "
     "vsqcus_delta dsqcus_delta\n";
  for (int i=0; i<10; i++) {
    double point=i/2.50;
    double derivative;
    double value;
    os<< point <<" : ";
    std::tie(value,derivative) = callplain::call(sw, point);
    os << value << " " << derivative <<" ";
    calculateAndSubtract<callSQR>(sw,point, value, derivative,os) << " ";
    //calculateAndSubtract<callplain>(swa,point, value, derivative,os) << " ";
    //calculateAndSubtract<callSQR>(swa,point, value, derivative,os) << " ";
    calculateAndSubtract<callplain>(swc,point, value, derivative,os) << " ";
    calculateAndSubtract<callSQR>(swc,point, value, derivative,os) << "\n";
  }
}


int main() {
//the test will output the value of the switch along with de difference of it and
//the same calculation under different settings (with squared/non squared input,
//with the "accelerable version and with custom)
  for(const auto &x: {
  std::tuple<std::string,std::string,std::string> {
  "COSINUS R_0=2.6",
  "CUSTOM FUNC=0.5*(cos(x*pi)+1) R_0=2.6 D_MAX=2.6",
  "cosinus"
}, {
  "EXP R_0=0.8 D_0=0.5 D_MAX=2.6",
  "CUSTOM FUNC=exp(-x) R_0=0.8 D_0=0.5 D_MAX=2.6",
  "exp"
}, {
  "GAUSSIAN R_0=1.0 D_0=0.0 D_MAX=2.6",
  "CUSTOM FUNC=exp(-x^2/2) R_0=1.0 D_0=0.0 D_MAX=2.6",
  "fastgaussian"
}, {
  "GAUSSIAN R_0=1.0 D_0=0.3 D_MAX=2.6",
  "CUSTOM FUNC=exp(-x^2/2) R_0=1.0 D_0=0.3 D_MAX=2.6",
  "gaussian"
}, {
  "RATIONAL R_0=1.3 NN=6 MM=10 D_MAX=2.6",
  "CUSTOM FUNC=(1-x^6)/(1-x^10) D_MAX=2.6 R_0=1.3",
  "fastrational"
}, {
  "RATIONAL R_0=1.3 D_MAX=2.6",
  "CUSTOM FUNC=(1-x^6)/(1-x^12) D_MAX=2.6 R_0=1.3",
  "fastrational_NNeq2MM"
}, {
  "RATIONAL R_0=1.3 NN=5 MM=11 D_MAX=2.6",
  "CUSTOM FUNC=(1-x^5)/(1-x^11) D_MAX=2.6 R_0=1.3",
  "rational"
}, {
  "RATIONAL R_0=1.3 NN=5 D_MAX=2.6",
  "CUSTOM FUNC=(1-x^5)/(1-x^10) D_MAX=2.6 R_0=1.3",
  "rational_NNeq2MM"
}, {
  "Q R_0=1.0 D_0=0.3 BETA=5.0 LAMBDA=1.0 REF=1.3 D_MAX=2.6",
  "CUSTOM FUNC=1/(1+exp(5.0*((x+0.3)-1.3))) R_0=1.0 D_0=0.3 D_MAX=2.6",
  "q"
}, {
  "TANH R_0=1.3 D_MAX=2.6",
  "CUSTOM FUNC=1-tanh(x) R_0=1.3 D_MAX=2.6",
  "tanh"
}, {
  "SMAP R_0=1.3 A=3 B=2 D_MAX=2.6",
  "CUSTOM FUNC=(1+(2^(3/2)-1)*x^3)^(-2/3) R_0=1.3 D_MAX=2.6",
  "smap"
}, {
  "CUBIC D_MAX=2.6 D_0=0.6",
//Here R_O = DMAX-D_0
  "CUSTOM FUNC=(x-1)^2*(1+2*x) R_0=2.0 D_0=0.6 D_MAX=2.6",
  "cubic"
}
    }) {
    //auto [definition, custom, name] = x;
    std::string definition, custom, name;
    std::tie(definition, custom, name) = x;
    test (name,definition,custom);
    test (name,definition,custom,false);

  }

  return 0;
}
