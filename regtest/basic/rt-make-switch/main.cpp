
#include "plumed/tools/SwitchingFunction.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <utility>


using PLMD::SwitchingFunction;
using PLMD::SwitchingFunctionAccelerable;


struct callSQR {
  constexpr static auto kind="sqr_";
  template<typename T>
  static std::pair <double,double> call(const T& sw, double point) {
    double deriv;
    point *= point;
    double value = sw.calculateSqr(point,deriv);
    return std::make_pair(value,deriv);
  }
};

struct callplain {
  constexpr static auto kind="";
  template<typename T>
  static std::pair <double,double> call(const T& sw, double point) {
    double deriv;
    double value = sw.calculate(point,deriv);
    return std::make_pair(value,deriv);
  }
};

template<typename caller>
void test(const std::string& name, std::string definition,bool stretch=true) {
  if (!stretch) {
    definition += " NOSTRETCH ";
  }

  std::ofstream os("out_"+std::string(caller::kind)+name+(stretch ? "" : "_nostretch"));
  os <<std::fixed<<std::setprecision(6);
  SwitchingFunction sw;
  std::string error;
  sw.set(definition,error);
  if (!error.empty()) {
    std::cerr<<error<<"\n";
  }

  error.clear();
  SwitchingFunctionAccelerable swa;
  swa.set(definition,error);
  if (!error.empty()) {
    std::cerr<<error<<"\n";
  }
  os << "point :\tvalue deriv\tvalue_acc deriv_acc\tvalue_delta deriv_delta\n";
  for (int i=0; i<10; i++) {
    double point=i/2.50;
    double deriv;
    double value;
    os<< point <<" :\t";
    std::tie(value,deriv) = caller::call(sw, point);
    os << value << " " << deriv <<"\t";
    double deriv_acc;
    double value_acc;
    std::tie(value_acc,deriv_acc) = caller::call(swa, point);
    os << value_acc << " " << deriv_acc <<"\t";
    os << value-value_acc << " " << deriv-deriv_acc <<"\n";
  }
}


int main() {
  for(const auto &x: {
  std::pair<std::string,std::string> {"COSINUS R_0=2.6","cosinus"},
{"EXP R_0=0.8 D_0=0.5 D_MAX=2.6","exp"},
{"GAUSSIAN R_0=1.0 D_0=0.0 D_MAX=2.6","fastgaussian"},
{"GAUSSIAN R_0=1.0 D_0=0.3 D_MAX=2.6","gaussian"},
{"RATIONAL R_0=1.3 NN=6 MM=10 D_MAX=2.6","fastrational"},
{"RATIONAL R_0=1.3 D_MAX=2.6","fastrational_MMeq2MM"},
{"RATIONAL R_0=1.3 NN=5 MM=11 D_MAX=2.6","rational"},
{"RATIONAL R_0=1.3 NN=5 D_MAX=2.6","rational_MMeq2MM"},
{"Q R_0=1.0 D_0=0.3 BETA=1.0 LAMBDA=1.0 REF=1.3 D_MAX=2.6","q"},
{"SMAP R_0=1.3 A=3 B=2 D_MAX=2.6","smap"},
{"TANH R_0=1.3 D_MAX=2.6","tanh"}
    }) {
    auto [definition, name] = x;
    test<callplain> (name,definition);
    test<callplain> (name,definition,false);
    test<callSQR> (name,definition);
    test<callSQR> (name,definition,false);

  }

  return 0;
}
