#include "plumed/tools/Tools.h"
#include "plumed/tools/AtomNumber.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

using namespace PLMD;

template<typename T>
void test(std::ofstream & out,const std::string & type,const std::string & str) {
  T t{};
  const bool ok=Tools::convertNoexcept(str,t);
  out<<type<<" \""<<str<<"\": ";
  if(!ok) {
    out<<"not converted\n";
    return;
  }
  if constexpr (std::is_floating_point_v<T>) {
    if(std::isnan(t)) {
      out<<"nan\n";
    } else if(std::isinf(t)) {
      out<<(t>0?"+inf":"-inf")<<"\n";
    } else {
      out<<t<<"\n";
    }
  } else {
    out<<t<<"\n";
  }
}

int main() {
  std::ofstream out("output");
  out.precision(17);
  for(auto s : {
        "1","-2.5"," 3.25 ","1e-3","+4","0.1","nan","-nan","NAN","NaN","inf","-inf","+Inf","infinity",
        "pi","2pi","-PI","1+1","exp(1)","1.5x","","nanx","1,5"
      }) {
    test<double>(out,"double",s);
  }
  for(auto s : {
        "0.1","nan","-inf","1e40"
      }) {
    test<float>(out,"float",s);
  }
  for(auto s : {
        "12"," 12 ","+12","-12","2*3","1.0","1.5","2147483647","2147483648","-2147483648"
      }) {
    test<int>(out,"int",s);
  }
  // negative numbers must not be converted to unsigned types
  for(auto s : {
        "12","+12","-1"," -1","-0","4294967295","4294967296","3-5","2*3"
      }) {
    test<unsigned>(out,"unsigned",s);
  }
  for(auto s : {
        "-1","18446744073709551615"
      }) {
    test<unsigned long>(out,"unsigned long",s);
  }
  {
    AtomNumber a;
    out<<"AtomNumber \"-1\": "<<(Tools::convertNoexcept("-1",a)?"converted":"not converted")<<"\n";
    out<<"AtomNumber \"7\": "<<(Tools::convertNoexcept("7",a)?std::to_string(a.serial()):"not converted")<<"\n";
  }
  return 0;
}
