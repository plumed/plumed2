#include "plumed/tools/Tools.h"

#include "testUtils.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace PLMD;

template<typename T>
using FSUM=PLMD::Tools::FastStringUnorderedMap<T>;

int main() {
  try {
    tee out("output");
    FSUM<int> mymap;
    //insertions
    std::string addr;
    for (int i=0; i<90; i+=10) {
      addr=std::to_string(i);
      mymap[addr]= i;
    }
    {
      //This tests that conv will actually insert a  "90\0"
      //and not a"900" with not null char at the end
      std::string tmp = "90123456";
      std::string_view k {tmp.c_str(),2};
      mymap[k] = 90;
    }
    for (const auto & x: {
           "90", "0", "10","20","30","40","50","60","70","80"
         } ) {
      auto kk = mymap.find(x);
      plumed_assert(kk != mymap.end())
          <<"\""<< x << "\" should have been inserted";
      out <<std::left<< std::setw(3) << x << mymap[x] <<"\n";
    }
  } catch(PLMD::Exception &e) {
    std::cerr << "Exception:" << e.what() << "\n";
  }
  return 0;
}
