#include "plumed/tools/Tools.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace PLMD;

class tee {
//TODO:: put this tee into a common test utilities lib
  std::ofstream ofs;
public:
  tee(std::string filename) : ofs(filename) {}
  template<typename T>
  tee& operator<<(const T& t) {
    ofs<<t;
    std::cout <<t;
    return *this;
  }
};

int main() {
  try {
//NOTE:: We are simply checking that gw has the expected behaviour
    tee out("output");
    for (auto test : {
           "{abcd},{a b d}","a,b,c,d","{1,2,3},4,5,{6}",
           "@ndx:{file1 second}","@ndx:file,@ndx:{file 1 second}",
           "@replicas:{1,2,3} {0,0,0} {2,4,6}",
           "{1 @replicas:{2 0 4} 3}"
         }) {
      out << std::quoted(test)<< '\n';
      gch::small_vector<std::string_view> gws ;
      Tools::getWordsSimple(gws,test,", \t\n");
      auto gw  = Tools::getWords(test,", \t\n");
      for(auto x : gws ) {

        out <<"\tgws:"<< std::quoted(x) << "\n";
      }
      for(auto x : gw ) {
        out <<"\tgw :"<< std::quoted(x) << "\n";
      }
      out <<'\n';
    }
    out << "\n\nunravelReplicas\n";
    for(auto test: {
          "@replicas:{1,2,3} {0,0,0} {2,4,6}",
          "@replicas:{{1,2,3} {0,0,0} {2,4,6}}",
          "@replicas:1 2 3"
        }) {
      out << std::quoted(test) << "\n";
      for (int rep : {
             0,1,2
           }) {
        auto res=Tools::unravelReplicas(test, rep);
        out << "Replica: "<< rep << "\n";
        out << "\t" << res << " ";
        out << "\n";
      }
    }
    out << "\n\nreadVector\n";
    for(auto test: {
          //this will be passed as intended by the LineToken class, the output of this first line is wrong!
          //"@replicas:{{1,2,3} {0,0,0} {2,4,6}}",
          //"@replicas:1 2 3",
          "1,@replicas:{1,2,3},31,999,3",
          "1,999,3",
          //this works as intended, sadly it is not a possible input, currently
          "@replicas:{1,2,3},999,3",
          "@replicas:{1,2,3},@replicas:{3,5,6},3"
        }) {
      out << std::quoted(test) << "\n";
      for (int rep : {
             0,1,2
           }) {
        std::vector<int> res;
        Tools::parseVector(test,res,rep);
        out << "Replica: "<< rep << "\n";
        out << "\t" ;
        for (auto x:res) {
          out << x << " ";
        }
        out << "\n";
      }
    }
  } catch(PLMD::Exception &e) {
    std::cerr << "Exception:" << e.what() << "\n";
  }
  return 0;
}
