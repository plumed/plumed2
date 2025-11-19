#ifndef __PLUMED_TEST_MACROS
#define __PLUMED_TEST_MACROS
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <string>

///A very barebone tee implementation
///Useful for runnig make test on the fly with
class tee {
  std::ofstream ofs;
public:
  tee(std::string filename) : ofs(filename) {}
  template<typename T>
  tee& operator<<(const T& t) {
    ofs<<t;
    std::cout <<t;
    return *this;
  }
  template<typename Iterable>
  tee&  dump(const Iterable& v) {
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(ofs," "));
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(std::cout," "));
    return *this;
  }
};

class etee {
  std::ofstream ofs;
public:
  etee(std::string filename) : ofs(filename) {}
  template<typename T>
  etee& operator<<(const T& t) {
    ofs<<t;
    std::cerr <<t;
    return *this;
  }
  template<typename Iterable>
  etee&  dump(const Iterable& v) {
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(ofs," "));
    std::copy(v.begin(),v.end(),std::ostream_iterator<double>(std::cerr," "));
    return *this;
  }
};

namespace PLMDTests {

template <typename T>
struct namedVal {
  std::string name;
  T val;
};

template <typename T>
using testcases=std::vector<namedVal<T>>;

}
#endif //__PLUMED_TEST_MACROS
