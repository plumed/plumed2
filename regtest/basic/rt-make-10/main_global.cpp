/*
Test with gcmd
*/
#include <vector>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <type_traits>

struct vec3d {
  double x;
  double y;
  double z;
};  
      
struct tens3d3d {
  double xx;
  double xy;
  double xz;
  double yx;
  double yy;
  double yz;
  double zx;
  double zy;
  double zz;
  double & operator[](unsigned i) {
    return xx;
  }   
  const double & operator[](unsigned i) const {
    return xx;
  } 
};    
      
namespace PLMD {
namespace wrapper {

template<typename T> struct is_custom_array;

template<>
struct is_custom_array<vec3d> : std::true_type {
  using value_type = double;
};

template<>
struct is_custom_array<tens3d3d> : std::true_type {
  using value_type = vec3d;
};
}
}

#include "plumed/wrapper/Plumed.h"
using namespace PLMD;


void test_global(){
  std::ofstream output("output_global");

#define P_CREATE Plumed::gcreate()
#define P_CMD Plumed::gcmd
#define P_FINALIZE Plumed::gfinalize();

#include "main.inc"

}

