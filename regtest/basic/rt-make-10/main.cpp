#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <fstream>
#include <cstdio>
#include <iostream>

using namespace PLMD;

struct vec3d {
  double x;
  double y;
  double z;
};

template<>
struct PLMD::wrapper::is_custom_array<vec3d> : std::true_type {
  using value_type = double;
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

template<>
struct PLMD::wrapper::is_custom_array<tens3d3d> : std::true_type {
  using value_type = vec3d;
};



// preliminary test on pass by value
void preliminary(){
  std::ofstream output("output_preliminary");

  {

    Plumed p;

    // ok
    unsigned natoms=9;

    p.cmd("setNatoms",natoms);
    output<<"ok natoms\n";
    p.cmd("setNatoms",std::move(natoms));
    output<<"ok std::move(natoms)\n";
    p.cmd("setNatoms",&natoms);
    output<<"ok &natoms\n";

    // preliminary tests on setNatoms:

    long long unsigned natoms_long=natoms;
    try {
      p.cmd("setNatoms",natoms_long);
      output<<"ko natoms_long\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught natoms_long\n";
    }

    try {
      p.cmd("setNatoms",std::move(natoms_long));
      output<<"ko std::move(natoms_long)\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught std::move(natoms_long)\n";
    }

    try {
      p.cmd("setNatoms",&natoms_long);
      output<<"ko &natoms_long\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught &natoms_long\n";
    }

    long long unsigned natoms_float=natoms;
    try {
      p.cmd("setNatoms",natoms_float);
      output<<"ko natoms_float\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught natoms_float\n";
    }

    try {
      p.cmd("setNatoms",std::move(natoms_float));
      output<<"ko std::move(natoms_float)\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught std::move(natoms_float)\n";
    }

    try {
      p.cmd("setNatoms",&natoms_float);
      output<<"ko &natoms_float\n";
    } catch(Plumed::ExceptionTypeError & e) {
      output<<"caught &natoms_float\n";
    }
  }
}

void test1(){
  std::ofstream output("output");

#include "main.inc"

}

void test_no_shapes(); // this is in main_no_shapes.cpp

int main() {
  preliminary();
  test1();
  test_no_shapes();
  return 0;
}
