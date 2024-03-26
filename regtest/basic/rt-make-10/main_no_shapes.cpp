/*
This is a separate version where shape detection is disabled.
It is here to test what would happen with an old compiler that supports C++11
but does not support proper disambiguation of overloads in Plumed.h,
(gcc<6 and icc < 17), or to test what happens if a user explicitly disable
shape detection. Notice that some of the tests are then switched off
using #ifdefs in main.inc
*/
#define __PLUMED_WRAPPER_CXX_DETECT_SHAPES 0
#include "plumed/wrapper/Plumed.h"
#include <vector>
#include <fstream>
#include <cstdio>
#include <iostream>

using namespace PLMD;

void test_no_shapes(){
  std::ofstream output("output_no_shapes");

#include "main.inc"

}
