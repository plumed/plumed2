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
