#include "plumed/tools/Vector.h"
#include "plumed/tools/Angle.h"
#include "plumed/tools/Stopwatch.h"
#include <fstream>
#include <iostream>

using namespace PLMD;


int testCompute90Degrees(std::ofstream& ofs)
{
  Vector vector1(1.0,0.0,0.0);
  Vector vector2(0.0,1.0,0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs<<result<<"\n";
  return 0;
}

int testCompute45Degrees(std::ofstream& ofs)
{
  Vector vector1(1.0,0.0,0.0);
  Vector vector2(1.0,1.0,0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs<<result<<"\n";
  return 0;
}

int main(){
  Stopwatch sw;
  sw.start();
  std::ofstream ofs("output");

  testCompute90Degrees(ofs);
  testCompute45Degrees(ofs);

  sw.stop();
  std::cout<<sw;
  return 0;
}
