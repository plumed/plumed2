#include "plumed/tools/Vector.h"
#include "plumed/tools/Angle.h"
#include "plumed/tools/Stopwatch.h"
#include <fstream>
#include <iostream>

using namespace PLMD;


int testCompute2Vectors90Degrees(std::ofstream& ofs)
{
  Vector vector1(1.0,0.0,0.0);
  Vector vector2(0.0,1.0,0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs<<result<<"\n";
  return 0;
}

int testCompute2Vectors45Degrees(std::ofstream& ofs)
{
  Vector vector1(1.0,0.0,0.0);
  Vector vector2(1.0,1.0,0.0);
  Angle angle;
  double result = angle.compute(vector1, vector2);
  ofs<<result<<"\n";
  return 0;
}

int testCompute4Vectors45Degrees(std::ofstream& ofs)
{
  Vector vector1(1.0,0.0,0.0);
  Vector vector2(1.0,1.0,0.0);
  Vector d1,d2; Angle angle;
  double result = angle.compute(vector1, vector2, d1, d2);
  ofs<<result<<"\n";
  ofs<< d1[0] << " " << d1[1] << " " << d1[2] << "\n";
  ofs<< d2[0] << " " << d2[1] << " " << d2[2] << "\n";
  return 0;
}



int main(){
  Stopwatch sw;
  sw.start();
  std::ofstream ofs("output");

  testCompute2Vectors90Degrees(ofs);
  testCompute2Vectors45Degrees(ofs);

  testCompute4Vectors45Degrees(ofs);

  sw.stop();
  std::cout<<sw;
  return 0;
}
