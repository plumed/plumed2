#include "plumed/tools/Angle.h"
#include "plumed/tools/Vector.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

template<typename T>
using myVector=PLMD::VectorTyped<T,3>;
//by removing the \n the diff becames clearer
template<typename T=double>
int testCompute2Vectors90Degrees(std::ofstream &ofs) {
  ofs << "testCompute2Vectors90Degrees ";
  myVector<T> vector1(1.0, 0.0, 0.0);
  myVector<T> vector2(0.0, 1.0, 0.0);
  Angle angle;
  auto result = angle.compute(vector1, vector2);
  ofs << result << std::endl;
  return 0;
}

template<typename T=double>
int testCompute2Vectors45Degrees(std::ofstream &ofs) {
  ofs << "testCompute2Vectors45Degrees ";
  myVector<T> vector1(1.0, 0.0, 0.0);
  myVector<T> vector2(1.0, 1.0, 0.0);
  Angle angle;
  auto result = angle.compute(vector1, vector2);
  ofs << result << std::endl;
  return 0;
}

template<typename T=double>
int testCompute4Vectors45Degrees(std::ofstream &ofs) {
  ofs << "testCompute4Vectors45Degrees ";
  myVector<T> vector1(1.0, 1.0, 0.0);
  myVector<T> vector2(1.0, 0.0, 0.0);
  myVector<T> d1, d2;
  Angle angle;
  auto result = angle.compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << "testCompute4Vectors45Degrees d1= " << d1 << std::endl;
  ofs << "testCompute4Vectors45Degrees d2= " << d2 << std::endl;
  return 0;
}

template<typename T=double>
int testCompute4VectorsParallel(std::ofstream &ofs) {
  ofs << "testCompute4VectorsParallel ";
  myVector<T> vector1(1.0, 1.0, 0.0);
  myVector<T> vector2(1.0, 1.0, 0.0);
  myVector<T> d1, d2;
  auto result = Angle::compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << "testCompute4VectorsParallel d1= " << d1 << std::endl;
  ofs << "testCompute4VectorsParallel d2= " << d2 << std::endl;
  return 0;
}
template<typename T=double>
int testCompute4VectorsAntiParallel(std::ofstream &ofs) {
  ofs << "testCompute4VectorsAntiParallel ";
  myVector<T> vector1(1.0, 1.0, 0.0);
  myVector<T> vector2(-1.0, -1.0, 0.0);
  myVector<T> d1, d2;
  Angle angle;
  auto result = angle.compute(vector1, vector2, d1, d2);
  ofs << result << std::endl;
  ofs << "testCompute4VectorsAntiParallel d1= " << d1 << std::endl;
  ofs << "testCompute4VectorsAntiParallel d2= " << d2 << std::endl;
  return 0;
}

int main() {
  std::ofstream ofs("output");
  testCompute2Vectors90Degrees<double>(ofs);
  testCompute2Vectors45Degrees<double>(ofs);
  testCompute4Vectors45Degrees<double>(ofs);
  testCompute4VectorsParallel<double>(ofs);
  testCompute4VectorsAntiParallel<double>(ofs);
  ofs << "float\n";
  testCompute2Vectors90Degrees<float>(ofs);
  testCompute2Vectors45Degrees<float>(ofs);
  testCompute4Vectors45Degrees<float>(ofs);
  testCompute4VectorsParallel<float>(ofs);
  testCompute4VectorsAntiParallel<float>(ofs);
  return 0;
}
