#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/View2D.h"

#include "testUtils.h"

#include <numeric>
#include <fstream>
#include <array>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <typeinfo>

using PLMD::View;
using PLMD::View2D;

#define displaycall(x) #x << " = " << (x)

void basics(tee& out);
void nonspanlikeInteractions(tee& out);

int main() {
  tee out("output");
  basics(out);

  return 0;
}

template<typename T, size_t N, size_t M>
void matrixDisplay (tee& out, View2D<T,N,M> v) {
  for(std::size_t i=0; i<v.size(); ++i) {
    for(std::size_t j=0; j<v[i].size(); ++j) {
      out << std::setw(2)<< v[i][j] << " ";
    }
    out << "\n";
  }
}

void basics(tee& out) {
  // View2D should be used as a way to passing data to function o to not owning classes
  // View2D can be used to "view" the data witha a 2D shape, or to access data as in a matrix

  std::array<double,24> data;
  std::iota(data.begin(), data.end(), 0.0);
  out << "Original data:\n";
  for (const auto& x: data)  {
    out << x <<" ";
  }
  out << "\n";

  //the three way of initializing a 2dView
  View2D<double,8,3> v(data.data());
  View2D<double> vv(data.data(),4,6);
  View2D<double,PLMD::helpers::dynamic_extent,8> vvv(data.data(),3);
  out << "8x3:\n";
  matrixDisplay(out,v);
  out << "4x6:\n";
  matrixDisplay(out,vv);
  out << "3x8:\n";
  matrixDisplay(out,vvv);

  //operator[] returns a View with same type of data, and the second dimension
  //For example you can use the View features to modify the data

  out << "v<4*,8> means that the view i flexible on the first dimension and fixed on the second\n";
  v[5]=PLMD::Vector{1,1,1};
  out << "data: (assigned {1,1,1} to v<8,3>[5])\n";
  for (const auto& x: data)  {
    out << x <<" ";
  }
  out << "\n";

  vv[3][4]=-15;
  out << "data: (assigned -15 to v<4*,6*>[3][5])\n";
  for (const auto& x: data)  {
    out << x <<" ";
  }
  out << "\n";

  vv[1][4]=-12;
  out << "data: (assigned -12 to v<3*,8>[1][4])\n";
  for (const auto& x: data)  {
    out << x <<" ";
  }
  out << "\n";
  //as no there is no iterator to test, yet
}

