#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include "plumed/tools/Stopwatch.h"
#include <fstream>
#include <iostream>

using namespace PLMD;

int main(){
  Stopwatch sw;
  sw.start();
  Vector a(1.0,2.0,3.0);
  Tensor A(Tensor::identity());
  std::ofstream ofs("output");
  Vector b=(matmul(A,a));

  Tensor B(1.0,2.0,3.0,5.0,4.0,3.0,10.0,8.0,2.0);
  Vector c=(matmul(a,matmul(B,inverse(B))));

  ofs<<a[0]<<" "<<a[1]<<" "<<a[2]<<"\n";
  ofs<<b[0]<<" "<<b[1]<<" "<<b[2]<<"\n";
  ofs<<c[0]<<" "<<c[1]<<" "<<c[2]<<"\n";

  B-=A;
  ofs<<determinant(B)<<"\n";

  TensorGeneric<3,2> D(Vector(0,1,2),Vector2d(3,4));
  ofs<<D[0][0]<<" "<<D[0][1]<<"\n";
  ofs<<D[1][0]<<" "<<D[1][1]<<"\n";
  ofs<<D[2][0]<<" "<<D[2][1]<<"\n";

  TensorGeneric<2,3> E=transpose(D);
  ofs<<E[0][0]<<" "<<E[0][1]<<" "<<E[0][2]<<"\n";
  ofs<<E[1][0]<<" "<<E[1][1]<<" "<<E[1][2]<<"\n";

  sw.stop();
  std::cout<<sw;
  
  return 0;
}
