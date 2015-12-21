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
  Vector b=matmul(A,a);

  Tensor B(1.0,2.0,3.0,5.0,4.0,3.0,10.0,8.0,2.0);
  Vector c=matmul(a,B,inverse(B));

  ofs<<a<<"\n";
  ofs<<b<<"\n";
  ofs<<c<<"\n";

  B-=A;
  ofs<<determinant(B)<<"\n";

  TensorGeneric<3,2> D(Vector(0,1,2),Vector2d(3,4));
  ofs<<D<<"\n";

  TensorGeneric<2,3> E=transpose(D);
  ofs<<E<<"\n";

  double f(matmul(a,B,c));
  double f1(dotProduct(a,matmul(B,c)));
  double f2(matmul(a,matmul(B,c)));
  ofs<<f<<" "<<f1<<" "<<f2<<"\n";

  sw.stop();
  std::cout<<sw;
  
  return 0;
}
