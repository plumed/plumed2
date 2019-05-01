#include "plumed/blas/blas.h"
#include "plumed/lapack/lapack.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace PLMD;

/*
  This test is designed to check whether the interface we
  provide for blas/lapack is compatible with the actually linked
  blas/lapack. Notice that the interface we provide is taken from
  gromacs and is compatible with fortran blas/lapack.
*/

int main(){
  std::ofstream ofs("output");
  int size=5;
  std::vector<double> ad(size),bd(size);
  std::vector<float>  af(size),bf(size);
  for(unsigned i=0;i<size;i++){
    af[i]=ad[i]=size;
    bf[i]=bd[i]=size-i;
  }
  int inca=1;
  int incb=1;
  ofs<<plumed_blas_ddot(&size,&ad[0],&inca,&bd[0],&incb)<<"\n";
  ofs<<plumed_blas_sdot(&size,&af[0],&inca,&bf[0],&incb)<<"\n";
  ofs<<plumed_lapack_dlapy2(&ad[0],&bd[0])<<"\n";
  ofs<<plumed_lapack_slapy2(&af[0],&bf[0])<<"\n";
  return 0;
}
