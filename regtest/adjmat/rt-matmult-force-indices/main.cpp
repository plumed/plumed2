#include "plumed/matrixtools/MatrixTimesMatrix.h"
#include <cstdio>

namespace PLMD {
namespace sparse_matrix_force_test {

struct CV {};

bool check(bool empty_column) {
  using Product=matrixtools::MatrixTimesMatrix<CV>;
  Product::input_type data;
  data.secondMatrixIsSparse=true;
  data.outmat.ncols=1;
  data.outmat.resize(2);
  data.outmat[0]=1;
  data.outmat[1]=0;
  data.locations_v= {2,0,4,0,0,0,0,0,0};
  data.update();

  // The first matrix row has three entries; the second matrix column
  // intersects only two of them. Its stored entries map to offsets 9 and 13.
  std::size_t ranks[]= {2,2}, ncols[]= {3,2}, starts[]= {0,9};
  std::size_t shape[]= {3,3,3,3}, shape_starts[]= {0,2};
  std::size_t book_starts[]= {0,12}, book_sizes[]= {12,9};
  std::size_t book[]= {3,0,1,2,0,0,0,0,0,0,0,0,2,0,2,0,0,0,0,0,0};
  if(empty_column) {
    book[12]=0;
  }
  Product::ParallelActionsInput input;
  input.ncomponents=1;
  input.threadunsafe_forces_start=9;
  input.ranks=ranks;
  input.ncols=ncols;
  input.argstarts=starts;
  input.shapedata=shape;
  input.shapestarts=shape_starts;
  input.bookstarts=book_starts;
  input.booksizes=book_sizes;
  input.bookeeping=book;

  std::size_t indices[8]= {};
  auto holder=ForceIndexHolder::create(1,6,indices);
  Product::getForceIndices(0,0,15,data,input,holder);
  // Validate indices before using them: the old boundary classified unused
  // zero indices as second-operand indices, causing unsigned underflow.
  for(std::size_t i=0; i<holder.threadsafe_derivatives_end[0]; ++i) {
    if(holder.indices[0][i]>=9) {
      return false;
    }
  }
  for(std::size_t i=holder.threadsafe_derivatives_end[0]; i<holder.tot_indices[0]; ++i) {
    if(holder.indices[0][i]<9 || holder.indices[0][i]>=15) {
      return false;
    }
  }

  double output_force[]= {-1};
  double derivatives[]= {2,3,0,5,7,0};
  double first[15]= {}, second[6]= {};
  auto forces=ForcesInput<double>::create(1,output_force,6,derivatives);
  ParallelTaskManager<Product>::gatherThreadSafeForces(input,holder,forces,View<double>(first,15));
  ParallelTaskManager<Product>::gatherThreadUnsafeForces(input,holder,forces,View<double>(second,6));
  for(std::size_t i=0; i<15; ++i) {
    const double expected=empty_column ? 0 : (i==0 ? -2 : (i==2 ? -3 : 0));
    if(first[i]!=expected) {
      return false;
    }
  }
  for(std::size_t i=0; i<6; ++i) {
    const double expected=empty_column ? 0 : (i==0 ? -5 : (i==4 ? -7 : 0));
    if(second[i]!=expected) {
      return false;
    }
  }
  return true;
}

int run() {
  if(!check(false) || !check(true)) {
    std::fprintf(stderr,"invalid sparse matrix force gathering\n");
    return 1;
  }
  std::FILE* output=std::fopen("output","w");
  if(!output) {
    return 1;
  }
  std::fprintf(output,"sparse matrix force gathering passed\n");
  return std::fclose(output)!=0;
}

}
}

int main() {
  return PLMD::sparse_matrix_force_test::run();
}
