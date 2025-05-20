#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Random.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <limits>

using rng=PLMD::Random;
using namespace PLMD;

namespace {

static Tensor4d makeMatrix(const Tensor & rr01) {
  Tensor4d m;
  m[0][0]= +rr01[0][0]+rr01[1][1]+rr01[2][2];
  m[1][1]= +rr01[0][0]-rr01[1][1]-rr01[2][2];
  m[2][2]= -rr01[0][0]+rr01[1][1]-rr01[2][2];
  m[3][3]= -rr01[0][0]-rr01[1][1]+rr01[2][2];
  m[0][1]= +rr01[1][2]-rr01[2][1];
  m[0][2]= -rr01[0][2]+rr01[2][0];
  m[0][3]= +rr01[0][1]-rr01[1][0];
  m[1][2]= +rr01[0][1]+rr01[1][0];
  m[1][3]= +rr01[0][2]+rr01[2][0];
  m[2][3]= +rr01[1][2]+rr01[2][1];
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];

  return m;
}

void compute_quaternion_from_K(const Tensor4d& K, Vector4d& q, double& lambda_min) {

  // empirical:
  constexpr unsigned nsquare=32;
  constexpr unsigned niter=11;

  // Estimate upper bound for largest eigenvalue using Gershgorin disks (branch-free)
  double lambda_shift = 0.0;
  for (int i = 0; i < 4; ++i) {
    double row_sum = 0.0;
    for (int j = 0; j < 4; ++j) {
      double off_diag = (i == j) ? 0.0 : std::fabs(K[i][j]);
      row_sum += off_diag;
    }
    double bound = std::fabs(K[i][i]) + row_sum;
    lambda_shift = (bound > lambda_shift) ? bound : lambda_shift;
  }

  Tensor4d A=-K;

  A[0][0]+=lambda_shift;
  A[1][1]+=lambda_shift;
  A[2][2]+=lambda_shift;
  A[3][3]+=lambda_shift;

  for(unsigned i=0; i<nsquare; i++) {
    A=matmul(A,A);
    // Compute Frobenius norm of A
    double frob = 0.0;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        frob += A[i][j] * A[i][j];
      }
    frob = std::sqrt(frob);
    frob = (frob > 0.0) ? frob : 1.0;

    // Scale A in-place
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        A[i][j] /= frob;
      }
  }

  // Inverse iteration to find smallest eigenvalue of K
  Vector4d v = {1.0, 1.0, 1.0, 1.0};
  Vector4d z;

  for (int iter = 0; iter < niter; ++iter) {
    for (int i = 0; i < 4; ++i) {
      double sum = 0.0, c = 0.0;
      for (int j = 0; j < 4; ++j) {
        double y = A[i][j] * v[j] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
      }
      z[i] = sum;
    }
    double norm = std::sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2] + z[3]*z[3]);
    norm = (norm > 0.0) ? norm : 1.0;
    for (int i = 0; i < 4; ++i) {
      v[i] = z[i] / norm;
    }
  }

  // Normalize vector
  double norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  norm = (norm > 0.0) ? norm : 1.0;
  for (int i = 0; i < 4; ++i) {
    v[i] /= norm;
  }

  // Compute Rayleigh quotient for lambda_min
  lambda_min = 0.0;
  for (int i = 0; i < 4; ++i) {
    double sum = 0.0, c = 0.0;
    for (int j = 0; j < 4; ++j) {
      double y = K[i][j] * v[j] - c;
      double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    lambda_min += v[i] * sum;
  }

  // Store eigenvector
  norm = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);
  norm = (norm > 0.0) ? norm : 1.0;
  q = v/ norm;
}

} // namespace

//from google benchmark - to avoid compiler optimizations
template <class T>
inline void DoNotOptimize(T& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}



PLMD::Tensor randomRMSDmat(rng & gen, const double referencedistance=1.0) {
  PLMD::Tensor toret;
  constexpr size_t n=30;

  std::vector<double>  align(n);
  std::vector<Vector> positions(n);
  std::vector<Vector> reference(n);
  double alignSum=0.0;
  for (unsigned i = 0; i < n; i++) {
    align[i]=gen.RandU01();
    alignSum+=align[i];
    reference[i]=10.0*Vector{
      gen.RandU01(),
      gen.RandU01(),
      gen.RandU01()
    };
    positions[i]=2.0*referencedistance*Vector{
      gen.RandU01()-0.5,
      gen.RandU01()-0.5,
      gen.RandU01()-0.5
    } + reference[i];
  }

  Vector cpositions;
  cpositions.zero();
  for(unsigned iat=0; iat<n; iat++) {
    align[iat]/=alignSum;
    const double w=align[iat];
    cpositions+=positions[iat]*w;
  }

// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0; iat<n; iat++) {
    double w=align[iat];
    toret+=Tensor(positions[iat]-cpositions,reference[iat])*w;
  }
  return toret;
}


int main() {
  std::cerr <<std::setprecision(9) << std::fixed;

  std::vector<PLMD::Tensor> tests= {
    // {1,1,1,1,1,1,1,1,1}, // 'wrong' eigenvector selected (due to the degeneracy)
    {1,2,3,4,5,6,7,8,9},
    {0,1,1,1,1,1,1,1,1},
    {
      1.0, 1.0, 1.0,
      1.0, .0, 1.0,
      1.0, 1.0, 1.0
    },
    // {
    //   0.0, 1.0, 1.0,
    //   1.0, 0.0, 1.0,
    //   1.0, 1.0, 0.0
    // }, // 'wrong' eigenvector selected (due to the degeneracy)
    {
      1.0, 1.0, 0.0,
      1.0, 0.0, 1.0,
      0.0, 1.0, 1.0
    },
    {
      1.0, 5.0, 1.0,
      1.0, 5.0, 3.0,
      1.0, 5.0, 1.0
    }
  };

  rng gen{"MatrixGenerator"};
  gen.setSeed(123456789);
  double error=100.0;
  for(int i=0; i<400; i++) {
    tests.push_back(randomRMSDmat(gen,error));
    if(i%50) {
      error/=10.0;
    }
  }
  std::ofstream matOut ("matrices");
  std::ofstream out ("eigenValues");
  std::ofstream outv ("eigenVectors");
// setting up an header for clarity

  matOut<<"# id:";
  for (unsigned i=0; i<16; ++i) {
    matOut<<" "<<std::setw(9)<<"mat"+std::to_string(i);
    // outv<<" "<<std::setw(9)<<"mat"+std::to_string(i);
  }
  matOut << "\n";

  out <<"# id:\t";
  for (unsigned i=0; i<1; ++i) {
    out<<" "<<std::setw(10) <<"eigv"+std::to_string(i);
  }
  out << "\n";


  outv <<"# id:\t";
  for (unsigned i=0; i<4; ++i) {
    outv<<" "<<std::setw(10) <<"eigvct"+std::to_string(0)+"_"+std::to_string(i);
  }
  outv << "\n";

  out <<std::setprecision(5) << std::fixed;
  outv <<std::setprecision(5) << std::fixed;

  for (unsigned index=0; index<tests.size(); ++index) {
    const auto& mat = tests[index];
    auto m=makeMatrix(mat);

    Vector1d eigenval_ref;
    TensorGeneric<1,4> eigenvec_ref;
    diagMatSym(m, eigenval_ref, eigenvec_ref );

    double eigenval;
    Vector4d eigenvec;
    compute_quaternion_from_K(m,eigenvec,eigenval);

    matOut <<std::setw(4)<<index <<":\t";
    for (unsigned i=0; i<16; ++i) {
      //it would be lovely if `&mat[0][0]` became `mat.data()`
      matOut<<" "<<std::setw(9) <<(&m[0][0])[i];
    }
    matOut<<"\n";

    out <<std::setw(4)<<index <<":\t";
    out<<" "<<std::setw(10) <<eigenval_ref[0] - eigenval;
    out<<"\n";

    //The eigenvectors sometimes are returned with the direction inverted
    if(eigenvec_ref[0][0] * eigenvec[0] < 0.0 ) {
      eigenvec *=-1;
      //first element is 0
    } else if (eigenvec_ref[0][1] * eigenvec[1] < 0.0) {
      eigenvec *=-1;
      //also second element is 0
    } else if (eigenvec_ref[0][2] * eigenvec[2] < 0.0) {
      eigenvec *=-1;
      //also third element is 0
    } else if (eigenvec_ref[0][3] * eigenvec[3] < 0.0) {
      eigenvec *=-1;
    }

    outv <<std::setw(4)<<index <<":\t";
    for (unsigned i=0; i<4; ++i) {
      outv<<" "<<std::setw(10) <<eigenvec_ref[0][i] - eigenvec[i];
    }
    outv<<"\n";
  }

  //Timings
  {
    std::vector<PLMD::Tensor4d> matrices(tests.size());
    for(size_t i=0; i< tests.size(); ++i) {
      matrices[i] = makeMatrix(tests[i]);
    }

    Stopwatch sw;
    for(int i=0; i<10; ++i) {
      {
        auto sww = sw.startStop("diagMatSym_all");
        for (const auto & mat : matrices) {
          auto sww = sw.startStop("diagMatSym");

          Vector1d eigenvals;
          TensorGeneric<1,4> eigenvecs;

          diagMatSym(mat, eigenvals, eigenvecs );
          DoNotOptimize(eigenvals);
          DoNotOptimize(eigenvecs);
        }
      }

      {
        auto sww = sw.startStop("compute_quaternion_from_K_all");
        for (const auto & mat : matrices) {
          auto sww = sw.startStop("compute_quaternion_from_K");

          double eigenvals;
          Vector4d eigenvecs;

          compute_quaternion_from_K(mat, eigenvecs, eigenvals );
          DoNotOptimize(eigenvals);
          DoNotOptimize(eigenvecs);
        }
      }
    }
    std::cerr<<sw;
  }
  return 0;
}
