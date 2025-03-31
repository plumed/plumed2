#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Random.h"
#include "plumed/tools/EigenSolverForRMSD.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>

using rng=PLMD::Random;

//from google benchmark - to avoid compiler optimizations
template <class T>
inline void DoNotOptimize(T& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}

using namespace PLMD;

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
    {1,1,1,1,1,1,1,1,1},
    {1,2,3,4,5,6,7,8,9},
    {0,1,1,1,1,1,1,1,1},
//problems
    { 1.0, 1.0, 1.0,
      1.0, .0, 1.0,
      1.0, 1.0, 1.0
    },
//problems
    { 0.0, 1.0, 1.0,
      1.0, 0.0, 1.0,
      1.0, 1.0, 0.0
    },
//problems in eigenvalues
    { 1.0, 1.0, 0.0, \
      1.0, 0.0, 1.0, \
      0.0, 1.0, 1.0
    },
//sign in some eivgenvectors
    { 1.0, 5.0, 1.0, \
      1.0, 5.0, 3.0, \
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

  std::ofstream out ("eigenValues");
  std::ofstream outv ("eigenVectors");
// setting up an header for clarity
  out <<"# id:";
  outv <<"# id:";
  for (unsigned i=0; i<9; ++i) {
    out<<" "<<std::setw(9)<<"mat"+std::to_string(i);
    outv<<" "<<std::setw(9)<<"mat"+std::to_string(i);
  }

  out <<"\t";
  for (unsigned i=0; i<4; ++i) {
    out<<" "<<std::setw(10) <<"eigv"+std::to_string(i);
  }
  outv <<"\t";
  for (unsigned r=0; r<4; ++r) {
    outv<<"\t";
    for (unsigned i=0; i<4; ++i) {
      outv<<" "<<std::setw(10) <<"eigvct"+std::to_string(r)+"_"+std::to_string(i);
    }
  }

  out << "\n";
  outv << "\n";

  out <<std::setprecision(5) << std::fixed;
  outv <<std::setprecision(5) << std::fixed;

  for (unsigned index=0; index<tests.size(); ++index) {

    const auto& mat = tests[index];
    out <<std::setw(4)<<index <<":";
    outv <<std::setw(4)<<index <<":";
    for (unsigned i=0; i<9; ++i) {
      out<<" "<<std::setw(9) <<mat.data()[i];
      outv<<" "<<std::setw(9) <<mat.data()[i];
    }
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    //The reference is created by calling the internal blas/lapack interface
    // Vector4d eigenvals;
    // Tensor4d eigenvecs;
    // diagMatSym(m, eigenvals, eigenvecs );

    auto [eigenvals,eigenvecs] = EigenSolverForRMSD::calculate(Tensor(mat));

    out <<"\t";
    for (unsigned i=0; i<4; ++i) {
      out<<" "<<std::setw(10) <<eigenvals[i];
    }
    outv <<"\t";
    for (unsigned r=0; r<4; ++r) {
      outv<<"\t";
      for (unsigned i=0; i<4; ++i) {
        outv<<" "<<std::setw(10) <<eigenvecs.getRow(r)[i];
      }
    }
    out<<"\n";
    outv<<"\n";

  }


  //timing
  Stopwatch sw;
  {
    auto sww = sw.startStop("diagMatSym_all");
    for (const auto & mat : tests) {
      auto sww = sw.startStop("diagMatSym");

      auto m=EigenSolverForRMSD::makeMatrix(mat);
      Vector4d eigenvals;
      Tensor4d eigenvecs;

      diagMatSym(m, eigenvals, eigenvecs );
      DoNotOptimize(eigenvals);
      DoNotOptimize(eigenvecs);
    }
  }

  {
    auto sww = sw.startStop("mysolver_all");
    for (const auto & mat : tests) {
      auto sww = sw.startStop("mysolver");

      auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));
      DoNotOptimize(ev);
      DoNotOptimize(evv);
    }
  }

//output
  {
    const auto& mat = tests[6];
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    std::cerr<<"->"<<eigenvals<<"\n";
    // out <<"->"<<eigenvals<<"\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<eigenvecs.getRow(i)<<"\n";
      // out<<eigenvecs.getRow(i)<<"\n";
    }
  }

  {

    const auto& mat = tests[6];
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));
    std::cerr<<"->"<<ev<<"\n";
    // out<<"->"<<ev<<"\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<evv.getRow(i)<<"\n";
      // out<<evv.getRow(i)<<"\n";
    }
  }
  std::cerr<<"\n";
  { int index=1;

    const auto& mat = tests[index];
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));

    std::cerr<<"->"<<eigenvals<<"\n"
             <<"-<"<<ev<<"\n\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<">"<<eigenvecs.getRow(i)<<"\n"
               <<"<"<<evv.getRow(i)<<"\n"
               <<"\n";
    }
  }
  std::cerr<<"\n";
  {

    auto mymat =randomRMSDmat(gen);
    const auto& mat = mymat;
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));

    std::cerr<<"->"<<eigenvals<<"\n"
             <<"-<"<<ev<<"\n\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<">"<<eigenvecs.getRow(i)<<"\n"
               <<"<"<<evv.getRow(i)<<"\n"
               <<"\n";
    }
  }

  {

    auto mymat =randomRMSDmat(gen,1000.0);
    const auto& mat = mymat;
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));

    std::cerr<<"->"<<eigenvals<<"\n"
             <<"-<"<<ev<<"\n\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<">"<<eigenvecs.getRow(i)<<"\n"
               <<"<"<<evv.getRow(i)<<"\n"
               <<"\n";
    }
  }

  {

    auto mymat =randomRMSDmat(gen,1e-8);
    const auto& mat = mymat;
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = EigenSolverForRMSD::calculate(Tensor(mat));

    std::cerr<<"->"<<eigenvals<<"\n"
             <<"-<"<<ev<<"\n\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<">"<<eigenvecs.getRow(i)<<"\n"
               <<"<"<<evv.getRow(i)<<"\n"
               <<"\n";
    }
  }


  {
std::cerr <<"prefactor -2\n";
    auto mymat =randomRMSDmat(gen,1000.0);
    const auto& mat = mymat;
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=-2.0*EigenSolverForRMSD::makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = EigenSolverForRMSD::calculate(-2.0*Tensor(mat));

    std::cerr<<"->"<<eigenvals<<"\n"
             <<"-<"<<ev<<"\n\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<">"<<eigenvecs.getRow(i)<<"\n"
               <<"<"<<evv.getRow(i)<<"\n"
               <<"\n";
    }
  }

  std::cerr<<sw;

  return 0;
}
