#include "plumed/tools/Vector.h"
#include "plumed/tools/Tensor.h"
#include "plumed/tools/Stopwatch.h"
#include "plumed/tools/Random.h"

#include <complex>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>

using rng=PLMD::Random;

//form google benchmark
template <class T>
inline void DoNotOptimize(T& value) {
#if defined(__clang__)
  asm volatile("" : "+r,m"(value) : : "memory");
#else
  asm volatile("" : "+m,r"(value) : : "memory");
#endif
}

#define vdbg(x) std::cerr <<__LINE__ <<": " #x << " = " << (x) << std::endl

using namespace PLMD;


PLMD::Vector4d eigen(double d2, double d1, double d0)
{
  using std::complex;
  // CForm[Simplify[Solve[l^4 + d2 l^2 + d1 l + d0 == 0, l]]]
  const complex<double> twoTothird = pow(2, 1.0 / 3.0);
  // const double twoTotwothird = twoTothird * twoTothird; // pow(2, 2.0 / 3.0);
  const complex<double> twoSqrt6 = 2. * sqrt(6);
  const complex<double> d1sq = d1 * d1;
  const complex<double> d2sq = d2 * d2;
  const complex<double> d2cb = d2sq * d2;
  const complex<double> P1 = 12.0 * d0 + d2sq;
  const complex<double> P1cb = P1 * P1 * P1;
  const complex<double> H = 27.0 * d1sq - 72.0 * d0 * d2 + 2.0 * d2cb;
  // vdbg(pow(H, 2) - 4.0 * P1cb);
  const complex<double> F1 = H + sqrt(pow(H, 2) - 4.0 * P1cb);
  const complex<double> D1 = pow(F1, 1.0 / 3.0);
  const complex<double> T1 = twoTothird * (2.0 * P1 / D1 + twoTothird * D1);
  // vdbg(T1 - 4.0 * d2);
  const complex<double> r1 = sqrt(T1 - 4.0 * d2);
  // vdbg(-8.0 * d2 - T1 - (6.0 * twoSqrt6 * d1) / r1);
  const complex<double> r2m = sqrt(-8.0 * d2 - T1 - (6.0 * twoSqrt6 * d1) / r1);
  const complex<double> r2p = sqrt(-8.0 * d2 - T1 + (6.0 * twoSqrt6 * d1) / r1);

  vdbg(-(r1 + r2p)); //
  vdbg((-r1 + r2p)); //
  vdbg((r1 - r2m) );  //
  vdbg((r1 + r2m) );
  vdbg(r1);
  PLMD::Vector4d toret{
    std::real(-(r1 + r2p) / (twoSqrt6)), //
    std::real((-r1 + r2p) / (twoSqrt6)), //
    std::real((r1 - r2m) / (twoSqrt6)),  //
    std::real((r1 + r2m) / (twoSqrt6))

    //
  };
  // std::sort(&toret[0], &toret[0]+4);
  return toret;
}

//in the code the matrix used is -2m of this:
Tensor4d makeMatrix(const Tensor & rr01) {
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

PLMD::Tensor randomatSmall(rng & gen) {
  PLMD::Tensor toret;
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      toret[i][j] = gen.RandU01()*2.0-1.0;
    }
  }
  return toret;
}

PLMD::Tensor randomat(rng & gen) {
  PLMD::Tensor toret;
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      toret[i][j] = gen.RandU01()*20.0-10.0;
    }
  }
  return toret;
}
class EigenVectors
{
  // double v4_0;
  // double v4_1;
  // double v4_2;

  double v1_0;
  double v1_1;
  double v1_2;

  // v1:
  // -pow(l,3) + pow(r00,3) + pow(r01,2)*r11 - pow(r02,2)*r11 + 2*r01*r10*r11 + pow(r10,2)*r11 +
  // pow(r11,3) + 2*r01*r02*r12 + 2*r02*r10*r12 + r11*pow(r12,2) - 2*r02*r11*r20 + 2*r01*r12*r20 +
  // 2*r10*r12*r20 - r11*pow(r20,2) + 2*r01*r02*r21 + 2*r02*r10*r21 + 2*r11*r12*r21 + 2*r01*r20*r21 +
  // 2*r10*r20*r21 + r11*pow(r21,2) - pow(r01,2)*r22 + pow(r02,2)*r22 - 2*r01*r10*r22 - pow(r10,2)*r22 -
  // pow(r11,2)*r22 + pow(r12,2)*r22 + 2*r02*r20*r22 + pow(r20,2)*r22 + 2*r12*r21*r22 + pow(r21,2)*r22 -
  // r11*pow(r22,2) + pow(r22,3) - pow(r00,2)*(r11 + r22) - pow(l,2)*(r00 + r11 + r22) +
  // r00*(pow(r01,2) + pow(r02,2) + 2*r01*r10 + pow(r10,2) - pow(r11,2) - pow(r12,2) + 2*r02*r20 +
  // pow(r20,2) - 2*r12*r21 - pow(r21,2) + 2*r11*r22 - pow(r22,2)) +
  // l*(pow(r00,2) + pow(r01,2) + pow(r02,2) + 2*r01*r10 + pow(r10,2) + pow(r11,2) + pow(r12,2) +
  // 2*r02*r20 + pow(r20,2) + 2*r12*r21 + pow(r21,2) - 2*r11*r22 + pow(r22,2) - 2*r00*(r11 + r22))

  double v2_0;
  double v2_1;
  double v2_2;

  // v2:
  // 2*r01*r02*r11 - pow(r01,2)*r12 + pow(r02,2)*r12 + pow(r10,2)*r12 + pow(r11,2)*r12 + pow(r12,3) -
  // 2*r10*r11*r20 - r12*pow(r20,2) + 2*r00*(r02*r10 - r01*r20) - pow(r01,2)*r21 + pow(r02,2)*r21 +
  // pow(r10,2)*r21 - pow(r11,2)*r21 + pow(r12,2)*r21 - pow(r20,2)*r21 - r12*pow(r21,2) - pow(r21,3) +
  // pow(l,2)*(-r12 + r21) + pow(r00,2)*(-r12 + r21) + 2*l*(r02*r10 - r00*r12 - r01*r20 + r00*r21) -
  // 2*r01*r02*r22 - 2*r11*r12*r22 + 2*r10*r20*r22 + 2*r11*r21*r22 + r12*pow(r22,2) - r21*pow(r22,2)

  double v3_0;
  double v3_1;
  double v3_2;

  //v3:
  // -(pow(r01,2)*r02) - pow(r02,3) + r02*pow(r10,2) + r02*pow(r11,2) - 2*r01*r11*r12 - r02*pow(r12,2) +
  // pow(l,2)*(r02 - r20) - pow(r01,2)*r20 - pow(r02,2)*r20 + pow(r10,2)*r20 - pow(r11,2)*r20 -
  // pow(r12,2)*r20 + r02*pow(r20,2) + pow(r20,3) + pow(r00,2)*(-r02 + r20) + 2*r10*r11*r21 +
  // r02*pow(r21,2) + r20*pow(r21,2) + 2*l*(r02*r11 - r01*r12 - r11*r20 + r10*r21) + 2*r10*r12*r22 -
  // 2*r01*r21*r22 - r02*pow(r22,2) + r20*pow(r22,2) - 2*r00*(r10*r12 - r01*r21 - r02*r22 + r20*r22)



  double v4_0;
  double v4_1;
  double v4_2;
  //v4:
  // pow(r01,3) + r01*pow(r02,2) + pow(r00,2)*(r01 - r10) + pow(r01,2)*r10 + pow(r02,2)*r10 -
  // r01*pow(r10,2) - pow(r10,3) + pow(l,2)*(-r01 + r10) + r01*pow(r11,2) - r10*pow(r11,2) +
  // 2*r02*r11*r12 - r01*pow(r12,2) - r10*pow(r12,2) - r01*pow(r20,2) - r10*pow(r20,2) - 2*r11*r20*r21 +
  // r01*pow(r21,2) + r10*pow(r21,2) - 2*r00*(r01*r11 - r10*r11 + r02*r12 - r20*r21) - 2*r12*r20*r22 +
  // 2*r02*r21*r22 - r01*pow(r22,2) + r10*pow(r22,2) - 2*l*(r12*r20 - r02*r21 + r01*r22 - r10*r22)

public:
  EigenVectors(double r00, double r01, double r02,  //
               double r10, double r11, double r12, //
               double r20, double r21, double r22)
    :
    // v1_0(
    //   pow(r00,3) + pow(r01,2)*r11 - pow(r02,2)*r11 + 2*r01*r10*r11 + pow(r10,2)*r11 +
    //   pow(r11,3) + 2*r01*r02*r12 + 2*r02*r10*r12 + r11*pow(r12,2) - 2*r02*r11*r20 + 2*r01*r12*r20 +
    //   2*r10*r12*r20 - r11*pow(r20,2) + 2*r01*r02*r21 + 2*r02*r10*r21 + 2*r11*r12*r21 + 2*r01*r20*r21 +
    //   2*r10*r20*r21 + r11*pow(r21,2) - pow(r01,2)*r22 + pow(r02,2)*r22 - 2*r01*r10*r22 - pow(r10,2)*r22 -
    //   pow(r11,2)*r22 + pow(r12,2)*r22 + 2*r02*r20*r22 + pow(r20,2)*r22 + 2*r12*r21*r22 + pow(r21,2)*r22 -
    //   r11*pow(r22,2) + pow(r22,3) - pow(r00,2)*(r11 + r22) +
    //   r00*(pow(r01,2) + pow(r02,2) + 2*r01*r10 + pow(r10,2) - pow(r11,2) - pow(r12,2) + 2*r02*r20 +
    //        pow(r20,2) - 2*r12*r21 - pow(r21,2) + 2*r11*r22 - pow(r22,2))
    // ),
    // v1_1 (pow(r00,2) + pow(r01,2) + pow(r02,2) + 2*r01*r10 + pow(r10,2) + pow(r11,2) + pow(r12,2) +
    //       2*r02*r20 + pow(r20,2) + 2*r12*r21 + pow(r21,2) - 2*r11*r22 + pow(r22,2) - 2*r00*(r11 + r22)),
    // v1_2(- (r00 + r11 + r22)),
    v2_0(
      2*r01*r02*r11 - pow(r01,2)*r12 + pow(r02,2)*r12 + pow(r10,2)*r12 + pow(r11,2)*r12 + pow(r12,3) -
      2*r10*r11*r20 - r12*pow(r20,2) + 2*r00*(r02*r10 - r01*r20) - pow(r01,2)*r21 + pow(r02,2)*r21 +
      pow(r10,2)*r21 - pow(r11,2)*r21 + pow(r12,2)*r21 - pow(r20,2)*r21 - r12*pow(r21,2) - pow(r21,3) + pow(r00,2)*(-r12 + r21) -
      2*r01*r02*r22 - 2*r11*r12*r22 + 2*r10*r20*r22 + 2*r11*r21*r22 + r12*pow(r22,2) - r21*pow(r22,2)
    ),
    v2_1(2*(r02*r10 - r00*r12 - r01*r20 + r00*r21)
        ),
    v2_2(
      (-r12 + r21)
    ),
    v3_0(-(pow(r01,2)*r02) - pow(r02,3) + r02*pow(r10,2) + r02*pow(r11,2) - 2*r01*r11*r12 - r02*pow(r12,2)
         - pow(r01,2)*r20 - pow(r02,2)*r20 + pow(r10,2)*r20 - pow(r11,2)*r20 -
         pow(r12,2)*r20 + r02*pow(r20,2) + pow(r20,3) + pow(r00,2)*(-r02 + r20) + 2*r10*r11*r21 +
         r02*pow(r21,2) + r20*pow(r21,2)  + 2*r10*r12*r22 -
         2*r01*r21*r22 - r02*pow(r22,2) + r20*pow(r22,2) - 2*r00*(r10*r12 - r01*r21 - r02*r22 + r20*r22)),
    v3_1(2*(r02*r11 - r01*r12 - r11*r20 + r10*r21) ),
    v3_2(r02 - r20),

    v4_0(
      pow(r01,3) + r01*pow(r02,2) + pow(r00,2)*(r01 - r10) + pow(r01,2)*r10 + pow(r02,2)*r10 -
      r01*pow(r10,2) - pow(r10,3)  + r01*pow(r11,2) - r10*pow(r11,2) +
      2*r02*r11*r12 - r01*pow(r12,2) - r10*pow(r12,2) - r01*pow(r20,2) - r10*pow(r20,2) - 2*r11*r20*r21 +
      r01*pow(r21,2) + r10*pow(r21,2) - 2*r00*(r01*r11 - r10*r11 + r02*r12 - r20*r21) - 2*r12*r20*r22 +
      2*r02*r21*r22 - r01*pow(r22,2) + r10*pow(r22,2)  ),
    v4_1(-2.0*(r12*r20 - r02*r21 + r01*r22 - r10*r22)),
    v4_2(-r01 + r10)
  {
    v1_0=
      (r00)*(r00)*(r00) + (r11)*(r11)*(r11) + 2*(r01 + r10)*(r02 + r20)*(r12 + r21) +
      r00*((r01 + r10)*(r01 + r10) + (r02 + r20)*(r02 + r20) - (r12 + r21)*(r12 + r21)) +
      r11*((r01 + r10)*(r01 + r10) - (r02 + r20)*(r02 + r20) + (r12 + r21)*(r12 + r21)) + (r00)*(r00)*(-r11 - r22) -
      r00*(r11 - r22)*(r11 - r22) + (-(r01 + r10)*(r01 + r10) + (r02 + r20)*(r02 + r20) + (r12 + r21)*(r12 + r21))*r22 + (r22)*(r22)*(r22) -
      r11*r22*(r11 + r22)
      ;
    v1_1 =pow(r00,2) + pow(r01,2) + pow(r02,2) + 2*r01*r10 + pow(r10,2) + pow(r11,2) + pow(r12,2) +
          2*r02*r20 + pow(r20,2) + 2*r12*r21 + pow(r21,2) - 2*r11*r22 + pow(r22,2) - 2*r00*(r11 + r22);
    v1_2=- (r00 + r11 + r22);
    vdbg(v1_0);
    vdbg(v1_1);
    vdbg(v1_2);
    vdbg(v2_0);
    vdbg(v2_1);
    vdbg(v2_2);
    vdbg(v3_0);
    vdbg(v3_1);
    vdbg(v3_2);
    vdbg(v4_0);
    vdbg(v4_1);
    vdbg(v4_2);
  }
  PLMD::Vector4d calculate(double l)
  {
    vdbg(l);
    vdbg(v1_0);
    vdbg(l*v1_1);
    vdbg(l* l*v1_2);
    vdbg(l* l*l);
    vdbg(v1_0 + l*v1_1 +l* l*v1_2-l* l*l);
    vdbg(v2_0);
    vdbg(l*v2_1);
    vdbg(l*l*v2_2);
    vdbg(v2_0 +l*v2_1 + l*l*v2_2);
    return {
      v1_0 +l*(v1_1 + l*(v1_2 -l)),
      v2_0 +l*(v2_1 + l*v2_2),
      v3_0 +l*(v3_1 + l*v3_2),
      v4_0 +l*(v4_1 + l*v4_2)

    };
  }
};

std::pair<PLMD::Vector4d, PLMD::Tensor4d> mysolver(PLMD::Tensor in)
{

  vdbg(in);
  auto toEigen = makeMatrix(in);
  double d2 = -2 * std::accumulate(in.data(), in.data() + 9, 0.0,
  [](double a, double b) { return a + b * b; });
  double d1 = -8 * in.determinant();
  double d0 = toEigen.determinant();
  vdbg(d2);
  vdbg(d1);
  vdbg(d0);
  auto ev = eigen(d2, d1, d0);
  auto evv = PLMD::Tensor4d{};
  auto evs = EigenVectors(in(0, 0), in(0, 1), in(0, 2), //
                          in(1, 0), in(1, 1), in(1, 2), //
                          in(2, 0), in(2, 1), in(2, 2));
  for (int i = 0; i < 4; i++)
  {
    auto t = evs.calculate(ev[i]);
    if (t[0] < 0) {
      t *= -1;
    }
    vdbg(t);
    evv.setRow(i, t
               / t.modulo()
              );
  }
  return {ev, evv};
}

int main() {
  std::cerr <<std::setprecision(9) << std::fixed;
  std::ofstream out ("result");
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
  for(int i=0; i<200; i++) {
    tests.push_back(randomatSmall(gen));
  }
  for(int i=0; i<200; i++) {
    tests.push_back(randomat(gen));
  }

  //timing
  Stopwatch sw;
  // {
  //   // auto sww = sw.startStop("diagMatSym");
  //   for (const auto & mat : tests) {
  //     auto sww = sw.startStop("diagMatSym");

  //     auto m=makeMatrix(mat);
  //     Vector4d eigenvals;
  //     Tensor4d eigenvecs;

  //     diagMatSym(m, eigenvals, eigenvecs );
  //     DoNotOptimize(eigenvals);
  //     DoNotOptimize(eigenvecs);
  //   }
  // }

  // {
  //   // auto sww = sw.startStop("mysolver");
  //   for (const auto & mat : tests) {
  //     auto sww = sw.startStop("mysolver");

  //     auto [ev,evv] = mysolver(Tensor(mat));
  //     DoNotOptimize(ev);
  //     DoNotOptimize(evv);
  //   }
  // }

//output
  {
    const auto& mat = tests[6];
    auto m=makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    std::cerr<<"->"<<eigenvals<<"\n";
    out <<"->"<<eigenvals<<"\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<eigenvecs.getRow(i)<<"\n";
      out<<eigenvecs.getRow(i)<<"\n";
    }
  }

  {

    const auto& mat = tests[6];
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto [ev,evv] = mysolver(Tensor(mat));
    std::cerr<<"->"<<ev<<"\n";
    out<<"->"<<ev<<"\n";
    for (unsigned i=0; i<4; ++i) {
      std::cerr<<evv.getRow(i)<<"\n";
      out<<evv.getRow(i)<<"\n";
    }
  }
  std::cerr<<"\n";
  { int index=1;

    const auto& mat = tests[index];
    std::cerr<<"mat:\n"<<mat<<"\n";
    auto m=makeMatrix(mat);
    Vector4d eigenvals;
    Tensor4d eigenvecs;
    diagMatSym(m, eigenvals, eigenvecs );
    auto [ev,evv] = mysolver(Tensor(mat));

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
