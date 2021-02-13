/*

The MIT License (MIT)

Copyright (c) 2015 Jacques-Henri Jourdan <jourgun@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef SIMD_MATH_PRIMS_H
#define SIMD_MATH_PRIMS_H

#include<math.h>
#include<stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Workaround a lack of optimization in gcc */
float exp_cst1_f = 2139095040.f;
float exp_cst2_f = 0.f;

/* Relative error bounded by 1e-5 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error */
inline float expapprox(float val) {
  union { int32_t i; float f; } xu, xu2;
  float val2, val3, val4, b;
  int32_t val4i;
  val2 = 12102203.1615614f*val+1065353216.f;
  val3 = val2 < exp_cst1_f ? val2 : exp_cst1_f;
  val4 = val3 > exp_cst2_f ? val3 : exp_cst2_f;
  val4i = (int32_t) val4;
  xu.i = val4i & 0x7F800000;
  xu2.i = (val4i & 0x7FFFFF) | 0x3F800000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x|],
               [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.509871020f + b * (0.312146713f + b * (0.166617139f + b *
              (-2.190619930e-3f + b * 1.3555747234e-2f))));
}

/* Same code with better precision, less performance and using double */

double exp_cst1_d = 9218868437227405312.;
double exp_cst2_d = 0.;

/* Relative error bounded by 3e-9 for normalized outputs
   Returns invalid outputs for nan inputs
   Continuous error

   Vectorizable only with AVX512dq extensions because of the
   double->int64 cast. On GCC, use option -mavx512dq. */
inline double expapprox_d(double val) {
  union { int64_t i; double f; } xu, xu2;
  double val2, val3, val4, b;
  int64_t val4i;
  val2 = 6497320848556798.092*val+4607182418800017408.;
  val3 = val2 < exp_cst1_d ? val2 : exp_cst1_d;
  val4 = val3 > exp_cst2_d ? val3 : exp_cst2_d;
  val4i = (int64_t) val4;
  xu.i = val4i & 0x7FF0000000000000;
  xu2.i = (val4i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
  b = xu2.f;

  /* Generated in Sollya with:
     > f=remez(1-x*exp(-(x-1)*log(2)),
               [|(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
                 (x-1)*(x-2)*x*x*x, (x-1)*(x-2)*x*x*x*x|],
                 [1.000001,1.999999], exp(-(x-1)*log(2)));
     > plot(exp((x-1)*log(2))/(f+x)-1, [1,2]);
     > f+x;
  */
  return
    xu.f * (0.5002494548377929861 + b * (0.3453457447382168695 + b *
            (0.1226618159032020501 + b * (2.4869768911673215212e-2 + b *
             (6.7148791796145116483e-3 + b * (-5.8813825553864185693e-5 + b *
              2.17150255054231565039e-4))))));
}

/* Absolute error bounded by 1e-5 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error. */
inline float logapprox(float val) {
  union { float f; int32_t i; } valu;
  float exp, addcst, x;
  valu.f = val;
  exp = valu.i >> 23;
  /* -89.970756366f = -127 * log(2) + constant term of polynomial bellow. */
  addcst = val > 0 ? -89.970756366f : -(float)INFINITY;
  valu.i = (valu.i & 0x7FFFFF) | 0x3F800000;
  x = valu.f;


  /* Generated in Sollya using:
    > f = remez(log(x)-(x-1)*log(2),
            [|1,(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
              (x-1)*(x-2)*x*x*x|], [1,2], 1, 1e-8);
    > plot(f+(x-1)*log(2)-log(x), [1,2]);
    > f+(x-1)*log(2)
 */
  return
    x * (3.529304993f + x * (-2.461222105f + x * (1.130626167f +
      x * (-0.288739945f + x * 3.110401639e-2f))))
    + (addcst + 0.6931471805f*exp);
}

/* Absolute error bounded by 5e-9 for normalized inputs
   Returns a finite number for +inf input
   Returns -inf for nan and <= 0 inputs.
   Continuous error.

   Vectorizable only with AVX512dq extensions because of the
   int64 left shift. On GCC, use option -mavx512dq. */
inline double logapprox_d(double val) {
  union { double f; int64_t i; } valu;
  double exp, addcst, x;
  valu.f = val;
  exp = valu.i >> 52;
  /* -711.57470503155781879 = -1023 * log(2) + constant term of polynomial bellow. */
  addcst = val > 0 ? -711.57470503155781879 : -(double)INFINITY;
  valu.i = (valu.i & 0xFFFFFFFFFFFFF) | 0x3FF0000000000000;
  x = valu.f;


  /* Generated in Sollya using:
    > f = remez(log(x)-(x-1)*log(2),
            [|1,(x-1)*(x-2), (x-1)*(x-2)*x, (x-1)*(x-2)*x*x,
              (x-1)*(x-2)*x*x*x|], [1,2], 1, 1e-8);
    > plot(f+(x-1)*log(2)-log(x), [1,2]);
    > f+(x-1)*log(2)
 */
  return x * (6.3599808385417079442 + x * (-8.9283961901995856725 +
          x * (9.6834305233551590662 + x * (-7.5448441898968122224 +
           x * (4.1526922366643242018 + x * (-1.5768742873140365359 +
            x * (0.3933681393710522449 + x * (-5.8063046913513248806e-2 +
             x * 3.8452996567427771004e-3))))))))
    + (addcst + 0.693147180559945309*exp);
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 5e-5
   Continuous error */
inline float cosapprox(float val) {
  float val2 = val*val;
  /* Generated in Sollya using:
     > f = remez(cos(x)-1, [|x*x, x*x*x*x, x*x*x*x*x*x, x*x*x*x*x*x*x*x|],
                            [0.000001, pi], 1, 1e-8);
     > plot(f-cos(x)+1, [0, pi]);
     > f+1
  */
  return
    1.f + val2 * (-0.4998515820f + val2 * (4.1518035216e-2f + val2 *
      (-1.3422947025e-3f + val2 * 1.8929864824e-5f)));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 2e-10
   Continuous error */
inline double cosapprox_d(double val) {
  double val2 = val*val;
  /* Generated in Sollya using:
     > f = remez(cos(x)-1, [|x*x, x*x*x*x, x*x*x*x*x*x, x*x*x*x*x*x*x*x,
                             x*x*x*x*x*x*x*x*x*x, x*x*x*x*x*x*x*x*x*x*x*x,
                             x*x*x*x*x*x*x*x*x*x*x*x*x*x|], [1e-10, pi], 1, 1e-10);
     > plot(f-cos(x)+1, [0, pi]);
     > f+1
 */
  return 1. + val2 * (-0.4999999989854595 + val2 * (4.1666664031946494e-2 +
               val2 * (-1.3888865732154298e-3 + val2 * (2.4800621164404954e-5 +
                val2 * (-2.7535709220609759e-7 + val2 * (2.0609775337948301e-9 +
                 val2 * (-9.7393350917248365e-12)))))));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 6e-6
   Continuous error */
inline float sinapprox(float val) {
  float val2 = val*val;
  return
    val * (0.9999793767f + val2 * (-0.1666243672f + val2 *
          (8.3089787513e-3f + val2 * (-1.9264918228e-4f + val2 *
           2.1478401777e-6f))));
}

/* Correct only in [-pi, pi]
   Absolute error bounded by 2e-9
   Continuous error */
inline double sinapprox_d(double val) {
  double val2 = val*val;
  return
    val * (0.9999999945159759653 + val2 * (-0.1666666458182987439 +
     val2 * (8.3333103922589284663e-3 + val2 * (-1.9840155355055654144e-4 +
      val2 * (2.7529454331962521774e-6 + val2 * (-2.4676970823046321831e-8 +
       val2 * 1.3451481340051383601e-10))))));
}


#ifdef __cplusplus
}
#endif

#endif
