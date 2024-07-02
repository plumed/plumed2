#ifndef __PLUMED_ttsketch_BasisFunc_h
#define __PLUMED_ttsketch_BasisFunc_h
#include "tools/Matrix.h"
#include <utility>
#include <vector>

namespace PLMD {
namespace ttsketch {

// Represents a 1D basis expansion for one CV dimension, used to express the
// bias potential as a tensor train over basis function indices. Two types are
// supported:
//
//   Fourier (kernel_=false): orthonormal sinusoidal functions on [a, b].
//     phi_1(x) = 1/sqrt(2L)
//     phi_{2k}(x) = 1/sqrt(L) * cos(pi*(x-shift)*k/L)   (k = 1, 2, ...)
//     phi_{2k+1}(x) = 1/sqrt(L) * sin(pi*(x-shift)*k/L) (k = 1, 2, ...)
//
//   Kernel (kernel_=true): periodized Gaussian kernels on a uniform grid,
//     plus a constant function phi_1(x) = 1.
//
// "Conv" (convolution) mode pre-smooths basis functions by convolving with a
// Gaussian of width w_. For Fourier bases this multiplies each harmonic by
// exp(-(pi*w*k)^2 / (2L^2)); for kernel bases it widens each kernel
// analytically. This yields smooth bias forces at no extra runtime cost.
class BasisFunc {

private:
  std::pair<double, double> dom_; // CV domain [a, b]
  int nbasis_;                    // number of basis functions n (including constant phi_1)
  double L_;                      // domain half-length: L = (b-a)/2
  double shift_;                  // domain midpoint: shift = (a+b)/2
  double w_;                      // smoothing width for conv mode (Fourier: sigma in x-units; kernel: Gaussian hill width)
  bool kernel_;                   // true = kernel basis, false = Fourier basis
  double dx_;                     // kernel basis function width (std dev of each periodized Gaussian)
  std::vector<double> centers_;   // centers of kernel basis functions (size nbasis-1, uniformly spaced)
  Matrix<double> ginv_;           // pseudo-inverse of the Gram matrix (kernel basis only)
  Matrix<double> gram_;           // Gram matrix: gram_(i,j) = int phi_i(x)*phi_j(x) dx (kernel basis only)

public:
  BasisFunc();
  BasisFunc(std::pair<double, double> dom, int nbasis, double w, bool kernel, double dx);
  // Evaluate Fourier basis function phi_pos(x) (1-indexed, not convolved).
  double fourier(double x, int pos) const;
  // Evaluate kernel basis function phi_pos(x) (1-indexed, not convolved).
  // Uses periodic images (k=-1,0,1) so the function wraps smoothly across domain boundaries.
  double gaussian(double x, int pos) const;
  // Derivative of fourier(x, pos) w.r.t. x.
  double fourierd(double x, int pos) const;
  // Derivative of gaussian(x, pos) w.r.t. x.
  double gaussiand(double x, int pos) const;
  // Evaluate basis function phi_pos(x). If conv=true, applies Gaussian smoothing of width w_.
  double operator()(double x, int pos, bool conv) const;
  // Derivative of operator()(x, pos, conv) w.r.t. x.
  double grad(double x, int pos, bool conv) const;
  int nbasis() const { return this->nbasis_; }
  const std::pair<double, double>& dom() const { return this->dom_; }
  double w() const { return this->w_; }
  // Integral of phi_pos over the domain: int_a^b phi_pos(x) dx.
  // Used to compute the partition function Z in covMat().
  double int0(int pos) const;
  // Integral of x*phi_pos over the domain: int_a^b x*phi_pos(x) dx.
  // Used to compute the marginal mean E[x] in covMat().
  double int1(int pos) const;
  // Integral of x^2*phi_pos over the domain: int_a^b x^2*phi_pos(x) dx.
  // Used to compute E[x^2] for the marginal variance in covMat().
  double int2(int pos) const;
  bool kernel() const { return this->kernel_; }
  const Matrix<double>& ginv() const { return this->ginv_; }
  const Matrix<double>& gram() const { return this->gram_; }
  double center(int pos) const { return this->kernel_ ? this->centers_[pos - 1] : 0.0; }
  double dx() const { return this->dx_; }
};

}
}

#endif
