#include "BasisFunc.h"

using namespace std;

namespace PLMD {
namespace ttsketch {

BasisFunc::BasisFunc()
  : dom_(make_pair(0.0, 0.0)), nbasis_(0), L_(0.0), shift_(0.0), w_(0.0), kernel_(false), dx_(0.0) {}

BasisFunc::BasisFunc(pair<double, double> dom, int nbasis, double w, bool kernel, double dx)
  : dom_(dom), nbasis_(nbasis), L_((dom.second - dom.first) / 2), shift_((dom.second + dom.first) / 2), w_(w), kernel_(kernel), dx_(dx)
{
  if(kernel) {
    double spacing = (dom.second - dom.first) / (nbasis - 1);
    if(dx == 0) {
      // default kernel width equals the grid spacing
      this->dx_ = spacing;
    }
    // nbasis-1 kernel centers, uniformly spaced starting at dom.first
    this->centers_ = vector<double>(nbasis - 1);
    for(int i = 0; i < nbasis - 1; ++i) {
      this->centers_[i] = dom.first + i * spacing;
    }
    // Compute the Gram matrix G where G(i,j) = int phi_i(x)*phi_j(x) dx,
    // integrating the periodized Gaussian kernels analytically over the domain.
    // phi_1 = 1 (constant), phi_{i+1} = sum_{k=-1,0,1} exp(-(x-c_i+2kL)^2/(2*dx^2)).
    // The pseudo-inverse of G is stored in ginv_ and applied after sketching to
    // recover dual basis coefficients from the Gram representation.
    this->gram_ = Matrix<double>(nbasis, nbasis);
    // G(0,0) = int 1*1 dx = domain length
    this->gram_(0, 0) = this->dom_.second - this->dom_.first;
    for(int i = 1; i < nbasis; ++i) {
      // G(0,i) = G(i,0) = int phi_{i+1}(x) dx (integral of a periodized Gaussian)
      this->gram_(i, 0) = this->gram_(0, i) = this->dx_ * sqrt(M_PI / 2) *
                                (erf((this->dom_.second -
                                2 * this->dom_.first + this->centers_[i - 1]) /
                                (sqrt(2) * this->dx_)) -
                                erf((this->dom_.first - 2 * this->dom_.second +
                                this->centers_[i - 1]) /
                                (sqrt(2) * this->dx_)));
      // G(i,j) = int phi_{i+1}(x)*phi_{j+1}(x) dx; the double sum over k,l
      // accounts for cross-terms between periodic images of each kernel.
      for(int j = i; j < nbasis; ++j) {
        double result = 0.0;
        for(int k = -1; k <= 1; ++k) {
          for(int l = -1; l <= 1; ++l) {
            result += this->dx_ / 2 * exp(-pow((this->dom_.first -
                      this->dom_.second) * (k - l) + this->centers_[i - 1] -
                      this->centers_[j - 1], 2) / (4 * pow(this->dx_, 2))) *
                      sqrt(M_PI) * (erf((this->dom_.first * (k + l - 2) -
                      this->dom_.second * (k + l) + this->centers_[i - 1] +
                      this->centers_[j - 1]) / (2 * this->dx_)) -
                      erf((this->dom_.first * (k + l) - this->dom_.second *
                      (k + l + 2) + this->centers_[i - 1] +
                      this->centers_[j - 1]) / (2 * this->dx_)));
          }
        }
        this->gram_(i, j) = this->gram_(j, i) = result;
      }
    }
    pseudoInvert(this->gram_, this->ginv_);
  }
}

// Orthonormal Fourier basis (not convolved):
//   pos=1: phi_1(x) = 1/sqrt(2L)
//   pos even: phi_pos(x) = 1/sqrt(L) * cos(pi*(x-shift)*(pos/2)/L)
//   pos odd:  phi_pos(x) = 1/sqrt(L) * sin(pi*(x-shift)*(pos/2)/L)
double BasisFunc::fourier(double x, int pos) const {
  if(pos == 1) {
    return 1 / sqrt(2 * this->L_);
  } else if(pos % 2 == 0) {
    return sqrt(1 / this->L_) * cos(M_PI * (x - this->shift_) * (pos / 2) / this->L_);
  } else {
    return sqrt(1 / this->L_) * sin(M_PI * (x - this->shift_) * (pos / 2) / this->L_);
  }
}

// Periodized Gaussian kernel basis:
//   pos=1: phi_1(x) = 1 (constant)
//   pos>1: phi_pos(x) = sum_{k=-1,0,1} exp(-(x - c_{pos-1} + 2kL)^2 / (2*dx^2))
// The k=-1,0,1 sum enforces periodicity across domain boundaries.
double BasisFunc::gaussian(double x, int pos) const {
  if(pos == 1) {
    return 1.0;
  } else {
    double result = 0.0;
    for(int k = -1; k <= 1; ++k) {
      result += exp(-pow(x - this->centers_[pos - 2] + 2 * k * this->L_, 2) / (2 * pow(this->dx_, 2)));
    }
    return result;
  }
}

double BasisFunc::fourierd(double x, int pos) const {
  if(pos == 1) {
    return 0.0;
  } else if(pos % 2 == 0) {
    return -pow(1 / this->L_, 1.5) * M_PI * (pos / 2) * sin(M_PI * (x - this->shift_) * (pos / 2) / this->L_);
  } else {
    return pow(1 / this->L_, 1.5) * M_PI * (pos / 2) * cos(M_PI * (x - this->shift_) * (pos / 2) / this->L_);
  }
}

double BasisFunc::gaussiand(double x, int pos) const {
  if(pos == 1) {
    return 0.0;
  } else {
    double result = 0.0;
    for(int k = -1; k <= 1; ++k) {
      result += (this->centers_[pos - 2] - x - 2 * k * this->L_) /
                pow(this->dx_, 2) * exp(-pow(x - this->centers_[pos - 2] +
                2 * k * this->L_, 2) / (2 * pow(this->dx_, 2)));
    }
    return result;
  }
}

// Evaluate basis function phi_pos(x), optionally in convolution (conv) mode.
// Conv mode for Fourier: multiplies each harmonic by exp(-(pi*w*k)^2/(2L^2)),
//   which is equivalent to convolving the basis function with a Gaussian of width w.
// Conv mode for kernel: widens each kernel analytically to std dev sqrt(dx^2+w^2),
//   equivalent to convolving the Gaussian kernel with a Gaussian of width w.
// In both cases the constant basis function (pos=1) is unchanged by convolution.
double BasisFunc::operator()(double x, int pos, bool conv) const {
  if(this->kernel_) {
    if(conv) {
      if(pos == 1) {
        return 1.0;
      } else {
        // Gaussian kernel convolved with Gaussian of width w_: std dev becomes sqrt(dx^2+w^2)
        double result = 0.0;
        for(int k = -1; k <= 1; ++k) {
          result += exp(-pow(x - this->centers_[pos - 2] +
                    2 * k * this->L_, 2) / (2 * (pow(this->dx_, 2) +
                    pow(this->w_, 2)))) / (sqrt(1 / pow(this->dx_, 2) + 1 /
                    pow(this->w_, 2)) * this->w_);
        }
        return result;
      }
    } else {
      return gaussian(x, pos);
    }
  } else {
    if(conv) {
      if(pos == 1) {
        return 1 / sqrt(2 * this->L_);
      } else {
        // Fourier convolution theorem: multiplying by exp(-(pi*w*k)^2/(2L^2))
        // is equivalent to convolving phi_pos with a Gaussian of width w_
        return exp(-pow(M_PI * this->w_ * (pos / 2), 2) / (2 * pow(this->L_, 2))) * fourier(x, pos);
      }
    } else {
      return fourier(x, pos);
    }
  }
}

double BasisFunc::grad(double x, int pos, bool conv) const {
  if(this->kernel_) {
    if(conv) {
      if(pos == 1) {
        return 0.0;
      } else {
        double result = 0.0;
        for(int k = -1; k <= 1; ++k) {
          result += pow(this->dx_, 2) * exp(-pow(x - this->centers_[pos - 2] +
                    2 * k * this->L_, 2) / (2 * (pow(this->dx_, 2) +
                    pow(this->w_, 2)))) * (this->centers_[pos - 2] -
                    2 * k * this->L_ - x) * sqrt(1 / pow(this->dx_, 2) +
                    1 / pow(this->w_, 2)) * this->w_ / pow(pow(this->dx_, 2) +
                    pow(this->w_, 2), 2);
        }
        return result;
      }
    } else {
      return gaussiand(x, pos);
    }
  } else {
    if(conv) {
      if(pos == 1) {
        return 0.0;
      } else {
        return exp(-pow(M_PI * this->w_ * (pos / 2), 2) / (2 * pow(this->L_, 2))) * fourierd(x, pos);
      }
    } else {
      return fourierd(x, pos);
    }
  }
}

// Compute int_a^b phi_pos(x) dx analytically.
// Used by covMat() to evaluate the partition function Z = <1> under the TT distribution.
double BasisFunc::int0(int pos) const {
  if(this->kernel_) {
    if(pos == 1) {
      return this->dom_.second - this->dom_.first;
    } else {
      return -this->dx_ * sqrt(M_PI / 2) * (erf((2 * this->dom_.first -
             this->dom_.second - this->centers_[pos - 2]) /
             (sqrt(2) * this->dx_)) + erf((this->dom_.first -
             2 * this->dom_.second + this->centers_[pos - 2]) /
             (sqrt(2) * this->dx_)));
    }
  } else {
    if(pos == 1) {
      return sqrt(2 * this->L_);
    } else if(pos % 2 == 0) {
      return 2 * sqrt(this->L_) / (M_PI * (pos / 2)) * sin(M_PI * (pos / 2));
    } else {
      return 0.0;
    }
  }
}

// Compute int_a^b x*phi_pos(x) dx analytically.
// Used by covMat() to evaluate E[x] = (1/Z) * <x> under the TT distribution.
double BasisFunc::int1(int pos) const {
  if(this->kernel_) {
    if(pos == 1) {
      return (pow(this->dom_.second, 2) - pow(this->dom_.first, 2)) / 2;
    } else {
      double a = this->dom_.first;
      double b = this->dom_.second;
      double c = this->centers_[pos - 2];
      double dx_sq = pow(this->dx_, 2);
      double sqrt2pi = sqrt(2 * M_PI);
      double term1 = exp(-pow(b - 2 * a + c, 2) / (2 * dx_sq)) -
                    exp(-pow(a - 2 * b + c, 2) / (2 * dx_sq));
      double term2 = (b - a + c) * sqrt2pi * erf((a - c) / (sqrt(2) * this->dx_));
      double term3 = (a - b - c) * sqrt2pi * erf((2 * a - b - c) / (sqrt(2) * this->dx_));
      double term4 = sqrt2pi * (c * erf((-a + c) / (sqrt(2) * this->dx_)) -
                      (a - b + c) * erf((a - 2 * b + c) / (sqrt(2) * this->dx_)) +
                      (a - b) * erf((-b + c) / (sqrt(2) * this->dx_)));
      return this->dx_ / 2 * (2 * this->dx_ * term1 + term2 + term3 + term4);
    }
  } else {
    if(pos == 1) {
      return this->shift_ * sqrt(2 * this->L_);
    } else if(pos % 2 == 0) {
      return 2 * this->shift_ * sqrt(this->L_) / (M_PI * (pos / 2)) * sin(M_PI * (pos / 2));
    } else {
      return 2 * pow(this->L_, 1.5) / pow(M_PI * (pos / 2), 2) * (sin(M_PI * (pos / 2)) - M_PI * (pos / 2) * cos(M_PI * (pos / 2)));
    }
  }
}

// Compute int_a^b x^2*phi_pos(x) dx analytically.
// Used by covMat() to evaluate E[x^2] for variance: Var[x] = E[x^2] - E[x]^2.
double BasisFunc::int2(int pos) const {
  if(this->kernel_) {
    if(pos == 1) {
      return (pow(this->dom_.second, 3) - pow(this->dom_.first, 3)) / 3;
    } else {
      double a = this->dom_.first;
      double b = this->dom_.second;
      double c = this->centers_[pos - 2];
      double dx_sq = pow(this->dx_, 2);
      double sqrt2pi = sqrt(2 * M_PI);
      double sqrt2_dx = sqrt(2) * this->dx_;
      double exp1 = exp(-pow(a - c, 2) / (2 * dx_sq));
      double exp2 = exp(-pow(b - c, 2) / (2 * dx_sq));
      double exp3 = exp(-pow(a - 2 * b + c, 2) / (2 * dx_sq));
      double exp4 = exp(-pow(b - 2 * a + c, 2) / (2 * dx_sq));
      double term1 = 2 * this->dx_ * (
          a * (2 * exp1 + 2 * exp2 - exp3) +
          b * (exp4 - 2 * exp1 - 2 * exp2) +
          c * (exp4 - exp3)
      );
      double diff1_sq = pow(b - a + c, 2) + dx_sq;
      double diff2_sq = pow(a - b + c, 2) + dx_sq;
      double c_sq = pow(c, 2) + dx_sq;
      double erf1 = erf((a - c) / sqrt2_dx);
      double erf2 = erf((2 * a - b - c) / sqrt2_dx);
      double erf3 = erf((c - a) / sqrt2_dx);
      double erf4 = erf((a - 2 * b + c) / sqrt2_dx);
      double erf5 = erf((c - b) / sqrt2_dx);
      double term2 = diff1_sq * sqrt2pi * erf1 - diff1_sq * sqrt2pi * erf2;
      double term3 = sqrt2pi * (
          c_sq * erf3 - diff2_sq * erf4 +
          (a - b) * (a - b + 2 * c) * erf5
      );
      return this->dx_ / 2 * (term1 + term2 + term3);
    }
  } else {
    if(pos == 1) {
      return sqrt(2 * this->L_) / 3 * (3 * pow(this->shift_, 2) + pow(this->L_, 2));
    } else if(pos % 2 == 0) {
      return 2 * sqrt(this->L_) / pow(M_PI * (pos / 2), 3) *
            (2 * pow(this->L_, 2) * M_PI * (pos / 2) * cos(M_PI * (pos / 2)) +
            (pow(M_PI * (pos / 2), 2) * (pow(this->shift_, 2) +
            pow(this->L_, 2)) - 2 * pow(this->L_, 2)) * sin(M_PI * (pos / 2)));
    } else {
      return 4 * this->shift_ * pow(this->L_, 1.5) / pow(M_PI * (pos / 2), 2) * (sin(M_PI * (pos / 2)) - M_PI * (pos / 2) * cos(M_PI * (pos / 2)));
    }
  }
}

}
}
