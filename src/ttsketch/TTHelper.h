#ifndef __PLUMED_ttsketch_TTHelper_h
#define __PLUMED_ttsketch_TTHelper_h
#include "BasisFunc.h"
#include "tools/Matrix.h"
#include "itensor/all.h"

namespace PLMD {
namespace ttsketch {

// Write TT to HDF5 file under dataset "tt_{count-1}".
// count=2 opens the file for writing (creating/truncating); subsequent calls append.
void ttWrite(const std::string& filename, const itensor::MPS& tt, unsigned count);
// Read TT from HDF5 file at dataset "tt_{count-1}".
itensor::MPS ttRead(const std::string& filename, unsigned count);
// Append the cumulative bias TT (vb) to HDF5 file under dataset "vb_{count-1}".
void ttSumWrite(const std::string& filename, const itensor::MPS& tt, unsigned count);
// Read cumulative bias TT from HDF5 file at dataset "vb_{count-1}".
itensor::MPS ttSumRead(const std::string& filename, unsigned count);
// Evaluate the TT bias at point `elements` by contracting each core with its
// 1D basis vector. If conv=true, basis functions are evaluated in convolution mode.
double ttEval(const itensor::MPS& tt, const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv);
// Gradient of ttEval w.r.t. `elements`. Computed by replacing the k-th dimension's
// basis vector with its derivative while keeping all others unchanged (chain rule).
std::vector<double> ttGrad(const itensor::MPS& tt, const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv);
// Compute the covariance matrix, mean vector, and partition function of the TT distribution.
// Normalizes tt by Z = int tt(x) dx, then returns:
//   sigma(k,l) = E[x_k * x_l] - E[x_k]*E[x_l]  (covariance matrix)
//   mu[k]      = E[x_k]                           (marginal means)
//   Z          = int tt(x) dx                     (partition function)
// All expectations use int0/int1/int2 evaluated analytically via BasisFunc.
std::tuple<Matrix<double>, std::vector<double>, double> covMat(const itensor::MPS& tt, const std::vector<BasisFunc>& basis);
// Fill `grid` (bins x bins) with the 2D marginal density of `tt` for dimensions
// pos1 and pos2, obtained by integrating out all other dimensions analytically.
void marginal2d(const itensor::MPS& tt, const std::vector<BasisFunc>& basis, int pos1, int pos2, std::vector<std::vector<double>>& grid, bool conv);

}
}

#endif
