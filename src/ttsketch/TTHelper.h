#ifndef __PLUMED_ttsketch_TTHelper_h
#define __PLUMED_ttsketch_TTHelper_h
#include "BasisFunc.h"
#include "tools/Matrix.h"
#ifdef __PLUMED_HAS_LIBITENSOR
#ifdef __PLUMED_HAS_LIBHDF5
#include <itensor/all.h>
#include <utility>

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
// Number of TTs stored in the file, i.e. the largest k such that "tt_1".."tt_k" all exist
// (0 if the file does not exist).
unsigned ttCount(const std::string& filename);
// Append the accumulated tensor moments of the NONINTRUSIVE variant (B and the cross-gram
// matrices A) under "bprev_{count-1}" and "aprev_{count-1}_{k}", k = 1..A.size().
void ttMomentsWrite(const std::string& filename, const itensor::MPS& B, const std::vector<itensor::ITensor>& A, unsigned count);
// Read what ttMomentsWrite() wrote; returns false if the moments are not in the file.
bool ttMomentsRead(const std::string& filename, unsigned count, unsigned nA, itensor::MPS& B, std::vector<itensor::ITensor>& A);
// Evaluate the TT bias at point `elements` by contracting each core with its
// 1D basis vector. If conv=true, basis functions are evaluated in convolution mode.
double ttEval(const itensor::MPS& tt, const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv);
// Compute the covariance matrix, mean vector, and partition function of the TT distribution.
// Normalizes tt by Z = int tt(x) dx, then returns:
//   sigma(k,l) = E[x_k * x_l] - E[x_k]*E[x_l]  (covariance matrix)
//   mu[k]      = E[x_k]                           (marginal means)
//   Z          = int tt(x) dx                     (partition function)
// All expectations use int0/int1/int2 evaluated analytically via BasisFunc.
std::tuple<Matrix<double>, std::vector<double>, double> covMat(const itensor::MPS& tt, const std::vector<BasisFunc>& basis);
// Fill `grid` (bins x bins) with the 2D marginal density of `tt` for dimensions
// pos1 and pos2, obtained by integrating out all other dimensions analytically.
void marginal2d(const itensor::MPS& tt, const std::vector<BasisFunc>& basis, int pos1, int pos2, Matrix<double>& grid, bool conv);
// Plain-array copy of a TT, used to evaluate the bias and its gradient at every step.
// Core i is stored as core_[i][(a*n_i + j)*r_i + b], with a/b the left/right bond
// indices (r_0 = r_d = 1) and j the basis index. Building it costs O(sum r n r) ITensor
// element accesses, so it is rebuilt only when the TT changes.
class DenseTT {
  std::vector<int> rank_;                 // bond dimensions r_0..r_d
  std::vector<int> nbasis_;               // site dimensions n_1..n_d
  std::vector<std::vector<double>> core_; // the cores
public:
  DenseTT() = default;
  explicit DenseTT(const itensor::MPS& tt);
  bool empty() const {
    return core_.empty();
  }
  unsigned size() const {
    return core_.size();
  }
  int nbasis(unsigned i) const {
    return nbasis_[i];
  }
  // bond dimension r_i, i = 0..d (r_0 = r_d = 1)
  int rank(unsigned i) const {
    return rank_[i];
  }
  // core i (0-based), element (a,j,b) at (a*nbasis(i) + j)*rank(i+1) + b
  const std::vector<double>& core(unsigned i) const {
    return core_[i];
  }
  // Value of the TT at `elements`; if grad is not null it is resized to d and filled
  // with the gradient. Same result as ttEvalAndGrad(), at O(d n r^2) cost.
  double eval(const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv, std::vector<double>* grad) const;
};
std::pair<double, std::vector<double>>
                                    ttEvalAndGrad(const itensor::MPS& tt,
                                        const std::vector<BasisFunc>& basis,
                                        const std::vector<double>& elements,
                                        bool conv);
}
}
#endif // __PLUMED_HAS_LIBHDF5
#endif // __PLUMED_HAS_LIBITENSOR
#endif
