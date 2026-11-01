/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Nils E. Strand, Siyao Yang, Yuehaw Khoo, Aaron R. Dinner

   See the COPYRIGHT file distributed with this software for license details.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "TTHelper.h"
#include "tools/OpenMP.h"
#include <cstdio>
#include <stdexcept>
#ifdef __PLUMED_HAS_LIBITENSOR
#ifdef __PLUMED_HAS_LIBHDF5

using namespace itensor;

namespace PLMD {
namespace ttsketch {

void ttWrite(const std::string& filename, const MPS& tt, unsigned count) {
  // count starts at 2 on the first sketch call; open for writing on first call, append thereafter
  auto f = count == 2 ? h5_open(filename, 'w') : h5_open(filename, 'a');
  h5_write(f, "tt_" + std::to_string(count - 1), tt);
}

MPS ttRead(const std::string& filename, unsigned count) {
  auto f = h5_open(filename, 'r');
  return h5_read<MPS>(f, "tt_" + std::to_string(count - 1));
}

void ttSumWrite(const std::string& filename, const MPS& tt, unsigned count) {
  auto f = h5_open(filename, 'a');
  h5_write(f, "vb_" + std::to_string(count - 1), tt);
}

MPS ttSumRead(const std::string& filename, unsigned count) {
  auto f = h5_open(filename, 'r');
  return h5_read<MPS>(f, "vb_" + std::to_string(count - 1));
}

unsigned ttCount(const std::string& filename) {
  {
    // h5::file throws (and HDF5 prints an error stack) on a missing file
    FILE* fp = std::fopen(filename.c_str(), "r");
    if(!fp) {
      return 0;
    }
    std::fclose(fp);
  }
  auto f = h5_open(filename, 'r');
  h5::group g(f);
  unsigned n = 0;
  while(g.has_subgroup("tt_" + std::to_string(n + 1))) {
    ++n;
  }
  return n;
}

void ttMomentsWrite(const std::string& filename, const MPS& B, const std::vector<ITensor>& A, unsigned count) {
  auto f = h5_open(filename, 'a');
  const std::string tag = std::to_string(count - 1);
  h5_write(f, "bprev_" + tag, B);
  for(unsigned k = 0; k < A.size(); ++k) {
    h5_write(f, "aprev_" + tag + "_" + std::to_string(k + 1), A[k]);
  }
}

bool ttMomentsRead(const std::string& filename, unsigned count, unsigned nA, MPS& B, std::vector<ITensor>& A) {
  auto f = h5_open(filename, 'r');
  h5::group g(f);
  const std::string tag = std::to_string(count - 1);
  if(!g.has_subgroup("bprev_" + tag)) {
    return false;
  }
  for(unsigned k = 0; k < nA; ++k) {
    if(!g.has_subgroup("aprev_" + tag + "_" + std::to_string(k + 1))) {
      return false;
    }
  }
  B = h5_read<MPS>(f, "bprev_" + tag);
  A.resize(nA);
  for(unsigned k = 0; k < nA; ++k) {
    A[k] = h5_read<ITensor>(f, "aprev_" + tag + "_" + std::to_string(k + 1));
  }
  return true;
}

// Evaluate TT and gradient together using a prefix-suffix sweep: O(d) contractions.
//
// For a TT of length d:
//   L[i] = (G_1 * phi_1) * ... * (G_i * phi_i)   [left partial product up to i]
//   R[i] = (G_i * phi_i) * ... * (G_d * phi_d)   [right partial product from i]
//
// val       = elt(L[d])
// grad[k]   = L[k-1] * (G_k * dphi_k) * R[k+1]
std::pair<double, std::vector<double>>
                                    ttEvalAndGrad(const itensor::MPS& tt,
                                        const std::vector<BasisFunc>& basis,
                                        const std::vector<double>& elements,
bool conv) {
  int d = length(tt);
  auto s = siteInds(tt);
  const unsigned nt = OpenMP::getNumThreads();

  // Cache cores and build phi/dphi — all independent across dimensions
  std::vector<ITensor> cores(d), phi(d), dphi(d);
  #pragma omp parallel for num_threads(nt) schedule(static)
  for (int i = 0; i < d; ++i) {
    cores[i] = tt(i + 1);
    phi[i] = dphi[i] = ITensor(s(i + 1));
    for (int j = 1; j <= dim(s(i + 1)); ++j) {
      phi[i].set(s(i + 1) = j,  basis[i](elements[i], j, conv));
      dphi[i].set(s(i + 1) = j, basis[i].grad(elements[i], j, conv));
    }
  }

  // Prefix-suffix sweep
  std::vector<ITensor> L(d + 1), R(d + 1);
  L[0] = ITensor(1.0);
  for (int i = 0; i < d; ++i) {
    L[i + 1] = L[i] * (cores[i] * phi[i]);
  }
  R[d] = ITensor(1.0);
  for (int i = d - 1; i >= 0; --i) {
    R[i] = (cores[i] * phi[i]) * R[i + 1];
  }

  double val = elt(L[d]);

  // Gradient combination — independent across k
  std::vector<double> grad(d);
  #pragma omp parallel for num_threads(nt) schedule(static)
  for (int k = 0; k < d; ++k) {
    grad[k] = elt(L[k] * (cores[k] * dphi[k]) * R[k + 1]);
  }

  return {val, grad};
}

DenseTT::DenseTT(const MPS& tt) {
  const int d = length(tt);
  rank_.assign(d + 1, 1);
  nbasis_.resize(d);
  core_.resize(d);
  for(int i = 1; i <= d; ++i) {
    const Index s = siteIndex(tt, i);
    const Index l = i > 1 ? linkIndex(tt, i - 1) : Index();
    const Index r = i < d ? linkIndex(tt, i) : Index();
    const int rl = i > 1 ? dim(l) : 1;
    const int rr = i < d ? dim(r) : 1;
    const int n = dim(s);
    rank_[i] = rr;
    nbasis_[i - 1] = n;
    // the core has exactly the indices (l,s,r) without the missing boundary bonds
    if(order(tt(i)) != 1 + (i > 1) + (i < d)) {
      throw std::runtime_error("unexpected index structure in TT core " + std::to_string(i));
    }
    auto& c = core_[i - 1];
    c.resize(static_cast<size_t>(rl) * n * rr);
    for(int a = 1; a <= rl; ++a) {
      for(int j = 1; j <= n; ++j) {
        for(int b = 1; b <= rr; ++b) {
          double v;
          if(i > 1 && i < d) {
            v = tt(i).elt(l = a, s = j, r = b);
          } else if(i > 1) {
            v = tt(i).elt(l = a, s = j);
          } else if(i < d) {
            v = tt(i).elt(s = j, r = b);
          } else {
            v = tt(i).elt(s = j);
          }
          c[(static_cast<size_t>(a - 1) * n + (j - 1)) * rr + (b - 1)] = v;
        }
      }
    }
  }
}

double DenseTT::eval(const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv, std::vector<double>* grad) const {
  const unsigned d = core_.size();
  // contracted cores: M[i] = sum_j core_i[:,j,:] phi_j(x_i), D[i] the same with dphi
  thread_local std::vector<std::vector<double>> M, D;
  thread_local std::vector<double> phi, dphi;
  M.resize(d);
  D.resize(d);
  for(unsigned i = 0; i < d; ++i) {
    const int rl = rank_[i], rr = rank_[i + 1], n = nbasis_[i];
    phi.resize(n);
    dphi.resize(n);
    basis[i].evalAll(elements[i], conv, phi.data(), dphi.data());
    M[i].assign(static_cast<size_t>(rl) * rr, 0.0);
    if(grad) {
      D[i].assign(static_cast<size_t>(rl) * rr, 0.0);
    }
    const double* c = core_[i].data();
    for(int a = 0; a < rl; ++a) {
      double* m = M[i].data() + static_cast<size_t>(a) * rr;
      double* g = grad ? D[i].data() + static_cast<size_t>(a) * rr : nullptr;
      for(int j = 0; j < n; ++j) {
        const double* cj = c + (static_cast<size_t>(a) * n + j) * rr;
        const double p = phi[j];
        for(int b = 0; b < rr; ++b) {
          m[b] += cj[b] * p;
        }
        if(g) {
          const double dp = dphi[j];
          for(int b = 0; b < rr; ++b) {
            g[b] += cj[b] * dp;
          }
        }
      }
    }
  }
  // left sweep: left[i] = M[0]...M[i-1] (row vector of size r_i)
  thread_local std::vector<std::vector<double>> left, right;
  left.resize(d + 1);
  left[0].assign(1, 1.0);
  for(unsigned i = 0; i < d; ++i) {
    const int rl = rank_[i], rr = rank_[i + 1];
    left[i + 1].assign(rr, 0.0);
    for(int a = 0; a < rl; ++a) {
      const double la = left[i][a];
      const double* m = M[i].data() + static_cast<size_t>(a) * rr;
      for(int b = 0; b < rr; ++b) {
        left[i + 1][b] += la * m[b];
      }
    }
  }
  const double val = left[d][0];
  if(grad) {
    // right sweep: right[i] = M[i]...M[d-1] (column vector of size r_i)
    right.resize(d + 1);
    right[d].assign(1, 1.0);
    for(int i = d - 1; i >= 0; --i) {
      const int rl = rank_[i], rr = rank_[i + 1];
      right[i].assign(rl, 0.0);
      for(int a = 0; a < rl; ++a) {
        const double* m = M[i].data() + static_cast<size_t>(a) * rr;
        double sum = 0.0;
        for(int b = 0; b < rr; ++b) {
          sum += m[b] * right[i + 1][b];
        }
        right[i][a] = sum;
      }
    }
    // grad_i = left[i] * D[i] * right[i+1]
    grad->resize(d);
    for(unsigned i = 0; i < d; ++i) {
      const int rl = rank_[i], rr = rank_[i + 1];
      double sum = 0.0;
      for(int a = 0; a < rl; ++a) {
        const double* g = D[i].data() + static_cast<size_t>(a) * rr;
        double t = 0.0;
        for(int b = 0; b < rr; ++b) {
          t += g[b] * right[i + 1][b];
        }
        sum += left[i][a] * t;
      }
      (*grad)[i] = sum;
    }
  }
  return val;
}

// ttEval delegates to ttEvalAndGrad and discards the gradient.
// Called infrequently (error measurement, well-tempered height) so the
// cost of computing dphi is negligible relative to avoiding code duplication.
double ttEval(const MPS& tt,
              const std::vector<BasisFunc>& basis,
              const std::vector<double>& elements,
              bool conv) {
  return ttEvalAndGrad(tt, basis, elements, conv).first;
}

// Compute covariance matrix, marginal means, and partition function of the TT distribution.
// Precomputes ITensors for the three moment integrals int0/int1/int2 per dimension,
// then evaluates expectations as TT contractions with these integral vectors.
std::tuple<Matrix<double>, std::vector<double>, double> covMat(const MPS& tt, const std::vector<BasisFunc>& basis) {
  int d = length(tt);
  auto s = siteInds(tt);
  // integral vectors: int0[i][j] = int phi_j(x_i) dx_i, etc.
  std::vector<ITensor> basis_int0(d), basis_int1(d), basis_int2(d);
  for(int i = 1; i <= d; ++i) {
    basis_int0[i - 1] = basis_int1[i - 1] = basis_int2[i - 1] = ITensor(s(i));
    for(int j = 1; j <= dim(s(i)); ++j) {
      basis_int0[i - 1].set(s(i) = j, basis[i - 1].int0(j));
      basis_int1[i - 1].set(s(i) = j, basis[i - 1].int1(j));
      basis_int2[i - 1].set(s(i) = j, basis[i - 1].int2(j));
    }
  }
  // Z = int tt(x) dx (partition function); normalize to get probability density rho
  auto Z = tt(1) * basis_int0[0];
  for(int i = 2; i <= d; ++i) {
    Z *= tt(i) * basis_int0[i - 1];
  }
  auto rho = tt;
  rho /= elt(Z);
  // ei[k]  = E[x_k],  eii[k] = E[x_k^2],  eij[k][l] = E[x_k * x_l] for k < l
  std::vector<double> ei(d), eii(d);
  Matrix<double> eij(d, d);
  for(int k = 1; k <= d; ++k) {
    // replace dimension k's integral vector with int1 (resp. int2) to get E[x_k] (resp. E[x_k^2])
    auto eival = rho(1) * (k == 1 ? basis_int1[0] : basis_int0[0]);
    auto eiival = rho(1) * (k == 1 ? basis_int2[0] : basis_int0[0]);
    for(int i = 2; i <= d; ++i) {
      eival *= rho(i) * (k == i ? basis_int1[i - 1] : basis_int0[i - 1]);
      eiival *= rho(i) * (k == i ? basis_int2[i - 1] : basis_int0[i - 1]);
    }
    ei[k - 1] = elt(eival);
    eii[k - 1] = elt(eiival);
    for(int l = k + 1; l <= d; ++l) {
      // replace both dimensions k and l with int1 to get E[x_k * x_l]
      auto eijval = rho(1) * (k == 1 ? basis_int1[0] : basis_int0[0]);
      for(int i = 2; i <= d; ++i) {
        eijval *= rho(i) * (k == i || l == i ? basis_int1[i - 1] : basis_int0[i - 1]);
      }
      eij(k - 1, l - 1) = elt(eijval);
    }
  }
  // sigma(k,l) = Cov(x_k, x_l) = E[x_k*x_l] - E[x_k]*E[x_l]
  Matrix<double> sigma(d, d);
  for(int k = 1; k <= d; ++k) {
    for(int l = k; l <= d; ++l) {
      sigma(k - 1, l - 1) = sigma(l - 1, k - 1) = k == l ? eii[k - 1] - pow(ei[k - 1], 2) : eij(k - 1, l - 1) - ei[k - 1] * ei[l - 1];
    }
  }
  return std::make_tuple(sigma, ei, elt(Z));
}

// Compute the 2D marginal density of the normalized TT distribution for dimensions
// pos1 and pos2 on a (bins x bins) grid. All other dimensions are integrated out
// analytically using int0, which collapses those TT cores to scalar factors.
void marginal2d(const MPS& tt, const std::vector<BasisFunc>& basis, int pos1, int pos2, Matrix<double>& grid, bool conv) {
  int bins = grid.nrows();
  int d = length(tt);
  auto s = siteInds(tt);
  std::vector<ITensor> basis_int0(d);
  for(int i = 1; i <= d; ++i) {
    basis_int0[i - 1] = ITensor(s(i));
    for(int j = 1; j <= dim(s(i)); ++j) {
      basis_int0[i - 1].set(s(i) = j, basis[i - 1].int0(j));
    }
  }
  // normalize to probability density
  auto Z = tt(1) * basis_int0[0];
  for(int i = 2; i <= d; ++i) {
    Z *= tt(i) * basis_int0[i - 1];
  }
  auto rho = tt;
  rho /= elt(Z);
  for(int i = 0; i < bins; ++i) {
    for(int j = 0; j < bins; ++j) {
      double x = basis[pos1 - 1].dom().first + i * (basis[pos1 - 1].dom().second - basis[pos1 - 1].dom().first) / bins;
      double y = basis[pos2 - 1].dom().first + j * (basis[pos2 - 1].dom().second - basis[pos2 - 1].dom().first) / bins;
      ITensor xevals(s(pos1)), yevals(s(pos2));
      for(int k = 1; k <= dim(s(pos1)); ++k) {
        xevals.set(s(pos1) = k, basis[pos1 - 1](x, k, conv));
        yevals.set(s(pos2) = k, basis[pos2 - 1](y, k, conv));
      }
      auto val = rho(1) * (pos1 == 1 ? xevals : basis_int0[0]);
      for(int k = 2; k <= d; ++k) {
        if(pos1 == k) {
          val *= rho(k) * xevals;
        } else if(pos2 == k) {
          val *= rho(k) * yevals;
        } else {
          val *= rho(k) * basis_int0[k - 1];
        }
      }
      grid(i, j) = elt(val);
    }
  }
}

}
}
#endif // __PLUMED_HAS_LIBHDF5
#endif // __PLUMED_HAS_LIBITENSOR
