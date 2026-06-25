/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Nils E. Strand, Siyao Yang, Yuehaw Khoo, Aaron R. Dinner

   See the COPYRIGHT file distributed with this software for license details.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "TTHelper.h"
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
  auto tt = h5_read<MPS>(f, "tt_" + std::to_string(count - 1));
  return tt;
}

void ttSumWrite(const std::string& filename, const MPS& tt, unsigned count) {
  auto f = h5_open(filename, 'a');
  h5_write(f, "vb_" + std::to_string(count - 1), tt);
}

MPS ttSumRead(const std::string& filename, unsigned count) {
  auto f = h5_open(filename, 'r');
  auto tt = h5_read<MPS>(f, "vb_" + std::to_string(count - 1));
  return tt;
}

// Precomputed basis evaluations, reusable across ttEval and ttGrad
struct BasisCache {
  std::vector<ITensor> phi;   // regular basis vectors
  std::vector<ITensor> dphi;  // derivative basis vectors (empty if not needed)
};

BasisCache buildBasisCache(const MPS& tt,
                           const std::vector<BasisFunc>& basis,
                           const std::vector<double>& elements,
                           bool conv,
                           bool with_grad = false) {
  int d = length(tt);
  auto s = siteInds(tt);
  BasisCache cache;
  cache.phi.resize(d);
  if (with_grad) {
    cache.dphi.resize(d);
  }

  // This loop is embarrassingly parallel and cheap; OMP overhead likely not worth it
  // unless d is large (d > ~32). Profile before adding #pragma omp parallel for here.
  for (int i = 0; i < d; ++i) {
    cache.phi[i] = ITensor(s(i + 1));
    if (with_grad) {
      cache.dphi[i] = ITensor(s(i + 1));
    }
    for (int j = 1; j <= dim(s(i + 1)); ++j) {
      cache.phi[i].set(s(i + 1) = j, basis[i](elements[i], j, conv));
      if (with_grad) {
        cache.dphi[i].set(s(i + 1) = j, basis[i].grad(elements[i], j, conv));
      }
    }
  }
  return cache;
}

// Evaluate TT using a prebuilt BasisCache (avoids recomputing phi when called
// immediately before ttGrad).
double ttEvalCached(const MPS& tt, const BasisCache& cache) {
  int d = length(tt);
  auto result = tt(1) * cache.phi[0];
  for (int i = 2; i <= d; ++i) {
    result *= tt(i) * cache.phi[i - 1];
  }
  return elt(result);
}

// Original ttEval kept for standalone use (no redundant work when called alone)
double ttEval(const MPS& tt,
              const std::vector<BasisFunc>& basis,
              const std::vector<double>& elements,
              bool conv) {
  auto cache = buildBasisCache(tt, basis, elements, conv, /*with_grad=*/false);
  return ttEvalCached(tt, cache);
}

// Gradient using prefix-suffix sweep: O(d) contractions instead of O(d^2).
//
// For a TT of length d, define:
//   L[k] = (G_1 * phi_1) * ... * (G_k * phi_k)   [left partial product up to k]
//   R[k] = (G_k * phi_k) * ... * (G_d * phi_d)   [right partial product from k]
//
// Then grad[k] = L[k-1] * (G_k * dphi_k) * R[k+1]
//
// We precompute all L[] (left-to-right pass) and all R[] (right-to-left pass),
// then combine them with a single derivative insertion per dimension.
std::vector<double> ttGradCached(const MPS& tt, const BasisCache& cache) {
  int d = length(tt);
  std::vector<double> grad(d, 0.0);

  // --- Left prefix products: L[i] = contraction of sites 1..i with phi ---
  // L[0] is a scalar 1 (identity for the left boundary)
  // L[i] has the bond index between site i and i+1 dangling (or is scalar for i==d)
  std::vector<ITensor> L(d + 1), R(d + 2);

  L[0] = ITensor(1.0);  // scalar identity
  for (int i = 1; i <= d; ++i) {
    L[i] = L[i - 1] * (tt(i) * cache.phi[i - 1]);
  }

  // --- Right suffix products: R[i] = contraction of sites i..d with phi ---
  R[d + 1] = ITensor(1.0);
  for (int i = d; i >= 1; --i) {
    R[i] = (tt(i) * cache.phi[i - 1]) * R[i + 1];
  }

  // --- Combine: grad[k] = L[k-1] * (G_k * dphi_k) * R[k+1] ---
  // Each combination is O(r^2 * n) instead of the full O(d * r^2 * n)
  for (int k = 1; k <= d; ++k) {
    grad[k - 1] = elt(L[k - 1] * (tt(k) * cache.dphi[k - 1]) * R[k + 1]);
  }

  return grad;
}

std::pair<double, std::vector<double>>
                                    ttEvalAndGrad(const itensor::MPS& tt,
                                        const std::vector<BasisFunc>& basis,
                                        const std::vector<double>& elements,
bool conv) {
  int d = length(tt);
  auto s = siteInds(tt);

  // Cache cores once — tt(i) does index lookup every call
  std::vector<ITensor> cores(d);
  for (int i = 0; i < d; ++i) {
    cores[i] = tt(i + 1);
  }

  // Build phi and dphi
  std::vector<ITensor> phi(d), dphi(d);
  for (int i = 0; i < d; ++i) {
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

  std::vector<double> grad(d);
  for (int k = 0; k < d; ++k) {
    grad[k] = elt(L[k] * (cores[k] * dphi[k]) * R[k + 1]);
  }

  return {val, grad};
}

// Standalone gradient (for places that don't need the value)
std::vector<double> ttGrad(const MPS& tt,
                           const std::vector<BasisFunc>& basis,
                           const std::vector<double>& elements,
                           bool conv) {
  auto cache = buildBasisCache(tt, basis, elements, conv, /*with_grad=*/true);
  return ttGradCached(tt, cache);
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
