#include "TTHelper.h"

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

// Evaluate TT by contracting each core with its 1D basis vector phi(x_i):
//   result = sum_{i_1,...,i_d} G_1[i_1] * ... * G_d[i_d] * phi_1(x_1,i_1) * ... * phi_d(x_d,i_d)
// The contraction is performed left-to-right so intermediate results stay rank-1 scalars.
double ttEval(const MPS& tt, const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv) {
  int d = length(tt);
  auto s = siteInds(tt);
  std::vector<ITensor> basis_evals(d);
  for(int i = 1; i <= d; ++i) {
    basis_evals[i - 1] = ITensor(s(i));
    for(int j = 1; j <= dim(s(i)); ++j) {
      basis_evals[i - 1].set(s(i) = j, basis[i - 1](elements[i - 1], j, conv));
    }
  }
  auto result = tt(1) * basis_evals[0];
  for(int i = 2; i <= d; ++i) {
    result *= tt(i) * basis_evals[i - 1];
  }
  return elt(result);
}

// Gradient of ttEval w.r.t. elements. For each dimension k, one TT contraction is performed
// with the k-th basis vector replaced by its derivative d phi_k/dx_k (chain rule).
// All basis evaluations are precomputed to avoid redundant work across the d contractions.
std::vector<double> ttGrad(const MPS& tt, const std::vector<BasisFunc>& basis, const std::vector<double>& elements, bool conv) {
  int d = length(tt);
  auto s = siteInds(tt);
  std::vector<double> grad(d, 0.0);
  std::vector<ITensor> basis_evals(d), basisd_evals(d);
  for(int i = 1; i <= d; ++i) {
    basis_evals[i - 1] = basisd_evals[i - 1] = ITensor(s(i));
    for(int j = 1; j <= dim(s(i)); ++j) {
      basis_evals[i - 1].set(s(i) = j, basis[i - 1](elements[i - 1], j, conv));
      basisd_evals[i - 1].set(s(i) = j, basis[i - 1].grad(elements[i - 1], j, conv));
    }
  }
  for(int k = 1; k <= d; ++k) {
    // replace dimension k's basis vector with its derivative, keep all others unchanged
    auto result = tt(1) * (k == 1 ? basisd_evals[0] : basis_evals[0]);
    for(int i = 2; i <= d; ++i) {
      result *= tt(i) * (k == i ? basisd_evals[i - 1] : basis_evals[i - 1]);
    }
    grad[k - 1] = elt(result);
  }
  return grad;
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
  std::vector<std::vector<double>> eij(d, std::vector<double>(d));
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
      eij[k - 1][l - 1] = elt(eijval);
    }
  }
  // sigma(k,l) = Cov(x_k, x_l) = E[x_k*x_l] - E[x_k]*E[x_l]
  Matrix<double> sigma(d, d);
  for(int k = 1; k <= d; ++k) {
    for(int l = k; l <= d; ++l) {
      sigma(k - 1, l - 1) = sigma(l - 1, k - 1) = k == l ? eii[k - 1] - pow(ei[k - 1], 2) : eij[k - 1][l - 1] - ei[k - 1] * ei[l - 1];
    }
  }
  return std::make_tuple(sigma, ei, elt(Z));
}

// Compute the 2D marginal density of the normalized TT distribution for dimensions
// pos1 and pos2 on a (bins x bins) grid. All other dimensions are integrated out
// analytically using int0, which collapses those TT cores to scalar factors.
void marginal2d(const MPS& tt, const std::vector<BasisFunc>& basis, int pos1, int pos2, std::vector<std::vector<double>>& grid, bool conv) {
  int bins = grid.size();
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
      grid[i][j] = elt(val);
    }
  }
}

}
}
