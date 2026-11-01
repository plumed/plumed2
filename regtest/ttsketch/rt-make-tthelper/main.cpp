// Unit tests of the ttsketch helpers:
// - BasisFunc::evalAll() against operator() and grad(), and grad() against finite differences
// - the analytical integrals int0/int1/int2 and the Gram matrix against a quadrature
// - DenseTT::eval() against the ITensor implementation ttEvalAndGrad()
// - the HDF5 helpers (ttWrite/ttRead/ttCount/ttMomentsWrite/ttMomentsRead)
// Only pass/fail is printed, so that the reference does not depend on rounding.
#include "plumed/ttsketch/TTHelper.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <random>

using namespace PLMD;
using namespace PLMD::ttsketch;

static double quad(const std::function<double(double)>& f, double a, double b) {
  // composite Simpson rule
  const int n = 20000;
  const double h = (b - a) / n;
  double s = f(a) + f(b);
  for(int i = 1; i < n; ++i) {
    s += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
  }
  return s * h / 3.0;
}

static const char* ok(bool b) {
  return b ? "ok" : "FAILED";
}

int main() {
  std::ofstream ofs("output");
  const std::pair<double, double> dom(-1.3, 2.1);
  struct Case {
    const char* name;
    int nbasis;
    double w;
    bool kernel;
    double dx;
    bool conv;
  };
  const Case cases[] = {{"fourier", 11, 0.0, false, 0.0, false},
    {"fourier-conv", 11, 0.3, false, 0.0, true},
    {"kernel", 9, 0.0, true, 0.35, false},
    {"kernel-default-dx", 9, 0.0, true, 0.0, false}
  };
  for(const auto& c : cases) {
    BasisFunc b(dom, c.nbasis, c.w, c.kernel, c.dx);
    const int n = b.nbasis();
    std::vector<double> phi(n), dphi(n);
    double eall = 0.0, efd = 0.0;
    // also outside the domain
    for(double x = -3.0; x < 4.0; x += 0.0731) {
      b.evalAll(x, c.conv, phi.data(), dphi.data());
      for(int p = 1; p <= n; ++p) {
        eall = std::max(eall, std::abs(phi[p - 1] - b(x, p, c.conv)));
        eall = std::max(eall, std::abs(dphi[p - 1] - b.grad(x, p, c.conv)));
        const double h = 1e-6;
        const double fd = (b(x + h, p, c.conv) - b(x - h, p, c.conv)) / (2 * h);
        efd = std::max(efd, std::abs(fd - b.grad(x, p, c.conv)) / (1.0 + std::abs(fd)));
      }
    }
    double eint = 0.0, egram = 0.0;
    for(int p = 1; p <= n; ++p) {
      eint = std::max(eint, std::abs(quad([&](double x) {
        return b(x, p, false);
      }, dom.first, dom.second) - b.int0(p)));
      eint = std::max(eint, std::abs(quad([&](double x) {
        return x * b(x, p, false);
      }, dom.first, dom.second) - b.int1(p)));
      eint = std::max(eint, std::abs(quad([&](double x) {
        return x * x * b(x, p, false);
      }, dom.first, dom.second) - b.int2(p)));
      if(c.kernel) {
        for(int q = 1; q <= n; ++q) {
          egram = std::max(egram, std::abs(quad([&](double x) {
            return b(x, p, false) * b(x, q, false);
          }, dom.first, dom.second) - b.gram()(p - 1, q - 1)));
        }
      }
    }
    ofs << c.name << " evalAll " << ok(eall < 1e-10) << "\n";
    ofs << c.name << " grad " << ok(efd < 1e-6) << "\n";
    ofs << c.name << " integrals " << ok(eint < 1e-8) << "\n";
    ofs << c.name << " gram " << ok(egram < 1e-8) << "\n";
  }

  // Gram pseudo-inverse of a nearly dependent kernel basis (default width = spacing: the
  // constant function is almost a combination of the kernels): one direction is left out,
  // and what is left is a proper pseudo-inverse, G ginv G = G
  {
    BasisFunc b(dom, 9, 0.0, true, 0.0);
    const int n = b.nbasis();
    Matrix<double> t1, t2;
    mult(b.gram(), b.ginv(), t1);
    mult(t1, b.gram(), t2);
    double e = 0.0;
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        e = std::max(e, std::abs(t2(i, j) - b.gram()(i, j)));
      }
    }
    ofs << "kernel gram pseudo-inverse: dropped " << ok(b.gramDropped() == 1) << " condition " << ok(b.gramCondition() > 1e12)
        << " G*ginv*G=G " << ok(e < 1e-10) << "\n";
  }

  // a random TT with three dimensions and different ranks
  std::mt19937 gen(7);
  std::normal_distribution<double> dist(0.0, 1.0);
  const int n = 7;
  auto sites = itensor::SiteSet(3, n);
  const itensor::Index l1(3, "Link,l=1"), l2(4, "Link,l=2");
  itensor::MPS tt(3);
  tt.ref(1) = itensor::ITensor(sites(1), l1);
  tt.ref(2) = itensor::ITensor(l1, sites(2), l2);
  tt.ref(3) = itensor::ITensor(l2, sites(3));
  for(int j = 1; j <= n; ++j) {
    for(int a = 1; a <= 3; ++a) {
      tt.ref(1).set(sites(1) = j, l1 = a, dist(gen));
      for(int b = 1; b <= 4; ++b) {
        tt.ref(2).set(l1 = a, sites(2) = j, l2 = b, dist(gen));
      }
    }
    for(int b = 1; b <= 4; ++b) {
      tt.ref(3).set(l2 = b, sites(3) = j, dist(gen));
    }
  }
  for(int kernel = 0; kernel < 2; ++kernel) {
    for(int conv = 0; conv < 2 - kernel; ++conv) {
      std::vector<BasisFunc> basis;
      basis.emplace_back(std::make_pair(-M_PI, M_PI), n, conv ? 0.2 : 0.0, kernel, 0.0);
      basis.emplace_back(std::make_pair(0.0, 3.0), n, conv ? 0.1 : 0.0, kernel, 0.0);
      basis.emplace_back(std::make_pair(-1.0, 1.0), n, conv ? 0.05 : 0.0, kernel, 0.0);
      const DenseTT dense(tt);
      double ev = 0.0, eg = 0.0;
      for(int k = 0; k < 20; ++k) {
        std::vector<double> x = {-3.0 + 0.3 * k, 0.15 * k, -1.0 + 0.1 * k};
        auto [val, grad] = ttEvalAndGrad(tt, basis, x, conv);
        std::vector<double> g;
        const double v = dense.eval(basis, x, conv, &g);
        ev = std::max(ev, std::abs(v - val) / (1.0 + std::abs(val)));
        for(int i = 0; i < 3; ++i) {
          eg = std::max(eg, std::abs(g[i] - grad[i]) / (1.0 + std::abs(grad[i])));
        }
        ev = std::max(ev, std::abs(dense.eval(basis, x, conv, nullptr) - val) / (1.0 + std::abs(val)));
      }
      ofs << "DenseTT kernel=" << kernel << " conv=" << conv << " value " << ok(ev < 1e-12) << " gradient " << ok(eg < 1e-12) << "\n";
    }
  }

  // HDF5 helpers
  const std::string fname = "tt.h5";
  ofs << "ttCount missing file " << ok(ttCount(fname) == 0) << "\n";
  auto tt2 = tt;
  tt2.ref(2) *= 2.0;
  ttWrite(fname, tt, 2);
  ttWrite(fname, tt2, 3);
  ofs << "ttCount " << ok(ttCount(fname) == 2) << "\n";
  auto back = ttRead(fname, 3);
  double eread = 0.0;
  for(int i = 1; i <= 3; ++i) {
    eread = std::max(eread, itensor::norm(back(i) - tt2(i)));
  }
  ofs << "ttRead " << ok(eread == 0.0) << "\n";
  std::vector<itensor::ITensor> A = {tt(1) * itensor::prime(tt(1), l1), tt(3) * itensor::prime(tt(3), l2)};
  itensor::MPS B;
  std::vector<itensor::ITensor> Aback;
  ofs << "ttMomentsRead missing " << ok(!ttMomentsRead(fname, 3, 2, B, Aback)) << "\n";
  ttMomentsWrite(fname, tt2, A, 3);
  const bool found = ttMomentsRead(fname, 3, 2, B, Aback);
  double emom = 0.0;
  for(int i = 1; found && i <= 3; ++i) {
    emom = std::max(emom, itensor::norm(B(i) - tt2(i)));
  }
  for(unsigned k = 0; found && k < A.size(); ++k) {
    emom = std::max(emom, itensor::norm(Aback[k] - A[k]));
  }
  ofs << "ttMomentsRead " << ok(found && emom == 0.0) << "\n";
  return 0;
}
