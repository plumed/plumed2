// One retained pair has weight one and zero derivatives.
#include "colvar/ReactiveVoronoiBackend.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

class TestVoronoiBackend final : public PLMD::colvar::ReactiveVoronoiBackend {
  static double separation(const double* x) {
    return std::hypot(x[3]-x[0],x[4]-x[1],x[5]-x[2]);
  }

public:
  TestVoronoiBackend(int centers, int assigned, int device) {
    if(centers!=1 || assigned!=1 || device>0) {
      throw std::runtime_error("test backend context creation failed");
    }
  }

  int device() const override {
    return 0;
  }

  void forward(const double* x, const double*, bool, double, double* occupancy) override {
    if(separation(x)<=0.0) {
      throw std::runtime_error("test backend forward failed");
    }
    occupancy[0]=1.0;
  }

  void forwardNeighbors(const double* x, const double* box, bool pbc, double kappa,
                        double* occupancy, double cutoff,
                        const std::vector<int>& offsets, const std::vector<int>& centers,
                        unsigned long) override {
    if(offsets.size()!=2 || offsets[0]!=0 || offsets[1]!=static_cast<int>(centers.size()) ||
        centers.size()>1 || (!centers.empty() && centers[0]!=0)) {
      throw std::runtime_error("invalid candidate list");
    }
    if(centers.empty() || separation(x)>cutoff) {
      throw std::runtime_error("no CENTER inside the candidate list");
    }
    forward(x,box,pbc,kappa,occupancy);
  }

  void backward(const double*, double* gradient, double* virial) override {
    std::fill_n(gradient,6,0.0);
    std::fill_n(virial,9,0.0);
  }
};

PLMD::colvar::ReactiveVoronoiBackendRegistration<TestVoronoiBackend> registerBackend("CUDA");

}
