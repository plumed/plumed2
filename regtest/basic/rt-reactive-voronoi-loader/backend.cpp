// One retained pair has weight one and zero derivatives. As a GPU backend, the occupancy
// is given only by finishForward(), after forward() or forwardNeighbors().
#include "colvar/ReactiveVoronoiBackend.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

class TestVoronoiBackend final : public PLMD::colvar::ReactiveVoronoiBackend {
  bool queued=false;
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

  void forward(const double* x, const PLMD::Pbc&, bool, double) override {
    if(separation(x)<=0.0) {
      throw InvalidPair(0,0,"test backend forward failed");
    }
    queued=true;
  }

  void forwardNeighbors(const double* x, const PLMD::Pbc& box, bool pbc, double kappa,
                        double cutoff, const std::vector<int>& offsets, const std::vector<int>& centers,
                        unsigned long) override {
    if(offsets.size()!=2 || offsets[0]!=0 || offsets[1]!=static_cast<int>(centers.size()) ||
        centers.size()>1 || (!centers.empty() && centers[0]!=0)) {
      throw std::runtime_error("invalid candidate list");
    }
    if(centers.empty() || separation(x)>cutoff) {
      throw std::runtime_error("no CENTER inside the candidate list");
    }
    forward(x,box,pbc,kappa);
  }

  void finishForward(double* occupancy) override {
    if(!queued) {
      throw std::runtime_error("finishForward without forward");
    }
    queued=false;
    occupancy[0]=1.0;
  }

  void backward(const double*, double* gradient, double* virial) override {
    std::fill_n(gradient,6,0.0);
    std::fill_n(virial,9,0.0);
  }
};

PLMD::colvar::ReactiveVoronoiBackendRegistration<TestVoronoiBackend> registerBackend("CUDA");

}
