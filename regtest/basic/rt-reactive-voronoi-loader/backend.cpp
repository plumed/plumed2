// Test fixture for one center and one assigned atom. The soft assignment is
// exactly one and its derivatives vanish whenever the pair is retained.
#include <algorithm>
#include <cmath>

struct Context {
  bool candidate=false;
};

extern "C" void* reactive_voronoi_create_v1(int centers, int assigned) {
  return centers==1 && assigned==1 ? new Context : nullptr;
}

extern "C" void reactive_voronoi_destroy_v1(void* state) {
  delete static_cast<Context*>(state);
}

extern "C" int reactive_voronoi_forward_v1(void*, const double* x, const double*,
    int, double, double* occupancy) {
  if(std::hypot(x[3]-x[0],x[4]-x[1],x[5]-x[2])<=0.0) {
    return 1;
  }
  occupancy[0]=1.0;
  return 0;
}

extern "C" int reactive_voronoi_backward_v1(void*, const double*, double* gradient, double* virial) {
  std::fill_n(gradient,6,0.0);
  std::fill_n(virial,9,0.0);
  return 0;
}

#ifndef LEGACY_VORONOI_ABI
extern "C" void* reactive_voronoi_create_v2(int centers, int assigned, int device) {
  return device<=0 ? reactive_voronoi_create_v1(centers,assigned) : nullptr;
}

extern "C" int reactive_voronoi_device_count_v1() {
  return 1;
}

extern "C" int reactive_voronoi_device_v1(void*) {
  return 0;
}

extern "C" int reactive_voronoi_forward_nlist_v1(void* state, const double* x,
    const double* box, int pbc, double kappa, double* occupancy, double cutoff,
    const int* offsets, const int* centers, int count) {
  auto& context=*static_cast<Context*>(state);
  if(offsets) {
    const bool valid=offsets[0]==0 && offsets[1]==count &&
                     (count==0 || (count==1 && centers && centers[0]==0));
    if(!valid) {
      return 1;
    }
    context.candidate=count==1;
  } else if(centers || count!=0) {
    return 1;
  }
  if(!context.candidate || std::hypot(x[3]-x[0],x[4]-x[1],x[5]-x[2])>cutoff) {
    return 2;
  }
  return reactive_voronoi_forward_v1(state,x,box,pbc,kappa,occupancy);
}
#endif
