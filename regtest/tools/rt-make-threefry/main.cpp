#include "plumed/tools/Threefry.h"
#include <cinttypes>
#include <cstdio>

using namespace PLMD;

// Threefry2x64-20 must reproduce the known-answer vectors of Random123
// (tests/kat_vectors: counter, key, expected output), and the stream must
// return the blocks (c0,0), (c0,1), ... in order
int main() {
  FILE* out=std::fopen("output","w");
  const Threefry2x64::Block kat[3][2]= {
    {{UINT64_C(0x0000000000000000),UINT64_C(0x0000000000000000)},{UINT64_C(0x0000000000000000),UINT64_C(0x0000000000000000)}},
    {{UINT64_C(0xffffffffffffffff),UINT64_C(0xffffffffffffffff)},{UINT64_C(0xffffffffffffffff),UINT64_C(0xffffffffffffffff)}},
    {{UINT64_C(0x243f6a8885a308d3),UINT64_C(0x13198a2e03707344)},{UINT64_C(0xa4093822299f31d0),UINT64_C(0x082efa98ec4e6c89)}}
  };
  for(const auto & v : kat) {
    const auto r=Threefry2x64::block(v[0],v[1]);
    std::fprintf(out,"%016" PRIx64 " %016" PRIx64 "\n",r[0],r[1]);
  }
  Threefry2x64 rng(7,3);
  rng.restart(42);
  bool ok=true;
  for(std::uint64_t c1=0; c1<5; c1++) {
    const auto r=Threefry2x64::block({42,c1}, {7,3});
    ok = ok && rng.next()==r[0] && rng.next()==r[1];
  }
  rng.restart(42);
  ok = ok && rng.next()==Threefry2x64::block({42,0}, {7,3})[0];
  double umin=1.0, umax=0.0;
  for(unsigned i=0; i<100000; i++) {
    const double u=rng.U01();
    umin = u<umin ? u : umin;
    umax = u>umax ? u : umax;
  }
  ok = ok && umin>0.0 && umax<1.0;
  std::fprintf(out,"stream %s\n",ok ? "OK" : "WRONG");
  std::fclose(out);
  return 0;
}
