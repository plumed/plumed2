#include "plumed/tools/Random.h"
#include <cmath>
#include <fstream>
#include <map>
#include <vector>

using namespace PLMD;

// Random::Shuffle must produce all permutations with the same probability
int main() {
  std::ofstream out("output");
  Random rnd;
  rnd.setSeed(-12345);
  for(unsigned n : {
        2u,3u,4u
      }) {
    std::map<std::vector<unsigned>,unsigned> counts;
    const unsigned ntrials=120000;
    for(unsigned t=0; t<ntrials; t++) {
      std::vector<unsigned> v(n);
      for(unsigned i=0; i<n; i++) {
        v[i]=i;
      }
      rnd.Shuffle(v);
      counts[v]++;
    }
    unsigned nperm=1;
    for(unsigned i=2; i<=n; i++) {
      nperm*=i;
    }
    double maxdev=0.0;
    for(const auto & c : counts) {
      maxdev=std::max(maxdev,std::fabs(c.second/double(ntrials)*nperm-1.0));
    }
    out<<"n="<<n<<" permutations found: "<<counts.size()<<" of "<<nperm
       <<" relative frequencies within 5%: "<<(maxdev<0.05?"yes":"no")<<"\n";
  }
  return 0;
}
