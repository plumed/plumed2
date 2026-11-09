#include "plumed/tools/Random.h"
#include "plumed/tools/Exception.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace PLMD;

// Tests of PLMD::Random (xoshiro256** engine)
int main() {
  std::ofstream out("output");
  const char* yes[2]= {"no","yes"};

  // 1. engine: reference outputs of xoshiro256** from the state {1,2,3,4}
  // (11520, 0, 1509978240, 1215971899390074240, ...), as U01=((x>>11)+0.5)/2^53
  {
    Random r;
    r.fromString("xoshiro256**|1|2|3|4");
    const double ref[8]= {6.106226635438361e-16,5.5511151231257827e-17,8.1856132982949248e-11,0.065917968750002165,
                          0.065928823519245616,0.0329591103084243,0.87673591144381624,0.45949417701161371
                         };
    bool ok=true;
    for(unsigned i=0; i<8; i++) {
      ok = ok && r.U01()==ref[i];
    }
    out<<"engine matches the reference xoshiro256** outputs: "<<yes[ok]<<"\n";
  }

  // 2. the stream of a seed (fixes the seeding, bit by bit)
  {
    Random r;
    r.setSeed(-2026);
    char buf[64];
    out<<"first numbers with seed -2026:";
    for(unsigned i=0; i<4; i++) {
      std::snprintf(buf,64," %.17g",r.U01());
      out<<buf;
    }
    out<<"\n";
    Random p;
    p.setSeed(2026);
    Random r2;
    r2.setSeed(-2026);
    bool same=true;
    for(unsigned i=0; i<100; i++) {
      same = same && p.U01()==r2.U01();
    }
    out<<"seeds 2026 and -2026 give the same stream: "<<yes[same]<<"\n";
    Random d;
    Random z;
    z.setSeed(0);
    same=true;
    for(unsigned i=0; i<100; i++) {
      same = same && d.U01()==z.U01();
    }
    out<<"default seed is 0: "<<yes[same]<<"\n";
  }

  // 3. moments of U01 and Gaussian
  {
    Random r;
    r.setSeed(7);
    const unsigned n=1000000;
    double s=0.0,s2=0.0,mn=1.0,mx=0.0;
    for(unsigned i=0; i<n; i++) {
      double u=r.U01();
      s+=u;
      s2+=u*u;
      mn=std::min(mn,u);
      mx=std::max(mx,u);
    }
    double m=s/n, v=s2/n-m*m;
    out<<"U01 in (0,1): "<<yes[mn>0.0 && mx<1.0]<<"\n";
    out<<"U01 mean and variance within 5 sigma: "
       <<yes[std::fabs(m-0.5)<5*std::sqrt(1.0/12/n) && std::fabs(v-1.0/12)<5*std::sqrt(1.0/180/n)]<<"\n";
    s=0.0;
    s2=0.0;
    double s4=0.0;
    for(unsigned i=0; i<n; i++) {
      double g=r.Gaussian();
      s+=g;
      s2+=g*g;
      s4+=g*g*g*g;
    }
    m=s/n;
    v=s2/n;
    double k=s4/n;
    out<<"Gaussian mean, variance and 4th moment within 5 sigma: "
       <<yes[std::fabs(m)<5*std::sqrt(1.0/n) && std::fabs(v-1.0)<5*std::sqrt(2.0/n) && std::fabs(k-3.0)<5*std::sqrt(96.0/n)]<<"\n";
  }

  // 4. RandInt: range and uniformity (chi square, 9 degrees of freedom, 99.99% quantile 33.7)
  {
    Random r;
    r.setSeed(11);
    for(int range : {
          1,7,10
        }) {
      std::vector<unsigned> c(range,0);
      const unsigned n=100000;
      bool inrange=true;
      for(unsigned i=0; i<n; i++) {
        int k=r.RandInt(range);
        if(k<0 || k>=range) {
          inrange=false;
        } else {
          c[k]++;
        }
      }
      double chi2=0.0;
      for(auto x : c) {
        chi2+=(x-double(n)/range)*(x-double(n)/range)/(double(n)/range);
      }
      out<<"RandInt("<<range<<") in range: "<<yes[inrange]<<", uniform: "<<yes[chi2<33.7]<<"\n";
    }
    bool thrown=false;
    try {
      r.RandInt(0);
    } catch(const Exception&) {
      thrown=true;
    }
    out<<"RandInt(0) is an error: "<<yes[thrown]<<"\n";
  }

  // 5. exact state round trip, also in the middle of a pair of Gaussian numbers
  {
    Random r;
    r.setSeed(3);
    r.Gaussian();
    std::string full=r.stateToString();
    std::string uniform;
    r.toString(uniform);
    std::vector<double> a;
    for(unsigned i=0; i<10; i++) {
      a.push_back(r.Gaussian());
    }
    Random q;
    q.stateFromString(full);
    bool same=true;
    for(unsigned i=0; i<10; i++) {
      same = same && q.Gaussian()==a[i];
    }
    out<<"stateToString/stateFromString continue the same stream: "<<yes[same]<<"\n";
    Random u;
    u.fromString(uniform);
    Random w;
    w.fromString(uniform);
    same=true;
    for(unsigned i=0; i<10; i++) {
      same = same && u.U01()==w.U01();
    }
    out<<"toString/fromString round trip: "<<yes[same]<<"\n";
    // a state written by the old engine (ran1) is refused
    bool thrown=false;
    try {
      Random old;
      old.fromString("-5|123|1|2|3");
    } catch(const Exception&) {
      thrown=true;
    }
    out<<"old ran1 state is an error: "<<yes[thrown]<<"\n";
  }
  return 0;
}
