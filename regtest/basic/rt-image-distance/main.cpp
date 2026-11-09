#include "plumed/wrapper/Plumed.h"
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <random>
#include <chrono>
#include <iostream>
#include <limits>
#include <algorithm>

using V=std::array<double,3>;
using B=std::array<double,9>;
using X=std::vector<V>;

void require(bool ok,const std::string& message) {
  if(!ok) {
    throw std::runtime_error(message);
  }
}
void close(double a,double b,double tol,const std::string& message) {
  if(std::abs(a-b)>tol || !std::isfinite(a)) {
    std::cerr<<message<<": "<<a<<" expected "<<b<<"\n";
    throw std::runtime_error(message);
  }
}

// Independent direct enumeration in the ORIGINAL basis. Fixtures and random
// cells below have their winning shifts inside this deliberately wide cube.
// This is an oracle for these fixtures, not a general arbitrary-cell algorithm.
double oracle(const X& x,const B& box,int width=12) {
  double best=1e100;
  for(const auto& p:x)
    for(const auto& q:x)
      for(int a=-width; a<=width; ++a)
        for(int b=-width; b<=width; ++b)
          for(int c=-width; c<=width; ++c) {
            if(a==0 && b==0 && c==0) {
              continue;
            }
            double d2=0;
            for(int k=0; k<3; ++k) {
              double d=q[k]-p[k]+a*box[k]+b*box[3+k]+c*box[6+k];
              d2+=d*d;
            }
            best=std::min(best,d2);
          }
  return std::sqrt(best);
}

struct Harness {
  PLMD::Plumed p;
  X forces;
  B virial{};
  std::vector<double> masses;
  int step=0;
  Harness(int n,const std::string& atoms="",bool whole=false,bool wall=false):
    forces(n),masses(n,1) {
    double dt=0.01;
    p.cmd("setNatoms",&n);
    p.cmd("setTimestep",&dt);
    p.cmd("setLogFile","test.log");
    p.cmd("init");
    if(whole) {
      p.cmd("readInputLine","WHOLEMOLECULES ENTITY0=1-3");
    }
    std::string selection=atoms;
    if(selection.empty()) {
      selection=n==1 ? "1" : "1-"+std::to_string(n);
    }
    p.cmd("readInputLine",("d: IMAGE_DISTANCE ATOMS="+selection).c_str());
    p.cmd("readInputLine",wall ? "w: LOWER_WALLS ARG=d AT=1.5 KAPPA=2 EXP=2" : "w: BIASVALUE ARG=d");
  }
  double sample(X x,B box) {
    forces.assign(x.size(),V{});
    virial.fill(0);
    p.cmd("setStep",&step);
    ++step;
    p.cmd("setPositions",x[0].data());
    p.cmd("setMasses",masses.data());
    p.cmd("setForces",forces[0].data());
    p.cmd("setBox",box.data());
    p.cmd("setVirial",virial.data());
    p.cmd("calc");
    double bias=0;
    p.cmd("getBias",&bias);
    return bias;
  }
};

void derivatives(const X& x,const B& box) {
  Harness h(x.size());
  h.sample(x,box);
  const auto f=h.forces;
  const auto v=h.virial;
  const double eps=1e-6;
  for(unsigned i=0; i<x.size(); ++i)
    for(int k=0; k<3; ++k) {
      X xp=x,xm=x;
      xp[i][k]+=eps;
      xm[i][k]-=eps;
      close(f[i][k],-(h.sample(xp,box)-h.sample(xm,box))/(2*eps),2e-7,"atomic derivative");
    }
  for(int k=0; k<3; ++k) {
    double sum=0;
    for(const auto& atom:f) {
      sum+=atom[k];
    }
    close(sum,0,1e-12,"net force");
  }
  // Deform BOTH positions and cell, r'_b=r_b+eps*r_a, including shear.
  for(int a=0; a<3; ++a)
    for(int b=0; b<3; ++b) {
      X xp=x,xm=x;
      B bp=box,bm=box;
      for(unsigned i=0; i<x.size(); ++i) {
        xp[i][b]+=eps*x[i][a];
        xm[i][b]-=eps*x[i][a];
      }
      for(int i=0; i<3; ++i) {
        bp[3*i+b]+=eps*box[3*i+a];
        bm[3*i+b]-=eps*box[3*i+a];
      }
      close(v[3*a+b],(h.sample(xp,bp)-h.sample(xm,bm))/(2*eps),2e-7,"box derivative");
    }
}

int main(int argc,char** argv) {
  std::ofstream out("output");
  const B cubic{{4,0,0,0,4,0,0,0,4}};
  const B rect{{4,0,0,0,5,0,0,0,6}};
  const B skew{{4,0,0,1.3,3.7,0,-0.8,0.7,5.1}};
  const X one{{{0.2,0.3,0.4}}};
  Harness h1(1);
  close(h1.sample(one,rect),4,1e-12,"one atom");
  close(h1.sample(one,cubic),4,1e-12,"cubic tie");
  for(const auto& f:h1.forces)
    for(double v:f) {
      close(v,0,0,"same atom cancellation");
    }
  derivatives(one,rect);
  out<<"self images and box-only derivatives PASS\n";

  const X pair{{{0,0,0}},{{1,0,0}}};
  Harness h2(2);
  close(h2.sample(pair,cubic),3,1e-12,"exclude all zero shifts");
  const X compact{{{0.1,0.2,0.3}},{{1.3,0.5,0.7}},{{0.5,1.1,0.2}}};
  for(const B& box: {
        rect,skew
      }) {
    Harness h(3),reverse(3,"3,1,2");
    const double d=h.sample(compact,box);
    close(d,oracle(compact,box),2e-12,"compact oracle");
    close(d,reverse.sample(compact,box),2e-12,"selection permutation");
    X translated=compact;
    for(auto& p:translated)
      for(int k=0; k<3; ++k) {
        p[k]+=2*box[k]-3*box[3+k]+0.27;
      }
    close(d,h.sample(translated,box),2e-12,"coherent translation");
    derivatives(compact,box);
  }
  out<<"values, permutations, translations, atomic and shear derivatives PASS\n";

  // Correct contact uses shift -3 in the original basis; fixed 26 fails.
  const X extended{{{0,0,0}},{{12.4,0.3,0.2}}};
  close(h2.sample(extended,rect),std::sqrt(0.29),2e-12,"extended group");
  require(oracle(extended,rect,1)>oracle(extended,rect,5)+1,"26-image counterexample");
  derivatives(extended,rect);
  // A shortest lattice vector itself needs a coefficient outside [-1,1].
  const B adversarial{{4,0,0,12,0.7,0,0,0,6}};
  close(h1.sample(one,adversarial),0.7,2e-12,"skew shortest vector");
  close(h1.sample(one,adversarial),oracle(one,adversarial),2e-12,"skew oracle");
  require(oracle(one,adversarial,1)>1,"skew 26-image counterexample");
  out<<"extended and skew counterexamples PASS\n";

  const X wrapped{{{3.8,0.2,0.1}},{{0.3,0.3,0.2}},{{0.9,0.4,0.3}}};
  X whole=wrapped;
  whole[1][0]+=4;
  whole[2][0]+=4;
  Harness hw(3,"",true),hp(3,"3,1,2",true),hd(3);
  const double expected=oracle(whole,rect);
  close(hw.sample(wrapped,rect),expected,2e-12,"WHOLEMOLECULES crossing");
  close(hp.sample(wrapped,rect),expected,2e-12,"whole selection permutation");
  close(hd.sample(whole,rect),expected,2e-12,"whole reference");
  X translated=wrapped;
  for(auto& p:translated) {
    p[0]+=8;
  }
  close(hw.sample(translated,rect),expected,2e-12,"whole translated");
  out<<"boundary reconstruction PASS\n";

  std::mt19937 rng(356);
  std::uniform_real_distribution<double> uni(-1,1);
  for(int trial=0; trial<20; ++trial) {
    B box{{4+uni(rng),0,0,3*uni(rng),4+uni(rng),0,3*uni(rng),3*uni(rng),4+uni(rng)}};
    X x(4);
    for(auto& p:x)
      for(double& v:p) {
        v=8*uni(rng);
      }
    Harness h(4);
    close(h.sample(x,box),oracle(x,box,34),2e-11,"random wide enumeration");
    if(trial<3) {
      derivatives(x,box);
    }
  }
  out<<"20 random triclinic wide-enumeration comparisons PASS\n";

  // Repeated exact ties select the same branch; no finite difference at a cusp.
  const X tie{{{0,0,0}},{{2,0,0}}};
  h2.sample(tie,cubic);
  const auto tieforce=h2.forces;
  for(int i=0; i<5; ++i) {
    close(h2.sample(tie,cubic),2,0,"tie value");
    require(h2.forces==tieforce,"deterministic tie branch");
  }
  X near{{{0,0,0}},{{4-1e-8,0,0}}};
  close(h2.sample(near,rect),1e-8,1e-14,"near contact");
  for(const auto& atom:h2.forces)
    for(double f:atom) {
      require(std::isfinite(f),"finite near-contact force");
    }
  out<<"ties and near contact PASS\n";

  // Actual force propagation through a lower wall: gradient descent expands
  // separation from an image and lowers the wall energy.
  Harness wall(2,"",false,true);
  X x{{{0,0,0}},{{3.1,0.2,0.1}}};
  const double initial=wall.sample(x,rect);
  for(int step=0; step<20; ++step) {
    wall.sample(x,rect);
    for(unsigned i=0; i<x.size(); ++i)
      for(int k=0; k<3; ++k) {
        x[i][k]+=0.01*wall.forces[i][k];
      }
  }
  require(wall.sample(x,rect)<initial*0.1,"wall relaxation");
  out<<"biased relaxation PASS\n";

  // Runtime error paths: use a valid invertible Pbc at initialization, then
  // exercise the action-level checks. Pbc itself owns singular triclinic input.
  for(int which=0; which<6; ++which) {
    bool rejected=false;
    try {
      Harness h(2);
      X bad=pair;
      B box=rect;
      if(which==0) {
        box.fill(0);
      }
      if(which==1)
        bad[1]= {{4,0,0}};
      if(which==2) {
        bad[1][0]=std::numeric_limits<double>::quiet_NaN();
      }
      if(which==3) {
        box[8]=0;
      }
      if(which==4) {
        box[0]=std::numeric_limits<double>::infinity();
      }
      if(which==5) {
        bad[1][0]=1e12;
      }
      h.sample(bad,box);
    } catch(const std::exception&) {
      rejected=true;
    }
    require(rejected,"invalid runtime input rejected");
  }
  bool rejected=false;
  try {
    Harness h(1," ");
  } catch(const std::exception&) {
    rejected=true;
  }
  require(rejected,"empty selection rejected");
  out<<"invalid inputs PASS\n";

  if(argc>1 && std::string(argv[1])=="--benchmark") {
    for(int n: {
          10,100,500,1000
        }) {
      X group(n);
      for(auto& p:group)
        for(double& c:p) {
          c=0.7*uni(rng);
        }
      Harness h(n);
      const auto start=std::chrono::steady_clock::now();
      for(int i=0; i<10; ++i) {
        h.sample(group,rect);
      }
      double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count()/10;
      std::cout<<n<<" atoms: "<<ms<<" ms/frame\n";
    }
  }
}
