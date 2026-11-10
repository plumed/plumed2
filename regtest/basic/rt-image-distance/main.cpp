#include "test-tools.h"

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
  const X tie{{{0,0,0}},{{1,1,0}}};
  h2.sample(tie,cubic);
  const auto tieforce=h2.forces;
  const auto tievirial=h2.virial;
  for(int i=0; i<5; ++i) {
    close(h2.sample(tie,cubic),std::sqrt(10.0),1e-12,"tie value");
    require(h2.forces==tieforce && h2.virial==tievirial,"deterministic tie branch");
  }
  X left=tie,right=tie;
  left[1][0]-=1e-6;
  right[1][0]+=1e-6;
  close(h2.sample(left,cubic),std::sqrt(10.0),2e-6,"left limit at switch");
  const auto leftforce=h2.forces;
  close(h2.sample(right,cubic),std::sqrt(10.0),2e-6,"right limit at switch");
  require(leftforce!=h2.forces,"distinct branches across the cusp");
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
