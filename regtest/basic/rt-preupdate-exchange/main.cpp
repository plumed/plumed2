#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Prescribe identical MD configurations in both orders. Trial evaluations must
// not change diffusion widths, ordinary forces, or acceleration statistics.
// Dense trials also exercise an exchange before the preceding one is consumed.
struct Options {
  bool pb=false;
  bool before=true;
  bool precal=false;
  bool auxiliary=false;
  bool acceleration=false;
  int stride=1;
  int pattern=0; // sparse accepted/rejected, consecutive trials, rejected only, none
};

static bool attempt(const Options& o,int step) {
  if(o.pattern==0) {
    return step==60 || step==120 || step==180;
  }
  if(o.pattern==1) {
    return step>=60 && step<160 && step%o.stride==0;
  }
  if(o.pattern==2) {
    return step>=60 && step<160 && step%(3*o.stride)==0;
  }
  return false;
}

static bool accepted(const Options& o,int step) {
  return attempt(o,step) && (o.pattern==0 ? step!=120 : o.pattern==1 && (step/o.stride)%3!=0);
}

static int configuration(const Options& o,int step) {
  int conf=0;
  for(int s=0; s<step; ++s)
    if(accepted(o,s)) {
      conf=1-conf;
    }
  return conf;
}

static std::array<double,9> positions(int conf,int step) {
  std::array<double,9> x{};
  x[3]=1.0+2.0*conf+0.05*std::sin(.3*step+conf)+.03*std::sin(1.7*step);
  x[7]=.8+1.3*conf+.04*std::cos(.4*step+conf)+.02*std::sin(.7*step);
  return x;
}

static std::string slurp(const std::string& name) {
  std::ifstream file(name);
  if(!file) {
    throw std::runtime_error("Missing "+name);
  }
  std::ostringstream text;
  text<<file.rdbuf();
  return text.str();
}

struct Result {
  std::vector<std::array<double,9>> forces;
  std::vector<std::string> hills;
  std::string acceleration;
};

static Result run(const std::string& name,const Options& o,int first=0,int last=240,bool restart=false,bool signal=false) {
  Result result;
  {
    PLMD::Plumed p;
    const int natoms=3;
    const double dt=.002;
    p.cmd("GREX init");
    p.cmd("setNatoms",natoms);
    p.cmd("setTimestep",dt);
    p.cmd("setLogFile",(name+(restart?"-restart.log":".log")).c_str());
    if(restart) {
      p.cmd("setRestart",1);
      p.cmd("setRestartStep",static_cast<long long>(first));
    }
    p.cmd("init");
    p.cmd("readInputLine",("CHECKPOINT FILE="+name+".cpt").c_str());
    p.cmd("readInputLine","d1: DISTANCE ATOMS=1,2 NOPBC");
    p.cmd("readInputLine","d2: DISTANCE ATOMS=1,3 NOPBC");
    double lambda=0,lambdaForce=0;
    p.cmd("setExtraCV lambda",&lambda);
    p.cmd("setExtraCVForce lambda",&lambdaForce);
    p.cmd("readInputLine","l: EXTRACV NAME=lambda");
    p.cmd("readInputLine","v: CUSTOM ARG=d1,l FUNC=x+0.1*y PERIODIC=NO");
    std::string line="m: "+std::string(o.pb?"PBMETAD":"METAD")+" ARG=v"+(o.pb?",d2":"")
                     +" HEIGHT=0.1 PACE=10 TEMP=300 BIASFACTOR=8 STRIDE="+std::to_string(o.stride)
                     +(o.acceleration?" SIGMA=0.2 ACCELERATION":" ADAPTIVE=DIFF SIGMA=20")
                     +(o.precal?" PRECAL_DIFF":"")
                     +" FILE="+name+".HA"+(o.pb?","+name+".HB":"")+" FMT=%14.9f";
    p.cmd("readInputLine",line.c_str());
    if(o.acceleration) {
      p.cmd("readInputLine",("PRINT ARG=m.acc FILE="+name+".acc STRIDE=1 FMT=%.12f").c_str());
    }
    std::array<double,9> x{},f{},virial{};
    const std::array<double,9> box{{10,0,0,0,10,0,0,0,10}};
    const std::array<double,3> masses{{1,1,1}};
    auto bind=[&](int step,const std::array<double,9>& coords) {
      x=coords;
      p.cmd("setStep",step);
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",f.data());
      p.cmd("setVirial",virial.data());
      p.cmd("setExtraCV lambda",&lambda);
      p.cmd("setExtraCVForce lambda",&lambdaForce);
    };
    auto foreign=[&](int step,int conf) {
      bind(step,positions(1-conf,step));
      p.cmd("setExchangeStep",1);
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      p.cmd("setExchangeStep",0);
    };
    for(int step=first; step<last; ++step) {
      const int conf=configuration(o,step);
      const auto local=positions(conf,step);
      f.fill(0);
      virial.fill(0);
      bind(step,local);
      p.cmd("doCheckPoint",int(signal && step==80));
      if(attempt(o,step) && o.before) {
        foreign(step,conf);
        bind(step,local);
      }
      p.cmd("prepareCalc");
      if(attempt(o,step) && !o.before) {
        p.cmd("GREX savePositions");
      }
      p.cmd("performCalcNoUpdate");
      if(attempt(o,step) && o.before) {
        p.cmd("GREX savePositions");
      }
      p.cmd("update");
      const auto normalForces=f;
      if(attempt(o,step) && !o.before) {
        foreign(step,conf);
      }
      if(o.auxiliary) {
        for(int trial : {
              1,2,0
            }) {
          lambda=trial;
          bind(step,local);
          // Full MTS dependency graph; local MCMC parameters, not foreign coordinates.
          p.cmd("setExchangeStep",1);
          p.cmd("prepareCalc");
          p.cmd("setExchangeStep",0);
          p.cmd("performCalcNoForces");
        }
      }
      if(f!=normalForces) {
        throw std::runtime_error("Calculation-only pass changed MD forces");
      }
      result.forces.push_back(f);
    }
  }
  result.hills.push_back(slurp(name+".0.HA"));
  if(o.pb) {
    result.hills.push_back(slurp(name+".0.HB"));
  }
  if(o.acceleration) {
    result.acceleration=slurp(name+".0.acc");
  }
  return result;
}

static void compare(const Result& a,const Result& b,const std::string& label,std::ofstream& output,int offset=0) {
  if(a.hills!=b.hills || a.acceleration!=b.acceleration) {
    throw std::runtime_error(label+": changed history/statistics");
  }
  if(a.forces.size()!=b.forces.size()+offset) {
    throw std::runtime_error(label+": incomplete forces");
  }
  for(unsigned step=0; step<b.forces.size(); ++step) {
    for(unsigned i=0; i<9; ++i) {
      if(std::abs(a.forces[step+offset][i]-b.forces[step][i])>1e-10) {
        throw std::runtime_error(label+": changed normal forces");
      }
    }
  }
  output<<label<<": pass\n";
}

int main(int argc,char** argv) {
  std::ofstream output("output");
  const std::string mode=argc>1 ? argv[1] : "all";
  if(mode=="--write-old-checkpoints" || mode=="--resume-old-checkpoints") {
    for(bool pb : {
          false,true
        }) {
      Options o;
      o.pb=pb;
      o.before=false;
      o.stride=4;
      o.pattern=1;
      const std::string name=pb?"legacy-pb":"legacy-metad";
      if(mode=="--write-old-checkpoints") {
        run(name,o,0,86,false,true);
      } else {
        compare(run(name+"-reference",o),run(name,o,80,240,true),name+"-checkpoint",output,80);
      }
    }
    return 0;
  }
  if(mode!="--acceleration-only")
    for(bool pb : {
          false,true
        }) {
      const std::string kind=pb?"PBMETAD":"METAD";
      for(int stride : {
            1,4,10
          }) {
        for(bool precal : {
              false,true
            }) {
          for(int pattern : {
                0,1
              }) {
            Options a;
            a.pb=pb;
            a.stride=stride;
            a.precal=precal;
            a.pattern=pattern;
            const std::string label=kind+"-stride"+std::to_string(stride)+"-precal"+std::to_string(precal)+"-pattern"+std::to_string(pattern);
            auto pre=run(label+"-pre",a);
            a.before=false;
            compare(pre,run(label+"-post",a),label,output);
          }
        }
      }
      for(int stride : {
            1,4
          }) {
        Options a;
        a.pb=pb;
        a.stride=stride;
        a.pattern=3;
        const std::string label=kind+"-stride"+std::to_string(stride);
        const auto base=run(label+"-base",a);
        a.pattern=2;
        compare(base,run(label+"-rejected",a),label+"-rejected",output);
        a.pattern=3;
        a.auxiliary=true;
        compare(base,run(label+"-aux",a),label+"-MCMC",output);
        a.auxiliary=false;
        a.pattern=1;
        const auto continuous=run(label+"-continuous",a);
        run(label+"-restart",a,0,86,false,true);
        compare(continuous,run(label+"-restart",a,80,240,true),label+"-restart",output,80);
      }
    }
  Options acc;
  acc.pattern=3;
  acc.acceleration=true;
  const auto baseline=run("acc-base",acc);
  acc.auxiliary=true;
  compare(baseline,run("acc-aux",acc),"METAD-acceleration-MCMC",output);
  run("acc-restart",acc,0,86,false,true);
  compare(baseline,run("acc-restart",acc,80,240,true),"METAD-acceleration-restart",output,80);
}
