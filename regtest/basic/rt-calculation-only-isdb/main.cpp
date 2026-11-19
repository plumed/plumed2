#include "plumed/wrapper/Plumed.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Result {
  std::vector<double> bias;
  std::vector<std::vector<double>> forces;
  std::string history;
};

static std::vector<double> structure() {
  std::ifstream file("structure.pdb");
  if(!file) {
    throw std::runtime_error("Missing structure.pdb");
  }
  std::vector<double> x;
  std::string line;
  while(std::getline(file,line))
    if(line.compare(0,4,"ATOM")==0) {
      const unsigned atom=std::stoul(line.substr(6,5))-1;
      x.resize(std::max(x.size(),3*static_cast<std::size_t>(atom+1)));
      for(unsigned k=0; k<3; ++k) {
        x[3*atom+k]=std::stod(line.substr(30+8*k,8))*.1;
      }
    }
  return x;
}

static Result run(const std::string& name,int kind,int mode,int stride) {
  Result out;
  {
    PLMD::Plumed p;
    const auto original=kind==2?structure():std::vector<double> {0,0,0,1.2,0,0,0,1,0,1.5,1,0};
    const int natoms=original.size()/3;
    p.cmd("GREX init");
    p.cmd("setNatoms",natoms);
    p.cmd("setTimestep",.002);
    p.cmd("setLogFile",(name+".log").c_str());
    p.cmd("init");
    if(kind<2) {
      p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
    }
    if(kind==0) {
      p.cmd("readInputLine",("m: METAINFERENCE ARG=d PARAMETERS=1.0 NOISETYPE=MGAUSS SIGMA0=0.5 SIGMA_MIN=0.1 SIGMA_MAX=1.0 DSIGMA=0.2 SIGMA_MEAN0=0.1 MC_STEPS=3 TEMP=300 SEED=11 WRITE_STRIDE=0 STRIDE="+std::to_string(stride)).c_str());
    } else if(kind==1) {
      p.cmd("readInputLine","r: RESTRAINT ARG=d AT=1 KAPPA=20");
      p.cmd("readInputLine","SELECTOR NAME=gamma VALUE=0");
      p.cmd("readInputLine",("s: RESCALE ARG=r.bias SELECTOR=gamma MAX_RESCALE=2 NBIN=4 W0=0.4 BIASFACTOR=5 BSTRIDE=1 BFILE="+name+".bias MC_STEPS=3 MC_STRIDE=1 PACE=1 SEED=1234 TEMP=300 STRIDE="+std::to_string(stride)).c_str());
    } else {
      p.cmd("readInputLine","MOLINFO STRUCTURE=structure.pdb");
      p.cmd("readInputLine","heavy: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H");
      p.cmd("readInputLine",("e: EMMI ATOMS=heavy GMM_FILE=map.dat ERR_FILE=err.dat NOPBC TEMP=310 NL_STRIDE=1 NL_CUTOFF=0.01 WRITE_STRIDE=10000 SEED=3 STATUS_FILE="+name+".status SIGMA0=0.5 SIGMA_MIN=0.2 DSIGMA=0.1 MC_STRIDE=1 RESOLUTION=1 NOISETYPE=OUTLIERS REGRESSION=4 REG_SCALE_MIN=100 REG_SCALE_MAX=400 REG_DSCALE=5").c_str());
      p.cmd("readInputLine",("BIASVALUE ARG=e.scoreb STRIDE="+std::to_string(stride)).c_str());
    }
    std::vector<double> x(original),forces(x.size()),masses(natoms,1.0);
    std::array<double,9> box{{10,0,0,0,10,0,0,0,10}},virial{};
    auto bind=[&](int step,bool foreign) {
      for(unsigned i=0; i<x.size(); ++i) {
        x[i]=original[i]+.002*std::sin(.3*step+i)+(foreign?.07*std::cos(i):0);
      }
      p.cmd("setStep",step);
      p.cmd("setPositions",x.data());
      p.cmd("setMasses",masses.data());
      p.cmd("setBox",box.data());
      p.cmd("setForces",forces.data());
      p.cmd("setVirial",virial.data());
    };
    auto trial=[&](int step,bool foreign) {
      bind(step,foreign);
      p.cmd("setExchangeStep",1);
      p.cmd("prepareCalc");
      if(!foreign) {
        p.cmd("setExchangeStep",0);
      }
      p.cmd("performCalcNoForces");
      p.cmd("setExchangeStep",0);
    };
    for(int step=0; step<24; ++step) {
      forces.assign(forces.size(),0);
      virial.fill(0);
      const bool exchange=mode==4 ? step>0 && step%stride==0 : step%3==1;
      if(mode==2 && exchange) {
        trial(step,true);
      }
      bind(step,false);
      p.cmd("calc");
      double bias=0;
      p.cmd("getBias",&bias);
      out.bias.push_back(bias);
      out.forces.push_back(forces);
      if(mode==1)
        for(int i=0; i<3; ++i) {
          trial(step,false);
        }
      if(mode==3 && exchange) {
        trial(step,true);
      }
      if(mode==4 && exchange) {
        p.cmd("GREX savePositions");
      }
      if(forces!=out.forces.back()) {
        throw std::runtime_error("Trial changed normal forces");
      }
    }
  }
  if(kind==1) {
    std::ifstream file(name+".0.bias");
    if(!file) {
      throw std::runtime_error("Missing RESCALE bias file");
    }
    std::ostringstream text;
    text<<file.rdbuf();
    out.history=text.str();
  }
  return out;
}

int main() {
  std::ofstream output("output");
  for(int kind=0; kind<3; ++kind) {
    const std::string label=kind==0?"METAINFERENCE":kind==1?"RESCALE":"EMMI";
    for(int mode=1; mode<=4; ++mode) {
      const int stride=mode==4?4:1;
      const std::string name=label+"-"+std::to_string(mode);
      const auto reference=run(name+"-reference",kind,0,stride);
      const auto actual=run(name,kind,mode,stride);
      double error=0;
      for(unsigned step=0; step<reference.bias.size(); ++step) {
        if(!std::isfinite(actual.bias[step])) {
          throw std::runtime_error(name+": nonfinite bias");
        }
        // With no active action the kernel leaves getBias at its last value;
        // only the scheduled normal calculation defines a fresh bias energy.
        if(step%stride==0) {
          error=std::max(error,std::abs(reference.bias[step]-actual.bias[step]));
        }
        for(unsigned i=0; i<reference.forces[step].size(); ++i) {
          if(!std::isfinite(actual.forces[step][i])) {
            throw std::runtime_error(name+": nonfinite force");
          }
          error=std::max(error,std::abs(reference.forces[step][i]-actual.forces[step][i]));
        }
      }
      if(error>1e-10 || reference.history!=actual.history) {
        throw std::runtime_error(name+": read-only trial changed normal sampling, energy, forces, or bias history");
      }
      output<<name<<": pass\n";
    }
  }
}
