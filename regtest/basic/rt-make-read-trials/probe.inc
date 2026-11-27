#define __PLUMED_WRAPPER_LIBCXX11 1
#define __PLUMED_WRAPPER_LIBCXX17 1
#include "plumed/wrapper/Plumed.h"
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
int runProbe(int argc,char**argv) {
  if(argc!=3) {
    return 2;
  }
  std::string action=argv[1];
  int mode=std::stoi(argv[2]);
  PLMD::Plumed p;
  p.cmd("setNatoms",3);
  p.cmd("setTimestep",.001);
  p.cmd("setKbT",2.5);
  p.cmd("setLogFile","plumed.log");
  p.cmd("init");
  p.cmd("readInputLine","d: DISTANCE ATOMS=1,2 NOPBC");
  if(action=="read") {
    std::ofstream input("input.dat");
    input<<"#! FIELDS time value\n";
    for(int i=0; i<40; ++i) {
      input<<i*.001<<" "<<i+1<<"\n";
    }
    input.close();
    p.cmd("readInputLine","r: READ FILE=input.dat VALUES=value IGNORE_FORCES");
    p.cmd("readInputLine","c: CUSTOM ARG=d,r FUNC=x*y PERIODIC=NO");
    p.cmd("readInputLine","a: RESTRAINT ARG=c AT=0 KAPPA=1");
    p.cmd("readInputLine","PRINT ARG=r,a.bias FILE=history FMT=%.17g");
  } else if(action=="rescale") {
    p.cmd("readInputLine","r: RESTRAINT ARG=d AT=0 KAPPA=1");
    p.cmd("readInputLine","SELECTOR NAME=gamma VALUE=0");
    p.cmd("readInputLine","a: RESCALE ARG=r.bias SELECTOR=gamma MAX_RESCALE=1.4 NBIN=4 W0=1 BIASFACTOR=5 BSTRIDE=1 BFILE=bias.dat MC_STEPS=5 MC_STRIDE=1 PACE=1 TEMP=300 SEED=1234");
    p.cmd("readInputLine","PRINT ARG=a.igamma,a.wtbias,a.bias FILE=history FMT=%.17g");
  } else if(action=="caliber") {
    std::ofstream out("target.dat");
    for(int i=0; i<30; ++i) {
      out<<i*.001<<" 0.8\n";
    }
    out.close();
    p.cmd("readInputLine","a: CALIBER ARG=d KAPPA=10 FILE=target.dat REGRES_ZERO=2 NOENSEMBLE");
    p.cmd("readInputLine","PRINT ARG=a.scale,a.bias FILE=history FMT=%.17g");
  } else if(action=="shadow") {
    p.cmd("readInputLine","s: SHADOW ATOMS=1-3 UPDATE=2 NOPBC REFERENCE");
    p.cmd("readInputLine","a: RESTRAINT ARG=s AT=0 KAPPA=1");
    p.cmd("readInputLine","PRINT ARG=s,a.bias FILE=history FMT=%.17g");
  } else if((action=="maze-bias" || action=="maze-no-print")) {
    p.cmd("readInputLine","MAZE_LOSS LABEL=l PARAMS=1,1,1");
    p.cmd("readInputLine","MAZE_STEERED_MD LABEL=opt LOSS=l PULLING=1,0,0 OPTIMIZER_STRIDE=2 LIGAND=2 PROTEIN=1,3");
    p.cmd("readInputLine","pos: POSITION ATOM=2 NOPBC");
    p.cmd("readInputLine","a: MAZE_OPTIMIZER_BIAS ARG=pos.x,pos.y,pos.z BIASING_RATE=0.02 ALPHA=7.2 OPTIMIZER=opt");
    if(action=="maze-no-print") {
      p.cmd("readInputLine","PRINT ARG=a.bias,a.tdist FILE=history FMT=%.17g");
    } else {
      p.cmd("readInputLine","PRINT ARG=opt.x,opt.y,opt.z,a.bias,a.tdist FILE=history FMT=%.17g");
    }
  } else if(action=="voronoi") {
    p.cmd("readInputLine","c: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3 KAPPA=2 REFERENCE=0.5 POWER=2 NLIST NL_CUTOFF=1.1 NL_STRIDE=2 NOPBC");
    p.cmd("readInputLine","a: BIASVALUE ARG=c");
  } else {
    return 3;
  }
  std::array<double,9>x{},f{},box{{10,0,0,0,10,0,0,0,10}},v{};
  std::array<double,3> mass{{1,1,1}},charge{};
  auto bind=[&](int step,bool foreign) {
    x.fill(0.);
    f.fill(0.);
    v.fill(0.);
    x[3]=foreign?1.:2.+.0003*step;
    x[7]=1.;
    if(action=="voronoi") {
      x[3]=1.;
      x[6]=foreign?1.2:.2;
      x[7]=0.;
    }
    p.cmd("setStep",step);
    p.cmd("setPositions",x.data());
    p.cmd("setMasses",mass.data());
    p.cmd("setCharges",charge.data());
    p.cmd("setForces",f.data());
    p.cmd("setBox",box.data());
    p.cmd("setVirial",v.data());
  };
  std::cout<<std::setprecision(17);
  for(int step=0; step<(action=="rescale"?200:12); ++step) {
    bind(step,mode==4);
    if(mode==6) {
      p.cmd("prepareDependencies");
    }
    p.cmd("calc");
    double b=0;
    p.cmd("getBias",&b);
    std::cout<<"normal "<<step<<" "<<b<<" "<<f[3]<<"\n";
    if(step==2 && mode>0 && mode!=4) {
      bind(step,mode!=2);
      if(mode!=3) {
        p.cmd("setExchangeStep",1);
      }
      p.cmd("prepareCalc");
      p.cmd("performCalcNoForces");
      p.cmd("setExchangeStep",0);
      p.cmd("getBias",&b);
      std::cout<<"trial "<<step<<" "<<b<<"\n";
    }
  }
  return 0;
}
#include <filesystem>
#include <sstream>
#include <vector>
#include <stdexcept>
struct ProbeResult {
  std::vector<std::string> normal;
  std::string history,biasFile;
};
ProbeResult evaluate(const std::string& action,int mode) {
  namespace fs=std::filesystem;
  const auto savedPath=fs::current_path();
  std::ostringstream capture;
  struct Restore {
    fs::path path;
    std::streambuf* output;
    ~Restore() {
      fs::current_path(path);
      std::cout.rdbuf(output);
    }
  } restore {savedPath,std::cout.rdbuf(capture.rdbuf())};
  const auto directory=action+"-"+std::to_string(mode);
  fs::create_directory(directory);
  fs::current_path(directory);
  std::string name=action,number=std::to_string(mode);
  char program[]="probe";
  char* argv[]= {program,name.data(),number.data()};
  if(runProbe(3,argv)!=0) {
    throw std::runtime_error("probe failed");
  }
  ProbeResult result;
  std::istringstream rows(capture.str());
  std::string line;
  while(std::getline(rows,line))
    if(line.rfind("normal ",0)==0) {
      result.normal.push_back(line);
    }
  auto read=[](const char* file) {
    std::ifstream in(file);
    return std::string(std::istreambuf_iterator<char>(in), {});
  };
  result.history=read("history");
  result.biasFile=read("bias.dat");
  if(result.normal.size()!=std::size_t(action=="rescale"?200:12)) {
    throw std::runtime_error("missing normal rows");
  }
  return result;
}
