#ifndef PLUMED_REGTEST_IMAGE_DISTANCE_TOOLS_H
#define PLUMED_REGTEST_IMAGE_DISTANCE_TOOLS_H
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
#include <functional>

using V=std::array<double,3>;
using B=std::array<double,9>;
using X=std::vector<V>;

inline void require(bool ok,const std::string& message) {
  if(!ok) {
    throw std::runtime_error(message);
  }
}
inline void close(double a,double b,double tol,const std::string& message) {
  if(std::abs(a-b)>tol || !std::isfinite(a)) {
    std::cerr<<message<<": "<<a<<" expected "<<b<<"\n";
    throw std::runtime_error(message);
  }
}

// Independent direct enumeration in the ORIGINAL basis. Fixtures and random
// cells below have their winning shifts inside this deliberately wide cube.
// This is an oracle for these fixtures, not a general arbitrary-cell algorithm.
inline double oracle(const X& x,const B& box,int width=12) {
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
  Harness(int n,const std::string& atoms="",bool whole=false,bool wall=false,
          const std::string& options="",
          const std::function<void(PLMD::Plumed&)>& setup= {},
          const std::string& logfile="test.log",
          const std::string& action="IMAGE_DISTANCE"):
    forces(n),masses(n,1) {
    double dt=0.01;
    p.cmd("setNatoms",&n);
    p.cmd("setTimestep",&dt);
    if(setup) {
      setup(p);
    }
    p.cmd("setLogFile",logfile.c_str());
    p.cmd("init");
    if(whole) {
      p.cmd("readInputLine","WHOLEMOLECULES ENTITY0=1-3");
    }
    std::string selection=atoms;
    if(selection.empty()) {
      selection=n==1 ? "1" : "1-"+std::to_string(n);
    }
    p.cmd("readInputLine",("d: "+action+" ATOMS="+selection+" "+options).c_str());
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

inline void derivatives(const X& x,const B& box) {
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


#endif
