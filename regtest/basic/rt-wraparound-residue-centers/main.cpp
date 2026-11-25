#include "plumed/core/PlumedMain.h"
#include "plumed/tools/Pbc.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace PLMD;
constexpr unsigned natoms=160;
using Coordinates=std::vector<Vector>;
using Groups=std::vector<std::vector<unsigned>>;
const Groups mixed= {{1,2,3},{4,5,6,7},{8,9,10},{11,12,13},{14}};
void require(bool condition,const std::string& message) {
  if(!condition) {
    throw std::runtime_error(message);
  }
}
void writeMetadata(const char* name,bool changed=false,bool reused=false) {
  FILE* file=std::fopen(name,"w");
  require(file,"metadata output");
  if(changed) {
    std::fprintf(file,"CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n");
  }
  for(unsigned i=0; i<14; ++i) {
    if(i==10 && !reused) {
      std::fprintf(file,"TER\n");
    }
    const char chain=i<7 || (reused && i>=10 && i<13) ? 'A' : 'B';
    const char insertion=i>=3 && i<7 ? 'A' : ' ';
    const char* residue=i<3 ? "ALA" : i==13 ? " NA" : "HOH";
    // Intentionally unrelated coordinates; no reference geometry should be used.
    const double x=changed ? 11*(int(i%5)-2) : 0;
    std::fprintf(file,"ATOM  %5u  C   %s %c%4d%c   %8.3f%8.3f%8.3f  1.00  1.00           C\n",
                 i+1,residue,chain,i==13 ? 2 : 1,insertion,x,changed ? -2*x : 0,0.0);
  }
  std::fprintf(file,"END\n");
  std::fclose(file);
}
struct Frame {
  Coordinates positions;
  Tensor box;
  int step;
};
struct Context {
  FILE* log=std::tmpfile();
  std::unique_ptr<PlumedMain> p=std::make_unique<PlumedMain>();
  Context(const std::string& metadata) {
    require(log,"temporary log");
    int precision=sizeof(double),n=natoms;
    double timestep=1;
    p->cmd("setRealPrecision",&precision);
    p->cmd("setNatoms",&n);
    p->cmd("setTimestep",&timestep);
    p->cmd("setLog",log);
    p->cmd("init");
    if(!metadata.empty()) {
      line("MOLINFO STRUCTURE="+metadata+" PYTHON_BIN=no");
    }
  }
  void line(const std::string& input) {
    std::istringstream lines(input);
    std::string entry;
    while(std::getline(lines,entry)) {
      p->cmd("readInputLine",entry.c_str());
    }
  }
  ~Context() {
    p.reset();
    std::fclose(log);
  }
};
std::vector<Coordinates> run(const std::string& input,const std::vector<Frame>& frames,
                             const std::string& metadata="metadata.pdb") {
  std::remove("wrapped.xyz");
  {
    Context context(metadata);
    context.line(input);
    context.line("DUMPATOMS ATOMS=1-160 FILE=wrapped.xyz PRECISION=12");
    std::array<double,natoms> masses{},charges{};
    masses.fill(1);
    for(const auto& frame:frames) {
      auto positions=frame.positions;
      auto box=frame.box;
      int step=frame.step;
      Coordinates forces(natoms);
      Tensor virial;
      context.p->cmd("setStep",&step);
      context.p->cmd("setBox",box.data());
      context.p->cmd("setPositions",positions.data()->data());
      context.p->cmd("setMasses",masses.data());
      context.p->cmd("setCharges",charges.data());
      context.p->cmd("setForces",forces.data()->data());
      context.p->cmd("setVirial",virial.data());
      context.p->cmd("calc");
    }
  }
  std::ifstream file("wrapped.xyz");
  std::vector<Coordinates> result;
  for(unsigned f=0; f<frames.size(); ++f) {
    std::string row,label;
    std::getline(file,row);
    require(row=="160","DUMPATOMS atom count");
    std::getline(file,row);
    Coordinates positions(natoms);
    for(auto& x:positions) {
      file>>label>>x[0]>>x[1]>>x[2];
      require(bool(file),"DUMPATOMS positions");
    }
    std::getline(file,row);
    result.push_back(positions);
  }
  return result;
}
Coordinates calculate(const std::string& input,const Frame& frame,
                      const std::string& metadata="metadata.pdb") {
  return run(input, {frame},metadata)[0];
}
void same(const Coordinates& actual,const Coordinates& expected,const std::string& label) {
  require(actual.size()==expected.size(),label+" size");
  for(unsigned i=0; i<actual.size(); ++i) {
    require((actual[i]-expected[i]).modulo()<2e-10,label+" atom "+std::to_string(i+1));
  }
}
// Independent oracle: the physical whole coordinates are known before imaging.
// No ordered reconstruction, production groups or reference coordinates are used.
Coordinates oracle(const Frame& frame,const Coordinates& physical,const Groups& groups,
                   const std::vector<Vector>& references,bool pair=false) {
  Pbc pbc;
  pbc.setBox(frame.box);
  Coordinates result=frame.positions;
  for(unsigned g=0; g<groups.size(); ++g) {
    Vector center;
    for(unsigned atom:groups[g]) {
      center+=physical[atom-1];
    }
    center/=groups[g].size();
    unsigned nearest=pair ? g : 0;
    if(!pair) {
      double best=modulo2(pbc.distance(references[0],center));
      for(unsigned j=1; j<references.size(); ++j) {
        const double distance=modulo2(pbc.distance(references[j],center));
        if(distance<best) {
          nearest=j;
          best=distance;
        }
      }
    }
    const Vector shift=references[nearest]+pbc.distance(references[nearest],center)-center;
    for(unsigned atom:groups[g]) {
      result[atom-1]=physical[atom-1]+shift;
    }
  }
  return result;
}
Coordinates physicalPositions() {
  Coordinates x(natoms);
  const std::vector<Vector> first= {
    {1.8,-0.05,0},{2.1,0,0},{2.2,0.05,0},
    {0.3,0.9,0},{0.4,1,0},{0.5,0.95,0},{0.6,0.85,0},
    {2.8,-0.2,0},{2.9,-0.1,0},{3.1,-0.15,0},
    {-0.6,0.1,0},{-0.5,0.2,0},{-0.4,0.15,0},{3.4,1.3,0}
  };
  std::copy(first.begin(),first.end(),x.begin());
  for(unsigned j=14; j<natoms; ++j) {
    x[j]=Vector(1.7*std::sin(1.31*j),1.5*std::cos(0.71*j),1.6*std::sin(0.43*j));
  }
  x[30]=Vector(0,0,0);
  x[31]=Vector(0.3,1.4,0.1);
  return x;
}
Frame image(const Coordinates& physical,const Tensor& box,int step=0) {
  Frame frame{physical,box,step};
  for(unsigned j=0; j<natoms; ++j) {
    frame.positions[j]+=matmul(Vector(int(j%5)-2,int(j%3)-1,int(j%7)-3),box);
  }
  return frame;
}
std::vector<Vector> refs(const Frame& frame,const std::vector<unsigned>& atoms) {
  std::vector<Vector> result;
  for(unsigned atom:atoms) {
    result.push_back(frame.positions[atom-1]);
  }
  return result;
}
void rejects(const std::string& input,const std::string& message,const std::string& metadata="metadata.pdb") {
  bool rejected=false;
  try {
    Context c(metadata);
    c.line(input);
  } catch(const std::exception& e) {
    rejected=std::string(e.what()).find(message)!=std::string::npos;
  }
  require(rejected,"expected diagnostic: "+message);
}
int main() {
  writeMetadata("metadata.pdb");
  writeMetadata("changed.pdb",true);
  writeMetadata("reused.pdb",false,true);
  std::ofstream output("output");
  for(bool skew: {
        false,true
      }) {
    const Tensor box=skew ? Tensor(4,0,0,0.8,4,0,-0.4,0.5,4) : Tensor(4,0,0,0,4,0,0,0,4);
    const auto physical=physicalPositions();
    const auto frame=image(physical,box);
    const std::string command="WRAPAROUND ATOMS=1-14 AROUND=31,32 GROUP_BY_RESIDUE";
    const auto actual=calculate(command,frame);
    same(actual,oracle(frame,physical,mixed,refs(frame, {31,32})),"mixed centroid oracle");
    same(calculate(command,frame,"changed.pdb"),actual,"reference coordinates and cell ignored");
    same(calculate(command,frame,"reused.pdb"),actual,"nonconsecutive reused residue identifiers");
    // Individual atom wrapping must give a different image for part of batch 1.
    Groups single;
    for(unsigned j=1; j<=14; ++j) {
      single.push_back({j});
    }
    same(calculate("WRAPAROUND ATOMS=1-14 AROUND=31,32",frame),
         oracle(frame,physical,single,refs(frame, {31,32})),"opt-in only");
    same(calculate("WRAPAROUND ATOMS=1,3,1,7,4 AROUND=31 GROUP_BY_RESIDUE GROUPBY=1",frame),
    oracle(frame,physical, {{1,3},{7,4}},refs(frame, {31})),"partial reversed groups and stable duplicates");
    same(calculate("WRAPAROUND ATOMS=1,4,2,5,3,7 AROUND=31,32 GROUP_BY_RESIDUE",frame),
    oracle(frame,physical, {{1},{4},{2},{5},{3},{7}},refs(frame, {31,32})),"separate selected runs");
    same(calculate("WRAPAROUND ATOMS=1,11 AROUND=31,32 GROUP_BY_RESIDUE PAIR",frame,"reused.pdb"),
    oracle(frame,physical, {{1},{11}},refs(frame, {31,32}),true),"reused identifiers in a sparse selection");
    same(calculate("WRAPAROUND ATOMS=1-14 AROUND=31,32,31,32,31 GROUP_BY_RESIDUE PAIR",frame),
         oracle(frame,physical,mixed,refs(frame, {31,32,31,32,31}),true),"PAIR batch order and repeated references");
    const auto overlap=calculate("WRAPAROUND ATOMS=1-14 AROUND=1-14 GROUP_BY_RESIDUE",frame);
    std::vector<unsigned> selected;
    for(unsigned j=1; j<=14; ++j) {
      selected.push_back(j);
    }
    same(overlap,oracle(frame,physical,mixed,refs(frame,selected)),"snapshot overlapping references");
    std::vector<unsigned> many;
    for(unsigned j=33; j<=160; ++j) {
      many.push_back(j);
    }
    same(calculate("WRAPAROUND ATOMS=1-14 AROUND=33-160 GROUP_BY_RESIDUE",frame),
         oracle(frame,physical,mixed,refs(frame,many)),"accelerated cell search versus direct oracle");
    Pbc pbc;
    pbc.setBox(box);
    const Vector com=frame.positions[30]+0.5*pbc.distance(frame.positions[30],frame.positions[31]);
    same(calculate("c: COM ATOMS=31,32\nWRAPAROUND ATOMS=1-14 AROUND=c GROUP_BY_RESIDUE",frame),
         oracle(frame,physical,mixed, {com}),"virtual COM reference without metadata");
    same(calculate("WRAPAROUND AROUND=31 GROUP_BY_RESIDUE",frame),frame.positions,"empty selection");
    auto nextPhysical=physical;
    nextPhysical[1][0]=1.9;
    nextPhysical[2][0]=2.1;
    const Tensor nextBox=box*1.1;
    const auto next=image(nextPhysical,nextBox,1);
    const auto frames=run(command, {frame,next});
    same(frames[0],actual,"first dynamic frame");
    same(frames[1],oracle(next,nextPhysical,mixed,refs(next, {31,32})),"new shape and box on next frame");
    const auto strided=run(command+" STRIDE=2", {frame,next});
    same(strided[0],actual,"active stride");
    same(strided[1],next.positions,"inactive stride");
  }
  output<<"current centroids, mixed sizes, metadata-only grouping and reference independence PASS\n";
  output<<"rectangular/skew cells, images, partial selections, duplicates and identity collisions PASS\n";
  output<<"PAIR, direct/cell reference searches, COM, overlap and unselected atoms PASS\n";
  output<<"dynamic shapes/boxes, STRIDE, empty selections and opt-in behavior PASS\n";
  rejects("WRAPAROUND ATOMS=1-3 AROUND=31 GROUP_BY_RESIDUE","requires a preceding MOLINFO","");
  rejects("WRAPAROUND ATOMS=15 AROUND=31 GROUP_BY_RESIDUE","atom 15 is missing");
  rejects("c: COM ATOMS=1-3\nWRAPAROUND ATOMS=c AROUND=31 GROUP_BY_RESIDUE","is missing");
  rejects("WRAPAROUND ATOMS=1-3 AROUND=31 GROUP_BY_RESIDUE GROUPBY=3","cannot be combined");
  rejects("WRAPAROUND ATOMS=1-14 AROUND=31 GROUP_BY_RESIDUE PAIR","one AROUND entry");
  rejects("WRAPAROUND ATOMS=1-3 AROUND=31,32 GROUPBY=3 PAIR","#ATOMS = #AROUND");
  output<<"invalid metadata, GROUPBY conflicts and PAIR counts PASS\n";
}
