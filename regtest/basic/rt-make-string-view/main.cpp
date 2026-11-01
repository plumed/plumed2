/*
Regression test for command keys passed as std::string_view (C++17).

In Plumed.h, a std::string_view key is copied to a local buffer so that it can be null terminated.
A bug in the buffer handling was erasing the first character of the key whenever it was short enough
to fit the local buffer, so that plumed received an empty key and silently ignored the command.

Keys are passed here in all the ways in which a std::string_view can be built, including
a view over a longer string, which is not null terminated, and a key longer than the local buffer.
*/
#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

using namespace PLMD;

int main() {
  using namespace std::literals;

  const unsigned natoms=2;
  std::vector<double> masses(natoms,1.0);
  // the distance between the two atoms is 3
  std::vector<double> positions= {0.0,0.0,0.0, 1.0,2.0,2.0};
  std::vector<double> forces(3*natoms,0.0);
  double box[3][3]= {{10.0,0.0,0.0},{0.0,10.0,0.0},{0.0,0.0,10.0}};
  double virial[3][3];
  double bias=0.0;
  int step=0;

  Plumed p;

  // string_view built from a literal
  p.cmd("setNatoms"sv,natoms);
  p.cmd("init"sv);

  // string_view built from a std::string
  const std::string readInputLine="readInputLine";
  p.cmd(std::string_view(readInputLine),"d: DISTANCE ATOMS=1,2");
  p.cmd(std::string_view(readInputLine),"RESTRAINT ARG=d AT=0 KAPPA=2");

  // string_view over a longer string: the viewed characters are not followed by a null terminator
  const std::string longer="setStepAndSomethingElse";
  p.cmd(std::string_view(longer).substr(0,7),step);

  p.cmd("setMasses"sv,masses);
  p.cmd("setPositions"sv,positions.data(), {natoms,3});
  p.cmd("setForces"sv,forces.data(), {natoms,3});
  p.cmd("setBox"sv,box);
  p.cmd("setVirial"sv,virial);
  p.cmd("calc"sv);
  p.cmd("getBias"sv,&bias);

  // a key longer than the local buffer: it must reach plumed unmodified (and thus be rejected as unknown)
  const std::string longKey(300,'x');
  bool failed=false;
  try {
    p.cmd(std::string_view(longKey));
  } catch(const Plumed::Exception & e) {
    failed=true;
  }

  std::ofstream os("ooo");
  os<<bias<<"\n";
  for(unsigned i=0; i<natoms; i++) {
    for(unsigned j=0; j<3; j++) {
      os<<" "<<forces[3*i+j];
    }
  }
  os<<"\n";
  os<<(failed?"unknown key rejected":"unknown key accepted")<<"\n";

  return 0;
}
