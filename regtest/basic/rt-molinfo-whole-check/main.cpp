#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <string>
#include <vector>

// MOLINFO WHOLE takes the reference structure as whole. Here each reference is read with
// MOLINFO ... WHOLE and used by WHOLEMOLECULES, and we report whether PLUMED accepts it or
// reports it as split by the periodic boundary of its CRYST1 record.
// torn.pdb is whole.pdb with 15 atoms wrapped across a face of its (rhombic dodecahedron) cell.
// Its copies without a usable cell (no CRYST1, the 1 A placeholder, a zero cell, a crystallographic
// space group) cannot be checked and are accepted.
int main(int, char**) {
  std::ofstream out("output");
  unsigned natoms=304;
  std::vector<double> positions(3*natoms,0.0);
  const std::vector<std::pair<std::string,std::string>> cases{
    {"whole.pdb","1-304"},
    {"torn.pdb","1-304"},
    {"torn.pdb","1-150"},
    {"torn_nocryst.pdb","1-304"},
    {"torn_placeholder.pdb","1-304"},
    {"torn_zero.pdb","1-304"},
    {"torn_crystal.pdb","1-304"}
  };
  for(const auto & c : cases) {
    PLMD::Plumed plumed;
    plumed.cmd("setNatoms",&natoms);
    plumed.cmd("setLogFile","test.log");
    plumed.cmd("init");
    out<<c.first<<" ENTITY0="<<c.second<<" : ";
    try {
      plumed.cmd("readInputLine",("MOLINFO STRUCTURE="+c.first+" WHOLE").c_str());
      plumed.cmd("readInputLine",("WHOLEMOLECULES ENTITY0="+c.second).c_str());
      out<<"accepted\n";
    } catch(PLMD::Plumed::ExceptionError &e) {
      const std::string m(e.what());
      const auto a=m.find("is not whole: ");
      const auto b=m.find("Atoms:");
      const auto e2=m.find(". This happens");
      if(a==std::string::npos || b==std::string::npos || e2==std::string::npos) {
        out<<"unexpected error: "<<m<<"\n";
      } else {
        out<<"not whole, "<<m.substr(a+14,m.find(" atoms",a)-a-14)<<" atoms;"<<m.substr(b+6,e2-b-6)<<"\n";
      }
    }
  }
  return 0;
}
