#include "plumed/wrapper/Plumed.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// MOLINFO WHOLE takes the reference structure as whole. Here each reference is read with
// MOLINFO ... WHOLE and used by WHOLEMOLECULES, and we report whether PLUMED accepts it or
// writes a warning that a chain looks split by the periodic boundary of its CRYST1 record.
// torn.pdb is whole.pdb with 15 atoms wrapped across a face of its (rhombic dodecahedron) cell.
// Its copies without a usable cell (no CRYST1, the 1 A placeholder, a zero cell, a crystallographic
// space group) cannot be checked and are accepted.
// dna_torn.pdb is dna_whole.pdb (chain A of 1BNA) wrapped nucleotide by nucleotide: each nucleotide is
// whole, and the chain is split between them.
int main(int, char**) {
  std::ofstream out("output");
  unsigned natoms=304;
  std::vector<double> positions(3*natoms,0.0);
  const std::vector<std::pair<std::string,std::string>> cases{
    {"whole.pdb","1-304"},
    {"torn.pdb","1-304"},
    {"torn_nocryst.pdb","1-304"},
    {"torn_placeholder.pdb","1-304"},
    {"torn_zero.pdb","1-304"},
    {"torn_crystal.pdb","1-304"},
    {"dna_whole.pdb","1-243"},
    {"dna_torn.pdb","1-243"}
  };
  for(unsigned i=0; i<cases.size(); ++i) {
    const auto & c(cases[i]);
    const std::string logfile="test"+std::to_string(i)+".log";
    out<<c.first<<" ENTITY0="<<c.second<<" : ";
    try {
      PLMD::Plumed plumed;
      plumed.cmd("setNatoms",&natoms);
      plumed.cmd("setLogFile",logfile.c_str());
      plumed.cmd("init");
      plumed.cmd("readInputLine",("MOLINFO STRUCTURE="+c.first+" WHOLE").c_str());
      plumed.cmd("readInputLine",("WHOLEMOLECULES ENTITY0="+c.second).c_str());
    } catch(PLMD::Plumed::ExceptionError &e) {
      out<<"unexpected error: "<<e.what()<<"\n";
      continue;
    }
    // the log is complete once plumed is destroyed
    std::ifstream in(logfile);
    std::stringstream ss;
    ss<<in.rdbuf();
    const std::string m(ss.str());
    const auto a=m.find("may not be whole: ");
    if(a==std::string::npos) {
      out<<"accepted\n";
      continue;
    }
    const auto e=m.find(". WHOLEMOLECULES",a);
    if(e==std::string::npos) {
      out<<"unexpected warning: "<<m.substr(a)<<"\n";
    } else {
      out<<"warning, "<<m.substr(a+18,e-a-18)<<"\n";
    }
  }
  return 0;
}
