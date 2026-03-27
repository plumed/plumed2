#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/AtomDistribution.h"

#include "frameGenerator.h"
#include "nlTools.h"

#include <fstream>
#include <iostream>
#include <vector>

using PLMD::test::printNeighbors;

#define check(arg) (((arg)) ? "pass\n" : "not pass\n")

constexpr bool serial = true;
void testSingleList(bool do_pbc, std::ostream& ofs);
void testDoubleList(bool do_pbc, std::ostream& ofs);
void testPairList(bool do_pbc, std::ostream& ofs);

int main(int, char **) {
  {
    std::ofstream ofs("unitTest_SL");
    testSingleList(false,ofs);
  }
  {
    std::ofstream ofs("unitTest_SL_pbc");
    testSingleList(true,ofs);
  }
  {
    std::ofstream ofs("unitTest_DL");
    testDoubleList(false,ofs);
  }
  {
    std::ofstream ofs("unitTest_DL_pbc");
    testDoubleList(true,ofs);
  }
}

void testSingleList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
  frameGenerator md(5*5*5,"sc");
  std::vector<Vector> atoms(5*5*5);
  auto d = AtomDistribution::getAtomDistribution("sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  double cutoff=(mybox[0][0]/5)*1.999;
  std::vector<AtomNumber> indexes(atoms.size());
  std::generate(indexes.begin(),indexes.end(),
                [i=0]() mutable {return AtomNumber().setIndex(i++);});
  auto nl=  NeighborList(indexes,
                         serial,
                         do_pbc,
                         pbc,
                         cm,
                         cutoff,
                         1,
                         true);//doCells
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  auto title = std::string("Single list, pbc ") +((do_pbc)?"on":"off");
  printNeighbors(title,nl,md.size(),ofs);
}

void testDoubleList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
  frameGenerator md(5*5*5-1,"sc");
  std::vector<Vector> atoms(5*5*5-1);
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  double cutoff=(mybox[0][0]/5)*1.999;

  std::vector<AtomNumber> indexesA(atoms.size()/2);
  std::vector<AtomNumber> indexesB(atoms.size()/2);
  std::generate(indexesA.begin(),indexesA.end(),
                [i=0]() mutable {auto x= AtomNumber().setIndex(i); i+=2; return x;}
               );

  std::generate(indexesB.begin(),indexesB.end(),
                [i=1]() mutable {auto x= AtomNumber().setIndex(i); i+=2; return x;}
               );

  auto nl=  NeighborList(indexesA,indexesB,
                         serial,
                         false,//do_pair
                         do_pbc,
                         pbc,
                         cm,
                         cutoff,
                         1,
                         true);//doCells

  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  auto title = std::string("Two lists, pbc ") +((do_pbc)?"on":"off");
  printNeighbors(title,nl,md.size(),ofs);
}
