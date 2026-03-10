#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/Random.h"

#include "frameGenerator.h"
#include "nlTools.h"

#include <fstream>
#include <iostream>
#include <vector>

using PLMD::test::printNeighbors;

constexpr bool serial = true;
void testNoList(std::ostream& ofs);
void testNoList_partial(std::ostream& ofs);
void testPairList(bool do_pbc, std::ostream& ofs);

int main(int, char **) {
  {
    std::ofstream ofs("unitTest_noList");
    testNoList(ofs);
  }
  {
    std::ofstream ofs("unitTest_partial");
    testNoList_partial(ofs);
  }
  {
    std::ofstream ofs("unitTest_Pair");
    testPairList(false,ofs);
  }
  {
    std::ofstream ofs("unitTest_Pair_pbc");
    testPairList(true,ofs);
  }
}
void testNoList(std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
//only 10 atoms for not cluttering the output
  frameGenerator md(10,"sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  std::vector<AtomNumber> indexes(md.size());
  std::generate(indexes.begin(),indexes.end(),
                [i=0]() mutable {return AtomNumber().setIndex(i++);});
  auto nl=  NeighborList(indexes, serial, true, pbc, cm);
  // "step 1"
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  printNeighbors("Step 0",nl, md.size(),ofs);
}

void testNoList_partial(std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
//only 10 atoms for not cluttering the output
  frameGenerator md(10,"sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  std::vector<AtomNumber> indexes(5);
  std::generate(indexes.begin(),indexes.end(),
                [i=8]() mutable {auto x=AtomNumber::index((i)); i-=2; return x;});
  auto nl=  NeighborList(indexes, serial, true, pbc, cm);
  // "step 1"
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  printNeighbors("Step 0",nl, md.size(),ofs);
}

void testPairList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
  Communicator cm{};
  std::vector<double> box(9);
  frameGenerator md(26,"sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  double cutoff=(mybox[0][0]/5)*1.999;

  std::vector<AtomNumber> indexesA;
  std::vector<AtomNumber> indexesB;
//curating some interestin couples
//Note that 0 and 10 are intentially skipped
  for(const auto& [i,j]: {
  std::pair<unsigned,unsigned> {1,2}, {3,4},
{4,5}, {5,4}, {6,7}, {7,8},
{12,15}, {13,4}, {16,17},
{17,8}, {18,19}, {19,22},
{21,18}, {25,24},
// works only for pbcs
{2,8}, {5,3}, {8,6},
{9,15}, {11,17},
{14,12}, {15,17}, {20,2},
{22,4}, {23,21},
//this is always not connected
{24,11},
    }) {
    indexesA.push_back(AtomNumber::index(i));
    indexesB.push_back(AtomNumber::index(j));
  }
  auto nl=  NeighborList(indexesA,indexesB,serial,
                         true,//do_pair
                         do_pbc, pbc, cm, cutoff,1);
  //reordeing the atoms to respect the list of indexes passed
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  auto title = std::string("Step 0, pbc ") +((do_pbc)?"on":"off");
  printNeighbors(title,nl,md.size(),ofs);
}
