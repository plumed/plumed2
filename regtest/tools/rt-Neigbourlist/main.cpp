#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include "plumed/tools/AtomDistribution.h"
#include "plumed/tools/Random.h"
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>


#define check(arg) (((arg)) ? "pass\n" : "not pass\n")

constexpr bool serial = true;
void testSingleList(bool do_pbc, std::ostream& ofs);
void testDoubleList(bool do_pbc, std::ostream& ofs);
void testPairList(bool do_pbc, std::ostream& ofs);

void testResult(std::string name,
                const std::vector<PLMD::AtomNumber> &indexes,
                const bool do_pbc,
                const PLMD::NeighborList& nl,
                std::ostream& ofs) {
  name = "["+name+", pbc " +((do_pbc)?"on":"off")+"] atom ";
  for( unsigned i=0; i < indexes.size() ; ++i) {
    auto mynl=nl.getNeighbors(i);
    ofs <<name<< indexes[i].index() << ":";
    for (auto j: mynl) {
      ofs << " " <<indexes[j].index();
    }
    ofs << std::endl;
  }
}

int main(int, char **) {
  std::ofstream ofs("unitTest");
  testSingleList(false,ofs);
  testSingleList(true,ofs);
  testDoubleList(false,ofs);
  testDoubleList(true,ofs);
  testPairList(false,ofs);
  testPairList(true,ofs);
}

void testSingleList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
  Communicator cm{};
  std::vector<double> box(9);
  std::vector<Vector> atoms(5*5*5);
  auto d = AtomDistribution::getAtomDistribution("sc");
//getting some base informations
  d->frame(atoms,box,0,rng);
  Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
  pbc.setBox(mybox);
  double cutoff=(box[0]/5)*1.999;
  std::vector<AtomNumber> indexes(atoms.size());
  std::generate(indexes.begin(),indexes.end(),
                [i=0]() mutable {return AtomNumber().setIndex(i++);});
  auto nl=  NeighborList(indexes,
                         serial,
                         do_pbc,
                         pbc,
                         cm,
                         cutoff,
                         1);
  nl.update(atoms);
  testResult("Single list",indexes,do_pbc,nl,ofs);
}

void testDoubleList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
  Communicator cm{};
  std::vector<double> box(9);
  std::vector<Vector> atoms(5*5*5-1);
  auto d = AtomDistribution::getAtomDistribution("sc");
//getting some base informations
  d->frame(atoms,box,0,rng);
  Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
  pbc.setBox(mybox);
  double cutoff=(box[0]/5)*1.999;

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
                         1);

  //reordeing the atoms to respect the list of indexes passed
  std::vector<Vector> atoms_indexed=atoms;
  for(unsigned i=0; i< indexesA.size(); ++i) {
    atoms_indexed[i] = atoms[indexesA[i].index()];
    atoms_indexed[i+indexesA.size()] = atoms[indexesB[i].index()];
  }
  nl.update(atoms_indexed);
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  testResult("Two lists",indexes,do_pbc,nl,ofs);
}
void testPairList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
  Communicator cm{};
  std::vector<double> box(9);
  std::vector<Vector> atoms(5*5*5-1);
  auto d = AtomDistribution::getAtomDistribution("sc");
//getting some base informations
  d->frame(atoms,box,0,rng);
  Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
  pbc.setBox(mybox);
  double cutoff=(box[0]/5)*1.999;

  std::vector<AtomNumber> indexesA(atoms.size()/2);
  std::vector<AtomNumber> indexesB(indexesA.size());
  std::generate(indexesA.begin(),indexesA.end(),
                [i=0]() mutable {auto x= AtomNumber().setIndex(i); i+=2; return x;}
               );

  std::generate(indexesB.begin(),indexesB.end(),
                [i=1]() mutable {auto x= AtomNumber().setIndex(i); i+=2; return x;}
               );
  auto nl=  NeighborList(indexesA,indexesB,
                         serial,
                         true,//do_pair
                         do_pbc,
                         pbc,
                         cm,
                         cutoff,
                         1);
  //reordeing the atoms to respect the list of indexes passed
  std::vector<Vector> atoms_indexed=atoms;
  for(unsigned i=0; i< indexesA.size(); ++i) {
    atoms_indexed[i] = atoms[indexesA[i].index()];
    atoms_indexed[i+indexesA.size()] = atoms[indexesB[i].index()];
  }
  nl.update(atoms_indexed);
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  testResult("List of pairs",indexes,do_pbc,nl,ofs);
}
