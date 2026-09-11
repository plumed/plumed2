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

//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif //USE_MPI

#define check(arg) (((arg)) ? "pass\n" : "not pass\n")

constexpr bool serial = false;
void testSingleList(bool do_pbc, std::ostream& ofs, std::ostream& initFile, PLMD::Communicator& comm);
void testDoubleList(bool do_pbc, std::ostream& ofs, std::ostream& initFile, PLMD::Communicator& comm);
void testPairList(bool do_pbc,   std::ostream& ofs, std::ostream& initFile, PLMD::Communicator& comm);
void testNoNL(bool do_pbc,       std::ostream& ofs, std::ostream& initFile, PLMD::Communicator& comm);

void testResult(std::string name,
                const std::vector<PLMD::AtomNumber> &indexes,
                const bool do_pbc,
                const PLMD::NeighborList& nl,
                std::ostream& ofs) {
  name = "["+name+", pbc " +((do_pbc)?"on":"off")+"] atom ";
  for( unsigned i=0; i < indexes.size() ; ++i) {
    auto mynl=nl.getNeighbors(i);
    ofs <<name<< indexes[i].index() << ":";
    //sorting for test stability:
    std::transform(mynl.begin(),mynl.end(),mynl.begin(),[&](unsigned ii) {
      return indexes[ii].index();
    });
    std::sort(mynl.begin(),mynl.end());
    for (auto j: mynl) {
      ofs << " " <<j;
    }
    ofs << std::endl;
  }
}

int main(int argc, char **argv) {

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
#endif //USE_MPI
  {
    PLMD::Communicator comm;
#ifdef USE_MPI
    MPI_Comm c;
    MPI_Comm_dup(MPI_COMM_WORLD,&c);
    comm.Set_comm(&c);
#endif //USE_MPI
    std::string rank="";
    if(auto myrank = comm.Get_rank(); myrank!=0) {
      rank=std::to_string(myrank);
    }
    std::ofstream ofs("unitTest"+rank);
    std::ofstream initfile("initCheck"+rank);
    testSingleList(false,ofs,initfile,comm);
    testSingleList(true, ofs,initfile,comm);
    testDoubleList(false,ofs,initfile,comm);
    testDoubleList(true, ofs,initfile,comm);
    testPairList(false,ofs,initfile,comm);
    testPairList(true, ofs,initfile,comm);
    // an extra test to check for the NL with no NL
    std::ofstream ofsnl("testNoNL"+rank);
    testNoNL(false,ofsnl,initfile,comm);
    testNoNL(true, ofsnl,initfile,comm);
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif //USE_MPI
}

void testSingleList(const bool do_pbc, std::ostream& ofs, std::ostream& initfile, PLMD::Communicator& cm) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
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
  initfile << "Single list: ready() before update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  nl.update(atoms);
  initfile << "Single list: ready() after update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  testResult("Single list",indexes,do_pbc,nl,ofs);
}

void testDoubleList(const bool do_pbc, std::ostream& ofs, std::ostream& initfile, PLMD::Communicator& cm) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
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
  initfile << "Two lists: ready() before update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  nl.update(atoms_indexed);
  initfile << "Two lists: ready() after update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  testResult("Two lists",indexes,do_pbc,nl,ofs);
}
void testPairList(const bool do_pbc, std::ostream& ofs, std::ostream& initfile, PLMD::Communicator& cm) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
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
  initfile << "List of pairs: ready() before update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  nl.update(atoms_indexed);
  initfile << "List of pairs: ready() after update returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  testResult("List of pairs",indexes,do_pbc,nl,ofs);
}

void testNoNL(const bool do_pbc, std::ostream& ofs, std::ostream& initfile, PLMD::Communicator& cm) {
  using namespace PLMD;
  Pbc pbc{};
  Random rng;
  std::vector<double> box(9);
  // a smaller system for sanity
  std::vector<Vector> atoms(3*3*3);
  auto d = AtomDistribution::getAtomDistribution("sc");
//getting some base informations
  d->frame(atoms,box,0,rng);
  Tensor mybox{box[0], box[1], box[2], box[3], box[4], box[5], box[6], box[7], box[8]};
  pbc.setBox(mybox);
  std::vector<AtomNumber> indexes(atoms.size());
  std::generate(indexes.begin(),indexes.end(),
                [i=0]() mutable {return AtomNumber().setIndex(i++);});
  auto nl=  NeighborList(indexes,
                         serial,
                         do_pbc,
                         pbc,
                         cm);

  initfile << "NoNL: ready() (always) returns \""
           << (nl.ready() ? "true" : "false")  << "\"\n";
  if(nl.getStride()>0) {
    //you should not be here
    plumed_assert(false)<<"getstride returned > 0!!!";
    nl.update(atoms);
  }
  testResult("NoNL",indexes,do_pbc,nl,ofs);
}
