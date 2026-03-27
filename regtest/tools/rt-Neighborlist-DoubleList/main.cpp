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
void testDoubleList(bool do_pbc, std::ostream& ofs);
void testRemapping(const bool do_pbc, std::ostream& ofs);
int main(int, char **) {
  {
    std::ofstream ofs("unitTest_DL");
    testDoubleList(false,ofs);
  }
  {
    std::ofstream ofs("unitTest_DL_pbc");
    testDoubleList(true,ofs);
  }
  {
    std::ofstream ofs("unitTest_DL_partial");
    testRemapping(false, ofs);
  }
}

void testDoubleList(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
  frameGenerator md(5*5*5-1,"sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  double cutoff=(mybox[0][0]/5)*1.999;

  std::vector<AtomNumber> indexesA(md.size()/2);
  std::vector<AtomNumber> indexesB(md.size()/2);
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
                         10);

  //reordeing the atoms to respect the list of indexes passed
  //let's try to imitate the behaviour of requestAtoms:
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  auto indexes=indexesA;
  indexes.insert(indexes.end(),indexesB.begin(),indexesB.end());
  auto title = std::string("Two lists, pbc ") +((do_pbc)?"on":"off");
  printNeighbors(title,nl,md.size(),ofs);
}

void testRemapping(const bool do_pbc, std::ostream& ofs) {
  using namespace PLMD;
  Pbc pbc{};
  Communicator cm{};
  frameGenerator md(26,"sc");
//getting some base informations
  Tensor mybox=md.getBox();
  pbc.setBox(mybox);
  double cutoff=(mybox[0][0]/5)*1.999;
  std::vector<AtomNumber> indexes(md.size()/2);
  std::generate(indexes.begin(),indexes.end(),
                [i=0]() mutable {auto x= AtomNumber().setIndex(i); i+=2; return x;});
  std::vector<AtomNumber> indexesB;
  for(auto i: {
        1,2,3,5,8,23,25
      }) {
    indexesB.push_back(AtomNumber::index(i));
  }
  auto nl=  NeighborList(indexes,indexesB,
                         serial,
                         false,//do_pair
                         do_pbc, pbc, cm, cutoff, 10);
  // "step 1"
  std::vector<Vector> atoms_indexed = md.requestAtoms(nl.getFullAtomList());
  nl.update(atoms_indexed);
  auto title = std::string("Reference, pbc ") +((do_pbc)?"on":"off");
  printNeighbors(title,nl,md.size(),ofs);
//building the reference:
  std::map<unsigned,std::vector<unsigned>> ref;
  auto fal = nl.getFullAtomList();
  for (unsigned at =0; at<md.size(); ++at) {
    auto idx = PLMD::AtomNumber::index(at);
    std::vector<std::size_t> idxs;
    idxs.reserve(4);
    for (auto it = std::find(fal.begin(), fal.end(), idx);
         it != fal.end();
         it = std::find(std::next(it), fal.end(), idx)) {
      idxs.push_back(static_cast<std::size_t>(std::distance(fal.begin(), it)));
    }

    if (idxs.size()>0) {
      std::vector<unsigned> mynl;
      for( auto i:idxs) {
        auto mynl_=nl.getNeighbors(i);
        //converting from fal relative to the NL to real atom indexes
        std::transform(mynl_.begin(),mynl_.end(),mynl_.begin(),[&](unsigned ii) {
          return fal[ii].index();
        });
        mynl.insert(mynl.end(),mynl_.begin(),mynl_.end());
      }
      if (mynl.size() >0) {
        std::sort(mynl.begin(),mynl.end());
        ref[at]=mynl;
      }
    }
  }
//triggering the remapping:
  auto reducedList=nl.getReducedAtomList();
//checking that the list is not changed due to the remapping
  for (unsigned at =0; at<md.size(); ++at) {
    auto idx = PLMD::AtomNumber::index(at);
    std::vector<std::size_t> idxs;
    idxs.reserve(4);
    for (auto it = std::find(reducedList.begin(), reducedList.end(), idx);
         it != reducedList.end();
         it = std::find(std::next(it), reducedList.end(), idx)) {
      idxs.push_back(static_cast<std::size_t>(std::distance(reducedList.begin(), it)));
    }

    if (idxs.size()>0) {
      std::vector<unsigned> mynl;
      for( auto i:idxs) {
        auto mynl_=nl.getNeighbors(i);
        //converting from reducedList relative to the NL to real atom indexes
        std::transform(mynl_.begin(),mynl_.end(),mynl_.begin(),[&](unsigned ii) {
          return reducedList[ii].index();
        });
        mynl.insert(mynl.end(),mynl_.begin(),mynl_.end());
      }
      std::sort(mynl.begin(),mynl.end());
      if (ref.count(at) == 0 ) {
        if (mynl.size() >0) {
          ofs <<"Atom " << at << ": reducing the NL list added some neigbors to atom "<<at<< "\n";
        }
      } else {
        if (mynl.size() == 0) {
          ofs <<"Atom " << at << ": reducing the NL list removed all the neigbors from atom "<<at<< "\n";
        } else {
          if(mynl!=ref.at(at)) {
            ofs <<"Atom " << at << ": reducing the NL list changed the neigbors of atom "<<at<< "\n";
            /* this can be useful for debugging locally, but will clutter the CI
            ofs << "From : ";
            std::copy(ref.at(at).begin(), ref.at(at).end(), std::ostream_iterator<unsigned>(ofs, " "));
            ofs<<"\n";
            ofs <<"To : ";
            std::copy(mynl.begin(), mynl.end(), std::ostream_iterator<unsigned>(ofs, " "));
            ofs<<"\n";
            */
          }
        }
      }

    }
  }
}
