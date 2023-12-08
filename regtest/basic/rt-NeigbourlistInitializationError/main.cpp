#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/Pbc.h"
#include <fstream>
#include <iostream>

using PLMD::AtomNumber;
using PLMD::Communicator;
using PLMD::NeighborList;
using PLMD::Pbc;

// Testing that the Neigbour list must throw an explanatory PLMD::Exception
// when asked to allocate too many pairs


#define check(arg) (((arg)) ? "pass\n" : "not pass\n")

int main(int, char **) {
  std::ofstream report("unitTest");
  Pbc pbc{};
  pbc.setBox(PLMD::Tensor({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}));
  Communicator cm{};
  bool serial = true;
  bool do_pbc = false;
  // nat0 times nat1 shold be ludicrous big (nat0*nat1*8B~=30~40GB)
  const size_t nat0=100000;
  const size_t nat1 = 90000;
  std::vector<AtomNumber> list0(nat0);
  size_t i = 0;
  for (auto &an : list0) {
    an.setIndex(i);
    ++i;
  }
  {
    report << "Single list:\n";
    std::string prepend="["+std::to_string(nat0)+"]";
    try{
    size_t expected = ((nat0 - 1) * nat0) / 2;
    auto nl = NeighborList(list0, serial, do_pbc, pbc, cm);
    //I need this line to ensure that nl is not optimized away
    report << prepend << "Initial number:      "
            << check(nl.size() == expected);
    } catch ( PLMD::Exception & error ){
      report << prepend <<"Exception text: "
              << error.what();
    }
    report << "\n";
  }
  //doing the same thing with two lists
  
  std::vector<AtomNumber> list1(nat1);
  i = 0;
  for (auto &an : list1) {
    an.setIndex(i);
    ++i;
  }

  {
    report << "Double list, no pairs:\n";
    std::string prepend="["+std::to_string(nat0)
        + ", " + std::to_string(nat1)  +"]";
    bool do_pair = false;
    size_t expected = nat1 * nat0;
    try{
      auto nl = NeighborList(list0, list1, serial, do_pair, do_pbc, pbc, cm);
      report << prepend << "Initial number:      "
              << check(nl.size() == expected);
      
    } catch( PLMD::Exception & error) {
      report << prepend <<"Exception text: "
              << error.what();
    }
    report << "\n";   
  }
}
  
