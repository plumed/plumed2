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

// Testing that the Neigbour list will be intialized with the desired number of
// couples
// We are initializing with distance and stride not set to check the default
// parameters

#define check(arg) (((arg)) ? "pass\n" : "not pass\n")

int main(int, char **) {
  std::ofstream report("unitTest");
  Pbc pbc{};
  pbc.setBox(PLMD::Tensor({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}));
  Communicator cm{};
  bool serial = true;
  bool do_pbc = false;
  for (const size_t nat0 : {100, 500, 1000, 10000}) {
    std::vector<AtomNumber> list0(nat0);
    size_t i = 0;
    for (auto &an : list0) {
      an.setIndex(i);
      ++i;
    }
    {
      report << "Single list:\n";
      std::string prepend="["+std::to_string(nat0)+"]";
      size_t expected = ((nat0 - 1) * nat0) / 2;
      auto nl = NeighborList(list0, serial, do_pbc, pbc, cm);
            
      bool expectedcouples = true;
      {
        size_t cID = 0;
        for (size_t i0 = 0; i0 < nat0 && expectedcouples; ++i0) {
          for (size_t i1 = i0+1; i1 < nat0 && expectedcouples; ++i1) {  
            auto couple = nl.getClosePair(cID);
            expectedcouples &= couple.first == i0;
            expectedcouples &= couple.second == i1;
            ++cID;
          }
        }
      }
      report << prepend << "Initial number:      "
        << check(nl.size() == expected);
      report << prepend << "getIndexPair():      "
        << check(expectedcouples);
      report << prepend << "Lastupdate is 0:     "
        << check(nl.getLastUpdate() == 0);
      report << prepend << "Default stride is 0: "
          << check(nl.getStride() == 0);
      report << "\n";
    }
    for (const size_t nat1 : {100, 500, 1000, 10000}) {

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
        auto nl = NeighborList(list0, list1, serial, do_pair, do_pbc, pbc, cm);
                
        bool expectedcouples = true;
        {
          size_t cID = 0;
          for (size_t i0 = 0; i0 < nat0 && expectedcouples; ++i0) {
            for (size_t i1 = 0; i1 < nat1 && expectedcouples; ++i1) {
              auto couple = nl.getClosePair(cID);
              //The getIndexPair for non couple input must return this be this
              //(cID / nat1);
              expectedcouples &= couple.first == i0;
              //(cID % nat1 + nat0);
              expectedcouples &= couple.second == nat0+i1;
              ++cID;
            }
          }
        }
        report << prepend << "Initial number:      "
          << check(nl.size() == expected);
        report << prepend << "getIndexPair():      "
          << check(expectedcouples);
        report << prepend << "Lastupdate is 0:     "
          << check(nl.getLastUpdate() == 0);
        report << prepend << "Default stride is 0: "
          << check(nl.getStride() == 0);
        report << "\n";
      }

      if (nat1 == nat0) {
        report << "Double list, with pairs:\n";
        std::string prepend="["+std::to_string(nat0)
            + ", " + std::to_string(nat1) +"]";
        bool do_pair = true;
        size_t expected = nat0;
        auto nl = NeighborList(list0, list1, serial, do_pair, do_pbc, pbc, cm);
        
        bool expectedcouples = true;
        for (size_t cID = 0; cID < nat0 && expectedcouples; ++cID) {
          auto couple = nl.getClosePair(cID);
          expectedcouples &= couple.first == cID;
          expectedcouples &= couple.second == cID + nat0;
        }
        report << prepend << "Initial number:      " 
          << check(nl.size() == expected);
        report << prepend << "getIndexPair():      "
          << check(expectedcouples);
        report << prepend <<  "Lastupdate is 0:     " 
          << check(nl.getLastUpdate() == 0);
        report << prepend << "Default stride is 0: "
          << check(nl.getStride() == 0);
        report << "\n";
      }
    }
  }
}