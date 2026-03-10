#include "plumed/tools/AtomNumber.h"
#include "plumed/tools/NeighborList.h"

#include "nlTools.h"

#include <iostream>
#include <algorithm>
#include <vector>


namespace PLMD {
namespace test {

void printNeighbors( std::string prefix,
                     const PLMD::NeighborList& nl,
                     unsigned const nat,
                     std::ostream& ofs) {
  // Apparently getFullAtomList not being const method is a old thing, this workarond guarantees that this test is backportable :D
  const std::vector<PLMD::AtomNumber> &indexes = const_cast<PLMD::NeighborList&>(nl).getFullAtomList();
  prefix = "["+prefix+"] atom ";
  for (unsigned at =0; at<nat; ++at) {
    auto idx = PLMD::AtomNumber::index(at);
    std::vector<std::size_t> idxs;
    idxs.reserve(10);
    for (auto it = std::find(indexes.begin(), indexes.end(), idx);
         it != indexes.end();
         it = std::find(std::next(it), indexes.end(), idx)) {
      idxs.push_back(static_cast<std::size_t>(std::distance(indexes.begin(), it)));
    }

    if (idxs.size()>0) {
      std::vector<unsigned> mynl;
      for( auto i:idxs) {
        auto mynl_=nl.getNeighbors(i);
        //converting from indexes relative to the NL to real atom indexes
        std::transform(mynl_.begin(),mynl_.end(),mynl_.begin(),[&](unsigned ii) {
          return indexes[ii].index();
        });
        mynl.insert(mynl.end(),mynl_.begin(),mynl_.end());
      }
      if (mynl.size() >0) {
        //printing the list of neighbors, this list will show duplicates
        ofs <<prefix<< idx.index() << ":";
        //sorting for test stability:
        std::sort(mynl.begin(),mynl.end());
        for (auto j: mynl) {
          ofs << " " <<j;
        }
        ofs << std::endl;
      }
    }
  }
}
} // namespace test
} // namespace PLMD
