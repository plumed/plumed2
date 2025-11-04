/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_NeighborList_h //{
#define __PLUMED_tools_NeighborList_h

#include "Vector.h"
#include "AtomNumber.h"

#include <vector>

namespace PLMD {

class Pbc;
class Communicator;

/// \ingroup TOOLBOX
/// A class that implements neighbor lists from two lists or a single list of atoms
class NeighborList {
public:
  using pairIDs=std::pair<unsigned,unsigned>;
  using pairAtomNumbers=std::pair<PLMD::AtomNumber,PLMD::AtomNumber>;
private:
  enum class NNStyle {Pair,TwoList,SingleList};
  bool reduced=false;
  bool serial_;
  bool do_pbc_;
  NNStyle style_;
  const PLMD::Pbc* pbc_;
  Communicator& comm;
  std::vector<PLMD::AtomNumber> fullatomlist_{};
  std::vector<PLMD::AtomNumber> requestlist_{};
  std::vector<pairIDs > neighbors_{};
  double distance_;
  size_t nlist0_=0;
  size_t nlist1_=0;
  size_t nallpairs_;
  unsigned stride_=0;
  unsigned lastupdate_=0;
/// Initialize the neighbor list with all possible pairs
  void initialize();
/// Return the pair of indexes in the positions array
/// of the two atoms forming the i-th pair among all possible pairs
  pairIDs getIndexPair(unsigned i);
/// Extract the list of atoms from the current list of close pairs
  void setRequestList();
public:
  NeighborList(const std::vector<PLMD::AtomNumber>& list0,
               const std::vector<PLMD::AtomNumber>& list1,
               bool serial,
               bool do_pair,
               bool do_pbc,
               const PLMD::Pbc& pbc,
               Communicator &cm,
               double distance=1.0e+30,
               unsigned stride=0);
  NeighborList(const std::vector<PLMD::AtomNumber>& list0,
               bool serial,
               bool do_pbc,
               const PLMD::Pbc& pbc,
               Communicator &cm,
               double distance=1.0e+30,
               unsigned stride=0);
  ~NeighborList();
/// Return the list of all atoms. These are needed to rebuild the neighbor list.
  std::vector<PLMD::AtomNumber>& getFullAtomList();
/// Update the indexes in the neighbor list to match the
/// ordering in the new positions array
/// and return the new list of atoms that must be requested to the main code
  std::vector<PLMD::AtomNumber>& getReducedAtomList();
/// Update the neighbor list and prepare the new
/// list of atoms that will be requested to the main code
  void update(const std::vector<PLMD::Vector>& positions);
/// Get the update stride of the neighbor list
  unsigned getStride() const;
/// Get the last step in which the neighbor list was updated
  unsigned getLastUpdate() const;
/// Set the step of the last update
  void setLastUpdate(unsigned step);
/// Get the size of the neighbor list
  unsigned size() const;
/// Get the distance used to create the neighbor list
  double distance() const;
/// Get the i-th pair of the neighbor list
  pairIDs getClosePair(unsigned i) const;
/// Get the list of neighbors of the i-th atom
  std::vector<unsigned> getNeighbors(unsigned i);
/// Get the i-th pair of AtomNumbers from the neighbor list
  pairAtomNumbers getClosePairAtomNumber(unsigned i) const;
};

} // namespace PLMD

#endif //__PLUMED_tools_NeighborList_h }
