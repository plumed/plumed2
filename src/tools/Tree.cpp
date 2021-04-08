/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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

#include "Tree.h"
#include "Tools.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "OpenMP.h"
#include "core/GenericMolInfo.h"
#include <vector>
#include <limits>

namespace PLMD {

Tree::Tree(GenericMolInfo* moldat, bool nopbc) {
// initialize class
  moldat_ = moldat;
// check if molinfo present
  if(!moldat_) plumed_merror("MOLINFO DATA not found");
// pbc
  nopbc_ = nopbc;
}

std::vector<AtomNumber> Tree::getTree(std::vector<AtomNumber> atoms)
{
  // OpenMP stuff
  unsigned nth = OpenMP::getNumThreads();
  // list of AtomNumbers ordered by proximity in PDB file
  std::vector<AtomNumber> tree;
  // clear root_ vector
  root_.clear();
  // initialize tree
  tree.push_back(atoms[0]);
  // remove first entry in atoms
  atoms.erase(atoms.begin());
  // loop on remaining atoms
  while(atoms.size()>0) {
    // TODO
    // This can be easily parallelized with OpenMP
    // It is too slow, we calculate the same distances multiple times
    // Why not precalculating in parallel the distance matrix beforehand?
    // Too much memory?
    // reset minimum distance
    std::vector<double> mindist(nth, std::numeric_limits<double>::max());
    std::vector<unsigned> iat(nth), itr(nth);
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    // find closest pair of atoms
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get position in atoms
      Vector posi = moldat_->getPosition(atoms[i]);
      for(unsigned j=0; j<tree.size(); ++j) {
        // get position in tree
        Vector posj = moldat_->getPosition(tree[j]);
        // calculate distance
        double dist = delta(posi,posj).modulo();
        // check if minimum distance
        if(dist<mindist[OpenMP::getThreadNum()]) {
          mindist[OpenMP::getThreadNum()] = dist;
          iat[OpenMP::getThreadNum()] = i;
          itr[OpenMP::getThreadNum()] = j;
        }
      }
    }
    // index of minimum distance across threads
    unsigned imind = std::distance(mindist.begin(), std::min_element(mindist.begin(), mindist.end()));
    // now update tree, root_ and atoms vectors
    tree.push_back(atoms[iat[imind]]);
    root_.push_back(tree[itr[imind]]);
    atoms.erase(atoms.begin()+iat[imind]);
  }
  // return
  return tree;
}

std::vector<AtomNumber> Tree::getRoot() const
{
  return root_;
}

}
