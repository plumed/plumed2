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
  // Implementation inspired from:
  // https://mayanknatani.wordpress.com/2013/05/31/euclidean-minimummaximum-spanning-tree-emst/
  //
  // list of AtomNumbers ordered by proximity in PDB file
  std::vector<AtomNumber> tree;
  // clear root_ vector
  root_.clear();
  // useful vectors
  std::vector<bool> intree(atoms.size(), false);
  std::vector<double> mindist(atoms.size(), std::numeric_limits<double>::max());
  // initialize tree with first atom
  mindist[0] = 0.0;
  // loops
  for(unsigned i=0; i<atoms.size(); ++i) {
    int selected_vertex = -1;
    for(unsigned j=0; j<atoms.size(); ++j) {
      if( !intree[j] && (selected_vertex==-1 || mindist[j] < mindist[selected_vertex]) )
        selected_vertex = j;
    }
    // add to tree
    tree.push_back(atoms[selected_vertex]);
    intree[selected_vertex] = true;
    // update distances
    double minroot = std::numeric_limits<double>::max();
    int iroot = -1;
    for(unsigned j=0; j<atoms.size(); ++j) {
      double dist = delta(moldat_->getPosition(atoms[selected_vertex]), moldat_->getPosition(atoms[j])).modulo();
      if(dist < mindist[j]) mindist[j] = dist;
      if(dist < minroot && intree[j] && dist>0.0) {
        minroot = dist;
        iroot = j;
      }
    }
    // add to root vector
    if(iroot>=0) root_.push_back(atoms[iroot]);
  }
  // return
  return tree;
}

std::vector<AtomNumber> Tree::getRoot() const
{
  return root_;
}

}
