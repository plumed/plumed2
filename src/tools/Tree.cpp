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

Tree::Tree(GenericMolInfo* moldat) {
// initialize class 
   moldat_ = moldat;
// check if molinfo present
   if(!moldat_) plumed_merror("MOLINFO DATA not found");
}

std::vector<AtomNumber> Tree::buildTree(std::vector<AtomNumber> atoms)
{
  // list of AtomNumbers ordered by proximity in PDB file
  std::vector<AtomNumber> tree;
  // clear root_ vector
  root_.clear();
  // store position of first atom
  refpos_ = moldat_->getPosition(atoms[0]);
  // initialize tree and root vectors
  tree.push_back(atoms[0]);
  root_.push_back(atoms[0]);
  // remove first entry in atoms
  atoms.erase(atoms.begin());
  // loop on remaining atoms
  while(atoms.size()>0) {
   // reset minimum distance
   double mindist = std::numeric_limits<double>::max();
   unsigned iat, itr;
   // find closest pair of atoms
   for(unsigned i=0; i<atoms.size(); ++i){
    // get position in atoms
    Vector posi = moldat_->getPosition(atoms[i]);
    for(unsigned j=0; j<tree.size(); ++j){
      // get position in tree
      Vector posj = moldat_->getPosition(tree[j]);
      // calculate distance
      double dist = delta(posi,posj).modulo();
      // check if minimum distance
      if(dist<mindist){
         mindist = dist;
         iat = i;
         itr = j;
      }
    }
   }
   // now update tree, root_ and atoms vectors
   tree.push_back(atoms[iat]);
   root_.push_back(tree[itr]);
   atoms.erase(atoms.begin()+iat);
  }
  // return 
  return tree;
}

std::vector<AtomNumber> Tree::getRoot() const
{
  return root_;
}

Vector Tree::getFirstPosition() const
{
  return refpos_;
}

}
