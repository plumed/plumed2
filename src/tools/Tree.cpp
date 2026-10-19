/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2021-2023 The plumed team
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

Tree::Tree(GenericMolInfo* moldat) {
// initialize class
  moldat_ = moldat;
// check if molinfo present
  if(!moldat_) {
    plumed_merror("MOLINFO DATA not found");
  }
// check if reference structure is whole
  if(!moldat_->isWhole()) {
    plumed_merror("Check that reference structure in PDB file is not broken by pbc and set WHOLE in MOLINFO line");
  }
}

void Tree::buildTree(const std::vector<AtomNumber> & atoms) {
  // Implementation inspired from:
  // https://mayanknatani.wordpress.com/2013/05/31/euclidean-minimummaximum-spanning-tree-emst/
  //
  // clear root_ vector
  root_.clear();
  tree_.clear();
  root_indexes_.clear();
  tree_indexes_.clear();

  std::vector<Vector> positions;

  // remove atoms not in PDB file
  std::vector<unsigned> addtotree, addtoroot;
  std::vector<unsigned> newatoms;
  newatoms.reserve(atoms.size());
  positions.reserve(atoms.size());

  tree_.reserve(atoms.size());
  root_.reserve(atoms.size());
  tree_indexes_.reserve(atoms.size());
  root_indexes_.reserve(atoms.size());
  if(!moldat_->checkForAtom(atoms[0])) {
    plumed_merror("The first atom in the list should be present in the PDB file");
  }
  // store first atom
  newatoms.push_back(0);
  positions.push_back(moldat_->getPosition(atoms[0]));
  for(unsigned i=1; i<atoms.size(); ++i) {
    if(!moldat_->checkForAtom(atoms[i])) {
      // store this atom for later
      addtotree.push_back(i);
      // along with its root (the previous atom)
      addtoroot.push_back(i-1);
    } else {
      newatoms.push_back(i);
      positions.push_back(moldat_->getPosition(atoms[i]));
    }
  }

  // start EMST
  std::vector<bool> intree(newatoms.size(), false);
  std::vector<double> mindist(newatoms.size(), std::numeric_limits<double>::max());
  // initialize tree with first atom
  mindist[0] = 0.0;
  // loops
  for(unsigned i=0; i<newatoms.size(); ++i) {
    int selected_vertex = -1;
    for(unsigned j=0; j<newatoms.size(); ++j) {
      if( !intree[j] && (selected_vertex==-1 || mindist[j] < mindist[selected_vertex]) ) {
        selected_vertex = j;
      }
    }
    // add to tree
    plumed_assert(selected_vertex>=0);
    tree_indexes_.push_back(newatoms[selected_vertex]);
    tree_.push_back(atoms[newatoms[selected_vertex]]);
    intree[selected_vertex] = true;
    // update distances
    double minroot = std::numeric_limits<double>::max();
    int iroot = -1;
    for(unsigned j=0; j<newatoms.size(); ++j) {
      double dist = delta(positions[selected_vertex], positions[j]).modulo2();
      if(dist < mindist[j]) {
        mindist[j] = dist;
      }
      if(dist < minroot && intree[j] && dist>0.0) {
        minroot = dist;
        iroot = j;
      }
    }
    // add to root vector
    if(iroot>=0) {
      root_.push_back(atoms[newatoms[iroot]]);
      root_indexes_.push_back(newatoms[iroot]);
    }
  }

  // now re-add atoms not present in the PDB
  for(unsigned i=0; i<addtotree.size(); ++i) {
    tree_.push_back(atoms[addtotree[i]]);
    tree_indexes_.push_back(addtotree[i]);
    root_.push_back(atoms[addtoroot[i]]);
    root_indexes_.push_back(addtoroot[i]);
  }
}

const std::vector<AtomNumber> & Tree::getTree(const std::vector<AtomNumber> & atoms) {
  buildTree(atoms);
  return tree_;
}

const std::vector<AtomNumber> & Tree::getTree() const noexcept {
  return tree_;
}

const std::vector<AtomNumber> & Tree::getRoot() const noexcept {
  return root_;
}

const std::vector<unsigned> & Tree::getTreeIndexes() const noexcept {
  return tree_indexes_;
}

const std::vector<unsigned> & Tree::getRootIndexes() const noexcept {
  return root_indexes_;
}


}
