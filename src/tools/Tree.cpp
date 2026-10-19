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
#include "Pbc.h"
#include "Tensor.h"
#include "core/GenericMolInfo.h"
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>

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
  // longest edge of the tree (squared)
  double longest2=0.0;
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
    if(i>0) {
      longest2=std::max(longest2,mindist[selected_vertex]);
    }
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

  checkReferenceIsWhole(atoms, newatoms, positions, longest2);

  // now re-add atoms not present in the PDB
  for(unsigned i=0; i<addtotree.size(); ++i) {
    tree_.push_back(atoms[addtotree[i]]);
    tree_indexes_.push_back(addtotree[i]);
    root_.push_back(atoms[addtoroot[i]]);
    root_indexes_.push_back(addtoroot[i]);
  }
}

namespace {
// True if some atom has a periodic image of an atom (a nonzero lattice translation away) within
// sqrt(cut2) of it. When there is none, minimum-image and plain distances agree on every pair up to
// the longest edge of the plain tree, so the periodic tree is the plain tree and the reference is whole.
// A grid search, linear in the number of atoms. It answers true when the box is too small for that
// search to be cheap, so that the full check decides.
bool hasCloseImage(const std::vector<Vector> & positions, const Tensor & box, double cut2) {
  const double cut=std::sqrt(cut2);
  Vector lo(positions[0]), hi(positions[0]);
  for(const auto & p : positions) {
    for(unsigned k=0; k<3; ++k) {
      lo[k]=std::min(lo[k],p[k]);
      hi[k]=std::max(hi[k],p[k]);
    }
  }
  // lattice vectors that can bring an atom within cut of the structure
  const double reach=(hi-lo).modulo()+cut;
  const Tensor inv=inverse(box);
  int range[3];
  for(unsigned k=0; k<3; ++k) {
    const double col=std::sqrt(inv[0][k]*inv[0][k]+inv[1][k]*inv[1][k]+inv[2][k]*inv[2][k]);
    range[k]=static_cast<int>(std::ceil(reach*col));
  }
  if((2.0*range[0]+1)*(2.0*range[1]+1)*(2.0*range[2]+1)>4096) {
    return true;
  }
  std::vector<Vector> shifts;
  for(int i=-range[0]; i<=range[0]; ++i) {
    for(int j=-range[1]; j<=range[1]; ++j) {
      for(int k=-range[2]; k<=range[2]; ++k) {
        if(i==0 && j==0 && k==0) {
          continue;
        }
        const Vector n=matmul(Vector(i,j,k),box);
        if(n.modulo()<=reach) {
          shifts.push_back(n);
        }
      }
    }
  }
  // atoms in cells of side cut
  const double side=std::max(cut,1e-6*(hi-lo).modulo());
  long nc[3];
  for(unsigned k=0; k<3; ++k) {
    nc[k]=static_cast<long>((hi[k]-lo[k])/side)+1;
  }
  if(static_cast<double>(nc[0])*nc[1]*nc[2]>1e8) {
    return true;
  }
  auto cellOf=[&](const Vector & p,unsigned k) {
    return static_cast<long>(std::floor((p[k]-lo[k])/side));
  };
  std::vector<long> head(nc[0]*nc[1]*nc[2],-1), next(positions.size(),-1);
  for(unsigned a=0; a<positions.size(); ++a) {
    const long c=(cellOf(positions[a],0)*nc[1]+cellOf(positions[a],1))*nc[2]+cellOf(positions[a],2);
    next[a]=head[c];
    head[c]=a;
  }
  for(const auto & n : shifts) {
    for(const auto & x : positions) {
      const Vector p=x+n;
      if(p[0]<lo[0]-cut || p[0]>hi[0]+cut || p[1]<lo[1]-cut || p[1]>hi[1]+cut || p[2]<lo[2]-cut || p[2]>hi[2]+cut) {
        continue;
      }
      long c[3];
      for(unsigned k=0; k<3; ++k) {
        c[k]=cellOf(p,k);
      }
      for(long i=std::max(c[0]-1,0L); i<=std::min(c[0]+1,nc[0]-1); ++i) {
        for(long j=std::max(c[1]-1,0L); j<=std::min(c[1]+1,nc[1]-1); ++j) {
          for(long k=std::max(c[2]-1,0L); k<=std::min(c[2]+1,nc[2]-1); ++k) {
            for(long b=head[(i*nc[1]+j)*nc[2]+k]; b>=0; b=next[b]) {
              if(delta(p,positions[b]).modulo2()<=cut2) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}
}

void Tree::checkReferenceIsWhole(const std::vector<AtomNumber> & atoms,
                                 const std::vector<unsigned> & newatoms,
                                 const std::vector<Vector> & positions,
                                 double longest2) const {
  // The tree is built on plain distances in the reference, so WHOLE is taken on trust:
  // if part of the reference is wrapped across the periodic boundary, its link to the rest is
  // a long edge that is already its own minimum image, and the walk then moves the same atoms
  // across the boundary in every frame.
  // Check it with the CRYST1 box: build a second spanning tree on minimum-image distances and
  // walk it over the reference. A whole reference comes back unchanged; a wrapped piece comes
  // back moved by a lattice vector.
  const Tensor box=moldat_->getReferenceBox();
  if(box.determinant()==0.0 || positions.size()<2) {
    return;
  }
  // nearly always, no atom comes close to a periodic image of another: then the reference is whole
  if(!hasCloseImage(positions, box, longest2*(1.0+1e-6))) {
    return;
  }
  Pbc pbc;
  pbc.setBox(box);
  const unsigned n=positions.size();
  std::vector<bool> intree(n, false);
  std::vector<double> mindist(n, std::numeric_limits<double>::max());
  std::vector<unsigned> parent(n, 0);
  std::vector<Vector> rebuilt(positions);
  mindist[0]=0.0;
  for(unsigned i=0; i<n; ++i) {
    unsigned s=n;
    for(unsigned j=0; j<n; ++j) {
      if(!intree[j] && (s==n || mindist[j]<mindist[s])) {
        s=j;
      }
    }
    intree[s]=true;
    if(i>0) {
      rebuilt[s]=rebuilt[parent[s]]+pbc.distance(positions[parent[s]],positions[s]);
    }
    for(unsigned j=0; j<n; ++j) {
      if(!intree[j]) {
        const double d=pbc.distance(positions[s],positions[j]).modulo2();
        if(d<mindist[j]) {
          mindist[j]=d;
          parent[j]=s;
        }
      }
    }
  }
  // a wrapped atom moves by a lattice vector; rounding along the walk is many orders smaller
  double shortest=std::numeric_limits<double>::max();
  for(unsigned k=0; k<3; ++k) {
    shortest=std::min(shortest,box.getRow(k).modulo());
  }
  std::vector<unsigned> moved;
  for(unsigned j=0; j<n; ++j) {
    if((rebuilt[j]-positions[j]).modulo()>1e-3*shortest) {
      moved.push_back(j);
    }
  }
  if(moved.size()>0) {
    std::string list;
    for(unsigned k=0; k<moved.size() && k<10; ++k) {
      const AtomNumber a=atoms[newatoms[moved[k]]];
      list+=" "+std::to_string(a.serial())+" ("+moldat_->getResidueName(a)+std::to_string(moldat_->getResidueNumber(a))+")";
    }
    if(moved.size()>10) {
      list+=" and "+std::to_string(moved.size()-10)+" more";
    }
    // A warning, not an error: a group of several molecules that are not in contact (a ligand away
    // from the protein, a membrane) is reported too, and WHOLEMOLECULES can still handle it.
    std::ostringstream msg;
    msg<<"the reference structure in the MOLINFO pdb file may not be whole: "<<moved.size()
       <<" of its "<<n<<" atoms here are closer to a periodic image of the others,"
       <<" in the box of its CRYST1 record, than to the others themselves. Atoms:"<<list
       <<". If these atoms are one molecule, WHOLEMOLECULES will split it in every frame:"
       <<" the pdb was then written with molecules split by the periodic boundary"
       <<" (make it whole, for example with gmx trjconv -pbc whole),"
       <<" or its CRYST1 record is not the box the structure was written in."
       <<" This is also reported when the group holds several molecules that are not in contact,"
       <<" such as a ligand away from the protein or a membrane, or atoms too sparse for that box";
    moldat_->warning(msg.str());
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
