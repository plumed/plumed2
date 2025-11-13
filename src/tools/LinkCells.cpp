/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "LinkCells.h"
#include "Communicator.h"
#include "Tools.h"
#include "View.h"
#include <algorithm>
#include <functional>
#include <numeric>

namespace PLMD {

LinkCells::LinkCells( Communicator& cc ) :
  comm(cc) {
}

void LinkCells::setCutoff( const double lcut ) {
  cutoffwasset=true;
  link_cutoff=lcut;
}

double LinkCells::getCutoff() const {
  plumed_assert( cutoffwasset );
  return link_cutoff;
}

void LinkCells::setupCells(View<const Vector> pos,
                           const Pbc& pbc ) {
  auto box = pbc.getBox();
  if(box(0,0)==0.0 && box(0,1)==0.0 && box(0,2)==0.0 &&
      box(1,0)==0.0 && box(1,1)==0.0 && box(1,2)==0.0 &&
      box(2,0)==0.0 && box(2,1)==0.0 && box(2,2)==0.0) {
    // Create an orthorhombic box around the atomic positions that encompasses every atomic position if there are no pbc
    Vector minp, maxp;
    minp = maxp = pos[0];
    for(unsigned k=0; k<3; ++k) {
      for(unsigned i=1; i<pos.size(); ++i) {
        if( pos[i][k]>maxp[k] ) {
          maxp[k] = pos[i][k];
        }
        if( pos[i][k]<minp[k] ) {
          minp[k] = pos[i][k];
        }
      }
      if( link_cutoff<std::sqrt(std::numeric_limits<double>::max()) ) {
        box[k][k] = link_cutoff*( 1 + std::ceil( (maxp[k] - minp[k])/link_cutoff ) );
      } else {
        box[k][k] = maxp[k] - minp[k] + 1;
      }
      origin[k] = ( minp[k] + maxp[k] ) / 2;
    }
    nopbc=true;
    // Setup the pbc object by copying it from action if the Pbc were set
  } else {
    auto determinant = box.determinant();
    nopbc=false;
    plumed_assert(determinant > epsilon) <<"Cell lists cannot be built when passing a box with null volume. Volume is "<<determinant;
  }
  mypbc.setBox( box );
  {
// This is the reciprocal lattice
// notice that reciprocal.getRow(0) is a vector that is orthogonal to b and c
// This allows to use linked cells in non orthorhomic boxes
    Tensor reciprocal(transpose(mypbc.getInvBox()));
    ncells[0] = std::floor( 1.0/ reciprocal.getRow(0).modulo() / link_cutoff );
    if( ncells[0]==0 ) {
      ncells[0]=1;
    }
    ncells[1] = std::floor( 1.0/ reciprocal.getRow(1).modulo() / link_cutoff );
    if( ncells[1]==0 ) {
      ncells[1]=1;
    }
    ncells[2] = std::floor( 1.0/ reciprocal.getRow(2).modulo() / link_cutoff );
    if( ncells[2]==0 ) {
      ncells[2]=1;
    }
  }
  // Setup the strides
  nstride[0]=1;
  nstride[1]=ncells[0];
  nstride[2]=ncells[0]*ncells[1];
}

void LinkCells::buildCellLists( View<const Vector> pos,
                                View<const unsigned> indices,
                                const Pbc& pbc ) {
  setupCells(pos,pbc);
  resetCollection(innerCollection,pos, indices,pbc);
}

LinkCells::CellCollection LinkCells::getCollection( View<const Vector> pos,
    View<const unsigned> indices,
    const Pbc& pbc ) {
  CellCollection collection;
  resetCollection(collection,pos, indices,pbc);
  return collection;
}
void LinkCells::resetCollection(LinkCells::CellCollection &collection,
                                View<const Vector> pos,
                                View<const unsigned> indices,
                                const Pbc& pbc ) {
  plumed_assert( cutoffwasset && pos.size()==indices.size() );
  const auto nat=pos.size();
  // Resize and resets the lists
  allcells.assign( nat, 0 );
  collection.lcell_lists.resize( nat );
  // Setup the storage for link cells
  const unsigned ncellstot=ncells[0]*ncells[1]*ncells[2];
  collection.lcell_starts.resize( ncellstot );
  // Clear nlcells
  collection.lcell_tots.assign( ncellstot, 0 );

  // Find out what cell everyone is in
  const unsigned rank=comm.Get_rank();
  const unsigned size=comm.Get_size();
  const unsigned elementsPerRank = std::ceil(double(nat)/size);
  const unsigned int start= rank*elementsPerRank;
  const unsigned int end = ((start + elementsPerRank)< nat)?(start + elementsPerRank): nat;
  for(unsigned i=start; i<end; ++i) {
    allcells[i]=findCell( pos[i] );
    collection.lcell_tots[allcells[i]]++;
  }
  // And gather all this information on every node
  comm.Sum( allcells );
  comm.Sum( collection.lcell_tots );

  // Now prepare the link cell lists
  unsigned tot=0;
  for(unsigned i=0; i<collection.lcell_tots.size(); ++i) {
    collection.lcell_starts[i]=tot;
    tot+=collection.lcell_tots[i];
    collection.lcell_tots[i]=0;
  }
  plumed_assert( tot==nat ) <<"Total number of atoms found in link cells is "<<tot<<" number of atoms is "<<nat;

  // And setup the link cells properly
  for(unsigned j=0; j<nat; ++j) {
    unsigned myind = collection.lcell_starts[ allcells[j] ] + collection.lcell_tots[ allcells[j] ];
    collection.lcell_lists[ myind ] = indices[j];
    collection.lcell_tots[allcells[j]]++;
  }
}

#define LINKC_MIN(n) ((n<2)? 0 : -1)
#define LINKC_MAX(n) ((n<3)? 1 : 2)
#define LINKC_PBC(n,num) ((n<0)? num-1 : n%num )

void LinkCells::addRequiredCells( const std::array<unsigned,3>& celn,
                                  unsigned& ncells_required,
                                  std::vector<unsigned>& cells_required ) const {
  unsigned nnew_cells=0;
  for(int nx=LINKC_MIN(ncells[0]); nx<LINKC_MAX(ncells[0]); ++nx) {
    int xval = celn[0] + nx;
    xval=LINKC_PBC(xval,ncells[0])*nstride[0];
    for(int ny=LINKC_MIN(ncells[1]); ny<LINKC_MAX(ncells[1]); ++ny) {
      int yval = celn[1] + ny;
      yval=LINKC_PBC(yval,ncells[1])*nstride[1];
      for(int nz=LINKC_MIN(ncells[2]); nz<LINKC_MAX(ncells[2]); ++nz) {
        int zval = celn[2] + nz;
        zval=LINKC_PBC(zval,ncells[2])*nstride[2];

        unsigned mybox=xval+yval+zval;
        bool added=false;
        for(unsigned k=0; k<ncells_required; ++k) {
          if( mybox==cells_required[k] ) {
            added=true;
            break;
          }
        }
        if( !added ) {
          cells_required[ncells_required+nnew_cells]=mybox;
          nnew_cells++;
        }
      }
    }
  }
  ncells_required += nnew_cells;
}

void LinkCells::retrieveNeighboringAtoms( const Vector& pos,
    std::vector<unsigned>& cell_list,
    unsigned& natomsper,
    std::vector<unsigned>& atoms ) const {
  if( cell_list.size()!=getNumberOfCells() ) {
    cell_list.resize( getNumberOfCells() );
  }
  unsigned ncellt=0;
  addRequiredCells( findMyCell( pos ), ncellt, cell_list );
  retrieveAtomsInCells( ncellt, make_const_view(cell_list), natomsper, atoms );
}

void LinkCells::retrieveAtomsInCells( const unsigned ncells_required,
                                      View<const unsigned> cells_required,
                                      unsigned& natomsper,
                                      std::vector<unsigned>& atoms,
                                      const unsigned avoidIndex) const {
  for(unsigned i=0; i<ncells_required; ++i) {
    unsigned mybox=cells_required[i];
    auto boxList = innerCollection.getCellIndexes(mybox);
    if (avoidIndex!=std::numeric_limits<unsigned>::max()) {
      for(const unsigned myatom : boxList) {
        if( myatom!=avoidIndex ) { // Ideally would provide an option to not do this
          atoms[natomsper]=myatom;
          natomsper++;
        }
      }
    } else {
      for(const unsigned myatom : boxList) {
        atoms[natomsper]=myatom;
        natomsper++;
      }
    }
  }
}

std::array<unsigned,3> LinkCells::findMyCell( Vector mypos ) const {
  std::array<unsigned,3> celn;
  if( nopbc ) {
    mypos = mypos - origin;
  }
  Vector fpos=mypbc.realToScaled( mypos );
  for(unsigned j=0; j<3; ++j) {
    celn[j] = std::floor( ( Tools::pbc(fpos[j]) + 0.5 ) * ncells[j] );
    plumed_assert( celn[j]>=0 && celn[j]<ncells[j] ) <<"in link cell "<<celn[j]
        <<" but should be between 0 and "<<ncells[j]
        <<" link cell cutoff is "<<link_cutoff
        <<" position is "<<fpos[0]<<" "<<fpos[1]<<" "<<fpos[2]
        <<" box is "<<mypbc.getBox()(0,0)<<" "<<mypbc.getBox()(1,1)<<" "<<mypbc.getBox()(2,2);
  }
  return celn;
}

unsigned LinkCells::convertIndicesToIndex( const unsigned nx,
    const unsigned ny,
    const unsigned nz ) const {
  return nx*nstride[0] + ny*nstride[1] + nz*nstride[2];
}

unsigned LinkCells::findCell( const Vector& pos ) const {
  std::array<unsigned,3> celn( findMyCell(pos ) );
  return convertIndicesToIndex( celn[0], celn[1], celn[2] );
}

unsigned LinkCells::getMaxInCell() const {
  unsigned maxn = innerCollection.lcell_tots[0];
  for(unsigned i=1; i<innerCollection.lcell_tots.size(); ++i) {
    if( innerCollection.lcell_tots[i]>maxn ) {
      maxn=innerCollection.lcell_tots[i];
    }
  }
  return maxn;
}

///
/// @param nat
/// @param pos
/// @param ind
/// @param tind
/// @param neigh_pos
/// @param neigh_ind
/// @param pbc
/// @param natoms_per_list
/// @param nlist
void LinkCells::createNeighborList( unsigned nat,
                                    View<const Vector> pos,
                                    View<const unsigned> ind,
                                    View<const unsigned> tind,
                                    View<const Vector> neigh_pos,
                                    View<const unsigned> neigh_ind,
                                    const Pbc& pbc,
                                    unsigned& natoms_per_list,
                                    std::vector<std::size_t>& nlist ) {
  buildCellLists( neigh_pos, neigh_ind, pbc );
//  natoms_per_list = 27*getMaxInCell();
//this should save a little memory
  natoms_per_list=innerCollection.getMaximimumCombination(27);
  if( natoms_per_list>innerCollection.lcell_lists.size() ) {
    natoms_per_list = innerCollection.lcell_lists.size();
  }
  const unsigned nlist_sz = nat*( 2 + natoms_per_list );
  nlist.resize( nlist_sz );
  std::vector<unsigned> indices( 1+natoms_per_list );
  std::vector<unsigned> cells_required( getNumberOfCells() );
  for(unsigned i=0; i<pos.size(); ++i) {
    unsigned ncells_required=0;
    addRequiredCells( findMyCell( pos[i] ), ncells_required, cells_required );
    unsigned natoms=1;
    indices[0] = ind[i];
    retrieveAtomsInCells( ncells_required,
                          make_const_view(cells_required),
                          natoms,
                          indices,
                          ind[i]);
    nlist[tind[i]] = natoms;
    const std::size_t lstart = nat + tind[i]*(1+natoms_per_list);
    std::copy(indices.begin(),
              indices.begin()+natoms,
              nlist.begin()+lstart);
  }
}

unsigned LinkCells::CellCollection::getMaximimumCombination(const unsigned ncells) const {
  //this is not efficient (there is a copy), but in principle it should be called not much times:
  //nth_element order by partition until the array is partially sorted
  //(it stops when the pivot is at the asked point, so all the elements to the right of the nth element are higher/lower than it)
  if (ncells >= lcell_tots.size()) {
    return std::accumulate(lcell_tots.begin(),lcell_tots.end(),0);
  }
  auto tmpCopy=lcell_tots;
  auto nth_place=tmpCopy.begin() + ncells-1;
  std::nth_element(tmpCopy.begin(),nth_place,tmpCopy.end(),std::greater<>());

  return std::accumulate(tmpCopy.begin(),nth_place,0);
}

} //namespace PLMD
