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

void LinkCells::buildCellLists( View<const Vector> pos,
                                View<const unsigned> indices,
                                const Pbc& pbc ) {
  plumed_assert( cutoffwasset && pos.size()==indices.size() );

  // Create an orthorhombic box around the atomic positions that encompasses every atomic position if there are no pbc
  auto box = pbc.getBox();
  if(box(0,0)==0.0 && box(0,1)==0.0 && box(0,2)==0.0 &&
      box(1,0)==0.0 && box(1,1)==0.0 && box(1,2)==0.0 &&
      box(2,0)==0.0 && box(2,1)==0.0 && box(2,2)==0.0) {
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

  // Setup the lists
  if( pos.size()!=allcells.size() ) {
    allcells.resize( pos.size() );
    lcell_lists.resize( pos.size() );
  }

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

  // Setup the storage for link cells
  unsigned ncellstot=ncells[0]*ncells[1]*ncells[2];
  if( lcell_tots.size()!=ncellstot ) {
    lcell_tots.resize( ncellstot );
    lcell_starts.resize( ncellstot );
  }
  // Clear nlcells
  lcell_tots.assign( lcell_tots.size(), 0 );
  // Clear allcells
  allcells.assign( allcells.size(), 0 );

  // Find out what cell everyone is in
  unsigned rank=comm.Get_rank(), size=comm.Get_size();
  for(unsigned i=rank; i<pos.size(); i+=size) {
    allcells[i]=findCell( pos[i] );
    lcell_tots[allcells[i]]++;
  }
  // And gather all this information on every node
  comm.Sum( allcells );
  comm.Sum( lcell_tots );

  // Now prepare the link cell lists
  unsigned tot=0;
  for(unsigned i=0; i<lcell_tots.size(); ++i) {
    lcell_starts[i]=tot;
    tot+=lcell_tots[i];
    lcell_tots[i]=0;
  }
  plumed_assert( tot==pos.size() ) <<"Total number of atoms found in link cells is "<<tot<<" number of atoms is "<<pos.size();

  // And setup the link cells properly
  for(unsigned j=0; j<pos.size(); ++j) {
    unsigned myind = lcell_starts[ allcells[j] ] + lcell_tots[ allcells[j] ];
    lcell_lists[ myind ] = indices[j];
    lcell_tots[allcells[j]]++;
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

void LinkCells::retrieveNeighboringAtoms( const Vector& pos, std::vector<unsigned>& cell_list,
    unsigned& natomsper, std::vector<unsigned>& atoms ) const {
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
                                      std::vector<unsigned>& atoms ) const {
  for(unsigned i=0; i<ncells_required; ++i) {
    unsigned mybox=cells_required[i];
    for(unsigned k=0; k<lcell_tots[mybox]; ++k) {
      unsigned myatom = lcell_lists[lcell_starts[mybox]+k];
      if( myatom!=atoms[0] ) { // Ideally would provide an option to not do this
        atoms[natomsper]=myatom;
        natomsper++;
      }
    }
  }
}

std::array<unsigned,3> LinkCells::findMyCell( const Vector& pos ) const {
  std::array<unsigned,3> celn;
  Vector mypos = pos;
  if( nopbc ) {
    mypos = pos - origin;
  }
  Vector fpos=mypbc.realToScaled( mypos );
  for(unsigned j=0; j<3; ++j) {
    celn[j] = std::floor( ( Tools::pbc(fpos[j]) + 0.5 ) * ncells[j] );
    plumed_assert( celn[j]>=0 && celn[j]<ncells[j] ) <<"in link cell "<<celn[j]<<" but should be between 0 and "<<ncells[j]<<" link cell cutoff is "<<link_cutoff<<" position is "<<fpos[0]<<" "<<fpos[1]<<" "<<fpos[2]<<" box is "<<mypbc.getBox()(0,0)<<" "<<mypbc.getBox()(1,1)<<" "<<mypbc.getBox()(2,2);
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
  unsigned maxn = lcell_tots[0];
  for(unsigned i=1; i<lcell_tots.size(); ++i) {
    if( lcell_tots[i]>maxn ) {
      maxn=lcell_tots[i];
    }
  }
  return maxn;
}

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
  natoms_per_list = 27*getMaxInCell();
  if( natoms_per_list>allcells.size() ) {
    natoms_per_list = allcells.size();
  }

  unsigned nlist_sz = nat*( 2 + natoms_per_list );
  nlist.resize( nlist_sz );
  std::vector<unsigned> indices( 1+natoms_per_list );
  std::vector<unsigned> cells_required( getNumberOfCells() );
  for(unsigned i=0; i<pos.size(); ++i) {
    unsigned ncells_required=0;
    addRequiredCells( findMyCell( pos[i] ), ncells_required, cells_required );
    unsigned natoms=1;
    indices[0] = ind[i];
    retrieveAtomsInCells( ncells_required, make_const_view(cells_required), natoms, indices );
    nlist[tind[i]] = 0;
    std::size_t lstart = nat + tind[i]*(1+natoms_per_list);
    for(unsigned j=0; j<natoms; ++j) {
      nlist[ lstart + nlist[tind[i]] ] = indices[j];
      nlist[tind[i]]++;
    }
  }
}

}
