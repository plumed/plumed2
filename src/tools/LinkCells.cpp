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
#include "OpenMP.h"
#include "Exception.h"
#include "Tensor.h"
#include "Tools.h"
#include "View.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

namespace PLMD {

namespace {
/// The index of the cell of a position from its scaled coordinates fpos, which are
/// wrapped in the box (between -0.5 and 0.5) in place. Returns the number of cells if
/// the position is not in a cell. The same operations are in forEachCandidate (LinkCellsView.h)
inline unsigned cellOfScaled( Vector& fpos, const std::array<unsigned,3>& ncells,
                              const std::array<unsigned,3>& nstride ) noexcept {
  unsigned index=0;
  for(unsigned j=0; j<3; ++j) {
    fpos[j]=Tools::pbc(fpos[j]);
    const double x=( fpos[j] + 0.5 ) * ncells[j];
    // as floor(x) in [0,ncells) (also false for NaN); then floor(x) is the truncation
    if(!(x>=0.0 && x<ncells[j])) {
      return ncells[0]*ncells[1]*ncells[2];
    }
    index+=static_cast<unsigned>(x)*nstride[j];
  }
  return index;
}

/// matmul(a,b) for a 3x3 matrix, written out: the generic loops of matmul are not unrolled at -O2
inline Vector matmul3( const Tensor& a, const Vector& b ) noexcept {
  return Vector( a[0][0]*b[0]+a[0][1]*b[1]+a[0][2]*b[2],
                 a[1][0]*b[0]+a[1][1]*b[1]+a[1][2]*b[2],
                 a[2][0]*b[0]+a[2][1]*b[1]+a[2][2]*b[2] );
}
}

LinkCells::LinkCells( Communicator& ) {
}

void LinkCells::setSubdivision( const unsigned n ) {
  plumed_assert( n>=1 && n<=maxSubdivision );
  requestedSubdivision=n;
}

void LinkCells::setLimitCellsByAtoms( bool limit ) {
  limitCellsByAtoms=limit;
}

void LinkCells::setCutoff( const double lcut ) {
  cutoffwasset=true;
  link_cutoff=lcut;
}

double LinkCells::getCutoff() const {
  plumed_assert( cutoffwasset );
  return link_cutoff;
}

void LinkCells::setupCells(View<const Vector> pos) {

  nopbc=true;
  Tensor box;
  box.zero();
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
  cellLimitAtoms = limitCellsByAtoms ? pos.size() : 0;
  createCells(box);
  cellLimitAtoms = 0;
}

void LinkCells::setupCells(const Pbc& pbc ) {

  nopbc=false;

  auto box = pbc.getBox();
  auto determinant = box.determinant();
  plumed_assert(determinant > epsilon) <<"Cell lists cannot be built when passing a box with null volume. Volume is "<<determinant;
  createCells(box);
}

void LinkCells::setupCells(View<const Vector> pos,
                           const Pbc& pbc ) {
  auto box = pbc.getBox();
  if(box(0,0)==0.0 && box(0,1)==0.0 && box(0,2)==0.0 &&
      box(1,0)==0.0 && box(1,1)==0.0 && box(1,2)==0.0 &&
      box(2,0)==0.0 && box(2,1)==0.0 && box(2,2)==0.0) {
    // Create an orthorhombic box around the atomic positions that encompasses every atomic position if there are no pbc
    setupCells(pos);
  } else {
    // Setup the pbc object by copying it from action if the Pbc were set
    cellLimitAtoms = limitCellsByAtoms ? pos.size() : 0;
    setupCells(pbc);
    cellLimitAtoms = 0;
  }
}

void LinkCells::createCells(const PLMD::Tensor& box) {
  mypbc.setBox( box );
// This is the reciprocal lattice
// notice that reciprocal.getRow(0) is a vector that is orthogonal to b and c
// This allows to use linked cells in non orthorhomic boxes
  const auto reciprocal = transpose(mypbc.getInvBox());
  // the width of the cells along each direction must be at least the cutoff divided by
  // the number of cells per cutoff
  std::array<double,3> width, nc;
  for(unsigned k=0; k<3; ++k) {
    width[k]=1.0/reciprocal.getRow(k).modulo();
  }
  nsub=requestedSubdivision;
  while(true) {
    for(unsigned k=0; k<3; ++k) {
      nc[k]=std::floor( width[k] / link_cutoff * nsub );
      // also catches NaNs
      if( !(nc[k]>=1.0) ) {
        nc[k]=1.0;
      }
    }
    // Limit the total number of cells: with a small cutoff in a large box the number of
    // cells (and the memory and the time spent looping over empty cells) would grow
    // with the volume. Larger cells are still correct, they only contain more atoms.
    // With setLimitCellsByAtoms they are also limited by the number of atoms (see maxCellsPerAtom)
    double maxcells=maxNumberOfCells;
    if( cellLimitAtoms>0 ) {
      maxcells=std::min(maxcells,std::max(minCellsForAtomLimit,maxCellsPerAtom*cellLimitAtoms));
    }
    while(nc[0]*nc[1]*nc[2]>maxcells) {
      const unsigned k=std::distance(nc.begin(),std::max_element(nc.begin(),nc.end()));
      const double factor=std::cbrt(nc[0]*nc[1]*nc[2]/maxcells);
      nc[k]=std::max(1.0,std::min(nc[k]-1.0,std::floor(nc[k]/factor)));
    }
    // smaller cells are used only if the cells to be searched are all different
    if( nsub==1 || (nc[0]>=2*nsub+1 && nc[1]>=2*nsub+1 && nc[2]>=2*nsub+1) ) {
      break;
    }
    nsub=1;
  }
  for(unsigned k=0; k<3; ++k) {
    ncells[k]=static_cast<unsigned>(nc[k]);
    cellWidth[k]=width[k]/nc[k];
  }
  // see the comment of skewFactor
  Tensor cosines;
  bool orthogonal=true;
  for(unsigned k=0; k<3; ++k) {
    for(unsigned l=0; l<3; ++l) {
      cosines[k][l]=dotProduct(reciprocal.getRow(k),reciprocal.getRow(l))*width[k]*width[l];
      if(l!=k && std::fabs(cosines[k][l])>1e-12) {
        orthogonal=false;
      }
    }
  }
  skewFactor=1.0;
  if( !orthogonal ) {
    Vector eigenvalues;
    Tensor eigenvectors;
    diagMatSym(cosines,eigenvalues,eigenvectors);
    skewFactor=std::max(eigenvalues[0],std::max(eigenvalues[1],eigenvalues[2]));
    // it cannot be smaller than one (the trace is three): this only protects from a failed diagonalization
    if( !(skewFactor>=1.0) ) {
      skewFactor=3.0;
    }
  }
  for(unsigned k=0; k<3; ++k) {
    if( nopbc ) {
      lattice[k].zero();
    } else {
      lattice[k]=box.getRow(k);
    }
  }
  // Setup the strides
  nstride[0]=1;
  nstride[1]=ncells[0];
  nstride[2]=ncells[0]*ncells[1];
}

void LinkCells::buildCellLists( View<const Vector> pos,
                                View<const unsigned> indices,
                                const Pbc& pbc,
                                const bool withPositions ) {
  setupCells(pos,pbc);
  resetCollection(innerCollection,pos, indices, withPositions);
}

LinkCells::CellCollection LinkCells::getCollection( View<const Vector> pos,
    View<const unsigned> indices,
    const bool withPositions ) {
  CellCollection collection;
  resetCollection(collection,pos, indices, withPositions);
  return collection;
}

void LinkCells::resetCollection(LinkCells::CellCollection &collection,
                                View<const Vector> pos,
                                View<const unsigned> indices,
                                const bool withPositions ) {
  plumed_assert( cutoffwasset && pos.size()==indices.size() );
  const unsigned nat=pos.size();
  // Resize the lists (every element is overwritten below)
  allcells.resize( nat );
  collection.lcell_lists.resize( nat );
  if( withPositions ) {
    allwrapped.resize( nat );
    collection.lcell_pos.resize( nat );
  } else {
    collection.lcell_pos.clear();
  }
  const unsigned ncellstot=ncells[0]*ncells[1]*ncells[2];
  collection.lcell_starts.resize( ncellstot );
  collection.lcell_tots.resize( ncellstot );

  // This is done locally on every process: finding the cell of an atom costs a few
  // nanoseconds, much less than the two MPI_Allreduce (of the cell of each atom and of
  // the number of atoms in each cell) that were needed when the atoms were split among processes.
  // The cells are split in nranges contiguous ranges, each one handled by one thread: the
  // thread goes through all the atoms, and counts and then places those in its cells. The
  // ranges do not share the cells, nor the parts of the lists they write, and the atoms
  // of each cell are in increasing order, whatever the number of ranges and of threads.
  // Each range goes twice through the cells of all the atoms (about a nanosecond per atom),
  // which is cheap compared to finding the cells only with enough atoms per range.
  // All the threads, with at least minAtomsPerThread atoms each, or one (see OpenMP.h),
  // and at most one range per thread
  const unsigned nt=OpenMP::getGoodNumThreads( nat, minAtomsPerThread );
  const unsigned nranges=std::max( 1u, std::min( { nt, ncellstot, nat/minAtomsPerRange } ) );
  rangetots.resize( nranges );
  // the transposed matrices used by Pbc::realToScaled and Pbc::scaledToReal
  const Tensor invBoxT=transpose( mypbc.getInvBox() );
  const Tensor boxT=transpose( mypbc.getBox() );
  unsigned invalid=0;
  // the three steps share one parallel region, separated by the implicit barriers of the
  // loops. No exceptions inside the parallel region: invalid cells are checked after it
  // the loops over the atoms are shared by the team, the loops over the ranges work
  // with any number of threads
  #pragma omp parallel num_threads(nt)
  {
    // 1. the cell of each atom (and its position wrapped in the box)
    #pragma omp for schedule(static) reduction(+:invalid)
    for(unsigned i=0; i<nat; ++i) {
      // as findCellUnchecked and wrapPosition
      Vector fpos=matmul3( invBoxT, nopbc ? pos[i]-origin : pos[i] );
      const unsigned mycell=cellOfScaled( fpos, ncells, nstride );
      allcells[i]=mycell;
      if( withPositions ) {
        allwrapped[i]=( nopbc ? pos[i] : matmul3( boxT, fpos ) );
      }
      if( mycell>=ncellstot ) {
        invalid++;
      }
    }
    // the loops over the ranges work whatever the number of threads of the team
    if( invalid==0 ) {
      // 2. the number of atoms in each cell, and in each range
      #pragma omp for schedule(static,1)
      for(unsigned r=0; r<nranges; ++r) {
        const unsigned first=std::size_t(ncellstot)*r/nranges;
        const unsigned size=std::size_t(ncellstot)*(r+1)/nranges-first;
        unsigned* tots=collection.lcell_tots.data()+first;
        std::fill( tots, tots+size, 0u );
        for(unsigned i=0; i<nat; ++i) {
          // also false if allcells[i]<first (unsigned)
          const unsigned c=allcells[i]-first;
          if( c<size ) {
            tots[c]++;
          }
        }
        unsigned tot=0;
        for(unsigned c=0; c<size; ++c) {
          tot+=tots[c];
        }
        rangetots[r]=tot;
      }
      // 3. where each cell starts, and the atoms in the cells
      #pragma omp for schedule(static,1)
      for(unsigned r=0; r<nranges; ++r) {
        const unsigned first=std::size_t(ncellstot)*r/nranges;
        const unsigned size=std::size_t(ncellstot)*(r+1)/nranges-first;
        const unsigned* tots=collection.lcell_tots.data()+first;
        unsigned* starts=collection.lcell_starts.data()+first;
        unsigned tot=0;
        for(unsigned q=0; q<r; ++q) {
          tot+=rangetots[q];
        }
        for(unsigned c=0; c<size; ++c) {
          starts[c]=tot;
          tot+=tots[c];
        }
        // starts[c] is where the next atom of the cell goes...
        for(unsigned i=0; i<nat; ++i) {
          const unsigned c=allcells[i]-first;
          if( c<size ) {
            const unsigned myind=starts[c]++;
            collection.lcell_lists[myind]=indices[i];
            if( withPositions ) {
              collection.lcell_pos[myind]=allwrapped[i];
            }
          }
        }
        // ...and then it is moved back to the first atom
        for(unsigned c=0; c<size; ++c) {
          starts[c]-=tots[c];
        }
      }
    }
  }
  if( invalid>0 ) {
    for(unsigned i=0; i<nat; ++i) {
      if( allcells[i]>=ncellstot ) {
        // this raises an error with an explanation
        findMyCell( pos[i] );
        plumed_error()<<"atom "<<i<<" is not inside a link cell";
      }
    }
  }
}

void LinkCells::addRequiredCells( const std::array<unsigned,3>& celn,
                                  unsigned& ncells_required,
                                  std::vector<unsigned>& cells_required,
                                  const bool usePbc) const {
#define LINKC_MIN(n) LinkCellsView::firstCell(celn[n],ncells[n],usePbc,nsub)
#define LINKC_MAX(n) LinkCellsView::endCell(celn[n],ncells[n],usePbc,nsub)
#define LINKC_PBC(n,num) ((n<0)? n+num : n%num )
  unsigned nnew_cells=0;
  plumed_dbg_assert((celn[0] < ncells[0]) && (celn[1] < ncells[1]) && (celn[2] < ncells[2]));

  const int minX=LINKC_MIN(0);
  const int maxX=LINKC_MAX(0);
  const int minY=LINKC_MIN(1);
  const int maxY=LINKC_MAX(1);
  const int minZ=LINKC_MIN(2);
  const int maxZ=LINKC_MAX(2);

  for(int nx=minX; nx<maxX; ++nx) {
    int xval=LINKC_PBC(nx,ncells[0])*nstride[0];
    for(int ny=minY; ny<maxY; ++ny) {
      int yval=LINKC_PBC(ny,ncells[1])*nstride[1];
      for(int nz=minZ; nz<maxZ; ++nz) {
        int zval=LINKC_PBC(nz,ncells[2])*nstride[2];

        unsigned mybox=xval+yval+zval;
        plumed_assert(mybox < getNumberOfCells());
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
#undef LINKC_MIN
#undef LINKC_MAX
#undef LINKC_PBC
}

bool LinkCells::cellsFixImages() const {
  const unsigned needed=2*nsub+1;
  return nopbc || ( ncells[0]>=needed && ncells[1]>=needed && ncells[2]>=needed );
}

Vector LinkCells::wrapPosition( const Vector& pos ) const {
  if( nopbc ) {
    return pos;
  }
  // the same operations used in findCellUnchecked
  Vector fpos=matmul3( transpose(mypbc.getInvBox()), pos );
  for(unsigned j=0; j<3; ++j) {
    fpos[j]=Tools::pbc(fpos[j]);
  }
  return matmul3( transpose(mypbc.getBox()), fpos );
}

LinkCellsView LinkCells::getView( const CellCollection& collection, const bool wrapCells ) const {
  plumed_assert( cellsFixImages() ) << "the cells do not fix the periodic images: use addRequiredCells and Pbc::distance";
  plumed_assert( collection.lcell_pos.size()==collection.lcell_lists.size() )
      << "the collection was built without its positions";
  plumed_assert( collection.lcell_starts.size()==getNumberOfCells() );
  LinkCellsView view;
  view.boxT=mypbc.getBox().transpose();
  view.invBoxT=mypbc.getInvBox().transpose();
  view.origin=( nopbc ? origin : Vector(0.0,0.0,0.0) );
  for(unsigned k=0; k<3; ++k) {
    view.lattice[k]=lattice[k];
    view.ncells[k]=static_cast<int>(ncells[k]);
    view.nstride[k]=nstride[k];
    view.cellWidth[k]=cellWidth[k];
  }
  view.nsub=static_cast<int>(nsub);
  // cells farther than this cannot contain atoms within the cutoff (see skewFactor)
  view.reach2=link_cutoff*link_cutoff*skewFactor*(1.0+1e-10);
  view.cutoff2=link_cutoff*link_cutoff;
  view.nopbc=nopbc;
  view.wrapCells=wrapCells;
  view.ncellstot=getNumberOfCells();
  view.natoms=collection.lcell_lists.size();
  view.starts=collection.lcell_starts.data();
  view.tots=collection.lcell_tots.data();
  view.lists=collection.lcell_lists.data();
  view.positions=collection.lcell_pos.data();
  return view;
}

LinkCellsView LinkCells::getView( const bool wrapCells ) const {
  return getView( innerCollection, wrapCells );
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

unsigned LinkCells::findCellUnchecked( Vector mypos ) const noexcept {
  Vector fpos=matmul3( transpose(mypbc.getInvBox()), nopbc ? mypos-origin : mypos );
  return cellOfScaled( fpos, ncells, nstride );
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

std::array<unsigned,3> LinkCells::findMyCell( unsigned cellIndex ) const {
  std::array<unsigned,3> cell;
//std::div returns quotient and remainder of an integer division
  auto r=std::div(static_cast<int>(cellIndex),nstride[2]);
  cell[2]=r.quot;
  r=std::div(r.rem,nstride[1]);
  cell[1]=r.quot;
  //nstride[0] is hardcoded to be 1;
  cell[0]=r.rem;
  return cell;
}

unsigned LinkCells::convertIndicesToIndex( std::array<unsigned,3> cellCoord) const {
  //nstride[0] is hardcoded to be 1;
  return cellCoord[0]//*nstride[0]
         + cellCoord[1]*nstride[1]
         + cellCoord[2]*nstride[2];
}

unsigned LinkCells::findCell( const Vector& pos ) const {
  return convertIndicesToIndex(findMyCell(pos));
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

namespace {
/// The visitor of createNeighborList: appends the candidates within the cutoff,
/// except the central atom, to indices
struct KeepWithinCutoff {
  unsigned* indices;
  unsigned n;
  unsigned self;
  double cutoff2;
  const unsigned* lists;
  inline void operator()( const std::size_t j, const Vector& d ) {
    const unsigned atom=lists[j];
    // no branch: most of the candidates are beyond the cutoff, but not in a predictable way
    // (indices has one more element: a candidate is written before it is known if it is kept)
    indices[n]=atom;
    n+=( (d.modulo2()<=cutoff2) & (atom!=self) );
  }
};
}

void LinkCells::createNeighborList( unsigned nat,
                                    View<const Vector> pos,
                                    View<const unsigned> ind,
                                    View<const unsigned> tind,
                                    View<const Vector> neigh_pos,
                                    View<const unsigned> neigh_ind,
                                    const Pbc& pbc,
                                    unsigned& natoms_per_list,
                                    std::vector<std::size_t>& nlist,
                                    const bool applyCutoff ) {
  const bool useCutoff = applyCutoff && link_cutoff<std::sqrt(std::numeric_limits<double>::max());
  // with the distance test smaller cells give fewer candidates beyond the cutoff
  requestedSubdivision = useCutoff ? maxSubdivision : 1;
  // the lists are built at every step: no more cells than about twice the atoms
  const bool limitWasSet = limitCellsByAtoms;
  limitCellsByAtoms = true;
  // with the distance test the wrapped positions are stored with the cells (see getView)
  buildCellLists( neigh_pos, neigh_ind, pbc, useCutoff );
  limitCellsByAtoms = limitWasSet;
  // the largest number of atoms in the cells searched for one atom
  unsigned maxcandidates=innerCollection.getMaximimumCombination((2*nsub+1)*(2*nsub+1)*(2*nsub+1));
  if( maxcandidates>innerCollection.lcell_lists.size() ) {
    maxcandidates = innerCollection.lcell_lists.size();
  }
  const double cutoff2 = link_cutoff*link_cutoff;
  // When the cell of a candidate fixes its periodic image, the minimum image distance
  // is not needed for each candidate (they are many more than the atoms within the cutoff):
  // see getView
  const bool useShifts = useCutoff && cellsFixImages();
  const LinkCellsView view = useShifts ? getView() : LinkCellsView{};
  // otherwise the distance test needs the position of a candidate from its index: the
  // cell lists hold the indices of neigh_ind, the positions are in neigh_pos
  std::vector<unsigned> posOfIndex;
  if( useCutoff && !useShifts ) {
    unsigned maxind=0;
    for(unsigned j=0; j<neigh_ind.size(); ++j) {
      maxind = std::max( maxind, neigh_ind[j] );
    }
    posOfIndex.resize( maxind+1 );
    for(unsigned j=0; j<neigh_ind.size(); ++j) {
      posOfIndex[neigh_ind[j]] = j;
    }
  }
  // The lists are first stored one after the other by each thread, and copied to nlist
  // when the length of the longest one is known. That length is much smaller than the
  // number of atoms in the cells that are searched: sizing the slots of nlist on the
  // latter spreads the lists over a memory that is ten times larger, which is slow to
  // write here and to read for whoever uses the lists.
  std::vector<unsigned> nneigh( pos.size() );
  unsigned maxneigh=1;
#ifdef _OPENMP
  // without OpenMP nt would be unused (an error with -Werror)
  // the search of the neighbors of an atom costs much more than reading it (see OpenMP.h)
  const unsigned nt=OpenMP::getGoodNumThreads( pos.size(), minAtomsPerThreadNeighborList );
#endif
  #pragma omp parallel num_threads(nt)
  {
    // the lists found by this thread
    std::vector<unsigned> found;
    // one more element: a candidate is written before it is known if it is kept
    std::vector<unsigned> indices( 2+maxcandidates );
    std::vector<unsigned> cells_required( getNumberOfCells() );
    unsigned mymax=1;
    #pragma omp for schedule(static) nowait
    for(unsigned i=0; i<pos.size(); ++i) {
      unsigned natoms=1;
      indices[0] = ind[i];
      if( useShifts ) {
        KeepWithinCutoff keep{indices.data(),natoms,ind[i],cutoff2,view.lists};
        forEachCandidate( view, pos[i], keep );
        natoms = keep.n;
      } else {
        unsigned ncells_required=0;
        addRequiredCells( findMyCell( pos[i] ), ncells_required, cells_required );
        retrieveAtomsInCells( ncells_required,
                              make_const_view(cells_required),
                              natoms,
                              indices,
                              ind[i]);
        if( natoms>neigh_ind.size() ) {
          // the cells wrap around the box: the same atom was found more than once,
          // use the whole list of candidates instead
          natoms = 1 + neigh_ind.size();
          indices[0] = ind[i];
          std::copy( neigh_ind.begin(), neigh_ind.end(), indices.begin()+1 );
        }
        if( useCutoff ) {
          unsigned nkept=1;
          for(unsigned j=1; j<natoms; ++j) {
            const Vector d = pbc.distance( pos[i], neigh_pos[posOfIndex[indices[j]]] );
            if( d.modulo2()<=cutoff2 ) {
              indices[nkept] = indices[j];
              nkept++;
            }
          }
          natoms = nkept;
        }
      }
      nneigh[i] = natoms;
      mymax = std::max( mymax, natoms );
      found.insert( found.end(), indices.begin(), indices.begin()+natoms );
    }
    #pragma omp critical
    maxneigh = std::max( maxneigh, mymax );
    #pragma omp barrier
    #pragma omp single
    {
      natoms_per_list = maxneigh-1;
      nlist.resize( std::size_t(nat)*( 2 + std::size_t(natoms_per_list) ) );
    }
    // the same atoms as in the first loop (static schedule, same parallel region):
    // each thread copies its lists in the order in which it found them
    std::size_t next=0;
    #pragma omp for schedule(static)
    for(unsigned i=0; i<pos.size(); ++i) {
      const std::size_t lstart = nat + std::size_t(tind[i])*(1+std::size_t(natoms_per_list));
      nlist[tind[i]] = nneigh[i];
      std::copy(found.begin()+next,
                found.begin()+next+nneigh[i],
                nlist.begin()+lstart);
      next += nneigh[i];
    }
  }
}

unsigned LinkCells::CellCollection::getMaximimumCombination(const unsigned numCells) const {
  //this is not efficient (there is a copy), but in principle it should be called not much times:
  //nth_element order by partition until the array is partially sorted
  //(it stops when the pivot is at the asked point, so all the elements to the right of the nth element are higher/lower than it)
  // with no more atoms than cells, the most populated numCells cells contain all the atoms
  if (numCells >= lcell_tots.size() || lcell_lists.size() <= numCells) {
    return lcell_lists.size();
  }
  auto tmpCopy=lcell_tots;
  auto nth_place=tmpCopy.begin() + numCells-1;
  std::nth_element(tmpCopy.begin(),nth_place,tmpCopy.end(),std::greater<>());

  //the upper limit is not included, hence the +1
  return std::accumulate(tmpCopy.begin(),nth_place+1,0);
}

} //namespace PLMD
