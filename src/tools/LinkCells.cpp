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

LinkCells::LinkCells( Communicator& ) {
}

void LinkCells::setSubdivision( const unsigned n ) {
  plumed_assert( n>=1 && n<=maxSubdivision );
  requestedSubdivision=n;
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
  createCells(box);
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
    setupCells(pbc);
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
    while(nc[0]*nc[1]*nc[2]>maxNumberOfCells) {
      const unsigned k=std::distance(nc.begin(),std::max_element(nc.begin(),nc.end()));
      const double factor=std::cbrt(nc[0]*nc[1]*nc[2]/maxNumberOfCells);
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
                                const Pbc& pbc ) {
  setupCells(pos,pbc);
  resetCollection(innerCollection,pos, indices);
}

LinkCells::CellCollection LinkCells::getCollection( View<const Vector> pos,
    View<const unsigned> indices) {
  CellCollection collection;
  resetCollection(collection,pos, indices);
  return collection;
}

void LinkCells::resetCollection(LinkCells::CellCollection &collection,
                                View<const Vector> pos,
                                View<const unsigned> indices) {
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

  // Find out what cell everyone is in.
  // This is done locally on every process: finding the cell of an atom costs a few
  // nanoseconds, much less than the two MPI_Allreduce (of the cell of each atom and of
  // the number of atoms in each cell) that were needed when the atoms were split among processes.
#ifdef _OPENMP
  // without OpenMP nt would be unused (an error with -Werror)
  const unsigned nt=(nat>=minAtomsPerThread)?
                    std::min<unsigned>(OpenMP::getNumThreads(),nat/minAtomsPerThread) : 1;
#endif
  #pragma omp parallel for num_threads(nt) schedule(static)
  for(unsigned i=0; i<nat; ++i) {
    // no exceptions inside the parallel region: invalid cells are checked below
    allcells[i]=findCellUnchecked( pos[i] );
  }
  for(unsigned i=0; i<nat; ++i) {
    if(allcells[i]>=ncellstot) {
      // this raises an error with an explanation
      findMyCell( pos[i] );
      plumed_error()<<"atom "<<i<<" is not inside a link cell";
    }
    collection.lcell_tots[allcells[i]]++;
  }
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

// The first of the cells to be searched along one direction, around the cell celn
constexpr inline int min_cell(int celn, const int ncells, const bool usePBC, const int nsub=1) {
  if(nsub>1) {
    // only used with at least 2*nsub+1 cells
    return (usePBC) ? celn-nsub : std::max(celn-nsub,0);
  }
  celn+=(ncells<2)?0:-1;
  return (usePBC) ? celn : std::max(celn,0);
}

// One past the last of the cells to be searched along one direction, around the cell celn
constexpr inline int max_cell(int celn, const int ncells, const bool usePBC, const int nsub=1) {
  if(nsub>1) {
    return (usePBC) ? celn+nsub+1 : std::min(celn+nsub+1,ncells);
  }
  //if we are not using pbcs the minimum is 0 and not -1 (if there are only 2 cells)
  //so the maximum should be 2, so that cell 0 can check cell 1
  celn+=(ncells<3 && usePBC)?1:2;
  return (usePBC) ? celn : std::min(celn,ncells);
}

void LinkCells::addRequiredCells( const std::array<unsigned,3>& celn,
                                  unsigned& ncells_required,
                                  std::vector<unsigned>& cells_required,
                                  const bool usePbc) const {
#define LINKC_MIN(n) min_cell(celn[n],ncells[n],usePbc,nsub)
#define LINKC_MAX(n) max_cell(celn[n],ncells[n],usePbc,nsub)
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
  Vector fpos=mypbc.realToScaled( pos );
  for(unsigned j=0; j<3; ++j) {
    fpos[j]=Tools::pbc(fpos[j]);
  }
  return mypbc.scaledToReal( fpos );
}

unsigned LinkCells::getNeighborRuns( const Vector& pos,
                                     const CellCollection& collection,
                                     const bool usePbc,
                                     Vector& wrapped,
                                     std::array<CellRun,maxRuns>& runs ) const {
  plumed_dbg_assert( cellsFixImages() );
  constexpr int maxoff=maxSubdivision;
  // the cell of the position and where it is inside it (from 0 to 1), as in findCellUnchecked
  Vector fpos=mypbc.realToScaled( nopbc ? pos-origin : pos );
  const int nc[3]= {static_cast<int>(ncells[0]),static_cast<int>(ncells[1]),static_cast<int>(ncells[2])};
  int cell[3], first[3], last[3];
  // the squared distance of the position from the cells that are o cells away along each direction
  double gap2[3][2*maxoff+1];
  for(unsigned k=0; k<3; ++k) {
    fpos[k]=Tools::pbc(fpos[k]);
    const double x=( fpos[k] + 0.5 ) * nc[k];
    cell[k]=std::min(std::max(static_cast<int>(std::floor(x)),0),nc[k]-1);
    const double inside=x-cell[k];
    first[k]=min_cell(cell[k],nc[k],usePbc,nsub)-cell[k];
    last[k]=max_cell(cell[k],nc[k],usePbc,nsub)-1-cell[k];
    // with few cells the same cell is on both sides, and its distance is not known
    const bool known=(nc[k]>=static_cast<int>(2*nsub+1));
    for(int o=first[k]; o<=last[k]; ++o) {
      double g=0.0;
      if(known && o!=0) {
        g=(o>0 ? o-inside : inside-o-1.0)*cellWidth[k];
      }
      gap2[k][o+maxoff]=g*g;
    }
  }
  wrapped=(nopbc ? pos : mypbc.scaledToReal( fpos ));
  // cells farther than this cannot contain atoms within the cutoff (see skewFactor)
  const double reach2=link_cutoff*link_cutoff*skewFactor*(1.0+1e-10);
  const auto & starts=collection.lcell_starts;
  const auto & tots=collection.lcell_tots;
  unsigned n=0;
  for(int oz=first[2]; oz<=last[2]; ++oz) {
    const double gz2=gap2[2][oz+maxoff];
    if(gz2>reach2) {
      continue;
    }
    // the cell wrapped in the box, and the number of boxes it was moved by
    int iz=cell[2]+oz, sz=0;
    if(iz<0) {
      iz+=nc[2];
      sz=-1;
    } else if(iz>=nc[2]) {
      iz-=nc[2];
      sz=1;
    }
    for(int oy=first[1]; oy<=last[1]; ++oy) {
      const double gyz2=gz2+gap2[1][oy+maxoff];
      if(gyz2>reach2) {
        continue;
      }
      int iy=cell[1]+oy, sy=0;
      if(iy<0) {
        iy+=nc[1];
        sy=-1;
      } else if(iy>=nc[1]) {
        iy-=nc[1];
        sy=1;
      }
      // the cells of this row that are not too far
      int a=first[0], b=last[0];
      while(a<0 && gyz2+gap2[0][a+maxoff]>reach2) {
        ++a;
      }
      while(b>0 && gyz2+gap2[0][b+maxoff]>reach2) {
        --b;
      }
      a+=cell[0];
      b+=cell[0];
      const std::size_t row=iy*nstride[1]+iz*nstride[2];
      const Vector rowshift=double(sy)*lattice[1]+double(sz)*lattice[2];
      // The cells of a row are consecutive in the collection: they are one run,
      // or two if some of them are on the other side of the box
      if(a<0) {
        const std::size_t xa=a+nc[0], xb=std::min(b,-1)+nc[0];
        runs[n].begin=starts[row+xa];
        runs[n].end=starts[row+xb]+tots[row+xb];
        runs[n].shift=rowshift-lattice[0];
        n+=(runs[n].end>runs[n].begin);
        a=0;
      }
      const int bmid=std::min(b,nc[0]-1);
      if(a<=bmid) {
        runs[n].begin=starts[row+a];
        runs[n].end=starts[row+bmid]+tots[row+bmid];
        runs[n].shift=rowshift;
        n+=(runs[n].end>runs[n].begin);
      }
      if(b>=nc[0]) {
        const std::size_t xa=std::max(a,nc[0])-nc[0], xb=b-nc[0];
        runs[n].begin=starts[row+xa];
        runs[n].end=starts[row+xb]+tots[row+xb];
        runs[n].shift=rowshift+lattice[0];
        n+=(runs[n].end>runs[n].begin);
      }
    }
  }
  return n;
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
  if( nopbc ) {
    mypos = mypos - origin;
  }
  const Vector fpos=mypbc.realToScaled( mypos );
  const unsigned ncellstot=ncells[0]*ncells[1]*ncells[2];
  unsigned index=0;
  for(unsigned j=0; j<3; ++j) {
    const double c=std::floor( ( Tools::pbc(fpos[j]) + 0.5 ) * ncells[j] );
    if(!(c>=0.0 && c<ncells[j])) {
      return ncellstot;
    }
    index+=static_cast<unsigned>(c)*nstride[j];
  }
  return index;
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
  // The distance test needs the position of a candidate from its index: the
  // cell lists hold the indices of neigh_ind, the positions are in neigh_pos
  const bool useCutoff = applyCutoff && link_cutoff<std::sqrt(std::numeric_limits<double>::max());
  // with the distance test smaller cells give fewer candidates beyond the cutoff
  requestedSubdivision = useCutoff ? maxSubdivision : 1;
  buildCellLists( neigh_pos, neigh_ind, pbc );
  // the largest number of atoms in the cells searched for one atom
  unsigned maxcandidates=innerCollection.getMaximimumCombination((2*nsub+1)*(2*nsub+1)*(2*nsub+1));
  if( maxcandidates>innerCollection.lcell_lists.size() ) {
    maxcandidates = innerCollection.lcell_lists.size();
  }
  const double cutoff2 = link_cutoff*link_cutoff;
  std::vector<unsigned> posOfIndex;
  if( useCutoff ) {
    unsigned maxind=0;
    for(unsigned j=0; j<neigh_ind.size(); ++j) {
      maxind = std::max( maxind, neigh_ind[j] );
    }
    posOfIndex.resize( maxind+1 );
    for(unsigned j=0; j<neigh_ind.size(); ++j) {
      posOfIndex[neigh_ind[j]] = j;
    }
  }
  // When the cell of a candidate fixes its periodic image, the minimum image distance
  // is not needed for each candidate (they are many more than the atoms within the cutoff):
  // the candidates are wrapped in the box once, and stored in the order of the cells
  const bool useShifts = useCutoff && cellsFixImages();
  std::vector<Vector> cellpos;
  if( useShifts ) {
    cellpos.resize( innerCollection.lcell_lists.size() );
    for(std::size_t j=0; j<cellpos.size(); ++j) {
      cellpos[j] = wrapPosition( neigh_pos[posOfIndex[innerCollection.lcell_lists[j]]] );
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
  const unsigned nt=(pos.size()>=minAtomsPerThread/16)?
                    std::min<unsigned>(OpenMP::getNumThreads(),pos.size()/(minAtomsPerThread/16)) : 1;
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
        std::array<CellRun,maxRuns> runs;
        Vector mypos;
        const unsigned nruns = getNeighborRuns( pos[i], innerCollection, true, mypos, runs );
        const unsigned myind = ind[i];
        for(unsigned k=0; k<nruns; ++k) {
          const Vector shift = runs[k].shift - mypos;
          for(std::size_t j=runs[k].begin; j<runs[k].end; ++j) {
            const Vector d = cellpos[j] + shift;
            const unsigned myatom = innerCollection.lcell_lists[j];
            // no branch: most of the candidates are beyond the cutoff, but not in a predictable way
            indices[natoms] = myatom;
            natoms += ( (d.modulo2()<=cutoff2) & (myatom!=myind) );
          }
        }
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
  if (numCells >= lcell_tots.size()) {
    return std::accumulate(lcell_tots.begin(),lcell_tots.end(),0);
  }
  auto tmpCopy=lcell_tots;
  auto nth_place=tmpCopy.begin() + numCells-1;
  std::nth_element(tmpCopy.begin(),nth_place,tmpCopy.end(),std::greater<>());

  //the upper limit is not included, hence the +1
  return std::accumulate(tmpCopy.begin(),nth_place+1,0);
}

} //namespace PLMD
