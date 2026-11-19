/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#ifndef __PLUMED_tools_LinkCellsView_h
#define __PLUMED_tools_LinkCellsView_h

#include "Vector.h"
#include "Tensor.h"
#include "Tools.h"
#include <cmath>
#include <cstddef>
#include <limits>

namespace PLMD {

/// \ingroup TOOLBOX
/// The link cells built by LinkCells, as a plain struct that can be copied to a GPU:
/// the geometry of the cells and pointers to the arrays of one LinkCells::CellCollection
/// (which must have been built with its positions). It does not own the arrays: it is valid
/// until that collection is rebuilt. Get it with LinkCells::getView, and search it with
/// forEachCandidate or forEachNeighbor
struct LinkCellsView {
/// The largest number of cells per cutoff (see LinkCells::setSubdivision)
  static constexpr unsigned maxSubdivision=2;
/// The transposed box and its inverse: the scaled coordinates of x are matmul(invBoxT,x),
/// as in Pbc::realToScaled. Without periodic boundary conditions, the box around the atoms
  Tensor boxT, invBoxT;
/// The center of the box without periodic boundary conditions (zero with them)
  Vector origin;
/// The lattice vectors (rows of the box), zero without periodic boundary conditions
  Vector lattice[3];
/// The number of cells in each direction
  int ncells[3]= {0,0,0};
/// The strides of the cells in each direction (nstride[0] is 1)
  unsigned nstride[3]= {1,0,0};
/// The number of cells per cutoff in use (1 or 2)
  int nsub=1;
/// The width of the cells: the distance between their opposite faces
  double cellWidth[3]= {0.0,0.0,0.0};
/// The cells farther than sqrt(reach2) are not searched (cutoff^2 times skewFactor, see LinkCells)
  double reach2=0.0;
/// The square of the cutoff
  double cutoff2=0.0;
/// True without periodic boundary conditions
  bool nopbc=false;
/// If true the cells on the other side of the box are searched (the usePbc argument of
/// LinkCells::addRequiredCells)
  bool wrapCells=true;
/// The number of cells
  unsigned ncellstot=0;
/// The number of atoms in the cells
  std::size_t natoms=0;
/// The first element of lists of each cell
  const unsigned* starts=nullptr;
/// The number of atoms in each cell
  const unsigned* tots=nullptr;
/// The indices of the atoms, in the order of the cells
  const unsigned* lists=nullptr;
/// The positions of the atoms wrapped in the box, in the same order as lists
  const Vector* positions=nullptr;

/// The first of the cells to be searched along one direction, around the cell celn
  static constexpr int firstCell(int celn, const int ncells, const bool usePBC, const int nsub=1) {
    if(nsub>1) {
      // only used with at least 2*nsub+1 cells
      return (usePBC) ? celn-nsub : (celn-nsub>0 ? celn-nsub : 0);
    }
    celn+=(ncells<2)?0:-1;
    return (usePBC) ? celn : (celn>0 ? celn : 0);
  }
/// One past the last of the cells to be searched along one direction, around the cell celn
  static constexpr int endCell(int celn, const int ncells, const bool usePBC, const int nsub=1) {
    if(nsub>1) {
      return (usePBC) ? celn+nsub+1 : (celn+nsub+1<ncells ? celn+nsub+1 : ncells);
    }
    //if we are not using pbcs the minimum is 0 and not -1 (if there are only 2 cells)
    //so the maximum should be 2, so that cell 0 can check cell 1
    celn+=(ncells<3 && usePBC)?1:2;
    return (usePBC) ? celn : (celn<ncells ? celn : ncells);
  }

#ifdef _OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], starts[0:ncellstot], tots[0:ncellstot], lists[0:natoms], positions[0:natoms])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(positions[0:natoms], lists[0:natoms], tots[0:ncellstot], starts[0:ncellstot], this[0:1])
  }
#endif //_OPENACC
};

/// Calls visit(j,d) for each atom of the cells that can contain atoms within the cutoff
/// from the position pos, in the order of the cells. j is the index of the atom in
/// view.lists and view.positions, d the vector from pos to the image of the atom that is
/// in those cells: the minimum image if the atom is within the cutoff. The cutoff is not
/// tested, and the atom at pos (if it is in the cells) is visited too. visit is a functor.
/// Without periodic boundary conditions pos can be outside the cells: it is then moved to
/// the other side of the box, and its candidates are found through the cells on the other
/// side if wrapCells is true. This is correct within half a cutoff from the cells, that is
/// for any position within the cutoff from an atom (the cells extend at least half a
/// cutoff beyond the atoms).
/// No `acc routine` pragma: as the templated kernels of the actions, the device versions
/// of these templates are generated where they are called from a compute region
template<typename Visitor>
inline void forEachCandidate( const LinkCellsView& view, const Vector& pos, Visitor& visit ) {
  constexpr int maxoff=LinkCellsView::maxSubdivision;
  // the cell of the position and where it is inside it (from 0 to 1), as in LinkCells::findCellUnchecked
  Vector fpos=matmul( view.invBoxT, view.nopbc ? pos-view.origin : pos );
  const int* nc=view.ncells;
  int cell[3], first[3], last[3];
  // the squared distance of the position from the cells that are o cells away along each direction
  double gap2[3][2*maxoff+1];
  for(unsigned k=0; k<3; ++k) {
    fpos[k]=Tools::pbc(fpos[k]);
    const double x=( fpos[k] + 0.5 ) * nc[k];
    const int c=static_cast<int>(std::floor(x));
    cell[k]=(c<0 ? 0 : (c>nc[k]-1 ? nc[k]-1 : c));
    const double inside=x-cell[k];
    first[k]=LinkCellsView::firstCell(cell[k],nc[k],view.wrapCells,view.nsub)-cell[k];
    last[k]=LinkCellsView::endCell(cell[k],nc[k],view.wrapCells,view.nsub)-1-cell[k];
    // with few cells the same cell is on both sides, and its distance is not known
    const bool known=(nc[k]>=2*view.nsub+1);
    for(int o=first[k]; o<=last[k]; ++o) {
      double g=0.0;
      if(known && o!=0) {
        g=(o>0 ? o-inside : inside-o-1.0)*view.cellWidth[k];
      }
      gap2[k][o+maxoff]=g*g;
    }
  }
  const Vector wrapped=( view.nopbc ? pos : matmul( view.boxT, fpos ) );
  const double reach2=view.reach2;
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
      const std::size_t row=iy*view.nstride[1]+iz*view.nstride[2];
      const Vector rowshift=double(sy)*view.lattice[1]+double(sz)*view.lattice[2];
      // The cells of a row are consecutive in lists: their atoms are one run of
      // consecutive elements, or two or three if some cells are on the other side of the box.
      // Each run comes with the vector that moves its atoms next to pos
      std::size_t begin[3], end[3];
      Vector shift[3];
      int nruns=0;
      if(a<0) {
        const std::size_t xa=a+nc[0], xb=(b<-1 ? b : -1)+nc[0];
        begin[nruns]=view.starts[row+xa];
        end[nruns]=view.starts[row+xb]+view.tots[row+xb];
        shift[nruns]=(rowshift-view.lattice[0])-wrapped;
        ++nruns;
        a=0;
      }
      const int bmid=(b<nc[0]-1 ? b : nc[0]-1);
      if(a<=bmid) {
        begin[nruns]=view.starts[row+a];
        end[nruns]=view.starts[row+bmid]+view.tots[row+bmid];
        shift[nruns]=rowshift-wrapped;
        ++nruns;
      }
      if(b>=nc[0]) {
        const std::size_t xa=(a>nc[0] ? a : nc[0])-nc[0], xb=b-nc[0];
        begin[nruns]=view.starts[row+xa];
        end[nruns]=view.starts[row+xb]+view.tots[row+xb];
        shift[nruns]=(rowshift+view.lattice[0])-wrapped;
        ++nruns;
      }
      for(int r=0; r<nruns; ++r) {
        const Vector s=shift[r];
        for(std::size_t j=begin[r]; j<end[r]; ++j) {
          visit( j, view.positions[j]+s );
        }
      }
    }
  }
}

/// Calls visit(atom,d,d2) for each atom within the cutoff from the position pos, except
/// the atom self (use std::numeric_limits<unsigned>::max() to visit all), in the order
/// of the cells. atom is the index stored in the cells, d the minimum image vector from pos
/// to the atom and d2 its squared modulus. visit is a functor
template<typename Visitor>
struct LinkCellsNeighborFilter {
  const LinkCellsView& view;
  const unsigned self;
  Visitor& visit;
  inline void operator()( const std::size_t j, const Vector& d ) {
    const double d2=d.modulo2();
    const unsigned atom=view.lists[j];
    if( d2<=view.cutoff2 && atom!=self ) {
      visit( atom, d, d2 );
    }
  }
};

template<typename Visitor>
inline void forEachNeighbor( const LinkCellsView& view, const Vector& pos, const unsigned self, Visitor& visit ) {
  LinkCellsNeighborFilter<Visitor> filter{view,self,visit};
  forEachCandidate( view, pos, filter );
}

} // namespace PLMD

#endif
