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
#ifndef __PLUMED_tools_LinkCells_h
#define __PLUMED_tools_LinkCells_h

#include <limits>
#include <vector>
#include <array>
#include "View.h"
#include "Vector.h"
#include "Pbc.h"

namespace PLMD {

class Communicator;

/// \ingroup TOOLBOX
/// A class for doing link cells
class LinkCells {
public:
  ///This struct contains the configuration of cells
  struct CellCollection {
/// The start of each block corresponding to each link cell
    std::vector<unsigned> lcell_starts;
/// The number of atoms in each link cell
    std::vector<unsigned> lcell_tots;
/// The atoms ordered by link cells
    std::vector<unsigned> lcell_lists;
    ///return the sum of the number of atoms in the n cells with most atoms
    unsigned getMaximimumCombination(unsigned numCells=27) const;
    ///returns a const view of the on indexes the given cell
    inline View<const unsigned> getCellIndexes(const unsigned cellno) const {
      return {lcell_lists.data()+lcell_starts[cellno],lcell_tots[cellno]};
      // this will make the debug test fail in the case of
      // a series of empty cells in the end of the array,
      // because lcell_starts[cellno] will be==cell_lists.size()
      // the option abuve have the same problem, but it is masked
      // with pointer arithmetic
      // and should be ignored by the size passed being 0
      //return {&lcell_lists[lcell_starts[cellno]],lcell_tots[cellno]};
    }
  };
private:
/// Symbolic link to plumed communicator
  Communicator & comm;
/// Check that the link cells were set up correctly
  bool cutoffwasset=false;
/// Are there periodic boundary conditions setup
  bool nopbc=false;
/// The cutoff to use for the sizes of the cells
  double link_cutoff=0.0;
/// The pbc we are using for link cells
  Pbc mypbc;
/// The location of the origin if we are not using periodic boundary conditions
  Vector origin;
/// The number of cells in each direction
  std::array<unsigned,3> ncells{0,0,0};
/// The number of cells to stride through to get the link cells
  std::array<unsigned,3> nstride{1,0,0};
/// Work vector with the list of cells each atom is inside
  std::vector<unsigned> allcells;
///The collection of indexes per cell created by buildCellLists
  CellCollection innerCollection;
  void createCells(const PLMD::Tensor&);
public:
///
  explicit LinkCells( Communicator& comm );
/// Have the link cells been enabled
  bool enabled() const ;
/// Set the value of the cutoff
  void setCutoff( double lcut );
/// Get the value of the cutoff
  double getCutoff() const ;
/// Get the total number of link cells
  unsigned getNumberOfCells() const ;
///Get the number of cells per dimension
  const std::array<unsigned,3>& getCellLimits() const ;
/// Get the number of atoms in the cell that contains the most atoms
  unsigned getMaxInCell() const ;
/// Sets up the cells ignoring the pbcs
///
/// This creates an orthogonal box that encloses all the atoms
  void setupCells( View<const Vector> pos);
  template <typename T>
  inline void setupCells( const T& pos ) {
    setupCells(make_const_view(pos));
  }
/// Sets up the cells using the box stored in the Pbc object
  void setupCells( const Pbc& pbc );
///setups the cells
///
/// Usually it only uses the pbcs
/// if the pbc box is  null (three zero vectors) it builds an orthorombic box from the given positions
  void setupCells( View<const Vector> pos,
                   const Pbc& pbc );
  template <typename T>
  inline void setupCells( const T& pos,
                          const Pbc& pbc ) {
    setupCells(make_const_view(pos),pbc);
  }
/// Build the link cell lists
  void buildCellLists( View<const Vector> pos,
                       View<const unsigned> indices,
                       const Pbc& pbc );
  template <typename T,typename Y>
  inline void buildCellLists( const T& pos,
                              const Y& indices,
                              const Pbc& pbc ) {
    buildCellLists(make_const_view(pos), make_const_view(indices),pbc);
  }
// Reset the passed collection on the positions and the passed indexes, needs setupCells to be launched in advance
  void resetCollection( CellCollection& collection,
                        View<const Vector> pos,
                        View<const unsigned> indices );
// Returns a collection on the positions and the passed indexes, needs setupCells to be launched in advance
  CellCollection getCollection(View<const Vector> pos,
                               View<const unsigned> indices );
  template <typename T,typename Y>
  inline CellCollection getCollection(const T& pos,
                                      const Y& indices) {
    return getCollection(make_const_view(pos),
                         make_const_view(indices));
  }
/// Take three indices and return the index of the corresponding cell
  unsigned convertIndicesToIndex( std::array<unsigned,3> cellCoord) const;
/// Find the cell index in which this position is contained
  unsigned findCell( const Vector& pos ) const ;
/// Find the cell in which this position is contained
  std::array<unsigned,3> findMyCell( Vector pos ) const ;
/// Returns the "coordinate" of the cell index, assumes that cellIndex < getNumberOfCells()
  std::array<unsigned,3> findMyCell( unsigned cellIndex ) const ;
/// Get the list of cells we need to surround the a particular cell
///
/// usePbc toggles the use of pbcs for the input cell (celln)
  void addRequiredCells( const std::array<unsigned,3>& celn,
                         unsigned& ncells_required,
                         std::vector<unsigned>& cells_required,
                         bool usePbc=true) const ;
/// Retrieve the atoms in a list of cells
  void retrieveAtomsInCells( unsigned ncells_required,
                             View<const unsigned> cells_required,
                             unsigned& natomsper,
                             std::vector<unsigned>& atoms,
                             unsigned avoidIndex=std::numeric_limits<unsigned>::max()) const ;
/// Retrieve the atoms we need to consider
  void retrieveNeighboringAtoms( const Vector& pos,
                                 std::vector<unsigned>& cell_list,
                                 unsigned& natomsper,
                                 std::vector<unsigned>& atoms ) const ;
/// Create a neighbour list for the specified input atoms
///
/// @param nat number of atoms
/// @param pos position of the "central atoms"
/// @param ind indexes of the central atoms
/// @param tind temporary indexes
/// @param neigh_pos position of the candidate neighbors
/// @param neigh_ind indexes of the candidate neighbors
/// @param pbc
/// @param natoms_per_list *output* number of atom per list (mapped on tind) `nlist[tind[i]]=natoms[i]`
/// @param nlist *output* list of the central atoms plus their neigbors
  void createNeighborList( unsigned nat,
                           View<const Vector> pos,
                           View<const unsigned> ind,
                           View<const unsigned> tind,
                           View<const Vector> neigh_pos,
                           View<const unsigned> neigh_ind,
                           const Pbc& pbc,
                           unsigned& natoms_per_list,
                           std::vector<std::size_t>& nlist ) ;
};

inline
bool LinkCells::enabled() const {
  return cutoffwasset;
}

inline
unsigned LinkCells::getNumberOfCells() const {
  return ncells[0]*ncells[1]*ncells[2];
}

inline
const std::array<unsigned,3>& LinkCells::getCellLimits() const {
  return ncells;
}
}

#endif
