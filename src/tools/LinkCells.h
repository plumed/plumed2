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
#include "LinkCellsView.h"

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
/// The positions of the atoms wrapped in the box, in the same order as lcell_lists
/// (only if they were requested when the collection was built, see resetCollection)
    std::vector<Vector> lcell_pos;
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
/// The largest value accepted by setSubdivision
  static constexpr unsigned maxSubdivision=LinkCellsView::maxSubdivision;
/// Minimum number of atoms per thread when searching the neighbors of the atoms
/// (createNeighborList, NeighborList::update): a search costs about a microsecond
  static constexpr unsigned minAtomsPerThreadNeighborList=16;
private:
/// Minimum number of atoms per thread when building the lists of the atoms in the cells
/// (see resetCollection): finding the cell of an atom costs a few nanoseconds
  static constexpr unsigned minAtomsPerThread=64;
/// Minimum number of atoms per range of cells (see resetCollection): each range goes through all the atoms
  static constexpr unsigned minAtomsPerRange=1024;
/// Maximum total number of cells
  static constexpr double maxNumberOfCells=1048576;
/// With setLimitCellsByAtoms(true) the cells are at most this many times the atoms in them
/// (but at least minCellsForAtomLimit): building the lists costs time proportional to the
/// number of cells, and more cells than atoms only add empty cells
  static constexpr double maxCellsPerAtom=2;
  static constexpr double minCellsForAtomLimit=64;
/// Set with setLimitCellsByAtoms
  bool limitCellsByAtoms=false;
/// If not zero, the number of atoms that limits the number of cells in createCells
  std::size_t cellLimitAtoms=0;
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
/// The number of cells per cutoff that was asked for with setSubdivision
  unsigned requestedSubdivision=1;
/// The number of cells per cutoff in use: the cells within nsub cells from a cell, in each
/// direction, contain all the atoms within the cutoff from its atoms. It is the requested one
/// if there are enough cells (2*nsub+1 in each direction), otherwise 1
  unsigned nsub=1;
/// The width of the cells: the distance between their opposite faces
  std::array<double,3> cellWidth{0.0,0.0,0.0};
/// Two points whose distances along the normals to the three families of faces of the
/// cells are (at least) p_k are farther than sqrt(sum_k p_k^2/skewFactor). It is the
/// largest eigenvalue of the matrix of the cosines of the angles between the normals:
/// 1 for orthorhombic boxes, and at most 3
  double skewFactor=1.0;
/// The lattice vectors of the box (zero without periodic boundary conditions)
  std::array<Vector,3> lattice;
/// The number of cells in each direction
  std::array<unsigned,3> ncells{0,0,0};
/// The number of cells to stride through to get the link cells
  std::array<unsigned,3> nstride{1,0,0};
/// Work vector with the list of cells each atom is inside
  std::vector<unsigned> allcells;
/// Work vector with the positions of the atoms wrapped in the box
  std::vector<Vector> allwrapped;
/// Work vector with the number of atoms in each range of cells (see resetCollection)
  std::vector<unsigned> rangetots;
///The collection of indexes per cell created by buildCellLists
  CellCollection innerCollection;
  void createCells(const PLMD::Tensor& box);
/// Find the index of the cell of a position, without checks:
/// returns getNumberOfCells() if the position is not in a cell
  unsigned findCellUnchecked( Vector pos ) const noexcept;
public:
/// The communicator is not used: the cells are built locally on every process
/// (the argument is kept for backward compatibility)
  explicit LinkCells( Communicator& comm );
/// The cells are built locally: no communicator is needed
  LinkCells() = default;
/// Have the link cells been enabled
  bool enabled() const ;
/// Set the number of cells per cutoff (1 or 2, the default is 1) used by the next setupCells.
/// With 2 the cells are smaller and those within two cells are searched: fewer atoms
/// beyond the cutoff are visited. It is used only if there are at least five cells in
/// each direction. The functions that work on the neighboring cells (addRequiredCells,
/// retrieveNeighboringAtoms, and the searches on getView) take it into account
  void setSubdivision( unsigned n );
/// If true, the next setupCells with positions uses at most max(64, 2*natoms) cells
/// (larger cells, still correct). Useful when the lists are rebuilt often and the atoms
/// are few for the size of the box; createNeighborList always does it
  void setLimitCellsByAtoms( bool limit );
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
/// if the pbc box is  null (three zero vectors) it builds an orthorhombic box from the given positions
  void setupCells( View<const Vector> pos,
                   const Pbc& pbc );
  template <typename T>
  inline void setupCells( const T& pos,
                          const Pbc& pbc ) {
    setupCells(make_const_view(pos),pbc);
  }
/// Build the link cell lists (withPositions: see resetCollection)
  void buildCellLists( View<const Vector> pos,
                       View<const unsigned> indices,
                       const Pbc& pbc,
                       bool withPositions=false );
  template <typename T,typename Y>
  inline void buildCellLists( const T& pos,
                              const Y& indices,
                              const Pbc& pbc,
                              bool withPositions=false ) {
    buildCellLists(make_const_view(pos), make_const_view(indices),pbc,withPositions);
  }
/// Reset the passed collection on the positions and the passed indexes, needs setupCells to be launched in advance.
/// With withPositions the wrapped positions are also stored (lcell_pos), as needed by getView
  void resetCollection( CellCollection& collection,
                        View<const Vector> pos,
                        View<const unsigned> indices,
                        bool withPositions=false );
/// Returns a collection on the positions and the passed indexes, needs setupCells to be launched in advance
  CellCollection getCollection(View<const Vector> pos,
                               View<const unsigned> indices,
                               bool withPositions=false );
  template <typename T,typename Y>
  inline CellCollection getCollection(const T& pos,
                                      const Y& indices,
                                      bool withPositions=false) {
    return getCollection(make_const_view(pos),
                         make_const_view(indices),
                         withPositions);
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
/// Returns true if the periodic image of an atom that is within the cutoff from an atom
/// of another cell is fixed by the cells of the two atoms. This is the case when the cells
/// to be searched in each direction are all different (the box is at least three cutoffs
/// wide), and without periodic boundary conditions. The minimum image can then be found
/// without Pbc::distance: see getView
  bool cellsFixImages() const ;
/// The position moved inside the box, as it is done to find the cell of an atom
/// (the position itself without periodic boundary conditions)
  Vector wrapPosition( const Vector& pos ) const ;
/// The cells and a collection built with its positions, to be searched with forEachCandidate
/// or forEachNeighbor (also on a GPU). The cells farther than the cutoff from a position are
/// not searched, and the vectors to the candidates are found without Pbc::distance: it can
/// be used only if cellsFixImages() is true. wrapCells is the usePbc of addRequiredCells
  LinkCellsView getView( const CellCollection& collection, bool wrapCells=true ) const ;
/// The view of the collection built by buildCellLists
  LinkCellsView getView( bool wrapCells=true ) const ;
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
/// @param natoms_per_list *output* the number of neighbors in the longest list: each list has
///        1+natoms_per_list slots in nlist (the first one is for the central atom)
/// @param nlist *output* the number of atoms in each list, central atom included
///        (`nlist[tind[i]]`), followed by the lists: the one of atom i starts at
///        `nat+tind[i]*(1+natoms_per_list)`
/// @param applyCutoff if true only the candidates whose (minimum image) distance from the
///        central atom is not larger than the link cell cutoff are kept, and cells of half
///        the cutoff are used when the box is large enough (fewer candidates are visited).
///        Without it the cells are at least as wide as the cutoff and the lists hold all
///        the atoms of the 27 cells around the central atom, up to ~2-3 cutoffs away: with
///        a short cutoff in a dense system most of them are beyond it. Use it only when
///        the quantity computed from the list is exactly zero beyond the cutoff.
  void createNeighborList( unsigned nat,
                           View<const Vector> pos,
                           View<const unsigned> ind,
                           View<const unsigned> tind,
                           View<const Vector> neigh_pos,
                           View<const unsigned> neigh_ind,
                           const Pbc& pbc,
                           unsigned& natoms_per_list,
                           std::vector<std::size_t>& nlist,
                           const bool applyCutoff=false ) ;
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
