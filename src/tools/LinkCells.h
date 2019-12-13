/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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

#include <vector>
#include "Vector.h"
#include "Pbc.h"

namespace PLMD {

class Communicator;

/// \ingroup TOOLBOX
/// A class for doing link cells
class LinkCells {
private:
/// Symbolic link to plumed communicator
  Communicator & comm;
/// Check that the link cells were set up correctly
  bool cutoffwasset;
/// The cutoff to use for the sizes of the cells
  double link_cutoff;
/// The pbc we are using for link cells
  Pbc mypbc;
/// The number of cells in each direction
  std::vector<unsigned> ncells;
/// The number of cells to stride through to get the link cells
  std::vector<unsigned> nstride;
/// The list of cells each atom is inside
  std::vector<unsigned> allcells;
/// The start of each block corresponding to each link cell
  std::vector<unsigned> lcell_starts;
/// The number of atoms in each link cell
  std::vector<unsigned> lcell_tots;
/// The atoms ordered by link cells
  std::vector<unsigned> lcell_lists;
public:
///
  explicit LinkCells( Communicator& comm );
/// Have the link cells been enabled
  bool enabled() const ;
/// Set the value of the cutoff
  void setCutoff( const double& lcut );
/// Get the value of the cutoff
  double getCutoff() const ;
/// Get the total number of link cells
  unsigned getNumberOfCells() const ;
/// Build the link cell lists
  void buildCellLists( const std::vector<Vector>& pos, const std::vector<unsigned>& indices, const Pbc& pbc );
/// Take three indices and return the index of the corresponding cell
  unsigned convertIndicesToIndex( const unsigned& nx, const unsigned& ny, const unsigned& nz ) const ;
/// Find the cell index in which this position is contained
  unsigned findCell( const Vector& pos ) const ;
/// Find the cell in which this position is contained
  std::array<unsigned,3> findMyCell( const Vector& pos ) const ;
/// Get the list of cells we need to surround the a particular cell
  void addRequiredCells( const std::array<unsigned,3>& celn, unsigned& ncells_required,
                         std::vector<unsigned>& cells_required ) const ;
/// Retrieve the atoms in a list of cells
  void retrieveAtomsInCells( const unsigned& ncells_required,
                             const std::vector<unsigned>& cells_required,
                             unsigned& natomsper, std::vector<unsigned>& atoms ) const ;
/// Retrieve the atoms we need to consider
  void retrieveNeighboringAtoms( const Vector& pos, std::vector<unsigned>& cell_list, unsigned& natomsper, std::vector<unsigned>& atoms ) const ;
};

inline
bool LinkCells::enabled() const {
  return cutoffwasset;
}

inline
unsigned LinkCells::getNumberOfCells() const {
  return ncells[0]*ncells[1]*ncells[2];
}

}

#endif
