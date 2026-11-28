/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_tools_GridCoordinatesObject_h
#define __PLUMED_tools_GridCoordinatesObject_h

#include <string>
#include <cstring>
#include <vector>
#include "Exception.h"
#include "Tools.h"
#include "View.h"

namespace PLMD {

class CheckpointReader;
class CheckpointWriter;

class GridCoordinatesObject {
private:
/// Have the bounds been setup on the grid
  bool bounds_set;
/// The way that grid points are constructed
  enum {flat,fibonacci} gtype;
/// The number of points in the grid
  unsigned npoints;
/// Stuff for fibonacci grids
  double root5, golden, igolden, log_golden2;
/// Fib increment here is equal to 2*pi*(INVERSE GOLDEN RATIO)
  double fib_offset, fib_increment, fib_shift;
  std::vector<std::vector<unsigned> > fib_nlist;
/// The minimum and maximum of the grid stored as doubles
  std::vector<double> min, max;
/// The numerical distance between adjacent grid points
  std::vector<unsigned> stride;
/// The number of bins in each grid direction
  std::vector<unsigned> nbin;
/// Is this direction periodic
  std::vector<bool> pbc;
/// The minimum and maximum in the grid stored as strings
  std::vector<std::string> str_min, str_max;
/// The spacing between grid points
  std::vector<double> dx;
/// The dimensionality of the grid
  unsigned dimension;
/// The minimum when the lattice was set (setBounds, setLattice): extend does not change it
  std::vector<double> lattice_start;
/// Index on the lattice of the first point of the grid (zero when the lattice is set, changed by extend)
  std::vector<long long> lattice_origin;
/// Set the strides and the number of points from nbin
  void setStrides();
/// Get the index of the closest point on the fibonacci sphere
  unsigned getFibonacciIndex( View<const double> p ) const ;
/// Get the flat grid coordinates
  void getFlatGridCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const ;
/// Get the coordinates on the Fibonacci grid
  void getFibonacciCoordinates( const unsigned& ipoint, std::vector<double>& x ) const ;
public:
/// Setup the grid
  void setup( const std::string& geom, const std::vector<bool>& ipbc, const unsigned& np, const double& fib_cutoff );
/// Set the minimum and maximum of the grid
  void setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax, const std::vector<std::size_t>& nbins, std::vector<double>& spacing );
/// Set a flat grid from the minimum, the spacing and the number of bins along each direction, as doubles (no
/// parsing). The bounds (getMin, getMax) are smin and smax if given (e.g. -pi and pi of a periodic argument),
/// otherwise they are written with Tools::boundString; as in setBounds, the maximum of a direction that is not
/// periodic is the last point.
  void setLattice( const std::vector<double>& gmin, const std::vector<double>& spacing, const std::vector<std::size_t>& nbins,
                   const std::vector<std::string>& smin= {}, const std::vector<std::string>& smax= {} );
/// Change the range along the directions that are not periodic: lower[i] (upper[i]) points are added at the lower
/// (upper) end of direction i if positive, removed if negative. The spacing does not change and the points stay on
/// the same lattice. Same arithmetic and bounds as GridBase::resize.
  void extend( const std::vector<long long>& lower, const std::vector<long long>& upper );
/// Are this grid and other on the same lattice (same spacing and periodicity, and the same lattice set by
/// setBounds or setLattice, then possibly extended)?
  bool onSameLattice( const GridCoordinatesObject& other ) const ;
/// For two grids on the same lattice: point k of other along direction i is point k+shift[i] of this grid
  std::vector<long long> getShiftFrom( const GridCoordinatesObject& other ) const ;
/// Store the flat grid in a checkpoint (names starting with prefix)
  void writeCheckpoint( CheckpointWriter& w, const std::string& prefix ) const ;
/// Restore a flat grid stored by writeCheckpoint (the grid has been set up with the same periodicity)
  void readCheckpoint( const CheckpointReader& r, const std::string& prefix );
/// Convert an index into indices
  void convertIndexToIndices( const unsigned& index, const std::vector<unsigned>& nnbin, std::vector<unsigned>& indices ) const ;
/// Check if a point is within the grid boundaries
  bool inbounds( const std::vector<double>& point ) const ;
  bool inbounds( View<const double> point ) const ;
/// Convert a point in space the correspoinding grid point
  unsigned getIndex( const std::vector<double>& p ) const ;
  unsigned getIndex( View<const double> point ) const ;
///  Flatten the grid and get the grid index for a point
  unsigned getIndex( const std::vector<unsigned>& indices ) const ;
/// Get the indices fof a point
  void getIndices( const unsigned& index, std::vector<unsigned>& indices ) const ;
/// Get the indices of a particular point
  void getIndices( const std::vector<double>& point, std::vector<unsigned>& indices ) const ;
  void getIndices( View<const double> point, std::vector<unsigned>& indices ) const ;
/// Get the number of points in the grid
  unsigned getNumberOfPoints() const;
/// Get the coordinates for a point in the grid
  void getGridPointCoordinates( const unsigned&, std::vector<double>& ) const ;
  void getGridPointCoordinates( const unsigned&, std::vector<unsigned>&, std::vector<double>& ) const ;
/// Create a coordinate that has this value of the grid
  void putCoordinateAtValue( const unsigned&, const double&, std::vector<double>& ) const ;
/// Get the dimensionality of the function
  unsigned getDimension() const ;
/// Is the grid periodic in the ith direction
  bool isPeriodic( const unsigned& i ) const ;
/// Get the number of grid points for each dimension
  std::vector<std::size_t> getNbin( const bool& shape ) const ;
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// The minimum of the grid in each dimension as a double
  const std::vector<double>& getMinValues() const ;
/// The coordinate of x on the lattice of direction i, in units of the spacing from the first point
  double getLatticeCoordinate( const unsigned& i, const double& x ) const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Return the volume of one of the grid cells
  double getCellVolume() const ;
/// Get the set of points neighouring a particular location in space
  void getNeighbors( const std::vector<double>& pp, const std::vector<unsigned>& nneigh,
                     unsigned& num_neighbours, std::vector<unsigned>& neighbors ) const ;
  void getNeighbors( const View<const double> pp, const std::vector<unsigned>& nneigh,
                     unsigned& num_neighbours, std::vector<unsigned>& neighbors ) const ;
/// Get the neighbors for a set of indices of a point
  void getNeighbors( const std::vector<unsigned>& indices, const std::vector<unsigned>& nneigh,
                     unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const ;
/// Get the points neighboring a particular spline point
  void getSplineNeighbors( const unsigned& mybox, unsigned& nneighbors, std::vector<unsigned>& mysneigh ) const ;
/// Get the spacing between grid points
  const std::vector<double>& getGridSpacing() const ;
/// Get the stride (the distance between the grid points of an index)
  const std::vector<unsigned>& getStride() const ;
/// Get the type of the grid
  std::string getGridType() const ;
};

inline
unsigned GridCoordinatesObject::getNumberOfPoints() const {
  return npoints;
}

inline
const std::vector<double>& GridCoordinatesObject::getGridSpacing() const {
  if( gtype==flat ) {
    return dx;
  }
  plumed_merror("dont understand what spacing means for spherical grids");
  return dx;
}

inline
const std::vector<double>& GridCoordinatesObject::getMinValues() const {
  plumed_dbg_assert( gtype==flat );
  return min;
}

inline
double GridCoordinatesObject::getLatticeCoordinate( const unsigned& i, const double& x ) const {
  plumed_dbg_assert( gtype==flat );
  return (x-min[i])/dx[i];
}

inline
double GridCoordinatesObject::getCellVolume() const {
  if( gtype==flat ) {
    double myvol=1.0;
    for(unsigned i=0; i<dimension; ++i) {
      myvol *= dx[i];
    }
    return myvol;
  } else {
    return 4*pi / static_cast<double>( getNumberOfPoints() );
  }
}

inline
unsigned GridCoordinatesObject::getDimension() const {
  return dimension;
}

inline
bool GridCoordinatesObject::isPeriodic( const unsigned& i ) const {
  plumed_dbg_assert( gtype==flat );
  return pbc[i];
}

inline
const std::vector<unsigned>& GridCoordinatesObject::getStride() const {
  plumed_dbg_assert( gtype==flat );
  return stride;
}

inline
std::string GridCoordinatesObject::getGridType() const {
  if( gtype==flat ) {
    return "flat";
  } else if( gtype==fibonacci ) {
    return "fibonacci";
  }
  return "";
}

}
#endif
