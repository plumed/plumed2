/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#ifndef __PLUMED_gridtools_GridVessel_h
#define __PLUMED_gridtools_GridVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "vesselbase/AveragingVessel.h"

namespace PLMD {
namespace gridtools {

class GridVessel : public vesselbase::AveragingVessel {
  friend class ActionWithInputGrid;
  friend class DumpGrid;
private:
/// The way that grid points are constructed
  enum {flat,fibonacci} gtype;
/// Have the minimum and maximum for the grid been set
  bool bounds_set;
/// The number of points in the grid
  unsigned npoints;
/// Stuff for fibonacci grids
  double root5, golden, igolden, log_golden2;
/// Fib increment here is equal to 2*pi*(INVERSE GOLDEN RATIO)
  double fib_offset, fib_increment, fib_shift;
  std::vector<std::vector<unsigned> > fib_nlist;
/// Units for Gaussian Cube file
  double cube_units;
/// This flag is used to check if the user has created a valid input
  bool foundprint;
/// The minimum and maximum of the grid stored as doubles
  std::vector<double> min, max;
/// The numerical distance between adjacent grid points
  std::vector<unsigned> stride;
/// The number of bins in each grid direction
  std::vector<unsigned> nbin;
/// The grid point that was requested last by getGridPointCoordinates
  unsigned currentGridPoint;
/// The forces that will be output at the end of the calculation
  std::vector<double> finalForces;
protected:
/// Is forced
  bool wasforced;
/// Forces acting on grid elements
  std::vector<double> forces;
/// Do we have derivatives
  bool noderiv;
/// The names of the various columns in the grid file
  std::vector<std::string> arg_names;
/// The number of pieces of information we are storing for each
/// point in the grid
  unsigned nper;
/// Is this direction periodic
  std::vector<bool> pbc;
/// The minimum and maximum in the grid stored as strings
  std::vector<std::string> str_min, str_max;
/// The spacing between grid points
  std::vector<double> dx;
/// The dimensionality of the grid
  unsigned dimension;
/// Which grid points are we actively accumulating
  std::vector<bool> active;
/// Convert a point in space the the correspoinding grid point
  unsigned getIndex( const std::vector<double>& p ) const ;
/// Get the index of the closest point on the fibonacci sphere
  unsigned getFibonacciIndex( const std::vector<double>& p ) const ;
/// Get the flat grid coordinates
  void getFlatGridCoordinates( const unsigned& ipoint, std::vector<unsigned>& tindices, std::vector<double>& x ) const ;
/// Get the coordinates on the Fibonacci grid
  void getFibonacciCoordinates( const unsigned& ipoint, std::vector<double>& x ) const ;
public:
/// keywords
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit GridVessel( const vesselbase::VesselOptions& );
/// Remove the derivatives
  void setNoDerivatives();
/// Get the type of grid we are using
  std::string getType() const ;
/// Set the minimum and maximum of the grid
  virtual void setBounds( const std::vector<std::string>& smin, const std::vector<std::string>& smax, const std::vector<unsigned>& nbins, const std::vector<double>& spacing );
/// Get the cutoff to use for the Fibonacci spheres
  virtual double getFibonacciCutoff() const ;
/// Setup the grid if it is a fibonacci grid on the surface of a sphere
  void setupFibonacciGrid( const unsigned& np );
/// Get a description of the grid to output to the log
  std::string description();
/// Convert an index into indices
  void convertIndexToIndices( const unsigned& index, const std::vector<unsigned>& nnbin, std::vector<unsigned>& indices ) const ;
///  Flatten the grid and get the grid index for a point
  unsigned getIndex( const std::vector<unsigned>& indices ) const ;
/// Get the indices fof a point
  void getIndices( const unsigned& index, std::vector<unsigned>& indices ) const ;
/// Get the indices of a particular point
  void getIndices( const std::vector<double>& point, std::vector<unsigned>& indices ) const ;
/// Operations on one of the elements of grid point i
  void setGridElement( const unsigned&, const unsigned&, const double& );
/// Add data to an element of the grid
  void addToGridElement( const unsigned& ipoint, const unsigned& jelement, const double& value );
/// Operations on one of the elements of grid point specified by vector
  double getGridElement( const std::vector<unsigned>&, const unsigned& ) const ;
  void setGridElement( const std::vector<unsigned>&, const unsigned&, const double& );
/// Set the values and derivatives of a particular element
  void setValueAndDerivatives( const unsigned&, const unsigned&, const double&, const std::vector<double>& );
/// Set the size of the buffer equal to nper*npoints
  virtual void resize();
/// Get the number of points in the grid
  unsigned getNumberOfPoints() const;
/// Get the coordinates for a point in the grid
  void getGridPointCoordinates( const unsigned&, std::vector<double>& ) const ;
  void getGridPointCoordinates( const unsigned&, std::vector<unsigned>&, std::vector<double>& ) const ;
/// Get the dimensionality of the function
  unsigned getDimension() const ;
/// Get the number of components in the vector stored on each grid point
  virtual unsigned getNumberOfComponents() const ;
/// Is the grid periodic in the ith direction
  bool isPeriodic( const unsigned& i ) const ;
/// Get the number of quantities we have stored at each grid point
  unsigned getNumberOfQuantities() const ;
/// Get the number of grid points for each dimension
  std::vector<unsigned> getNbin() const ;
/// Get the name of the ith component
  std::string getComponentName( const unsigned& i ) const ;
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Get the number of points needed in the buffer
  virtual unsigned getNumberOfBufferPoints() const ;
/// Get the stride (the distance between the grid points of an index)
  const std::vector<unsigned>& getStride() const ;
/// Return the volume of one of the grid cells
  double getCellVolume() const ;
/// Get the value of the ith grid element
  virtual double getGridElement( const unsigned&, const unsigned& ) const ;
/// Get the set of points neighouring a particular location in space
  void getNeighbors( const std::vector<double>& pp, const std::vector<unsigned>& nneigh,
                     unsigned& num_neighbours, std::vector<unsigned>& neighbors ) const ;
/// Get the neighbors for a set of indices of a point
  void getNeighbors( const std::vector<unsigned>& indices, const std::vector<unsigned>& nneigh,
                     unsigned& num_neighbors, std::vector<unsigned>& neighbors ) const ;
/// Get the points neighboring a particular spline point
  void getSplineNeighbors( const unsigned& mybox, unsigned& nneighbors, std::vector<unsigned>& mysneigh ) const ;
/// Get the spacing between grid points
  const std::vector<double>& getGridSpacing() const ;
/// Get the extent of the grid in one of the axis
  double getGridExtent( const unsigned& i ) const ;
/// Copy data from the action into the grid
  virtual void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
/// Finish the calculation
  virtual void finish( const std::vector<double>& buffer );
/// This ensures that Gaussian cube fies are in correct units
  void setCubeUnits( const double& units );
/// This ensures that Gaussian cube files are in correct units
  double getCubeUnits() const ;
/// Return a string containing the input to the grid so we can clone it
  std::string getInputString() const ;
/// Does this have derivatives
  bool noDerivatives() const ;
/// Get the value and derivatives at a particular location using spline interpolation
  double getValueAndDerivatives( const std::vector<double>& x, const unsigned& ind, std::vector<double>& der ) const ;
/// Deactivate all the grid points
  void activateThesePoints( const std::vector<bool>& to_activate );
/// Is this point active
  bool inactive( const unsigned& ip ) const ;
/// This retrieves the final force
  virtual void getFinalForces( const std::vector<double>& buffer, std::vector<double>& finalForces ) { plumed_error(); }
/// Apply the forces
  void setForce( const std::vector<double>& inforces );
/// Was a force added to the grid
  bool wasForced() const ;
/// And retrieve the forces
  bool applyForce( std::vector<double>& fforces );
};

inline
unsigned GridVessel::getNumberOfQuantities() const {
  return nper;
}

inline
unsigned GridVessel::getNumberOfPoints() const {
  return npoints;
}

inline
const std::vector<double>& GridVessel::getGridSpacing() const {
  if( gtype==flat ) return dx;
  plumed_merror("dont understand what spacing means for spherical grids");
  return dx;
}

inline
double GridVessel::getCellVolume() const {
  if( gtype==flat ) {
    double myvol=1.0; for(unsigned i=0; i<dimension; ++i) myvol *= dx[i];
    return myvol;
  } else {
    return 4*pi / static_cast<double>( getNumberOfPoints() );
  }
}

inline
unsigned GridVessel::getDimension() const {
  return dimension;
}

inline
bool GridVessel::isPeriodic( const unsigned& i ) const {
  plumed_dbg_assert( gtype==flat );
  return pbc[i];
}

inline
std::string GridVessel::getComponentName( const unsigned& i ) const {
  return arg_names[i];
}

inline
unsigned GridVessel::getNumberOfComponents() const {
  if( noderiv ) return nper;
  return nper / ( dimension + 1 );
}

inline
double GridVessel::getGridExtent( const unsigned& i ) const {
  plumed_dbg_assert( gtype==flat );
  return max[i] - min[i];
}

inline
bool GridVessel::noDerivatives() const {
  return noderiv;
}

inline
bool GridVessel::inactive( const unsigned& ip ) const {
  plumed_dbg_assert( ip<npoints );
  return !active[ip];
}

inline
const std::vector<unsigned>& GridVessel::getStride() const {
  plumed_dbg_assert( gtype==flat );
  return stride;
}

inline
unsigned GridVessel::getNumberOfBufferPoints() const {
  return npoints;
}

inline
std::string GridVessel::getType() const {
  if( gtype==flat ) return "flat";
  else if( gtype==fibonacci ) return "fibonacci";
  plumed_error();
}

inline
double GridVessel::getFibonacciCutoff() const {
  return 0.0;
}

}
}
#endif
