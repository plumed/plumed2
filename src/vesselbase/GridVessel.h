/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_vesselbase_GridVessel_h
#define __PLUMED_vesselbase_GridVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"

namespace PLMD {
namespace vesselbase {

class GridVessel : public Vessel {
private:
/// These two variables are used to 
/// remember the box we were in on the previous call
 unsigned bold;
/// The number of points in the grid
 unsigned npoints;
/// Remember the neighbors that were used last time
 std::vector<unsigned> current_neigh; 
/// The names of the various columns in the grid file
 std::vector<std::string> arg_names;
/// The minimum and maximum in the grid stored as strings 
 std::vector<std::string> str_min, str_max;
/// The minimum and maximum of the grid stored as doubles
 std::vector<double> min, max;
/// The spacing between grid points
 std::vector<double> dx;
/// The numerical distance between adjacent grid points
 std::vector<unsigned> stride;
/// The number of bins in each grid direction
 std::vector<unsigned> nbin;
/// Is this direction periodic
 std::vector<bool> pbc;
/// A tempory array that can be used to store indices for a point
 std::vector<unsigned> tmp_indices;
/// The grid with all the data values on it
 std::vector<double> data;
///  Flatten the grid and get the grid index for a point
 unsigned getIndex( const std::vector<unsigned>& indices ) const ;
/// The number of pieces of information we are storing for each 
/// point in the grid
 unsigned nper;
/// The dimensionality of the grid
 unsigned dimension;
/// The grid point that was requested last by getGridPointCoordinates
 unsigned currentGridPoint;
public:
/// keywords
  static void registerKeywords( Keywords& keys );
/// Constructor
  GridVessel( const VesselOptions& );
/// Get a description of the grid to output to the log
 std::string getGridDescription() const ;
/// Get the indices fof a point
 void getIndices( const unsigned& index, std::vector<unsigned>& indices ) const ;

/// Operations on one of the elements of grid point i
 void setGridElement( const unsigned&, const unsigned&, const double& );
 void addToGridElement( const unsigned&, const unsigned&, const double& );

/// Operations on one of the elements of grid point specified by vector
 double getGridElement( const std::vector<unsigned>&, const unsigned& ) const ;
 void setGridElement( const std::vector<unsigned>&, const unsigned&, const double& );
 void addToGridElement( const std::vector<unsigned>&, const unsigned&, const double& );
/// Set the size of the buffer equal to nper*npoints
  virtual void resize();
/// Get the number of points in the grid
  unsigned getNumberOfPoints() const;
/// Get the coordinates for a point in the grid
  void getGridPointCoordinates( const unsigned& , std::vector<double>& );
/// Get the dimensionality of the function
  unsigned getDimension() const ;
/// Get the number of grid points for each dimension
  std::vector<unsigned> getNbin() const ;
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Return the volume of one of the grid cells
  double getCellVolume() const ;
/// Get the value of the ith grid element 
  double getGridElement( const unsigned&, const unsigned& ) const ;
/// Get the numerical index for the box that contains a particular point
 unsigned getLocationOnGrid( const std::vector<double>& x, std::vector<double>& dd );
/// Get the points neighboring a particular spline point
 void getSplineNeighbors( const unsigned& mybox, std::vector<unsigned>& mysneigh );
/// Get the spacing between grid points
  const std::vector<double>& getGridSpacing() const ;
};

inline
unsigned GridVessel::getNumberOfPoints() const {
  return npoints;
}

inline
const std::vector<double>& GridVessel::getGridSpacing() const {
  return dx;
}

inline
double GridVessel::getCellVolume() const {
  double myvol=1.0; for(unsigned i=0;i<dimension;++i) myvol *= dx[i];
  return myvol;
}

inline
unsigned GridVessel::getDimension() const {
  return dimension;
}

}
}
#endif
