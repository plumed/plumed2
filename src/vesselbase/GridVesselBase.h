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
#ifndef __PLUMED_vesselbase_GridVesselBase_h
#define __PLUMED_vesselbase_GridVesselBase_h

#include <string>
#include <cstring>
#include <vector>
#include "Vessel.h"

namespace PLMD {
namespace vesselbase {

class GridVesselBase : public Vessel {
friend class InterpolationBase;
friend class PrintGrid;
friend class ReduceGridDimension;
friend class FunctionOnGrid;
private:
/// Are we interpolating this grid
 bool interpolating;
/// These two variables are used to 
/// remember the box we were in on the previous call
 unsigned bold;
/// The number of points in the grid
 unsigned npoints;
/// Write the contents of the grid to the checkpoint file
  bool checkpoint;
/// This is if we are interpolating function to ensure that 
/// tables are updated if data is changed
 bool dataHasChangedSinceInterpol;
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
///  Flatten the grid and get the grid index for a point
 unsigned getIndex( const std::vector<unsigned>& indices ) const ;
protected:
/// The number of pieces of information we are storing for each 
/// point in the grid
 unsigned nper;
/// The dimensionality of the grid
 unsigned dimension;
/// The grid point that was requested last by getGridPointCoordinates
 unsigned currentGridPoint;
/// Get a description of the grid to output to the log
 std::string getGridDescription() const ;
/// Get the indices fof a point
 void getIndices( const unsigned& index, std::vector<unsigned>& indices ) const ;
/// Get the indices at which a particular point resides
// void getIndices(const std::vector<double>& x, std::vector<unsigned>& indices) const ; 

/// Operations on one of the elements of grid point i
 void setGridElement( const unsigned&, const unsigned&, const double& );
 void addToGridElement( const unsigned&, const unsigned&, const double& );

/// Operations on one of the elements of grid point specified by vector
 double getGridElement( const std::vector<unsigned>&, const unsigned& ) const ;
 void setGridElement( const std::vector<unsigned>&, const unsigned&, const double& );
 void addToGridElement( const std::vector<unsigned>&, const unsigned&, const double& );
/// Finish setup of grid object when it is storing the values of a more abstract quantity
  void finishSetup( const unsigned& nelem, const std::vector<bool>& ipbc, const std::vector<std::string>& names );
public:
/// keywords
  static void registerKeywords( Keywords& keys );
/// Constructor
  GridVesselBase( const VesselOptions& );
/// Set the size of the buffer equal to nper*npoints
  virtual void resize();
/// Get the number of points in the grid
  unsigned getNumberOfPoints() const;
/// Get the coordinates for a point in the grid
  void getGridPointCoordinates( const unsigned& , std::vector<double>& );
/// Return a description of the quantity stored in a particular column of the grid
 std::string getQuantityDescription(const unsigned& ) const;
/// Write the grid on a file
  void writeToFile( OFile& , const std::string& fmt );
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
/// Use grid checkpointing
  void storeInCheckpoint();
/// Write data to the checkpoint file
  void writeToCheckpoint( OFile& cfile );
/// Read data from checkpoint file
  void readFromCheckpoint( IFile& cifile );
/// Get the numerical index for the box that contains a particular point
 unsigned getLocationOnGrid( const std::vector<double>& x, std::vector<double>& dd );
/// Get the points neighboring a particular spline point
 void getSplineNeighbors( const unsigned& mybox, std::vector<unsigned>& mysneigh );
/// Calculate the vector from the grid point to point x then normalize by grid spacing
/// This is useful for interpolation
 void getFractionFromGridPoint( const unsigned& igrid, const std::vector<double>& x, std::vector<double>& dd );
/// Will the grid be interpolated
  bool gridWillBeInterpolated() const;
};

inline
unsigned GridVesselBase::getNumberOfPoints() const {
  return npoints;
}

inline
double GridVesselBase::getCellVolume() const {
  double myvol=1.0; for(unsigned i=0;i<dimension;++i) myvol *= dx[i];
  return myvol;
}

inline
unsigned GridVesselBase::getDimension() const {
  return dimension;
}

inline
bool GridVesselBase::gridWillBeInterpolated() const {
  return interpolating;
}

}
}
#endif
