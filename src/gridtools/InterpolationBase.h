/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_tools_InterpolationBase_h
#define __PLUMED_tools_InterpolationBase_h

#include <vector>
#include "GridVesselBase.h"

namespace PLMD {
namespace vesselbase {

/// Abstract base class for interpolation
class InterpolationBase {
private:
/// Tells us if we can use tables for interpolation
  bool tablesAreSet;
/// The start point for the data in the grid
  unsigned data;
/// The dimensionality of the interpolated function
  unsigned dimension;
/// The grid underlying the interpolation class
  GridVesselBase* mygrid;
/// A holder for the vector telling us where in the cell we are
  std::vector<double> ddx;
protected:
/// Get the dimension of the grid
  unsigned getDimension() const ;
/// Get number of points in the underlying grid
  unsigned getNumberOfSplinePoints() const ;
/// Get the spacings between adjacent grid points
  void getGridPointSpacing( std::vector<double>& spacing ) const ;
/// Get number of grid points along each axis
  void getNumberOfGridPoints( std::vector<unsigned>& bin ) const ;
/// Get the index for a particular numbered cell
  unsigned getBoxIndex( const std::vector<unsigned>& pn  ) const ;
/// Return the number of points between adjacent grid points
  unsigned getGridStride( const unsigned& i ) const ;
/// Get the value on the grid at a particular point labelled using numerical sequence
  double getValue( const unsigned& mybox ); 
/// Get the value at a particular grid point
  double getValue( const std::vector<unsigned>& pn );
/// Get the value and derivatives at a particular grid point
  double getValueAndDerivatives( const std::vector<unsigned>& boxno, std::vector<double>& der );
public:
  InterpolationBase( GridVesselBase* gg, const unsigned dstart=0 );
/// Interpolation tables were setup
  bool interpolationTablesWereSet() const ;
/// Setup interpolation tables
  void set_table();
/// Finish setup of interpolation functions
  virtual void setInterpolationTables()=0;
/// Get the value of the function using the interpolation tables
  double getFunctionValue( const std::vector<double>& pos );
/// This does an interpolation inside the box numbered mybox
  virtual double interpolateFunction( const unsigned& mybox, const std::vector<double>& dd )=0;
/// Eventually you could add something like this for metadynamics
//  double getValueNoTables( const std::vector<double>& pos );
/// Create an object to do cubic interpolation
  static InterpolationBase* createCubicInterpolator( GridVesselBase* gg, const unsigned dstart=0 );
};

inline
bool InterpolationBase::interpolationTablesWereSet() const {
  return ( tablesAreSet && !mygrid->dataHasChangedSinceInterpol );
}

inline
unsigned InterpolationBase::getDimension() const {
  return dimension;
}

inline
unsigned InterpolationBase::getNumberOfSplinePoints() const {
  return mygrid->npoints;
}

}
}

#endif
