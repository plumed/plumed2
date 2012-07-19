/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_Field_h
#define __PLUMED_Field_h

#include "PlumedException.h"
#include "Value.h"
#include "CubicInterpolation.h"

namespace PLMD{

class PlumedCommunicator;

// A class for storing the field
class Field {
private:
/// The type of field this is
  enum {identity,gaussian} fstyle;
/// Total number of points
  unsigned npoints;
/// Number of high dimensionality vectors
  unsigned ndX;
/// Number of low dimensionality vectors
  unsigned ndx;
/// Number of doubles for each point in the grid
  unsigned nper;
/// The value of sigma
  double sigma;
/// The sizes of all the base quantities
  std::vector<unsigned> baseq_nder;
/// The start in baseq_buffer for each function
  std::vector<unsigned> baseq_starts;
/// Storage space for the input data
  std::vector<double> baseq_buffer;
/// Storage space for the grid
  std::vector<double> grid_buffer;
/// The error message
  std::string errormsg;
/// The interpolator for the field
  CInterpolation* f_interpolator;
/// The interpolator for the derivatives
  std::vector<CInterpolation*> df_interpolators;
/// Everything to pass forces like Values 
  bool wasforced;
  std::vector<double> forces;
public:
  Field( const std::string ftype, const unsigned d );
  ~Field();
/// Write the documentation for the field
  static std::string documentation();
/// Write out the keywords for this action
  static void printKeywords( Log& log );
/// Read in the stuff for the field
  void read( const std::string& input, const unsigned nfunc, std::string& report );
/// Store an error message if something goes wrong during readin
  void error( const std::string & msg );
/// Return the value of the error message
  std::string errorMessage() const ;
/// Check if there was a problem on readin
  bool check() const ;
/// Set the number of high dimensional coordinates
  void resizeDerivatives( const unsigned D );
/// Resize the final number of derivatives
  void resizeOutDerivatives( const unsigned D );
/// Get the upper and lower boundaries of the grid
  void retrieveBoundaries( std::vector<double>& min, std::vector<double>& max );
/// Set up to recalculate everything in the field
  void clear();
/// Set the sizes of all the base quantity buffers
  void resizeBaseQuantityBuffers( const std::vector<unsigned>& cv_sizes );
/// Set the derivatives for one of the base quantities
  void setBaseQuantity( const unsigned nn, Value* val );
/// Gather the values and derivatives for all the base quantities
  void gatherBaseQuantities( PlumedCommunicator& comm );
/// Extract one of the base quantities
  void extractBaseQuantity( const unsigned nn, Value* val );
/// Get number of high dimensional derivatives 
  unsigned get_NdX() const ;
/// Get number of low dimensional derivatives
  unsigned get_Ndx() const ;
/// Return the number of base quantities
  unsigned get_nbase() const ;
/// Return the total number of derivatives
  unsigned get_nderivatives() const ;
/// Get number of splines
  void get_nspline( std::vector<unsigned>& nspline ) const ;
/// Get the value of sigma
  double get_sigma() const ;
/// Get the number of spline points
  unsigned getNumberOfSplinePoints() const ;
/// Get the coordinates of one of the spline points
  void getSplinePoint( const unsigned nn, std::vector<double>& pp ) const ;
/// Add some stress
  void addStress( const unsigned ii, const Value& x);
/// Add some derivative
  void addDerivative( const unsigned ii, const unsigned nn, const Value& dx );
/// Complete the field calculation by doing mpi_sum
  void gatherField( PlumedCommunicator& comm );
/// Setup the interpolation tables
  void set_tables();
/// Calculate the value of the field at a particular point in space
  double calculateField( const std::vector<double>& pp ) const ;
/// Calcualte the derivatives of the field at a particular point in space
  void calculateFieldDerivatives( const std::vector<double>& pp, std::vector<double>& tmpforce ) const ;
/// Add some forces to this field
  void addForces( std::vector<double>& inforces );
/// Apply the forces from this field
  bool applyForces( std::vector<double>& outforces ) const ;  
};

inline
unsigned Field::get_NdX() const {
  return ndX;
}

inline
unsigned Field::get_Ndx() const {
  return ndx;
}

inline
unsigned Field::get_nbase() const {
  return baseq_nder.size();
}

inline
double Field::get_sigma() const {
  return sigma;
}

inline
unsigned Field::getNumberOfSplinePoints() const {
  return f_interpolator->getNumberOfSplinePoints();
}

inline
void Field::addStress( const unsigned ii, const Value& x){
  plumed_assert( x.getNumberOfDerivatives()==ndx ); plumed_assert( ii<npoints );
  unsigned kk=ii*nper;
  grid_buffer[kk]+=x.get(); kk++;
  for(unsigned k=0;k<ndx;++k){ grid_buffer[kk]+=x.getDerivative(k); kk++; }
}

inline
void Field::addDerivative( const unsigned ii, const unsigned nn, const Value& dx ){
  plumed_assert( dx.getNumberOfDerivatives()==ndx );
  plumed_assert( ii<npoints && nn<ndX );
  unsigned kk=ii*nper + (nn+1)*(1+ndx);
  grid_buffer[kk]+=dx.get(); kk++;
  for(unsigned k=0;k<ndx;++k){ grid_buffer[kk]+=dx.getDerivative(k); kk++;}
}

}

#endif
