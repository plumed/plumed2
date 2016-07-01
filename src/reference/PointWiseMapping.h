/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#ifndef __PLUMED_reference_PointWiseMapping_h
#define __PLUMED_reference_PointWiseMapping_h

#include "MultiReferenceBase.h"

namespace PLMD {

class PointWiseMapping : public MultiReferenceBase {
private:
/// This is a path
  bool ispath;
/// The matrix of interframe distances
  Matrix<double> dmat;
/// The names of the projection coordinates
  std::vector<std::string> property;
/// These are where the reference configurations should be projected
  std::vector< std::vector<double> > low_dim;  
public:
  PointWiseMapping( const std::string& type, const bool& checksoff );
/// Set the names of the low dimensional properties
  void setPropertyNames( const std::vector<std::string>& prop, const bool isp );
/// Check if setup was completed
  bool mappingNeedsSetup() const;
/// Delete the low dimensional projections
  void clearRestOfData();
/// Read in the data from a file 
  void readRestOfFrame();
/// Resize everything else from a file
  void resizeRestOfFrame();
/// Make a second copy of the frame list 
  void duplicateFrameList();
/// Get the number of points we are mapping into the lower dimensional space
  unsigned getNumberOfMappedPoints() const ;
/// Get the number of properties
  unsigned getNumberOfProperties() const ;
/// Get the name of the ith property
  std::string getPropertyName( const unsigned& i ) const ;
/// Get the index of the property with name
  unsigned getPropertyIndex( const std::string& name ) const ;
/// Get the value of the ith property for th jth frame
  double getPropertyValue( const unsigned& iframe, const unsigned& jprop ) const ;
/// Get the derivatives wrt to the position of an atom
//  Vector getAtomDerivatives( const unsigned& iframe, const unsigned& jatom );
/// Get the derivatives wrt to the box
//  bool getVirial( const unsigned& iframe, Tensor& vir );
/// Ge the derivatives wrt to one of the arguments
//  double getArgumentDerivative( const unsigned& iframe, const unsigned& jarg );
/// Copy derivative information from frame number from to frame number to
  void copyFrameDerivatives( const unsigned& from, const unsigned& to );
/// Get a pointer to the matrix of pairwise distances
  Matrix<double>& modifyDmat();
/// Print out the low dimensional mapping
  void print( const std::string& method, const double & time, OFile& afile, 
              const std::string& fmt, const double& lunits );
/// Get the low dimensional embedding coordinate
  double getProjectionCoordinate( const unsigned& iframe, const unsigned& jcoord ) const ;
/// Set the value of the projection coordinate
  void setProjectionCoordinate( const unsigned& iframe, const unsigned& jcoord, const double& coord );
};

inline
bool PointWiseMapping::mappingNeedsSetup() const {
  bool didsetup=(frames.size()==2*low_dim.size());
  return !didsetup;
}

inline
void PointWiseMapping::copyFrameDerivatives( const unsigned& from, const unsigned& to ){
  plumed_dbg_assert( to>=frames.size()/2 && from<frames.size()/2 );
  frames[to]->copyDerivatives( frames[from] );
}

inline
unsigned PointWiseMapping::getNumberOfMappedPoints() const {
  return low_dim.size();
}

inline
unsigned PointWiseMapping::getNumberOfProperties() const {
  return property.size();
}

inline
std::string PointWiseMapping::getPropertyName( const unsigned& i ) const {
  plumed_dbg_assert( i<property.size() );
  return property[i];
}

inline
double PointWiseMapping::getPropertyValue( const unsigned& iframe, const unsigned& jprop ) const {
  plumed_dbg_assert( iframe<low_dim.size() && jprop<property.size() );
  return low_dim[iframe][jprop];
}

// inline
// Vector PointWiseMapping::getAtomDerivatives( const unsigned& iframe, const unsigned& jatom ){
//   return frames[iframe]->getAtomDerivative(jatom);
// }
// 
// inline
// bool PointWiseMapping::getVirial( const unsigned& iframe, Tensor& vir ){
//   return frames[iframe]->getVirial( vir );
// }

// inline
// double PointWiseMapping::getArgumentDerivative( const unsigned& iframe, const unsigned& jarg ){
//   return frames[iframe]->getArgumentDerivative(jarg);
// }

inline
Matrix<double>& PointWiseMapping::modifyDmat(){
  if( dmat.nrows()!=frames.size() || dmat.ncols()!=frames.size() ) dmat.resize( frames.size(), frames.size() );
  return dmat;
}

inline
double PointWiseMapping::getProjectionCoordinate( const unsigned& iframe, const unsigned& jcoord ) const {
  plumed_dbg_assert( iframe<frames.size() && jcoord<property.size() );
  return low_dim[iframe][jcoord];
}

inline
void PointWiseMapping::setProjectionCoordinate( const unsigned& iframe, const unsigned& jcoord, const double& coord ){
  plumed_dbg_assert( iframe<frames.size() && jcoord<property.size() );
  low_dim[iframe][jcoord]=coord;
} 

}
#endif
