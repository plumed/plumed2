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
#ifndef __PLUMED_FieldVessel_h
#define __PLUMED_FieldVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "VesselValueAccess.h"

namespace PLMD {

class CInterpolation;

class FieldVessel : public VesselStoreAllValues {
private:
/// When do we merge derivatives
  bool mergeBeforeInterpol;
/// Number of high dimensions
  unsigned ndX; 
/// Number of low dimensions
  unsigned ndx;
/// Number of values at each point on the grid
  unsigned nper;
/// Storage space for the grid
  std::vector<double> grid_buffer;
/// The interpolator for the field
  CInterpolation* f_interpolator;
/// The interpolator for the derivatives
  std::vector<CInterpolation*> df_interpolators;
/// Are there forces on this field
  bool wasForced;
/// The forces on each of the derivatives
  std::vector<double> forces;
//// A value used in calculations
  Value tmpvalue;
protected:
/// Return a bool that tells us when to interpolate
  bool mergeBeforeInterpolation() const ;
/// Return the field at the point
  double getFieldAt( const std::vector<double>& pp ); 
/// Return the derivatives of the field at the point
  void getDerivativesAt( const std::vector<double>& pp, std::vector<double>& tmpforce );
public:
/// This returns some documentation for fields (used to create docs)
  static std::string documentation();
/// This returns the keywords (used when code crashes)
  static void getKeywords( Keywords& );
  FieldVessel( const VesselOptions& da );
  virtual ~FieldVessel();
/// This does the resizing of the field
  void local_resizing();
/// This calculates the field on a grid
  void finish( const double& tolerance );
/// This calculates the energy at a point due to the kth member
  virtual void calculateEnergy( const unsigned&, const std::vector<double>&, const Value& , Value& , std::vector<Value>& )=0;
///  Get number of low dimensional derivatives
  unsigned get_Ndx() const ;
/// Get number of high dimensional derivatives 
  unsigned get_NdX() const ;
/// Get the number of spline points
  std::vector<unsigned> get_nspline() const ;
/// Get the minimums in each direction
  std::vector<double> get_min() const ;
/// Get the maximums in each direction
  std::vector<double> get_max() const ;
/// Calculate the value of the field at a particular point in space
  virtual double calculateField( const std::vector<double>& pp )=0;
/// Calcualte the derivatives of the field at a particular point in space
  virtual void calculateFieldDerivatives( const std::vector<double>& pp, std::vector<double>& tmpforce )=0;
/// Add some forces to this field
  void addForces( const std::vector<double>& inforces );
/// Apply the forces to the atoms
  bool applyForce( std::vector<double>& outforces );
/// Merge the derivatives onto a value for output
  void mergeFieldDerivatives( const std::vector<double>& inder, Value* outval );
};

inline
unsigned FieldVessel::get_Ndx() const {
  return ndx;
}

inline
unsigned FieldVessel::get_NdX() const {
  return ndX;
}  

inline
bool FieldVessel::mergeBeforeInterpolation() const {
  return mergeBeforeInterpol;
}

}
#endif
