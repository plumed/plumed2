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
#ifndef __PLUMED_vesselbase_FieldGridBase_h
#define __PLUMED_vesselbase_FieldGridBase_h

#include "GridVesselBase.h"

namespace PLMD {
namespace vesselbase{

class FieldGridBase : public GridVesselBase {
private:
/// The derivatives wrt to the low-dimensional properties
  std::vector<double> derlow;
/// Were forces applied on this object
  bool wasforced;
/// The forces that are acting on each of the derivatives in this object
  std::vector<double> forces;
protected:
/// Clear the forces from the previous step
  void clearForces();
/// Add value to the field of values, add low-dimensional derivatives and high-dimensional derivatives
  void accumulate( const double& , const double& , const double& , const double& , const unsigned& ); 
public:
/// Create the keywords
  static void registerKeywords( Keywords& keys );
/// The constructor
  FieldGridBase( const VesselOptions& );
/// Resize the field
  virtual void resize();
/// Apply some forces to the field
  bool applyForce(std::vector<double>&);
/// Set the forces on the quantities underlying the fields
  void setForces( const std::vector<double>& );
/// Get the number of base cvs
  unsigned getNumberOfBaseCVs() const ;
/// Return the name of the base CVs
  std::string getBaseCVName( const unsigned& ) const ;
/// Get the value of the field
  double getValue( const unsigned& ip ) const ;
/// Get the derivative
  double getDerivative( const unsigned& ip, const unsigned& jd ) const ;
};

inline
void FieldGridBase::clearForces(){
  wasforced=false; forces.assign( forces.size(), 0.0 );
}

inline
void FieldGridBase::setForces( const std::vector<double>& ff ){
  plumed_dbg_assert( ff.size()==forces.size() );
  wasforced=true; for(unsigned i=0;i<ff.size();++i) forces[i]=ff[i];
}

inline
unsigned FieldGridBase::getNumberOfBaseCVs() const {
  return nper/(1+dimension) - 1;
}

inline
double FieldGridBase::getValue( const unsigned& ip ) const {
  return getGridElement( ip, 0 );
}

inline
double FieldGridBase::getDerivative( const unsigned& ip, const unsigned& jd ) const {
  return getGridElement( ip, (jd+1)*(dimension+1) );
}

}
}
#endif

