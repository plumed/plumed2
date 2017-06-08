/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#ifndef __PLUMED_mapping_Mapping_h
#define __PLUMED_mapping_Mapping_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "vesselbase/ActionWithVessel.h"
#include "reference/PointWiseMapping.h"
#include <vector>

namespace PLMD {

class PDB;

namespace mapping {

class Mapping :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
{
  friend class TrigonometricPathVessel;
private:
//  The derivative wrt to the distance from the frame
  std::vector<double> dfframes;
/// This holds all the reference information
  PointWiseMapping* mymap;
/// The forces on each of the derivatives (used in apply)
  std::vector<double> forcesToApply;
protected:
/// The (transformed) distance from each frame
  std::vector<double> fframes;
/// Get the number of frames in the path
  unsigned getNumberOfReferencePoints() const ;
/// Finish the setup of the referenceValuePack by transfering atoms and args
  void finishPackSetup( const unsigned& ifunc, ReferenceValuePack& mypack ) const ;
/// Calculate the value of the distance from the ith frame
  double calculateDistanceFunction( const unsigned& ifunc, ReferenceValuePack& myder, const bool& squared ) const ;
/// Store the distance function
  void storeDistanceFunction( const unsigned& ifunc );
/// Get the value of the weight
  double getWeight( const unsigned& weight ) const ;
/// Return the vector of refernece configurations
  std::vector<ReferenceConfiguration*>& getAllReferenceConfigurations();
/// Return a pointer to one of the reference configurations
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& ifunc );
public:
  static void registerKeywords( Keywords& keys );
  explicit Mapping(const ActionOptions&);
  ~Mapping();
/// Overload the virtual functions that appear in both ActionAtomistic and ActionWithArguments
  void turnOnDerivatives();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void lockRequests();
  void unlockRequests();
/// Distance from a point is never periodic
  bool isPeriodic() { return false; }
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Get the value of lambda for paths and property maps
  virtual double getLambda();
/// This does the transformation of the distance by whatever function is required
  virtual double transformHD( const double& dist, double& df ) const=0;
/// Get the number of properties we are projecting onto
  unsigned getNumberOfProperties() const ;
/// Get the name of the ith property we are projecting
  std::string getPropertyName( const unsigned& iprop ) const ;
/// Get the index of a particular named property
  unsigned getPropertyIndex( const std::string& name ) const ;
/// Set the value of one of the projection coordinates
  void setPropertyValue( const unsigned& iframe, const unsigned& iprop, const double& property );
/// Get the name of the ith argument
  std::string getArgumentName( unsigned& iarg );
/// Get the value of the ith property for the current frame
  double getPropertyValue( const unsigned& current, const unsigned& iprop ) const ;
/// Stuff to do before we do the calculation
  void prepare();
/// Apply the forces
  void apply();
};

inline
unsigned Mapping::getNumberOfReferencePoints() const {
  return mymap->getNumberOfMappedPoints();
}

inline
unsigned Mapping::getNumberOfDerivatives() {
  unsigned nat=getNumberOfAtoms();
  if(nat>0) return 3*nat + 9 + getNumberOfArguments();
  return getNumberOfArguments();
}

inline
void Mapping::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

inline
void Mapping::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

inline
unsigned Mapping::getNumberOfProperties() const {
  return mymap->getNumberOfProperties();
}

inline
std::string Mapping::getPropertyName( const unsigned& iprop ) const {
  return mymap->getPropertyName(iprop);
}

inline
double Mapping::getPropertyValue( const unsigned& cur, const unsigned& iprop ) const {
  plumed_dbg_assert( iprop<getNumberOfProperties() );
  return mymap->getPropertyValue( cur, iprop );
}

inline
double Mapping::getWeight( const unsigned& current ) const {
  return mymap->getWeight( current );
}

inline
void Mapping::storeDistanceFunction( const unsigned& ifunc ) {
  plumed_dbg_assert( ifunc<getNumberOfReferencePoints() );
  unsigned storef=getNumberOfReferencePoints()+ifunc;
  fframes[storef]=fframes[ifunc]; dfframes[storef]=dfframes[ifunc];
  mymap->copyFrameDerivatives( ifunc, storef );
}

inline
std::vector<ReferenceConfiguration*>& Mapping::getAllReferenceConfigurations() {
  return mymap->getReferenceConfigurations();
}

}
}
#endif
