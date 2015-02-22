/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#ifndef __PLUMED_multicolvar_VolumeGradientBase_h
#define __PLUMED_multicolvar_VolumeGradientBase_h

#include "BridgedMultiColvarFunction.h"

namespace PLMD {
namespace multicolvar {

class VolumeGradientBase : public BridgedMultiColvarFunction {
friend class MultiColvarBase;   
private:
/// This is used for storing positions properly
  Vector tmp_p;
/// This is used to store forces temporarily in apply
  std::vector<double> tmpforces;
protected:
/// Get the cell box
  const Tensor & getBox() const;
/// Get reference to Pbc
  const Pbc & getPbc() const;
/// Calculate distance between two points
  Vector pbcDistance( const Vector& v1, const Vector& v2) const;
/// Get position of atom
  const Vector & getPosition( int iatom );
/// Request the atoms 
  void requestAtoms( const std::vector<AtomNumber>& atoms );
/// Set the number in the volume
  void setNumberInVolume( const unsigned& ivol, const double& weight, const Vector& wdf );
/// Add derivatinve to one of the reference atoms here
  void addReferenceAtomDerivatives( const unsigned& jvol, const unsigned& iatom, const Vector& der );
/// Add derivatives wrt to the virial
  void addBoxDerivatives( const unsigned& jvol, const Tensor& vir );
public:
  static void registerKeywords( Keywords& keys );
  VolumeGradientBase(const ActionOptions&);
/// Do jobs required before tasks are undertaken
  void doJobsRequiredBeforeTaskList();
/// Actually do what we are asked
  void completeTask();
/// Calculate what is in the volumes
  virtual void calculateAllVolumes()=0;
/// Setup the regions that this is based on 
  virtual void setupRegions()=0;
/// Forces here are applied through the bridge
  void addBridgeForces( const std::vector<double>& bb );
};

inline
const Tensor & VolumeGradientBase::getBox()const{
  return getPntrToMultiColvar()->getBox();
} 

inline
const Pbc & VolumeGradientBase::getPbc() const {
 return getPntrToMultiColvar()->getPbc();
}

inline
Vector VolumeGradientBase::pbcDistance( const Vector& v1, const Vector& v2) const {
 return getPntrToMultiColvar()->pbcDistance(v1,v2);
}

inline
const Vector & VolumeGradientBase::getPosition( int iatom ){
 if( !checkNumericalDerivatives() ) return ActionAtomistic::getPosition(iatom);
 // This is for numerical derivatives of quantity wrt to the local atoms
 tmp_p = ActionAtomistic::getPosition(iatom);
 if( bridgeVariable<3*getNumberOfAtoms() ){
    if( bridgeVariable>=3*iatom && bridgeVariable<(iatom+1)*3 ) tmp_p[bridgeVariable%3]+=sqrt(epsilon);
 }
 // This makes sure that numerical derivatives of virial are calculated correctly
 tmp_p = ActionAtomistic::getPbc().realToScaled( tmp_p ); 
 tmp_p = getPbc().scaledToReal( tmp_p );
 return tmp_p;
} 

inline
void VolumeGradientBase::addReferenceAtomDerivatives( const unsigned& jvol, const unsigned& iatom, const Vector& der ){
  // This is used for storing the derivatives wrt to the 
  // positions of any additional reference atoms
  double pref=getPntrToMultiColvar()->getElementValue(1); 
  unsigned nstart = jvol*getNumberOfDerivatives() + getPntrToMultiColvar()->getNumberOfDerivatives() + 3*iatom;
  addElementDerivative( nstart + 0, pref*der[0] );
  addElementDerivative( nstart + 1, pref*der[1] );
  addElementDerivative( nstart + 2, pref*der[2] );
}

inline
void VolumeGradientBase::addBoxDerivatives( const unsigned& jvol, const Tensor& vir ){
  double pref=getPntrToMultiColvar()->getElementValue(1);
  unsigned nstart = jvol*getNumberOfDerivatives() + getPntrToMultiColvar()->getNumberOfDerivatives() - 9;
  addElementDerivative( nstart + 0, pref*vir(0,0) );  
  addElementDerivative( nstart + 1, pref*vir(0,1) ); 
  addElementDerivative( nstart + 2, pref*vir(0,2) ); 
  addElementDerivative( nstart + 3, pref*vir(1,0) ); 
  addElementDerivative( nstart + 4, pref*vir(1,1) ); 
  addElementDerivative( nstart + 5, pref*vir(1,2) ); 
  addElementDerivative( nstart + 6, pref*vir(2,0) ); 
  addElementDerivative( nstart + 7, pref*vir(2,1) ); 
  addElementDerivative( nstart + 8, pref*vir(2,2) ); 
}

}
}
#endif
