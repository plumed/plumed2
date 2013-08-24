/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_multicolvar_ActionVolume_h
#define __PLUMED_multicolvar_ActionVolume_h

#include "core/ActionAtomistic.h"
#include "tools/HistogramBead.h"
#include "tools/Pbc.h"
#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include "vesselbase/BridgeVessel.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing a new way of definining a particular region of the simulation
box. You can use this to calculate the number of atoms inside that part or the average value of a quantity like the 
coordination number inside that part of the cell. 
*/

class ActionVolume :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
friend class Region;
private:
/// The value of sigma
  double sigma;
/// Are we interested in the area outside the colvar
  bool not_in;
/// This is used for storing positions properly
  Vector tmp_p;
/// The bead for the histogram
  HistogramBead bead;
/// The action that is calculating the colvars of interest
  MultiColvarBase* mycolv;
/// The vessel that bridges
  vesselbase::BridgeVessel* myBridgeVessel;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  unsigned lastUpdate;
/// Fast merging of derivatives (automatic skips of zero contributions)
  DynamicList<unsigned> activeAtoms;
/// This is used to store forces temporarily in apply
  std::vector<double> tmpforces;
/// This sets up array above
  void resizeLocalArrays();
protected:
  double getSigma() const ;
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
  MultiColvarBase* getPntrToMultiColvar();
/// Add derivatinve to one of the reference atoms here
  void addReferenceAtomDerivatives( const unsigned& iatom, const Vector& der );
/// Add derivatives wrt to the virial
  void addBoxDerivatives( const Tensor& vir );
public:
  static void registerKeywords( Keywords& keys );
  ActionVolume(const ActionOptions&);
/// Don't actually clear the derivatives when this is called from plumed main.  
/// They are calculated inside another action and clearing them would be bad  
  void clearDerivatives(){}
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Is the output quantity periodic
  bool isPeriodic();
/// Jobs to be done when the action is activated
  void prepare();
/// Do jobs required before tasks are undertaken
  void doJobsRequiredBeforeTaskList();
/// This calculates all the vessels and is called from within a bridge vessel
  void performTask();
/// Routines that have to be defined so as not to have problems with virtual methods 
  void deactivate_task();
  void calculate(){}
/// We need our own calculate numerical derivatives here
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  virtual void setupRegion()=0;
  virtual double calculateNumberInside( const Vector& cpos, HistogramBead& bead, Vector& derivatives )=0;
/// Forces here are applied through the bridge
  void applyBridgeForces( const std::vector<double>& bb );
  void apply(){};
/// These routines replace the virtual routines in ActionWithVessel for 
/// code optimization
  void mergeDerivatives( const unsigned& ider, const double& df );
  void clearDerivativesAfterTask( const unsigned& ider );
};

inline
const Tensor & ActionVolume::getBox()const{
  return mycolv->getBox();
} 

inline
const Pbc & ActionVolume::getPbc() const {
 return mycolv->getPbc();
}

inline
Vector ActionVolume::pbcDistance( const Vector& v1, const Vector& v2) const {
 return mycolv->pbcDistance(v1,v2);
}

inline
const Vector & ActionVolume::getPosition( int iatom ){
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
double ActionVolume::getSigma() const {
  return sigma;
}

inline
MultiColvarBase* ActionVolume::getPntrToMultiColvar(){
  return mycolv;
}

inline
unsigned ActionVolume::getNumberOfDerivatives(){
  return mycolv->getNumberOfDerivatives() + 3*getNumberOfAtoms();
}

inline
void ActionVolume::addReferenceAtomDerivatives( const unsigned& iatom, const Vector& der ){
  // This is used for storing the derivatives wrt to the 
  // positions of any additional reference atoms
  double pref=mycolv->getElementValue(1); 
  if( not_in ) pref*=-1;
  unsigned nstart = getNumberOfDerivatives() + mycolv->getNumberOfDerivatives() + 3*iatom;
  addElementDerivative( nstart + 0, pref*der[0] );
  addElementDerivative( nstart + 1, pref*der[1] );
  addElementDerivative( nstart + 2, pref*der[2] );
}

inline
void ActionVolume::addBoxDerivatives( const Tensor& vir ){
  double pref=mycolv->getElementValue(1);
  if( not_in ) pref*=-1;
  unsigned nstart = getNumberOfDerivatives() + mycolv->getNumberOfDerivatives() - 9;
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
