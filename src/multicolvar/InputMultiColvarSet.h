/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#ifndef __PLUMED_multicolvar_InputMultiColvarSet_h
#define __PLUMED_multicolvar_InputMultiColvarSet_h

#include "vesselbase/StoreDataVessel.h"
#include "CatomPack.h"
#include "MultiColvarBase.h"

namespace PLMD {

class ActionSet;

namespace multicolvar {

class InputMultiColvarSet {
private:
/// This is used to keep track of what is calculated where
  std::vector<unsigned> colvar_label;
/// The multicolvars from which we construct these quantities
  std::vector<MultiColvarBase*> mybasemulticolvars;
/// The vessels in these multicolvars in which the data is stored
  std::vector<vesselbase::StoreDataVessel*> mybasedata;
/// Convert an index in the global array to an index in the individual base colvars
  unsigned convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const ;
public:
/// This sets up all the base data
  void setup( const std::vector<std::string>& names, const ActionSet& myactions, const double& wtolerance, Action* action );
/// Get the number of base tasks
  unsigned getFullNumberOfBaseTasks() const ;
/// Get the number of colvars in the set
  unsigned getNumberOfBaseMultiColvars() const ;
/// Is the underlying multicolvar currently active
  bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
/// Get the position of the iatom th atom
  Vector getPosition( const unsigned& iatom ) const ;
/// Get the derivative of the position of the iatom th atom
  CatomPack getPositionDerivatives( const unsigned& ind ) const ;
/// Get the vector for a particular task
  void getVectorForTask( const unsigned& ind, const bool& normed, std::vector<double>& orient0 ) const ;
/// Get the derivatives for a particular task
  void getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder0 ) const ;
/// This ensures the derivatives are merged
  void mergeVectorDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end,
                               const unsigned& jatom, const std::vector<double>& der,
                               MultiValue& myder, AtomValuePack& myatoms ) const ;
/// Get the number of colvars calculated by the jth base task
  unsigned getNumberOfTasks( const unsigned& j ) const ;
/// Incorporate all dependencies into underlying action
  void makeDataRequests( Action* action );
/// Recalculate the base colvars (used for numerical derivatives)
  void recalculateBaseColvars( ActionAtomistic* action );
/// Get the total number of atoms involved
  unsigned getTotalNumberOfAtoms() const ;
/// Get the number of atoms in the ith base task
  unsigned getNumberOfAtoms( const unsigned& i ) const ;
/// Get the ith base colvar
  MultiColvarBase* getBaseColvar( const unsigned& i ) const ;
/// Get the colvar this atom is part of 
  unsigned getBaseColvarNumber( const unsigned& i ) const ;
/// Add derivatives to the centre of mass
  void addComDerivatives( const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const ;
};

inline
unsigned InputMultiColvarSet::getFullNumberOfBaseTasks() const {
  return colvar_label.size();
}

inline
unsigned InputMultiColvarSet::getNumberOfBaseMultiColvars() const {
  return mybasemulticolvars.size();
}

inline
unsigned InputMultiColvarSet::getNumberOfTasks( const unsigned& j ) const {
  plumed_dbg_assert( j<mybasemulticolvars.size() ); return mybasemulticolvars[j]->getFullNumberOfTasks();
}

inline
unsigned InputMultiColvarSet::convertToLocalIndex( const unsigned& index, const unsigned& mcv_code ) const {
  unsigned t1 = index;
  for(unsigned k=0;k<mcv_code;++k) t1 -= mybasemulticolvars[k]->getFullNumberOfTasks();
  return t1;
}

inline
bool InputMultiColvarSet::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  plumed_dbg_assert( code<getFullNumberOfBaseTasks() ); unsigned mmc=colvar_label[code];
  return mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(code,mmc) ); 
} 

inline
Vector InputMultiColvarSet::getPosition( const unsigned& iatom ) const {
  plumed_dbg_assert( iatom<getFullNumberOfBaseTasks() ); unsigned mmc=colvar_label[ iatom ];
  return mybasemulticolvars[mmc]->getCentralAtomPos( convertToLocalIndex(iatom,mmc) );
}

inline
CatomPack InputMultiColvarSet::getPositionDerivatives( const unsigned& ind ) const {
  plumed_dbg_assert( ind<getFullNumberOfBaseTasks() ); unsigned mmc=colvar_label[ind];
  unsigned basen=0;
  for(unsigned i=0;i<mmc;++i) basen+=mybasemulticolvars[i]->getNumberOfAtoms();
  return mybasemulticolvars[mmc]->getCentralAtomPack( basen, convertToLocalIndex(ind,mmc) );
}

inline
void InputMultiColvarSet::getVectorForTask( const unsigned& ind, const bool& normed, std::vector<double>& orient ) const {
  plumed_dbg_assert( ind<getFullNumberOfBaseTasks() ); unsigned mmc=colvar_label[ind];
  plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(ind,mmc) ) );
  mybasedata[mmc]->retrieveValue( convertToLocalIndex(ind,mmc), normed, orient );
}

inline
unsigned InputMultiColvarSet::getTotalNumberOfAtoms() const {
  unsigned nat=0; for(unsigned i=0;i<mybasemulticolvars.size();++i) nat += mybasemulticolvars[i]->getNumberOfAtoms();
  return nat;
}

inline
unsigned InputMultiColvarSet::getNumberOfAtoms( const unsigned& i ) const {
  plumed_dbg_assert( i<mybasemulticolvars.size() ); return mybasemulticolvars[i]->getNumberOfAtoms();
}

inline
MultiColvarBase* InputMultiColvarSet::getBaseColvar( const unsigned& i ) const {
  plumed_dbg_assert( i<mybasemulticolvars.size() ); return mybasemulticolvars[i]; 
}

inline
unsigned InputMultiColvarSet::getBaseColvarNumber( const unsigned& i ) const {
  plumed_dbg_assert( i<colvar_label.size() ); return colvar_label[i];
}

}
}
#endif
