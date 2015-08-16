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
#include "core/ActionSet.h"
#include "AtomValuePack.h" 
#include "BridgedMultiColvarFunction.h"
#include "InputMultiColvarSet.h"

namespace PLMD {
namespace multicolvar {

void InputMultiColvarSet::setup( const std::vector<std::string>& mlabs, const ActionSet& myactions, const double& wtolerance, Action* action ){
   std::string mname;
   for(unsigned i=0;i<mlabs.size();++i){
      MultiColvarBase* mycolv = myactions.selectWithLabel<MultiColvarBase*>(mlabs[i]);
      if(!mycolv) action->error("action labeled " + mlabs[i] + " does not exist or is not a multicolvar");
      // Check all base multicolvars are of same type
      if( i==0 ){ 
          mname = mycolv->getName();
          if( mycolv->isPeriodic() ) action->error("multicolvar functions don't work with this multicolvar");
      } else {
          if( mname!=mycolv->getName() ) action->error("All input multicolvars must be of same type"); 
      }
      // And track which variable stores each colvar
      for(unsigned j=0;j<mycolv->getFullNumberOfTasks();++j) colvar_label.push_back( mybasemulticolvars.size() );
      // And store the multicolvar base
      mybasemulticolvars.push_back( mycolv );
      // And create a basedata stash
      mybasedata.push_back( mybasemulticolvars[mybasemulticolvars.size()-1]->buildDataStashes( true, wtolerance ) );
      plumed_assert( mybasemulticolvars.size()==mybasedata.size() );
   }
}

void InputMultiColvarSet::makeDataRequests( const std::vector<AtomNumber>& atoms, Action* action ){

  // Copy lists of atoms involved from base multicolvars 
  std::vector<AtomNumber> tmp_atoms, all_atoms; 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      BridgedMultiColvarFunction* mybr=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[i] );
      if( mybr ) tmp_atoms=(mybr->getPntrToMultiColvar())->getAbsoluteIndexes();
      else tmp_atoms=mybasemulticolvars[i]->getAbsoluteIndexes();
      for(unsigned j=0;j<tmp_atoms.size();++j) all_atoms.push_back( tmp_atoms[j] );
  }  
  // Get additional atom requests
  for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );     
          
  // Now make sure we get all the atom positions 
  ActionAtomistic* myatoms = dynamic_cast<ActionAtomistic*>( action );
  plumed_assert( myatoms ); myatoms->requestAtoms( all_atoms );

  for(unsigned i=0;i<mybasemulticolvars.size();++i) action->addDependency( mybasemulticolvars[i] );
}

void InputMultiColvarSet::getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder ) const {
  plumed_dbg_assert( ind<colvar_label.size() ); unsigned mmc=colvar_label[ind];
  plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(ind,mmc) ) );
  if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
          myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  mybasedata[mmc]->retrieveDerivatives( convertToLocalIndex(ind,mmc), normed, myder );
}

void InputMultiColvarSet::mergeVectorDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end,
                                                  const unsigned& jatom, const std::vector<double>& der,
                                                  MultiValue& myder, AtomValuePack& myatoms ) const {
  plumed_dbg_assert( ival<myatoms.getUnderlyingMultiValue().getNumberOfValues() );
  plumed_dbg_assert( start<myder.getNumberOfValues() && end<=myder.getNumberOfValues() );
  plumed_dbg_assert( der.size()==myder.getNumberOfValues() && jatom<getFullNumberOfBaseTasks() );

  unsigned mmc=colvar_label[jatom]; plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(jatom,mmc) ) );

  // Get start of indices for this atom
  unsigned basen=0; for(unsigned i=0;i<mmc;++i) basen+=3*mybasemulticolvars[i]->getNumberOfAtoms();

  MultiValue& myvals=myatoms.getUnderlyingMultiValue();
  // Now get the start of the virial
  unsigned virbas = myvals.getNumberOfDerivatives()-9;
  for(unsigned j=0;j<myder.getNumberActive();++j){
     unsigned jder=myder.getActiveIndex(j);
     if( jder<3*mybasemulticolvars[mmc]->getNumberOfAtoms() ){
         unsigned kder=basen+jder;
         for(unsigned icomp=start;icomp<end;++icomp){
             myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
         }
     } else {
         unsigned kder=virbas + (jder - 3*mybasemulticolvars[mmc]->getNumberOfAtoms());
         for(unsigned icomp=start;icomp<end;++icomp){
             myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
         }
     }
  }
}

void InputMultiColvarSet::addComDerivatives( const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned mmc=colvar_label[iatom];
  unsigned basen=0; for(unsigned i=0;i<mmc;++i) basen+=mybasemulticolvars[i]->getNumberOfAtoms();
  multicolvar::CatomPack atom0=mybasemulticolvars[mmc]->getCentralAtomPack( basen, convertToLocalIndex(iatom,mmc) );
  myatoms.addComDerivatives( 1, der, atom0 );
}

void InputMultiColvarSet::recalculateBaseColvars( ActionAtomistic* action ){
     // Clear any derivatives in base colvar that were 
     // accumulated from previous calculations
     for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->clearDerivatives(); 
     // And recalculate
     for(unsigned i=0;i<mybasemulticolvars.size();++i){
        BridgedMultiColvarFunction* bb=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[i] );
        if( bb ) (bb->getPntrToMultiColvar())->calculate();
        else mybasemulticolvars[i]->calculate();
     }
     // Copy the box from the base multicolvar here
     unsigned maxb=mybasemulticolvars.size() - 1;
     BridgedMultiColvarFunction* bb=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[maxb] );
     if( bb ) action->changeBox( (bb->getPntrToMultiColvar())->getBox() ); 
     else action->changeBox( mybasemulticolvars[maxb]->getBox() );
}

}
}
