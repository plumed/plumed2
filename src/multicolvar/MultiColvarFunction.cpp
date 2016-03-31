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
#include "MultiColvarFunction.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Pbc.h"
#include "MultiColvarBase.h"
#include "BridgedMultiColvarFunction.h"

namespace PLMD {
namespace multicolvar { 

void MultiColvarFunction::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","the labels of the action that calculates the multicolvars we are interested in");
  keys.remove("NUMERICAL_DERIVATIVES");
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the arguments
  if( keywords.exists("DATA") ){ 
      std::vector<AtomNumber> fake_atoms; 
      if( !parseMultiColvarAtomList("DATA",-1,fake_atoms) ) error("missing DATA keyword");
      if( fake_atoms.size()>0 ) error("no atoms should appear in the specification for this object.  Input should be other multicolvars");
  }
}

void MultiColvarFunction::buildSets(){
  nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
     if( mybasemulticolvars[i]->getFullNumberOfTasks()!=nblock ){
          error("mismatch between numbers of tasks in various base multicolvars");
     }
  }
  nblock=0; ablocks.resize( mybasemulticolvars.size() );
  usespecies=false; // current_atoms.resize( mybasemulticolvars.size() );
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      ablocks[i].resize( nblock ); 
      for(unsigned j=0;j<nblock;++j) ablocks[i][j]=i*nblock+j;  
  }
  for(unsigned i=0;i<nblock;++i) addTaskToList( i );
  std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms ); 
  mybasedata[0]->resizeTemporyMultiValues( mybasemulticolvars.size() ); 
}

MultiValue& MultiColvarFunction::getVectorDerivatives( const unsigned& ind, const bool& normed ) const {
  plumed_dbg_assert( atom_lab[ind].first>0 ); unsigned mmc=atom_lab[ind].first-1;
  plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( atom_lab[ind].second ) );
  // Get a tempory multi value from the base class
  MultiValue& myder=mybasedata[0]->getTemporyMultiValue();

  if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
          myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  mybasedata[mmc]->retrieveDerivatives( atom_lab[ind].second, normed, myder );
  return myder;
}

void MultiColvarFunction::mergeVectorDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end, 
                                                  const unsigned& jatom, const std::vector<double>& der, 
                                                  MultiValue& myder, AtomValuePack& myatoms ) const {
  plumed_dbg_assert( ival<myatoms.getUnderlyingMultiValue().getNumberOfValues() );
  plumed_dbg_assert( start<myder.getNumberOfValues() && end<=myder.getNumberOfValues() );
  plumed_dbg_assert( der.size()==myder.getNumberOfValues() && jatom<getFullNumberOfBaseTasks() );

  unsigned mmc=atom_lab[jatom].first - 1; plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( atom_lab[ind].second ) );

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

}
}

