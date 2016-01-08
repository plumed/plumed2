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
  keys.reserve("compulsory","WTOL","if the base multicolvars have weights then you must define a hard cutoff on those you want to consider explicitally");
  keys.remove("NUMERICAL_DERIVATIVES");
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the arguments
  if( keywords.exists("DATA") ){ 
      std::vector<std::string> mlabs; parseVector("DATA",mlabs);

      if( keywords.exists("WTOL") ){
          double wtolerance; parse("WTOL",wtolerance); 
          log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
          bool found_acts=interpretInputMultiColvars(mlabs,wtolerance);
          if( !found_acts ) error("one or more items in input is not the label of a multicolvar");
      } else {
          bool found_acts=interpretInputMultiColvars(mlabs,0.0);
          if( !found_acts ) error("one or more items in input is not the label of a multicolvar");
      }
  }
}

void MultiColvarFunction::setupAtomLists(){
  // Make all atom requests and setup dependencies
  std::vector<AtomNumber> fake_atoms; 
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase( fake_atoms );
}

void MultiColvarFunction::buildSymmetryFunctionLists(){
  if( mybasemulticolvars.size()>2 ) error("Found too many multicolvars in DATA specification. You can use either 1 or 2");

  usespecies=true; ablocks.resize( 1 );
  for(unsigned i=0;i<mybasemulticolvars[0]->getFullNumberOfTasks();++i) addTaskToList( i );

  unsigned ntotal=0;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      ntotal += mybasemulticolvars[i]->getFullNumberOfTasks();
  }
  unsigned k=0, start=0;
  ablocks[0].resize( ntotal ); 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      for(unsigned j=0;j<mybasemulticolvars[i]->getFullNumberOfTasks();++j){
          ablocks[0][k]=start + j; k++;
      }
      start += mybasemulticolvars[i]->getFullNumberOfTasks();
  }  
  mybasedata[0]->resizeTemporyMultiValues(2); setupAtomLists();
}

void MultiColvarFunction::buildSets(){
  nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
     if( mybasemulticolvars[i]->getFullNumberOfTasks()!=nblock ){
          error("mismatch between numbers of tasks in various base multicolvars");
     }
  }
  ablocks.resize( mybasemulticolvars.size() );
  usespecies=false; 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      ablocks[i].resize( nblock ); 
      for(unsigned j=0;j<nblock;++j) ablocks[i][j]=i*nblock+j;  
  }
  for(unsigned i=0;i<nblock;++i){
      if( mybasemulticolvars.size()<4 ){
          unsigned cvcode=0, tmpc=1;
          for(unsigned j=0;j<ablocks.size();++j){ cvcode +=i*tmpc; tmpc *= nblock; }
          addTaskToList( cvcode );
      } else {
          addTaskToList( i );
      }
  }
  mybasedata[0]->resizeTemporyMultiValues( mybasemulticolvars.size() ); setupAtomLists();
}

void MultiColvarFunction::buildAtomListWithPairs( const bool& allow_intra_group ){
  if( !allow_intra_group && mybasemulticolvars.size()>2 ) error("only two input multicolvars allowed with this function"); 
  
  usespecies=false; ablocks.resize(2); 
  if( !allow_intra_group && mybasemulticolvars.size()==2 ){
     nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
     if( mybasemulticolvars[1]->getFullNumberOfTasks()>nblock ) nblock = mybasemulticolvars[1]->getFullNumberOfTasks();
    
     ablocks[0].resize( mybasemulticolvars[0]->getFullNumberOfTasks() );
     for(unsigned i=0;i<mybasemulticolvars[0]->getFullNumberOfTasks();++i) ablocks[0][i] = i;
     ablocks[1].resize( mybasemulticolvars[1]->getFullNumberOfTasks() ); unsigned istart = ablocks[0].size();
     for(unsigned i=0;i<mybasemulticolvars[1]->getFullNumberOfTasks();++i) ablocks[1][i] = istart + i;
     resizeBookeepingArray( ablocks[0].size(), ablocks[1].size() );
     for(unsigned i=0;i<ablocks[0].size();++i){
         for(unsigned j=0;j<ablocks[1].size();++j){
            bookeeping(i,j).first=getFullNumberOfTasks();
            addTaskToList( i*nblock + j );
            bookeeping(i,j).second=getFullNumberOfTasks();
         }
     }
  } else {
     nblock = 0; for(unsigned i=0;i<mybasemulticolvars.size();++i) nblock += mybasemulticolvars[i]->getFullNumberOfTasks(); 
     ablocks[0].resize( nblock ); ablocks[1].resize( nblock ); resizeBookeepingArray( nblock, nblock );
     for(unsigned i=0;i<nblock;++i){ ablocks[0][i] = i; ablocks[1][i] = i; }
     for(unsigned i=1;i<nblock;++i){
        for(unsigned j=0;j<i;++j){
           bookeeping(i,j).first=getFullNumberOfTasks();
           addTaskToList( i*nblock + j );
           bookeeping(i,j).second=getFullNumberOfTasks();
        }
     }
  }
  mybasedata[0]->resizeTemporyMultiValues(2); setupAtomLists(); 
}

MultiValue& MultiColvarFunction::getVectorDerivatives( const unsigned& ind, const bool& normed ) const {
  plumed_dbg_assert( ind<colvar_label.size() ); unsigned mmc=colvar_label[ind];
  plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(ind,mmc) ) );
  // Get a tempory multi value from the base class
  MultiValue& myder=mybasedata[0]->getTemporyMultiValue();

  if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
          myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  mybasedata[mmc]->retrieveDerivatives( convertToLocalIndex(ind,mmc), normed, myder );
  return myder;
}

void MultiColvarFunction::mergeVectorDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end, 
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

}
}

