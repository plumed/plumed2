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
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the arguments
  std::string mname; 
  std::vector<std::string> mlabs; parseVector("DATA",mlabs);

  if( keywords.exists("WTOL") ){
      double wtolerance; parse("WTOL",wtolerance); 
      log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
      myinputdata.setup( mlabs, plumed.getActionSet(), wtolerance, this );
  } else {
      myinputdata.setup( mlabs, plumed.getActionSet(), 0.0, this );
  }
  log.printf("  using colvars calculated by actions "); 
  for(unsigned i=0;i<mlabs.size();++i) log.printf("%s ",mlabs[i].c_str() );
  log.printf("\n");
}

void MultiColvarFunction::setupAtomLists(){
  // Make all atom requests and setup dependencies
  myinputdata.makeDataRequests( this );
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase();
}

void MultiColvarFunction::buildSymmetryFunctionLists(){
  if( myinputdata.getNumberOfBaseMultiColvars()>2 ) error("Found too many multicolvars in DATA specification. You can use either 1 or 2");

  usespecies=true; ablocks.resize( 1 );
  for(unsigned i=0;i<myinputdata.getNumberOfTasks(0);++i) addTaskToList( i );

  unsigned ntotal=0;
  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
      ntotal += myinputdata.getNumberOfTasks(i); // mybasemulticolvars[i]->getFullNumberOfTasks();
  }
  unsigned k=0, start=0;
  // current_atoms.resize( 1 + ntotal ); 
  ablocks[0].resize( ntotal ); 
  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
      for(unsigned j=0;j<myinputdata.getNumberOfTasks(i);++j){
          ablocks[0][k]=start + j; k++;
      }
      start += myinputdata.getNumberOfTasks(i); //mybasemulticolvars[i]->getFullNumberOfTasks();
  }  
  setupAtomLists();
}

void MultiColvarFunction::buildSets(){
  nblock = myinputdata.getNumberOfTasks(0); // mybasemulticolvars[0]->getFullNumberOfTasks();
  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
     if( myinputdata.getNumberOfTasks(i)!=nblock ){
          error("mismatch between numbers of tasks in various base multicolvars");
     }
  }
  ablocks.resize( myinputdata.getNumberOfBaseMultiColvars() );
  usespecies=false; // current_atoms.resize( mybasemulticolvars.size() );
  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
      ablocks[i].resize( nblock ); 
      for(unsigned j=0;j<nblock;++j) ablocks[i][j]=i*nblock+j;  
  }
  for(unsigned i=0;i<nblock;++i){
      if( myinputdata.getNumberOfBaseMultiColvars()<4 ){
          unsigned cvcode=0, tmpc=1;
          for(unsigned j=0;j<ablocks.size();++j){ cvcode +=i*tmpc; tmpc *= nblock; }
          addTaskToList( cvcode );
      } else {
          addTaskToList( i );
      }
  }
  setupAtomLists();
}

void MultiColvarFunction::buildAtomListWithPairs( const bool& allow_intra_group ){
  if( !allow_intra_group && myinputdata.getNumberOfBaseMultiColvars()>2 ) error("only two input multicolvars allowed with this function"); 
  
  usespecies=false; ablocks.resize(2); // current_atoms.resize( 2 );
  if( !allow_intra_group && myinputdata.getNumberOfBaseMultiColvars()==2 ){
     nblock = myinputdata.getNumberOfTasks(0);
     if( myinputdata.getNumberOfTasks(1)>nblock ) nblock = myinputdata.getNumberOfTasks(1);
    
     ablocks[0].resize(myinputdata.getNumberOfTasks(0) );
     for(unsigned i=0;i<myinputdata.getNumberOfTasks(0);++i) ablocks[0][i] = i;
     ablocks[1].resize( myinputdata.getNumberOfTasks(1) ); unsigned istart = ablocks[0].size();
     for(unsigned i=0;i<myinputdata.getNumberOfTasks(1);++i) ablocks[1][i] = istart + i;
     resizeBookeepingArray( ablocks[0].size(), ablocks[1].size() );
     for(unsigned i=0;i<ablocks[0].size();++i){
         for(unsigned j=0;j<ablocks[1].size();++j){
            bookeeping(i,j).first=getFullNumberOfTasks();
            addTaskToList( i*nblock + j );
            bookeeping(i,j).second=getFullNumberOfTasks();
         }
     }
  } else {
     nblock = 0; for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i) nblock += myinputdata.getNumberOfTasks(i); 
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
  setupAtomLists(); 
}

void MultiColvarFunction::calculate(){
  if( checkNumericalDerivatives() ) myinputdata.recalculateBaseColvars( this );
  // Setup the link cells
  setupLinkCells(); 
  // And run all tasks
  runAllTasks();
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  // Construct matrix to store numerical derivatives
  unsigned pstart=3*myinputdata.getTotalNumberOfAtoms(); 
//  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
//     BridgedMultiColvarFunction* bb=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[i] );
//     if( bb ){
//         BridgedMultiColvarFunction* bb2=dynamic_cast<BridgedMultiColvarFunction*>( bb->getPntrToMultiColvar() );
//         plumed_massert( !bb2, "double filtered multicolvars and NumericalDerivatives are not compatible" );
//         pstart+=3*(bb->getPntrToMultiColvar())->getNumberOfAtoms();
//     } else {
//        pstart+=3*mybasemulticolvars[i]->getNumberOfAtoms();
//     }
//  }
  Matrix<double> numder_store( getNumberOfComponents(), pstart + 9 );

  pstart=0; 
  for(unsigned i=0;i<myinputdata.getNumberOfBaseMultiColvars();++i){
     BridgedMultiColvarFunction* bb=dynamic_cast<BridgedMultiColvarFunction*>( myinputdata.getBaseColvar(i) );
     if( bb ){
        ( bb->getPntrToMultiColvar() )->calculateAtomicNumericalDerivatives( this, pstart );
        for(unsigned k=0;k<getNumberOfComponents();++k){
           Value* val=getPntrToComponent(k);
           for(unsigned j=0;j<3*(bb->getPntrToMultiColvar())->getNumberOfAtoms();++j){
              numder_store(k,pstart+j) = val->getDerivative(pstart + j);
           }
        }   
     } else {
        myinputdata.getBaseColvar(i)->calculateAtomicNumericalDerivatives( this, pstart );
        for(unsigned k=0;k<getNumberOfComponents();++k){
           Value* val=getPntrToComponent(k);
           for(unsigned j=0;j<3*myinputdata.getNumberOfAtoms(i);++j){
              numder_store(k,pstart+j) = val->getDerivative(pstart + j);
           }
        }
     }
     pstart += myinputdata.getNumberOfAtoms(i);
  }

  // Note numerical derivatives only work for virial if mybasemulticolvars.size()==1
  if( myinputdata.getNumberOfBaseMultiColvars()==1 ){
      for(unsigned k=0;k<getNumberOfComponents();++k){
         Value* val=getPntrToComponent(k);
         for(unsigned j=0;j<9;++j) numder_store(k,pstart+j) = val->getDerivative(pstart + j);
      }
  }

  // Now transfer to multicolvar
  clearDerivatives();
  for(unsigned j=0;j<getNumberOfComponents();++j){
     Value* val=getPntrToComponent(j);
     for(unsigned i=0;i<getNumberOfDerivatives();++i) val->addDerivative( i, numder_store(j,i) );
  }
}

void MultiColvarFunction::updateActiveAtoms( AtomValuePack& myatoms ) const {
  myatoms.updateDynamicList();
}

void MultiColvarFunction::getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder ) const {
  myinputdata.getVectorDerivatives( ind, normed, myder );
}

void MultiColvarFunction::mergeVectorDerivatives( const unsigned& ival, const unsigned& start, const unsigned& end, 
                                                  const unsigned& jatom, const std::vector<double>& der, 
                                                  MultiValue& myder, AtomValuePack& myatoms ) const {
  myinputdata.mergeVectorDerivatives( ival, start, end, jatom, der, myder, myatoms );
}

}
}

