/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
MultiColvarBase(ao),
wtolerance(0.0)
{
  // Read in the arguments
  std::string mname; std::vector<std::string> mlabs; parseVector("DATA",mlabs);
  log.printf("  using colvars calculated by actions "); bool useweights=false;
  for(unsigned i=0;i<mlabs.size();++i){
      log.printf("%s ",mlabs[i].c_str() );
      MultiColvarBase* mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlabs[i]); 
      if(!mycolv) error("action labeled " + mlabs[i] + " does not exist or is not a multicolvar");
      // Check all base multicolvars are of same type
      if( i==0 ){ 
          mname = mycolv->getName();
          tvals.resize( mycolv->getNumberOfQuantities()-4 );
          if( mycolv->isPeriodic() ) error("multicolvar functions don't work with this multicolvar");
      } else {
          if( mname!=mycolv->getName() ) error("All input multicolvars must be of same type"); 
      }
      // Make sure we use low memory option in base colvar
      mycolv->setLowMemOption( usingLowMem() );
      // Add the dependency
      addDependency(mycolv);
      // Ensure weights are considered
      if( mycolv->weightHasDerivatives ) useweights=true;
      // And store the multicolvar base
      mybasemulticolvars.push_back( mycolv );
      // And track which variable stores each colvar
      for(unsigned j=0;j<mycolv->getFullNumberOfTasks();++j) colvar_label.push_back( i );
  }
  log.printf("\n");
  if( keywords.exists("WTOL") && useweights ){
      parse("WTOL",wtolerance);
      log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
  }
}

void MultiColvarFunction::setupAtomLists(){
  // Copy lists of atoms involved from base multicolvars 
  std::vector<AtomNumber> all_atoms, tmp_atoms;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      BridgedMultiColvarFunction* mybr=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[i] );
      if( mybr ) tmp_atoms=(mybr->mycolv)->getAbsoluteIndexes();
      else tmp_atoms=mybasemulticolvars[i]->getAbsoluteIndexes();
      for(unsigned j=0;j<tmp_atoms.size();++j) all_atoms.push_back( tmp_atoms[j] );
  }   
  
  // Now make sure we get all the atom positions 
  ActionAtomistic::requestAtoms( all_atoms );
  for(unsigned i=0;i<mybasemulticolvars.size();++i) addDependency(mybasemulticolvars[i]);
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase();
}

void MultiColvarFunction::buildSymmetryFunctionLists(){
  if( mybasemulticolvars.size()>2 ) error("Found too many multicolvars in DATA specification. You can use either 1 or 2");

  // Make sure information is stored in the required multicolvars
  if( keywords.exists("WTOL") ){
      for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes( true, wtolerance );
  } else {
      for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes( false, 0.0 ); 
  }

  usespecies=true; ablocks.resize( 1 );
  for(unsigned i=0;i<mybasemulticolvars[0]->getFullNumberOfTasks();++i) addTaskToList( i );

  unsigned ntotal=0;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      ntotal += mybasemulticolvars[i]->getFullNumberOfTasks();
  }
  unsigned k=0, start=0;
  current_atoms.resize( 1 + ntotal ); ablocks[0].resize( ntotal ); 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      for(unsigned j=0;j<mybasemulticolvars[i]->getFullNumberOfTasks();++j){
          ablocks[0][k]=start + j; k++;
      }
      start += mybasemulticolvars[i]->getFullNumberOfTasks();
  }  
  setupAtomLists();
}

void MultiColvarFunction::buildSets(){
  nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
     if( mybasemulticolvars[i]->getFullNumberOfTasks()!=nblock ){
          error("mismatch between numbers of tasks in various base multicolvars");
     }
     mybasemulticolvars[i]->buildDataStashes( false, 0.0 );
  }
  ablocks.resize( mybasemulticolvars.size() );
  usespecies=false; current_atoms.resize( mybasemulticolvars.size() );
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
  setupAtomLists();
}

void MultiColvarFunction::buildAtomListWithPairs( const bool& allow_intra_group ){
  if( !allow_intra_group && mybasemulticolvars.size()>2 ) error("only two input multicolvars allowed with this function"); 

  // Make sure information is stored in the required multicolvars
  if( keywords.exists("WTOL") ){
      for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes( true, wtolerance );
  } else {
      for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes( false, 0.0 );
  }
  
  usespecies=false; ablocks.resize(2); current_atoms.resize( 2 );
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
  setupAtomLists(); 
}

void MultiColvarFunction::calculate(){
  if( checkNumericalDerivatives() ){
     // Clear any derivatives in base colvar that were 
     // accumulated from previous calculations
     for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->clearDerivatives(); 
     // And recalculate
     for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->calculate();
     // Copy the box from the base multicolvar here
     unsigned maxb=mybasemulticolvars.size() - 1;
     changeBox( mybasemulticolvars[maxb]->getBox() );
  }
  setupLinkCells(); 
  // And run all tasks
  runAllTasks();
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  // Construct matrix to store numerical derivatives
  unsigned pstart=0; 
  for(unsigned i=0;i<mybasemulticolvars.size();++i) pstart+=3*mybasemulticolvars[i]->getNumberOfAtoms();
  Matrix<double> numder_store( getNumberOfComponents(), pstart + 9 );

  pstart=0; 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
     mybasemulticolvars[i]->calculateAtomicNumericalDerivatives( this, pstart );
     for(unsigned k=0;k<getNumberOfComponents();++k){
        Value* val=getPntrToComponent(k);
        for(unsigned j=0;j<3*mybasemulticolvars[i]->getNumberOfAtoms();++j){
           numder_store(k,pstart+j) = val->getDerivative(pstart + j);
        }
     }
     pstart += 3*mybasemulticolvars[i]->getNumberOfAtoms(); 
  }

  // Note numerical derivatives only work for virial if mybasemulticolvars.size()==1
  if( mybasemulticolvars.size()==1 ){
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

void MultiColvarFunction::addStoredDerivative( const unsigned& jout, const unsigned& base_cv_no, const unsigned& base_index, const double& der ){
  plumed_dbg_assert( jout<getNumberOfQuantities() && base_cv_no<mybasemulticolvars.size() && base_index<mybasemulticolvars[base_cv_no]->getNumberOfDerivatives() );

  unsigned mbase = 3*mybasemulticolvars[base_cv_no]->getSizeOfAtomsWithDerivatives(), tbase = 3*getNumberOfAtoms();
  if( base_index>=mbase ){
      // Add virial element
      unsigned jindex = base_index - mbase + tbase;
      addElementDerivative( jout*getNumberOfDerivatives() + jindex, der ); 
  } else {
      // Add atomic element
      unsigned offset=0; for(unsigned i=0;i<base_cv_no;++i) offset += 3*mybasemulticolvars[i]->getNumberOfAtoms();
      unsigned jindex = offset + base_index; 
      plumed_dbg_assert( jindex<3*getNumberOfAtoms() );
      addElementDerivative( jout*getNumberOfDerivatives() + jindex, der );
      unsigned iatom = ( jindex / 3 );
      atoms_with_derivatives.activate( iatom );
  }
}

void MultiColvarFunction::updateActiveAtoms(){
  if( atoms_with_derivatives.updateComplete() ) return;
  atoms_with_derivatives.updateActiveMembers();
}

Vector MultiColvarFunction::calculateCentralAtomPosition(){
  Vector catom=getCentralAtom();
  atomsWithCatomDer.updateActiveMembers();
  return catom;
}

}
}

