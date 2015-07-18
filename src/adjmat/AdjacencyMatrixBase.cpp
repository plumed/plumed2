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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/BridgedMultiColvarFunction.h"
#include "multicolvar/AtomValuePack.h"
#include "multicolvar/CatomPack.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixBase::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","WTOL","0.0","if the base multicolvars have weights then you must define a hard cutoff on those you want to consider explicitally");
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Weight of this does have derivatives
  weightHasDerivatives=true;

  // Create the storeAdjacencyMatrixVessel
  std::string param; vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; AdjacencyMatrixVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  mat = new AdjacencyMatrixVessel(da2);
  // Set a cutoff for clustering
  mat->setHardCutoffOnWeight( getTolerance() );   
  // Add the vessel to the base
  addVessel( mat );
}

void AdjacencyMatrixBase::parseAtomList(const std::string& key, const int& num, const bool& isnodes, std::vector<AtomNumber>& t){
  std::string newkey; t.resize(0);
  if( num<0 ){  
     newkey=key;
  } else {
     std::string snum;
     Tools::convert( num, snum ); newkey=key+snum;
  }

  if( isnodes ){
     std::vector<std::string> mlabs; parseVector(newkey,mlabs);
     if( mlabs.size()==0 ) return;

     std::string mname; 
     for(unsigned i=0;i<mlabs.size();++i){
        multicolvar::MultiColvarBase* mycolv = plumed.getActionSet().selectWithLabel<multicolvar::MultiColvarBase*>(mlabs[i]);
        if( i==0 && !mycolv ) break;
        if( !mycolv ) error("action labeled " + mlabs[i] + " does not exist or is not a multicolvar");

        if( i==0 ){
           mname = mycolv->getName();
           if( mycolv->isPeriodic() ) error("multicolvar functions don't work with this multicolvar");
           log.printf("  using colvars calculated by actions %s",mlabs[i].c_str() );
        } else {
           if( mname!=mycolv->getName() ) error("All input multicolvars must be of same type"); 
           log.printf("%s ", mlabs[i].c_str() );
        }
        // And store the multicolvar base
        mat->mybasemulticolvars.push_back( mycolv );
        // And track which variable stores each colvar
        for(unsigned j=0;j<mycolv->getFullNumberOfTasks();++j) mat->colvar_label.push_back( i );
     } 
     if( mname.length()>0 ){
        // Build objects to calculate data in the underlying multicolvars
        log.printf("\n"); double wtolerance; parse("WTOL",wtolerance);
        log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
        mat->buildRemoteDataStashes( wtolerance );
        // Retrieve the atoms involved in the underlying multicolvars
        for(unsigned i=0;i<mat->mybasemulticolvars.size();++i){
            std::vector<AtomNumber> tmp_atoms; 
            multicolvar::BridgedMultiColvarFunction* mybr=dynamic_cast<multicolvar::BridgedMultiColvarFunction*>( mat->mybasemulticolvars[i] );
            if( mybr ) tmp_atoms=(mybr->getPntrToMultiColvar())->getAbsoluteIndexes();
            else tmp_atoms=(mat->mybasemulticolvars[i])->getAbsoluteIndexes();
            for(unsigned j=0;j<tmp_atoms.size();++j) t.push_back( tmp_atoms[j] );
        }
     } else {
        // This allows you just to read in atom numbers
        ActionAtomistic::interpretAtomList( mlabs, t );
        for(unsigned j=0;j<t.size();++j) mat->colvar_label.push_back( -1 );
     }
  } else {
     ActionAtomistic::parseAtomList( key, num, t );
  }
}

// Would like to get rid of this ideally and use the one in MultiColvarBase
void AdjacencyMatrixBase::addTaskToList( const unsigned& taskCode ){
  ActionWithVessel::addTaskToList( taskCode );
}

unsigned AdjacencyMatrixBase::getNumberOfNodeTypes() const {
  if( mat->mybasemulticolvars.size()==0 ) return 1;
  return mat->mybasemulticolvars.size();
}

unsigned AdjacencyMatrixBase::getNumberOfNodes() const {
  unsigned nb=0; for(unsigned i=0;i<(mat->mybasemulticolvars).size();++i) nb += (mat->mybasemulticolvars[i])->getFullNumberOfTasks();
  return nb;
}

void AdjacencyMatrixBase::requestAtoms( const std::vector<AtomNumber>& atoms ){
  ActionAtomistic::requestAtoms( atoms );
  for(unsigned i=0;i<mat->mybasemulticolvars.size();++i) addDependency(mat->mybasemulticolvars[i]);
  setupMultiColvarBase();
}

void AdjacencyMatrixBase::doJobsRequiredBeforeTaskList(){
  // Do jobs required by ActionWithVessel
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  // Dont calculate derivatives on first loop
  if( usingLowMem() ) dertime=false;
  else dertime=true;
}

void AdjacencyMatrixBase::calculate(){
  printf("IN ADJACENCY BASE CALCULATE \n");
  if( checkNumericalDerivatives() ) error("numerical derivatives currently broken");
  // Setup the linke cells
  setupLinkCells();
  // And run all tasks
  runAllTasks();
}

Vector AdjacencyMatrixBase::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  if( iatom>mat->colvar_label.size() ) return getPosition( iatom );
  if( mat->colvar_label[iatom]<0 ) return getPosition( iatom );
  return mat->getCentralAtomPos( iatom );
}

void AdjacencyMatrixBase::updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const {
  myatoms.updateDynamicList();
}

bool AdjacencyMatrixBase::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  if( code>mat->colvar_label.size() ) return true;
  return mat->isCurrentlyActive( code );
}

void AdjacencyMatrixBase::addAtomDerivatives( const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned jatom=myatoms.getIndex(iatom);

  if( jatom>mat->colvar_label.size() ){
      myatoms.addAtomsDerivatives( 1, iatom, der );
  } else if( mat->colvar_label[jatom]<0 ){
      myatoms.addAtomsDerivatives( 1, iatom, der );
  } else {
      unsigned mmc=mat->colvar_label[jatom];
      unsigned basen=0; for(unsigned i=0;i<mmc;++i) basen+=(mat->mybasemulticolvars[i])->getNumberOfAtoms();
      multicolvar::CatomPack atom0=(mat->mybasemulticolvars[mmc])->getCentralAtomPack( basen, mat->convertToLocalIndex(jatom,mmc) );
      myatoms.addComDerivatives( 1, der, atom0 );
  }
}

}
}
