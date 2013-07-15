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
#include "MultiColvarBase.h"
#include "vesselbase/Vessel.h"
#include "tools/Pbc.h"
#include <vector>
#include <string>

using namespace std;

namespace PLMD{
namespace multicolvar{

void MultiColvarBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  ActionWithVessel::registerKeywords( keys );
  keys.use("NL_TOL");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                "that contributed less than TOL at the previous neighbor list update step are ignored.");
} 

MultiColvarBase::MultiColvarBase(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
usepbc(false),
updateFreq(0),
lastUpdate(0)
{
  if( keywords.exists("NOPBC") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  } 
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0) log.printf("  Updating contributors every %d steps.\n",updateFreq);
  else log.printf("  Updating contributors every step.\n");
}

void MultiColvarBase::addColvar( const std::vector<unsigned>& newatoms ){
  DynamicList<unsigned> newlist; newlist.setupMPICommunication( comm );
  for(unsigned i=0;i<newatoms.size();++i) newlist.addIndexToList( newatoms[i] );
  taskList.addIndexToList( colvar_atoms.size() );
  colvar_atoms.push_back( newlist );
}

void MultiColvarBase::copyAtomListToFunction( MultiColvarBase* myfunction ){
  for(unsigned i=0;i<all_atoms.fullSize();++i) myfunction->all_atoms.addIndexToList( all_atoms(i) );
}

void MultiColvarBase::copyActiveAtomsToFunction( MultiColvarBase* myfunction ){
  plumed_dbg_assert( myfunction->all_atoms.fullSize()==all_atoms.fullSize() );
  myfunction->all_atoms.deactivateAll();
  for(unsigned i=0;i<all_atoms.getNumberActive();++i){
      unsigned iatom=all_atoms.linkIndex( i );
      myfunction->all_atoms.activate( iatom );
  }
  myfunction->all_atoms.updateActiveMembers();
}

void MultiColvarBase::setupMultiColvarBase(){
  // Activate everything
  taskList.activateAll();
  for(unsigned i=0;i<colvar_atoms.size();++i) colvar_atoms[i].activateAll();
  // Resize stuff in derived classes 
  resizeDynamicArrays();
  // Resize local arrays
  resizeLocalArrays();
  // And resize the local vessels
  resizeFunctions();
}

void MultiColvarBase::requestAtoms(){
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
} 

void MultiColvarBase::prepare(){
  bool updatetime=false;
  if( contributorsAreUnlocked ){
      taskList.mpi_gatherActiveMembers( comm );
      mpi_gatherActiveMembers( comm, colvar_atoms ); 
      lockContributors(); updatetime=true;
  }
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
      taskList.activateAll(); 
      for(unsigned i=0;i<taskList.getNumberActive();++i) colvar_atoms[i].activateAll();
      unlockContributors(); updatetime=true; lastUpdate=getStep();
  }
  if(updatetime){
     // Resize stuff in derived classes
     resizeDynamicArrays();
     // Resize local arrays
     resizeLocalArrays();
     // Resize vessels
     resizeFunctions(); 
  }
}

void MultiColvarBase::resizeLocalArrays(){
  atoms_with_derivatives.clear(); 
  for(unsigned i=0;i<getNumberOfAtoms();++i) atoms_with_derivatives.addIndexToList( i );
  atoms_with_derivatives.deactivateAll();
  // Set up stuff for central atoms
  atomsWithCatomDer.clear();
  for(unsigned i=0;i<getNumberOfAtoms();++i) atomsWithCatomDer.addIndexToList( i );
  atomsWithCatomDer.deactivateAll();
  // Resize tempory forces array
  forcesToApply.resize( getNumberOfDerivatives() );
}

void MultiColvarBase::performTask( const unsigned& j ){
  // Currently no atoms have derivatives so deactivate those that are active
  atoms_with_derivatives.deactivateAll();
  // Currently no central atoms have derivatives so deactive them all
  atomsWithCatomDer.deactivateAll();

  // Do nothing if there are no active atoms in the colvar
  if( colvar_atoms[current].getNumberActive()==0 ){  
     setElementValue(1,0.0);
     return;                      
  }
  // Do a quick check on the size of this contribution  
  calculateWeight();
  if( getElementValue(1)<getTolerance() ){
     updateActiveAtoms();
     return;   
  }

  // Compute everything
  double vv=doCalculation( j );
  // Set the value of this element in ActionWithVessel
  setElementValue( 0, vv );
  return;
}

double MultiColvarBase::doCalculation( const unsigned& j ){
  double val=compute(j); updateActiveAtoms();
  return val;
}

Vector MultiColvarBase::retrieveCentralAtomPos(){
  ibox=getPbc().getInvBox().transpose();

  if( atomsWithCatomDer.getNumberActive()==0 ){
      Vector cvec = calculateCentralAtomPosition();
      for(unsigned i=0;i<3;++i) setElementValue( 2+i, cvec[i] );
      return cvec;
  }
  Vector cvec; 
  for(unsigned i=0;i<3;++i) cvec[i]=getElementValue( 2+i );
  return cvec;
}

void MultiColvarBase::addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der ){
  plumed_dbg_assert( iatom<getNumberOfAtoms() );
  atomsWithCatomDer.activate(iatom);
  unsigned nder = 3*getNumberOfAtoms() + 9;
  for(unsigned i=0;i<3;++i){ 
    for(unsigned j=0;j<3;++j){
        addElementDerivative( (2+j)*nder + 3*iatom + i, der(j,i) );         
     }
  }
}

double MultiColvarBase::getCentralAtomDerivative( const unsigned& iatom, const unsigned& jcomp, const Vector& df ){
  plumed_dbg_assert( atomsWithCatomDer.isActive(iatom) && jcomp<3 );
  unsigned nder = 3*getNumberOfAtoms() + 9;
  return df[0]*getElementDerivative( 2*nder + 3*iatom + jcomp ) + 
         df[1]*getElementDerivative( 3*nder + 3*iatom + jcomp ) + 
         df[2]*getElementDerivative( 4*nder + 3*iatom + jcomp ); 
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvarBase::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
  indices[jstore]=3*atoms_with_derivatives.getNumberActive() + 9;
  if( indices[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");

  unsigned kder = ntotal + jstore*maxder;
  for(unsigned jder=0;jder<atoms_with_derivatives.getNumberActive();++jder){
     unsigned iatom = 3*atoms_with_derivatives[jder];
     for(unsigned icomp=0;icomp<3;++icomp){ indices[ kder ] = iatom+icomp; kder++; }
  }
  unsigned nbase = 3*getNumberOfAtoms(); 
  for(unsigned icomp=0;icomp<9;++icomp){ indices[ kder ] = nbase + icomp; kder++; }   
}   

void MultiColvarBase::getCentralAtomIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ) const {
  indices[jstore]=3*atomsWithCatomDer.getNumberActive();
  if( indices[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");

  unsigned kder  = ntotal + jstore*maxder;
  for(unsigned jder=0;jder<atomsWithCatomDer.getNumberActive();++jder){
     unsigned iatom = 3*atomsWithCatomDer[jder];
     for(unsigned icomp=0;icomp<3;++icomp){ indices[ kder ] = iatom+icomp; kder++; }
  }
}

void MultiColvarBase::activateIndexes( const unsigned& istart, const unsigned& number, const std::vector<unsigned>& indexes ){
  for(unsigned i=0;i<number-9;i+=3){
      plumed_dbg_assert( indexes[istart+i]%3==0 ); unsigned iatom=indexes[istart+i]/3; 
      atoms_with_derivatives.activate( iatom ); 
  }
}

void MultiColvarBase::quotientRule( const unsigned& uder, const unsigned& vder, const unsigned& iout ){
  unsigned ustart=uder*getNumberOfDerivatives();
  unsigned vstart=vder*getNumberOfDerivatives();
  unsigned istart=iout*getNumberOfDerivatives();
  double weight = getElementValue( vder ), pref = getElementValue( uder ) / (weight*weight);
  for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
      unsigned n=3*atoms_with_derivatives[i], nx=n, ny=n+1, nz=n+2;
      setElementDerivative( istart + nx, getElementDerivative(ustart+nx) / weight - pref*getElementDerivative(vstart+nx) );
      setElementDerivative( istart + ny, getElementDerivative(ustart+ny) / weight - pref*getElementDerivative(vstart+ny) );
      setElementDerivative( istart + nz, getElementDerivative(ustart+nz) / weight - pref*getElementDerivative(vstart+nz) );
  }
  unsigned vbase=3*getNumberOfAtoms();
  for(unsigned i=0;i<9;++i){ 
      setElementDerivative( istart + vbase + i, getElementDerivative(ustart+vbase+i) / weight - pref*getElementDerivative(vstart+vbase+i) );
  }
  thisval_wasset[iout]=false; setElementValue( iout, getElementValue(uder) / weight );
}

void MultiColvarBase::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
     unsigned iatom=3*atoms_with_derivatives[i];
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) );
  }
  unsigned nvir=3*getNumberOfAtoms();
  for(unsigned j=0;j<9;++j){
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
  }
}

void MultiColvarBase::clearDerivativesAfterTask( const unsigned& ider ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  thisval_wasset[ider]=false; setElementValue( ider, 0.0 );
  thisval_wasset[ider]=false;
  if( ider>1 && ider<5 ){
     for(unsigned i=0;i<atomsWithCatomDer.getNumberActive();++i){
        unsigned iatom=vstart+3*atomsWithCatomDer[i];
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 );
     }  
  } else {
     for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
        unsigned iatom=vstart+3*atoms_with_derivatives[i];
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 );
     }   
     unsigned nvir=vstart+3*getNumberOfAtoms();
     for(unsigned j=0;j<9;++j){
        setElementDerivative( nvir, 0.0 ); nvir++;
     }
  }
}

void MultiColvarBase::apply(){
  if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

StoreCentralAtomsVessel* MultiColvarBase::getCentralAtoms(){
  // Look to see if vectors have already been created
  StoreCentralAtomsVessel* mycatoms;
  for(unsigned i=0;i<getNumberOfVessels();++i){
     mycatoms=dynamic_cast<StoreCentralAtomsVessel*>( getPntrToVessel(i) );
     if( mycatoms ) return mycatoms;
  }

  // Create the vessel
  vesselbase::VesselOptions da("","",0,"",this); 
  StoreCentralAtomsVessel* sv=new StoreCentralAtomsVessel(da);
  addVessel(sv); resizeFunctions(); // This makes sure resizing of vessels is done
  return sv;
}
     
}
}
