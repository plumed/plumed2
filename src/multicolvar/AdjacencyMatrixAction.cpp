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
#include "AdjacencyMatrixAction.h"

namespace PLMD {
namespace multicolvar {

void AdjacencyMatrixAction::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys );
  keys.reserveFlag("USE_ORIENTATION",false,"When computing whether two atoms/molecules are adjacent also take their orientations into account");
  keys.add("numbered","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
}

AdjacencyMatrixAction::AdjacencyMatrixAction(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao),
tmpdf(1)
{
  if( keywords.exists("USE_ORIENTATION") ) parseFlag("USE_ORIENTATION",use_orient);
  else use_orient=true;
  // Weight of this does have derivatives
  weightHasDerivatives=true;
  // Read in the switching function
  if( getNumberOfBaseMultiColvars()==1 ){
      switchingFunction.resize(1,1);
      std::string sw, errors; parse("SWITCH",sw);
      if(sw.length()==0) error("missing SWITCH keyword");
      switchingFunction(0,0).set(sw,errors);
      log.printf("  constructing adjacency matrix between atoms that are within %s\n", ( switchingFunction(0,0).description() ).c_str() );
  } else {
      unsigned nfunc=getNumberOfBaseMultiColvars();
      switchingFunction.resize( nfunc,nfunc ); 
      for(unsigned i=0;i<nfunc;++i){
          // Retrieve the base number 
          unsigned ibase;
          if( nfunc<10 ){ 
             ibase=(i+1)*10;
          } else if ( nfunc<100 ){
             ibase=(i+1)*100;
          } else {
             error("wow this is an error I never would have expected");  
          }

          for(unsigned j=i;j<nfunc;++j){
             std::string sw, errors; parseNumbered("SWITCH",ibase+j+1,sw);
             if(sw.length()==0){
                std::string num; Tools::convert(ibase+j+1,num);
                error("could not find SWITCH" + num + " keyword. Need one SWITCH keyword for each distinct base-multicolvar-pair type");
             }
             switchingFunction(j,i).set(sw,errors);  
             if( j!=i) switchingFunction(i,j).set(sw,errors);
             log.printf("  %d th and %d th multicolvar groups must be within %s\n",i+1,j+1,(switchingFunction(i,j).description()).c_str() );      
          }
      }
  }

  // Build atom lists
  buildAtomListWithPairs( true );
  // Build active elements array
  for(unsigned i=0;i<getFullNumberOfTasks();++i) active_elements.addIndexToList( i );
  active_elements.setupMPICommunication( comm );

  // Find the largest sf cutoff
  double sfmax=switchingFunction(0,0).get_dmax();
  for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
      for(unsigned j=0;j<<getNumberOfBaseMultiColvars();++j){
          double tsf=switchingFunction(i,j).get_dmax(); 
          if( tsf>sfmax ) sfmax=tsf;
      }
  }
  // And set the link cell cutoff
  setLinkCellCutoff( sfmax );

  // Create the storeAdjacencyMatrixVessel
  std::string param; vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; AdjacencyMatrixVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  mat = new AdjacencyMatrixVessel(da2);
  // Add the vessel to the base
  addVessel( mat );
  // And resize everything
  resizeFunctions();

  // One component for regular multicolvar and nelements for vectormulticolvar
  unsigned ncomp=getBaseMultiColvar(0)->getNumberOfQuantities() - 5;;
  orient0.resize( ncomp ); orient1.resize( ncomp );
}

void AdjacencyMatrixAction::doJobsRequiredBeforeTaskList(){
  // Do jobs required by ActionWithVessel
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  // Make it possible only to go through loop once
  gathered=false; active_elements.deactivateAll();
  // Dont calculate derivatives on first loop
  if( usingLowMem() ) dertime=false;
  else dertime=true;
}

void AdjacencyMatrixAction::calculateWeight(){
  Vector distance = getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(1) );
  double dfunc, sw = switchingFunction( getBaseColvarNumber(0),getBaseColvarNumber(1) ).calculate( distance.modulo(), dfunc );
  setElementValue(1,sw);
}

double AdjacencyMatrixAction::compute(){
  active_elements.activate( getCurrentPositionInTaskList() );

  double f_dot, dot_df;
  if( use_orient ){
      getVectorForBaseTask( 0, orient0 );
      getVectorForBaseTask( 1, orient1 );   

      double dot; dot=0;
      for(unsigned k=0;k<orient0.size();++k) dot+=orient0[k]*orient1[k];
      f_dot=0.5*( 1 + dot ); dot_df=0.5;
      // Add smac stuff here if required
  } else {
      f_dot=1.0; dot_df=0.0;
  }

  // Retrieve the weight of the connection
  double weight = getElementValue(1); 

  if( dertime && !doNotCalculateDerivatives() ){
     // Add contribution due to separation between atoms
     Vector distance = getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(1) );
     double dfunc, sw = switchingFunction( getBaseColvarNumber(0), getBaseColvarNumber(1) ).calculate( distance.modulo(), dfunc ); 
     addCentralAtomsDerivatives( 0, 0, (-dfunc)*f_dot*distance );
     addCentralAtomsDerivatives( 1, 0, (dfunc)*f_dot*distance );
     MultiColvarBase::addBoxDerivatives( 0, (-dfunc)*f_dot*Tensor(distance,distance) );

     // And derivatives of orientation
     if( use_orient ){
        for(unsigned k=0;k<orient0.size();++k){
           orient0[k]*=sw*dot_df; orient1[k]*=sw*dot_df;
        }
        addOrientationDerivatives( 0, orient1 );
        addOrientationDerivatives( 1, orient0 );
     }
  }
  return weight*f_dot;
}

void AdjacencyMatrixAction::setMatrixIndexesForTask( const unsigned& ii ){
  unsigned icolv = active_elements[ii], tcode = getTaskCode( icolv );
  bool check = MultiColvarBase::setupCurrentAtomList( tcode );
  plumed_assert( check );
}

void AdjacencyMatrixAction::retrieveMatrix( Matrix<double>& mymatrix ){
  // Gather active elements in matrix
  if(!gathered) active_elements.mpi_gatherActiveMembers( comm ); 
  gathered=true;
  
  for(unsigned i=0;i<active_elements.getNumberActive();++i){
      setMatrixIndexesForTask( i );
      unsigned j = current_atoms[1], k = current_atoms[0];
      mymatrix(k,j)=mymatrix(j,k)=getMatrixElement( i );
  }
}

void AdjacencyMatrixAction::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ){
  plumed_dbg_assert( nneigh.size()==getFullNumberOfBaseTasks() && adj_list.nrows()==getFullNumberOfBaseTasks() && 
                       adj_list.ncols()==getFullNumberOfBaseTasks() );
  // Gather active elements in matrix
  if(!gathered) active_elements.mpi_gatherActiveMembers( comm );
  gathered=true;

  // Currently everything has zero neighbors
  for(unsigned i=0;i<nneigh.size();++i) nneigh[i]=0;

  // And set up the adjacency list
  for(unsigned i=0;i<active_elements.getNumberActive();++i){
      setMatrixIndexesForTask( i );
      unsigned j = current_atoms[1], k = current_atoms[0];
      adj_list(k,nneigh[k])=j; nneigh[k]++;
      adj_list(j,nneigh[j])=k; nneigh[j]++;
  } 
} 

void AdjacencyMatrixAction::addDerivativesOnMatrixElement( const unsigned& ielem, const unsigned& jrow, const double& df, Matrix<double>& der ){
  plumed_dbg_assert( ielem<active_elements.getNumberActive() );
  tmpdf[0]=df; unsigned jelem=active_elements[ielem];

  if( usingLowMem() ){
     mat->recompute( jelem, 0 ); mat->chainRule( 0, tmpdf );
     for(unsigned i=0;i<mat->getNumberOfDerivatives(0);++i) der( jrow, mat->getStoredIndex(0,i) ) += mat->getFinalDerivative(i);
  } else {
     mat->chainRule( jelem, tmpdf );
     for(unsigned i=0;i<mat->getNumberOfDerivatives(jelem);++i) der( jrow, mat->getStoredIndex(jelem,i) ) += mat->getFinalDerivative(i);
  }
}

}
}
