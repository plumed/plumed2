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
  use_orient=false;
  if( keywords.exists("USE_ORIENTATION") ) parseFlag("USE_ORIENTATION",use_orient);
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

  if( use_orient && getBaseMultiColvar(0)->getNumberOfQuantities()<3 ) error("using orientation but no orientations in base colvars"); 

  // Build active elements array
  for(unsigned i=0;i<getFullNumberOfTasks();++i) active_elements.addIndexToList( i );

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
  // Set a cutoff for clustering
  mat->setHardCutoffOnWeight( getTolerance() );
  // Add the vessel to the base
  addVessel( mat );
  // And resize everything
  resizeFunctions();
}

void AdjacencyMatrixAction::doJobsRequiredBeforeTaskList(){
  // Do jobs required by ActionWithVessel
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  // Make it possible only to go through loop once
  gathered=false; 
  // Dont calculate derivatives on first loop
  if( usingLowMem() ) dertime=false;
  else dertime=true;
}

void AdjacencyMatrixAction::calculateWeight( AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
  double dfunc, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ),getBaseColvarNumber( myatoms.getIndex(1) ) ).calculate( distance.modulo(), dfunc );
  myatoms.setValue(0,sw);
}

double AdjacencyMatrixAction::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
//  active_elements.activate( tindex );

  double f_dot, dot_df; 
  unsigned ncomp=getBaseMultiColvar(0)->getNumberOfQuantities();
  std::vector<double> orient0(ncomp), orient1(ncomp);
  if( use_orient ){
      getVectorForTask( myatoms.getIndex(0), true, orient0 );
      getVectorForTask( myatoms.getIndex(1), true, orient1 );      

      double dot; dot=0;
      for(unsigned k=2;k<orient0.size();++k) dot+=orient0[k]*orient1[k];
      f_dot=0.5*( 1 + dot ); dot_df=0.5;
      // Add smac stuff here if required
  } else {
      f_dot=1.0; dot_df=0.0;
  }

  // Retrieve the weight of the connection
  double weight = myatoms.getValue(0); 

  if( dertime && !doNotCalculateDerivatives() ){
     // Add contribution due to separation between atoms
     Vector distance = getSeparation( myatoms.getPosition(0), myatoms.getPosition(1) );
     double dfunc, sw = switchingFunction( getBaseColvarNumber( myatoms.getIndex(0) ), getBaseColvarNumber( myatoms.getIndex(0) ) ).calculate( distance.modulo(), dfunc ); 
     CatomPack atom0=getCentralAtomPackFromInput( myatoms.getIndex(0) );
     myatoms.addComDerivatives( 1, (-dfunc)*f_dot*distance, atom0 );
     CatomPack atom1=getCentralAtomPackFromInput( myatoms.getIndex(1) );

     myatoms.addComDerivatives( 1, (dfunc)*f_dot*distance, atom1 );
     myatoms.addBoxDerivatives( 1, (-dfunc)*f_dot*Tensor(distance,distance) );

     // And derivatives of orientation
     if( use_orient ){
        for(unsigned k=2;k<orient0.size();++k){
           orient0[k]*=sw*dot_df; orient1[k]*=sw*dot_df;
        }
        MultiValue myder0(0,0); getVectorDerivatives( myatoms.getIndex(0), true, myder0 );
        mergeVectorDerivatives( 1, 2, orient1.size(), myatoms.getIndex(0), orient1, myder0, myatoms );
        MultiValue myder1(0,0); getVectorDerivatives( myatoms.getIndex(1), true, myder1 );
        mergeVectorDerivatives( 1, 2, orient0.size(), myatoms.getIndex(1), orient0, myder1, myatoms );
     }
  }
  return weight*f_dot;
}

void AdjacencyMatrixAction::retrieveMatrix( Matrix<double>& mymatrix ){
  // Gather active elements in matrix
  if(!gathered){
     active_elements.deactivateAll();
     for(unsigned i=0;i<getFullNumberOfTasks();++i){
        if( mat->storedValueIsActive(i) ) active_elements.activate(i);
     }
     active_elements.updateActiveMembers(); gathered=true;
  }
 
  std::vector<unsigned> myatoms(2); std::vector<double> vals(2);
  for(unsigned i=0;i<active_elements.getNumberActive();++i){
      decodeIndexToAtoms( getTaskCode(active_elements[i]), myatoms ); 
      unsigned j = myatoms[1], k = myatoms[0];
      mat->retrieveValue( active_elements[i], false, vals );
      mymatrix(k,j)=mymatrix(j,k)=vals[1];   
  }
}

void AdjacencyMatrixAction::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ){
  plumed_dbg_assert( nneigh.size()==getFullNumberOfBaseTasks() && adj_list.nrows()==getFullNumberOfBaseTasks() && 
                       adj_list.ncols()==getFullNumberOfBaseTasks() );
  // Gather active elements in matrix
  if(!gathered){
     active_elements.deactivateAll();
     for(unsigned i=0;i<getFullNumberOfTasks();++i){
        if( mat->storedValueIsActive(i) ) active_elements.activate(i);
     }
     active_elements.updateActiveMembers(); gathered=true;
  }

  // Currently everything has zero neighbors
  for(unsigned i=0;i<nneigh.size();++i) nneigh[i]=0;

  // And set up the adjacency list
  std::vector<unsigned> myatoms(2);
  for(unsigned i=0;i<active_elements.getNumberActive();++i){
      decodeIndexToAtoms( getTaskCode(active_elements[i]), myatoms );
      unsigned j = myatoms[1], k = myatoms[0];
      adj_list(k,nneigh[k])=j; nneigh[k]++;
      adj_list(j,nneigh[j])=k; nneigh[j]++;
  } 
} 

}
}
