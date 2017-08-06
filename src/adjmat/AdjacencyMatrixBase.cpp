/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","");
  keys.add("atoms","GROUPB","");
  keys.reserve("atoms","GROUPC","");
  keys.addFlag("COMPONENTS",false,"also calculate the components of the vector connecting the atoms in the contact matrix");
  keys.addFlag("NOPBC",false,"don't use pbc");
  keys.addOutputComponent("w","default","the weight of the connection");
  keys.addOutputComponent("x","COMPONENT","the projection of the bond on the x axis");
  keys.addOutputComponent("y","COMPONENT","the projection of the bond on the y axis");
  keys.addOutputComponent("z","COMPONENT","the projection of the bond on the z axis");
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  linkcells(comm)
{
  std::vector<unsigned> shape(2); std::vector<AtomNumber> t; parseAtomList("GROUP", t );
  if( t.size()==0 ){
      std::vector<AtomNumber> ta; parseAtomList("GROUPA",ta); 
      std::vector<AtomNumber> tb; parseAtomList("GROUPB",tb);
      if( ta.size()==0 || tb.size()==0 ) error("no atoms have been specified in input");

      // Create list of tasks
      log.printf("  atoms in GROUPA ");
      for(unsigned i=0;i<ta.size();++i){ log.printf("%d ", ta[i].serial()); addTaskToList( i ); t.push_back(ta[i]); }
      log.printf("\n");
      log.printf("  atoms in GROUPB "); ablocks.resize( tb.size() ); unsigned n=0;
      for(unsigned i=0;i<tb.size();++i){ 
          log.printf("%d ", tb[i].serial()); 
          bool found=false; unsigned inum=0;
          for(unsigned j=0; j<ta.size(); ++j) {
              if( ta[j]==tb[i] ){ found=true; inum=j; break; }
          }
          if( found ){ ablocks[i]=inum; }
          else { ablocks[i]=ta.size()+n; t.push_back(tb[i]); n++; } 
      }
      log.printf("\n"); shape[0]=ta.size(); shape[1]=tb.size(); 
      // Create a group of the atoms in the rows 
      atoms.insertGroup( getLabel(), ta );
  } else {
      // Create list of tasks 
      log.printf("  atoms in GROUP "); ablocks.resize( t.size() );
      for(unsigned i=0;i<t.size();++i){ log.printf("%d ", t[i].serial()); addTaskToList( i ); ablocks[i]=i; }
      log.printf("\n"); shape[0]=shape[1]=t.size();
      // Create a group of the atoms in the rows
      atoms.insertGroup( getLabel(), t );
  }
  if( keywords.exists("GROUPC") ){
      std::vector<AtomNumber> tc; parseAtomList("GROUPC",tc);
      if( tc.size()==0 ) error("no atoms in GROUPC specified");
      log.printf("  atoms in GROUPC "); 
      for(unsigned i=0;i<tc.size();++i){ log.printf("%d ", tc[i].serial()); t.push_back(tc[i]); }
      log.printf("\n");
  }
  // Request the atoms from the ActionAtomistic
  requestAtoms( t ); parseFlag("COMPONENTS",components); parseFlag("NOPBC",nopbc);
  addComponentWithDerivatives( "w", shape ); componentIsNotPeriodic("w"); 
  if( components ){ 
     addComponentWithDerivatives( "x", shape ); componentIsNotPeriodic("x");
     addComponentWithDerivatives( "y", shape ); componentIsNotPeriodic("y");
     addComponentWithDerivatives( "z", shape ); componentIsNotPeriodic("z");
  }
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
}

void AdjacencyMatrixBase::setLinkCellCutoff( const double& lcut, double tcut ) {
  linkcells.setCutoff( lcut );
}

void AdjacencyMatrixBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  // Need to retrieve atoms here in case this is done in a stream
  if( actionInChain() ) retrieveAtoms();
  // Make sure all tasks are done
  tflags.assign(tflags.size(),1);
  // Build link cells here so that this is done in stream if it needed in stream
  std::vector<Vector> ltmp_pos( ablocks.size() );
  for(unsigned i=0; i<ablocks.size(); ++i) {
      ltmp_pos[i]=getPosition( ablocks[i] );
  }
  linkcells.buildCellLists( ltmp_pos, ablocks, getPbc() );
}

void AdjacencyMatrixBase::calculate(){
  // Now run all the tasks
  runAllTasks();
}

void AdjacencyMatrixBase::updateWeightDerivativeIndices( const unsigned& sno, const std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned w_ind = getPntrToOutput(0)->getPositionInStream();
  // Update dynamic list indices for central atom
  myvals.updateIndex( w_ind, 3*indices[0]+0 ); myvals.updateIndex( w_ind, 3*indices[0]+1 ); myvals.updateIndex( w_ind, 3*indices[0]+2 );
  // Update dynamic list indices for atom forming this bond
  myvals.updateIndex( w_ind, 3*indices[sno]+0 ); myvals.updateIndex( w_ind, 3*indices[sno]+1 ); myvals.updateIndex( w_ind, 3*indices[sno]+2 );
  // Update dynamic list indices for virial
  unsigned base = 3*getNumberOfAtoms(); for(unsigned j=0;j<9;++j) myvals.updateIndex( w_ind, base+j );
}

void AdjacencyMatrixBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Retrieve cells required from link cells
  std::vector<unsigned> cells_required( linkcells.getNumberOfCells() ); unsigned ncells_required=0;
  linkcells.addRequiredCells( linkcells.findMyCell( getPosition(current) ), ncells_required, cells_required );

  // Now retrieve bookeeping arrays
  std::vector<Vector> & atoms( myvals.getAtomVector() );
  std::vector<unsigned> & indices( myvals.getIndices() ); 
  if( indices.size()!=(1+ablocks.size()) ){
     indices.resize( 1+ablocks.size() ); atoms.resize( 1+ablocks.size() );
  }

  // Now get the positions
  unsigned natoms=1; indices[0]=current;
  linkcells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  myvals.setNumberOfIndices( natoms );

  // Apply periodic boundary conditions
  for(unsigned i=0; i<natoms; ++i) atoms[i] = getPosition(indices[i]) - getPosition(current);
  if( !nopbc ) pbcApply( atoms, natoms );

  // Now loop over all atoms in coordination sphere
  MatrixElementPack mypack( myvals, components, this ); 
  mypack.setIndex( 0, current ); Vector zero; zero.zero(); mypack.setPosition( 0, zero );
  for(unsigned i=1;i<natoms;++i){
      mypack.setIndex( 1, indices[i] ); mypack.setPosition( 1, atoms[i] );
      double weight = calculateWeight( mypack ); 
      if( weight<epsilon ){
          if( !doNotCalculateDerivatives() ) {
              updateWeightDerivativeIndices( i, indices, myvals ); clearMatrixElements( myvals );
          }
          continue;
      }
      // Now set the matrix weight and the vector if required  
      myvals.setValue( getPntrToOutput(0)->getPositionInStream(), weight ); 
      if( components ) {
          unsigned x_index = getPntrToOutput(1)->getPositionInStream();
          unsigned y_index = getPntrToOutput(2)->getPositionInStream();
          unsigned z_index = getPntrToOutput(3)->getPositionInStream();
          myvals.setValue( x_index, atoms[i][0] ); myvals.setValue( y_index, atoms[i][1] ); myvals.setValue( z_index, atoms[i][2] );
          if( !doNotCalculateDerivatives() ){
              myvals.addDerivative( x_index, 3*indices[0]+0, -1 ); myvals.addDerivative( x_index, 3*indices[i]+0, +1 );
              myvals.addDerivative( y_index, 3*indices[0]+1, -1 ); myvals.addDerivative( y_index, 3*indices[i]+1, +1 );
              myvals.addDerivative( z_index, 3*indices[0]+2, -1 ); myvals.addDerivative( z_index, 3*indices[i]+2, +1 );
              // Update dynamic lists for central atom
              myvals.updateIndex( x_index, 3*indices[0]+0 ); myvals.updateIndex( y_index, 3*indices[0]+1 ); myvals.updateIndex( z_index, 3*indices[0]+2 );
              // Update dynamic lists for bonded atom
              myvals.updateIndex( x_index, 3*indices[i]+0 ); myvals.updateIndex( y_index, 3*indices[i]+1 ); myvals.updateIndex( z_index, 3*indices[i]+2 );
              // Add derivatives of virial
              unsigned base = 3*getNumberOfAtoms();
              // Virial for x
              myvals.addDerivative( x_index, base+0, -atoms[i][0] ); myvals.addDerivative( x_index, base+3, -atoms[i][1] ); myvals.addDerivative( x_index, base+6, -atoms[i][2] ); 
              myvals.updateIndex( x_index, base+0 ); myvals.updateIndex( x_index, base+3 ); myvals.updateIndex( x_index, base+6 );
              // Virial for y
              myvals.addDerivative( y_index, base+1, -atoms[i][0] ); myvals.addDerivative( y_index, base+4, -atoms[i][1] ); myvals.addDerivative( y_index, base+7, -atoms[i][2] );  
              myvals.updateIndex( y_index, base+1 ); myvals.updateIndex( y_index, base+4 ); myvals.updateIndex( y_index, base+7 );
              // Virial for z
              myvals.addDerivative( z_index, base+2, -atoms[i][0] ); myvals.addDerivative( z_index, base+5, -atoms[i][1] ); myvals.addDerivative( z_index, base+8, -atoms[i][2] );  
              myvals.updateIndex( z_index, base+2 ); myvals.updateIndex( z_index, base+5 ); myvals.updateIndex( z_index, base+8 );
          }
      } 
      // Update derivatives
      if( !doNotCalculateDerivatives() ) updateWeightDerivativeIndices( i, indices, myvals );
      // This does everything in the stream that is done with single matrix elements 
      runTask( myvals.getTaskIndex(), current, indices[i], myvals );
      // Now clear only elements that are not accumulated over whole row
      clearMatrixElements( myvals );
  }
}

void AdjacencyMatrixBase::apply() {

}

}
}
