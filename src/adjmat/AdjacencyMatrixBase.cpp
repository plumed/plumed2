/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2018 The plumed team
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
  linkcells(comm),
  threecells(comm)
{
  std::vector<unsigned> shape(2); std::vector<AtomNumber> t; parseAtomList("GROUP", t );
  if( t.size()==0 ) {
    std::vector<AtomNumber> ta; parseAtomList("GROUPA",ta);
    std::vector<AtomNumber> tb; parseAtomList("GROUPB",tb);
    if( ta.size()==0 || tb.size()==0 ) error("no atoms have been specified in input");

    // Create list of tasks
    log.printf("  atoms in GROUPA ");
    for(unsigned i=0; i<ta.size(); ++i) { log.printf("%d ", ta[i].serial()); addTaskToList( i ); t.push_back(ta[i]); }
    log.printf("\n");
    log.printf("  atoms in GROUPB "); ablocks.resize( tb.size() ); unsigned n=0;
    for(unsigned i=0; i<tb.size(); ++i) {
      log.printf("%d ", tb[i].serial());
      ablocks[i]=ta.size()+n; t.push_back(tb[i]); n++;
    }
    log.printf("\n"); shape[0]=ta.size(); shape[1]=tb.size();
    // Create a group of the atoms in the rows
    atoms.insertGroup( getLabel(), ta );
  } else {
    // Create list of tasks
    log.printf("  atoms in GROUP "); ablocks.resize( t.size() );
    for(unsigned i=0; i<t.size(); ++i) { log.printf("%d ", t[i].serial()); addTaskToList( i ); ablocks[i]=i; }
    log.printf("\n"); shape[0]=shape[1]=t.size();
    // Create a group of the atoms in the rows
    atoms.insertGroup( getLabel(), t );
  }
  if( keywords.exists("GROUPC") ) {
    std::vector<AtomNumber> tc; parseAtomList("GROUPC",tc);
    if( tc.size()==0 ) error("no atoms in GROUPC specified");
    log.printf("  atoms in GROUPC "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("BRIDGING_ATOMS") ) {
    std::vector<AtomNumber> tc; parseAtomList("BRIDGING_ATOMS",tc);
    if( tc.size()==0 ) error("no BRIDGING_ATOMS specified");
    log.printf("  bridging atoms are "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("HYDROGENS") ) {
    std::vector<AtomNumber> tc; parseAtomList("HYDROGENS",tc);
    if( tc.size()==0 ) error("no HYDROGEN atoms specified");
    log.printf("  hydrogen atoms are "); setupThirdAtomBlock( tc, t );
  } else if( keywords.exists("ATOMS") ) {
    std::vector<AtomNumber> tc; parseAtomList("ATOMS",tc);
    if( tc.size()==0 ) error("no ATOMS atoms specified");
    log.printf("  atoms for background density are "); setupThirdAtomBlock( tc, t );
  }
  // Request the atoms from the ActionAtomistic
  requestAtoms( t ); forcesToApply.resize( getNumberOfDerivatives() );
  parseFlag("COMPONENTS",components); parseFlag("NOPBC",nopbc);
  addComponent( "w", shape ); componentIsNotPeriodic("w");
  if( components ) {
    addComponent( "x", shape ); componentIsNotPeriodic("x");
    addComponent( "y", shape ); componentIsNotPeriodic("y");
    addComponent( "z", shape ); componentIsNotPeriodic("z");
  }
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
}

void AdjacencyMatrixBase::setupThirdAtomBlock( const std::vector<AtomNumber>& tc, std::vector<AtomNumber>& t ) {
  threeblocks.resize( tc.size() ); unsigned base=t.size();
  for(unsigned i=0; i<tc.size(); ++i) { log.printf("%d ", tc[i].serial()); t.push_back(tc[i]); threeblocks[i]=base+i; }
  log.printf("\n");
}

void AdjacencyMatrixBase::setLinkCellCutoff( const double& lcut, double tcut ) {
  if( tcut<0 ) tcut=lcut;
  linkcells.setCutoff( lcut );
  threecells.setCutoff( tcut );
}

void AdjacencyMatrixBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  // Need to retrieve atoms here in case this is done in a stream
  if( actionInChain() ) retrieveAtoms();
  // Make sure all tasks are done
  tflags.assign(tflags.size(),1);
  // Build link cells here so that this is done in stream if it needed in stream
  std::vector<Vector> ltmp_pos( ablocks.size() );
  for(unsigned i=0; i<ablocks.size(); ++i) {
    ltmp_pos[i]=ActionAtomistic::getPosition( ablocks[i] );
  }
  linkcells.buildCellLists( ltmp_pos, ablocks, getPbc() );
  if( threeblocks.size()>0 ) {
    std::vector<Vector> ltmp_pos2( threeblocks.size() );
    for(unsigned i=0; i<ablocks.size(); ++i) {
      ltmp_pos[i]=ActionAtomistic::getPosition( threeblocks[i] );
    }
    threecells.buildCellLists( ltmp_pos2, threeblocks, getPbc() );
  }
}

void AdjacencyMatrixBase::calculate() {
  // Now run all the tasks
  runAllTasks();
}

void AdjacencyMatrixBase::updateWeightDerivativeIndices( const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned w_ind = getPntrToOutput(0)->getPositionInStream();
  // Update dynamic list indices for central atom
  myvals.updateIndex( w_ind, 3*index1+0 ); myvals.updateIndex( w_ind, 3*index1+1 ); myvals.updateIndex( w_ind, 3*index1+2 );
  // Update dynamic list indices for atom forming this bond
  myvals.updateIndex( w_ind, 3*index2+0 ); myvals.updateIndex( w_ind, 3*index2+1 ); myvals.updateIndex( w_ind, 3*index2+2 );
  // Now look after all the atoms in the third block
  std::vector<unsigned> & indices( myvals.getIndices() );
  for(unsigned i=myvals.getSplitIndex(); i<myvals.getNumberOfIndices(); ++i) {
    myvals.updateIndex( w_ind, 3*indices[i]+0 ); myvals.updateIndex( w_ind, 3*indices[i]+1 ); myvals.updateIndex( w_ind, 3*indices[i]+2 );
  }
  // Update dynamic list indices for virial
  unsigned base = 3*getNumberOfAtoms(); for(unsigned j=0; j<9; ++j) myvals.updateIndex( w_ind, base+j );
}

void AdjacencyMatrixBase::setupForTask( const unsigned& current, MultiValue& myvals, std::vector<unsigned> & indices, std::vector<Vector>& atoms ) const {
  // Retrieve cells required from link cells - for matrix blocks
  std::vector<unsigned> cells_required( linkcells.getNumberOfCells() ); unsigned ncells_required=0;
  linkcells.addRequiredCells( linkcells.findMyCell( ActionAtomistic::getPosition(current) ), ncells_required, cells_required );

  // Now retrieve bookeeping arrays
  if( indices.size()!=(1+ablocks.size()+threeblocks.size()) ) indices.resize( 1+ablocks.size()+threeblocks.size() );

  // Now get the positions
  unsigned natoms=1; indices[0]=current;
  linkcells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  unsigned ntwo_atoms=natoms; myvals.setSplitIndex( ntwo_atoms );

  // Now retrieve everything for the third atoms
  if( threeblocks.size()>0 ) {
    if( cells_required.size()!=threecells.getNumberOfCells() ) cells_required.resize( threecells.getNumberOfCells() );
    unsigned ncells_required=0;
    threecells.addRequiredCells( threecells.findMyCell( ActionAtomistic::getPosition(current) ), ncells_required, cells_required );
    threecells.retrieveAtomsInCells( ncells_required, cells_required, natoms, indices );
  }
  myvals.setNumberOfIndices( natoms );
  // Ensure that things that come later know if we have used GROUPA + GROUPB style symmetry function
  if( indices[1]>getFullNumberOfTasks() ) myvals.setNumberOfIndicesInFirstBlock( getFullNumberOfTasks() );
  else myvals.setNumberOfIndicesInFirstBlock( 0 );

// Apply periodic boundary conditions to atom positions
  std::vector<Vector> & t_atoms( myvals.getSecondAtomVector() );
  if( t_atoms.size()<getNumberOfAtoms() ) t_atoms.resize( getNumberOfAtoms() );
  for(unsigned i=0; i<natoms; ++i) t_atoms[i] = ActionAtomistic::getPosition(indices[i]) - ActionAtomistic::getPosition(current);
  if( !nopbc ) pbcApply( t_atoms, natoms );
  // And collect atom position data
  if( atoms.size()<getNumberOfAtoms() ) atoms.resize( getNumberOfAtoms() );
  for(unsigned i=0; i<natoms; ++i) atoms[ indices[i] ] = t_atoms[i];
}

void AdjacencyMatrixBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  // And collect atom position data
  std::vector<unsigned> & indices( myvals.getIndices() );
  std::vector<Vector> & atoms( myvals.getFirstAtomVector() );
  setupForTask( current, myvals, indices, atoms );

  // Now loop over all atoms in coordination sphere
  unsigned natoms = myvals.getNumberOfIndices();
  unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned i=1; i<ntwo_atoms; ++i) {
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), myvals.getTaskIndex(), current, indices[i], myvals );
    // Now clear only elements that are not accumulated over whole row
    clearMatrixElements( myvals );
  }
  // Now update the matrix indices
  if( !doNotCalculateDerivatives() ) updateMatrixIndices( indices, myvals );
}

void AdjacencyMatrixBase::updateMatrixIndices( const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  plumed_assert( !doNotCalculateDerivatives() );
  unsigned natoms = myvals.getNumberOfIndices(), nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( nmat ) );
  plumed_dbg_assert( mat_indices.size()>=(3*getNumberOfAtoms()+9) );
  myvals.setNumberOfMatrixIndices( nmat, 3*natoms+9 );
  for(unsigned i=0; i<natoms; ++i) {
    mat_indices[3*i+0] = 3*indices[i]; mat_indices[3*i+1] = 3*indices[i]+1; mat_indices[3*i+2]=3*indices[i]+2;
  }
  unsigned nbase=3*natoms, vbase=3*getNumberOfAtoms();
  for(unsigned i=0; i<9; ++i) mat_indices[nbase+i] = vbase + i;
}

bool AdjacencyMatrixBase::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  // This makes sure other AdjacencyMatrixBase actions in the stream don't get their matrix elements calculated here
  if( controller!=getLabel() ) return false;

  Vector zero; zero.zero(); plumed_dbg_assert( index2<myvals.getAtomVector().size() );
  double weight = calculateWeight( zero, myvals.getAtomVector()[index2], myvals.getNumberOfIndices()-myvals.getSplitIndex(), myvals );
  if( weight<epsilon ) {
    if( !doNotCalculateDerivatives() ) {
      updateWeightDerivativeIndices( index1, index2, myvals );
      clearMatrixElements( myvals );
    }
    return false;
  }
  // Now set the matrix weight and the vector if required
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), weight );
  if( components ) {
    unsigned x_index = getPntrToOutput(1)->getPositionInStream();
    unsigned y_index = getPntrToOutput(2)->getPositionInStream();
    unsigned z_index = getPntrToOutput(3)->getPositionInStream();
    Vector atom = myvals.getAtomVector()[index2];
    myvals.setValue( x_index, atom[0] ); myvals.setValue( y_index, atom[1] ); myvals.setValue( z_index, atom[2] );
    if( !doNotCalculateDerivatives() ) {
      myvals.addDerivative( x_index, 3*index1+0, -1 ); myvals.addDerivative( x_index, 3*index2+0, +1 );
      myvals.addDerivative( y_index, 3*index1+1, -1 ); myvals.addDerivative( y_index, 3*index2+1, +1 );
      myvals.addDerivative( z_index, 3*index1+2, -1 ); myvals.addDerivative( z_index, 3*index2+2, +1 );
      // Update dynamic lists for central atom
      myvals.updateIndex( x_index, 3*index1+0 ); myvals.updateIndex( y_index, 3*index1+1 ); myvals.updateIndex( z_index, 3*index1+2 );
      // Update dynamic lists for bonded atom
      myvals.updateIndex( x_index, 3*index2+0 ); myvals.updateIndex( y_index, 3*index2+1 ); myvals.updateIndex( z_index, 3*index2+2 );
      // Add derivatives of virial
      unsigned base = 3*getNumberOfAtoms();
      // Virial for x
      myvals.addDerivative( x_index, base+0, -atom[0] ); myvals.addDerivative( x_index, base+3, -atom[1] ); myvals.addDerivative( x_index, base+6, -atom[2] );
      myvals.updateIndex( x_index, base+0 ); myvals.updateIndex( x_index, base+3 ); myvals.updateIndex( x_index, base+6 );
      // Virial for y
      myvals.addDerivative( y_index, base+1, -atom[0] ); myvals.addDerivative( y_index, base+4, -atom[1] ); myvals.addDerivative( y_index, base+7, -atom[2] );
      myvals.updateIndex( y_index, base+1 ); myvals.updateIndex( y_index, base+4 ); myvals.updateIndex( y_index, base+7 );
      // Virial for z
      myvals.addDerivative( z_index, base+2, -atom[0] ); myvals.addDerivative( z_index, base+5, -atom[1] ); myvals.addDerivative( z_index, base+8, -atom[2] );
      myvals.updateIndex( z_index, base+2 ); myvals.updateIndex( z_index, base+5 ); myvals.updateIndex( z_index, base+8 );
    }
  }
  // Update derivatives
  if( !doNotCalculateDerivatives() ) updateWeightDerivativeIndices( index1, index2, myvals );
  return true;
}

void AdjacencyMatrixBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply, mm );
}

}
}
