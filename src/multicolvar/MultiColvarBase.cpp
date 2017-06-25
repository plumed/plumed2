/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"

namespace PLMD {
namespace multicolvar {

void MultiColvarBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the atoms involved in each of the colvars you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one or more scalars will be "
           "calculated for each ATOM keyword you specify");
  keys.reset_style("ATOMS","atoms");
}

MultiColvarBase::MultiColvarBase(const ActionOptions& ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
usepbc(true)
{
  if( keywords.exists("NOPBC") ) {
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  }
  if( usepbc ) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  std::vector<AtomNumber> all_atoms; parseAtomList( "ATOMS", all_atoms );
  if( all_atoms.size()>0 ){
      ablocks.resize(all_atoms.size());
      log.printf("  Colvar is calculated from atoms : ");
      for(unsigned j=0; j<ablocks.size(); ++j){ ablocks[j].push_back(j); log.printf("%d ",all_atoms[j].serial() ); }
      log.printf("\n");
  } else {
      std::vector<AtomNumber> t;
      for(int i=1;; ++i ) {
        parseAtomList("ATOMS", i, t );
        if( t.empty() ) break;

        log.printf("  Colvar %d is calculated from atoms : ", i);
        for(unsigned j=0; j<t.size(); ++j) log.printf("%d ",t[j].serial() );
        log.printf("\n");

        if( i==1 ) { ablocks.resize(t.size()); }
        if( t.size()!=ablocks.size() ) {
          std::string ss; Tools::convert(i,ss);
          error("ATOMS" + ss + " keyword has the wrong number of atoms");
        }
        for(unsigned j=0; j<ablocks.size(); ++j) {
          ablocks[j].push_back( ablocks.size()*(i-1)+j ); all_atoms.push_back( t[j] );
        }
        t.resize(0);
      }
  }
  requestAtoms(all_atoms);
  if( all_atoms.size()>0 ) {
      for(unsigned i=0; i<ablocks[0].size(); ++i) addTaskToList( i );
  }
}

void MultiColvarBase::addValue(){
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ){ shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addValue( shape );
}

void MultiColvarBase::addComponent( const std::string& name ){
  std::vector<unsigned> shape;
  if( getFullNumberOfTasks()>1 ){ shape.resize(1); shape[0]=getFullNumberOfTasks(); }
  ActionWithValue::addComponent( name, shape );
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc) { return pbcDistance( vec1, vec2 ); }
  else { return delta( vec1, vec2 ); }
}

void MultiColvarBase::buildCurrentTaskList( std::vector<unsigned>& tflags ) const {
  tflags.assign(tflags.size(),1);
}

void MultiColvarBase::calculate(){
  runAllTasks();
}

void MultiColvarBase::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Set the atoms pack up for the calculation 
  AtomValuePack myatoms( myvals, this ); myatoms.setNumberOfAtoms( ablocks.size() );
  for(unsigned i=0; i<ablocks.size(); ++i) myatoms.setAtom( i, ablocks[i][task_index] );
  // If we are using pbc make whole
  if(usepbc) myatoms.makeWhole();
  // And compute
  compute( task_index, myatoms ); 
  // Now update the active derivatives
  if( !doNotCalculateDerivatives() ) myatoms.updateUsingIndices();
}

void MultiColvarBase::apply(){
  if( getFullNumberOfTasks()>1 ) return;
  std::vector<Vector>&   f(modifyForces());
  Tensor&           v(modifyVirial());
  const unsigned    nat=getNumberOfAtoms();
  const unsigned    ncp=getNumberOfComponents();

  std::vector<double> forces(3*nat+9);
  for(unsigned i=0; i<ncp; ++i) {
    if(getPntrToComponent(i)->applyForce(forces)) {
      for(unsigned j=0; j<nat; ++j) {
        f[j][0]+=forces[3*j+0];
        f[j][1]+=forces[3*j+1];
        f[j][2]+=forces[3*j+2];
      }
      v(0,0)+=forces[3*nat+0];
      v(0,1)+=forces[3*nat+1];
      v(0,2)+=forces[3*nat+2];
      v(1,0)+=forces[3*nat+3];
      v(1,1)+=forces[3*nat+4];
      v(1,2)+=forces[3*nat+5];
      v(2,0)+=forces[3*nat+6];
      v(2,1)+=forces[3*nat+7];
      v(2,2)+=forces[3*nat+8];
    }
  }
}

}
}
