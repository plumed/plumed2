/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionVolume.h"
#include "core/ActionToPutData.h"
#include "core/PlumedMain.h"
#include "multicolvar/MultiColvarBase.h"

namespace PLMD {
namespace volumes {

void ActionVolume::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that you would like to investigate");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao)
{
  std::vector<AtomNumber> atoms; parseAtomList("ATOMS",atoms);
  if( atoms.size()==0 ) error("no atoms were specified");
  log.printf("  examining positions of atoms ");
  for(unsigned i=0; i<atoms.size(); ++i) { log.printf(" %d", atoms[i].serial() );  addTaskToList( i ); }
  log.printf("\n"); ActionAtomistic::requestAtoms( atoms );

  parseFlag("OUTSIDE",not_in); sigma=0.0;
  if( keywords.exists("SIGMA") ) parse("SIGMA",sigma);
  if( keywords.exists("KERNEL") ) parse("KERNEL",kerneltype);

  if( getFullNumberOfTasks()==1 ) { ActionWithValue::addValueWithDerivatives(); }
  else { std::vector<unsigned> shape(1); shape[0]=getFullNumberOfTasks(); ActionWithValue::addValue( shape ); }
  setNotPeriodic();
}

void ActionVolume::requestAtoms( const std::vector<AtomNumber> & a ) {
  std::vector<AtomNumber> all_atoms( getAbsoluteIndexes() );
  for(unsigned i=0; i<a.size(); ++i) all_atoms.push_back( a[i] );
  ActionAtomistic::requestAtoms( all_atoms ); forcesToApply.resize( 3*all_atoms.size()+9 );
}

void ActionVolume::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  setupRegions(); actionsThatSelectTasks.push_back( getLabel() );

  Vector wdf; Tensor vir; std::vector<Vector> refders( getNumberOfAtoms()-getFullNumberOfTasks() );
  for(unsigned i=0;i<tflags.size();++i) {
      // Calculate weight for this position
      double weight=calculateNumberInside( ActionAtomistic::getPosition(i), wdf, vir, refders );
      if( not_in ) weight = 1.0 - weight;
      // Now activate only those tasks that have a significant weight
      if( weight>epsilon ) tflags[i]=1;
  }
}

void ActionVolume::calculate() {
  if( actionInChain() ) return;
  runAllTasks();
}

void ActionVolume::performTask( const unsigned& curr, MultiValue& outvals ) const {
  Vector wdf; Tensor vir; std::vector<Vector> refders( getNumberOfAtoms()-getFullNumberOfTasks() );
  double weight=calculateNumberInside( ActionAtomistic::getPosition(curr), wdf, vir, refders );

  if( not_in ) {
    weight = 1.0 - weight; wdf *= -1.; vir *=-1;
    for(unsigned i=0; i<refders.size(); ++i) refders[i]*=-1;
  }
  outvals.setValue( getPntrToOutput(0)->getPositionInStream(), weight );
  if( !doNotCalculateDerivatives() ) {
    // Atom position
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+0, wdf[0] );
    outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+0 );
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+1, wdf[1] );
    outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+1 );
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), 3*curr+2, wdf[2] );
    outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), 3*curr+2 );
    // Add derivatives with respect to reference positions
    unsigned vbase = 3*getFullNumberOfTasks();
    for(unsigned i=0; i<refders.size(); ++i) {
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][0] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][1] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
      outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, refders[i][2] );
      outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase ); vbase++;
    }
    // Add virial
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,0) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,1) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(0,2) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,0) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,1) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(1,2) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,0) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,1) ); vbase++;
    outvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), vbase, vir(2,2) ); vbase++;
    vbase = 3*getNumberOfAtoms();
    for(unsigned i=0; i<9; ++i) outvals.updateIndex( getPntrToOutput(0)->getPositionInStream(), vbase + i );
  }
}

void ActionVolume::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply, mm );
}

}
}
