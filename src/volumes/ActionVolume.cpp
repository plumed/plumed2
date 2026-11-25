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

namespace PLMD {
namespace volumes {

void ActionVolume::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that you would like to investigate");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
  keys.setValueDescription("vector of numbers between 0 and 1 that measure the degree to which each atom is within the volume of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if( atoms.size()==0 ) {
    error("no atoms were specified");
  }
  log.printf("  examining positions of atoms ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf(" %d", atoms[i].serial() );
  }
  log.printf("\n");
  ActionAtomistic::requestAtoms( atoms );

  parseFlag("OUTSIDE",not_in);
  sigma=0.0;
  if( keywords.exists("SIGMA") ) {
    parse("SIGMA",sigma);
  }
  if( keywords.exists("KERNEL") ) {
    parse("KERNEL",kerneltype);
  }

  if( atoms.size()==1 ) {
    ActionWithValue::addValueWithDerivatives();
  } else {
    std::vector<unsigned> shape(1);
    shape[0]=atoms.size();
    ActionWithValue::addValue( shape );
  }
  setNotPeriodic();
  getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
}

bool ActionVolume::isInSubChain( unsigned& nder ) {
  nder = 0;
  getFirstActionInChain()->getNumberOfStreamedDerivatives( nder, getPntrToComponent(0) );
  nder = nder - getNumberOfDerivatives();
  return true;
}

void ActionVolume::requestAtoms( const std::vector<AtomNumber> & a ) {
  std::vector<AtomNumber> all_atoms( getAbsoluteIndexes() );
  for(unsigned i=0; i<a.size(); ++i) {
    all_atoms.push_back( a[i] );
  }
  ActionAtomistic::requestAtoms( all_atoms );
  if( getPntrToComponent(0)->getRank()==0 ) {
    getPntrToComponent(0)->resizeDerivatives( 3*getNumberOfAtoms()+9 );
  }
}

void ActionVolume::areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {
  task_reducing_actions.push_back(this);
}

void ActionVolume::getNumberOfTasks( unsigned& ntasks ) {
  setupRegions();
  ActionWithVector::getNumberOfTasks( ntasks );
}

int ActionVolume::checkTaskStatus( const unsigned& taskno, int& flag ) const {
  unsigned nref=getNumberOfAtoms()-getConstPntrToComponent(0)->getShape()[0];
  Vector wdf;
  Tensor vir;
  std::vector<Vector> refders( nref );
  double weight=calculateNumberInside( ActionAtomistic::getPosition(taskno), wdf, vir, refders );
  if( not_in ) {
    weight = 1.0 - weight;
  }
  if( weight>epsilon ) {
    return 1;
  }
  return 0;
}

void ActionVolume::calculate() {
  if( actionInChain() ) {
    return;
  }
  if( getPntrToComponent(0)->getRank()==0 ) {
    setupRegions();
    unsigned nref = getNumberOfAtoms() - 1;
    Vector wdf;
    Tensor vir;
    std::vector<Vector> refders( nref );
    double weight=calculateNumberInside( ActionAtomistic::getPosition(0), wdf, vir, refders );
    if( not_in ) {
      weight = 1.0 - weight;
      wdf *= -1.;
      vir *=-1;
      for(unsigned i=0; i<refders.size(); ++i) {
        refders[i]*=-1;
      }
    }
    // Atom position
    Value* v = getPntrToComponent(0);
    v->set( weight );
    for(unsigned i=0; i<3; ++i ) {
      v->addDerivative( i, wdf[i] );
    }
    // Add derivatives with respect to reference positions
    for(unsigned i=0; i<refders.size(); ++i) {
      for(unsigned j=0; j<3; ++j ) {
        v->addDerivative( 3 + 3*i + j, refders[i][j] );
      }
    }
    // Add virial
    unsigned vbase = 3*getNumberOfAtoms();
    for(unsigned i=0; i<3; ++i)
      for(unsigned j=0; j<3; ++j) {
        v->addDerivative( vbase + 3*i + j, vir(i,j) );
      }
  } else {
    runAllTasks();
  }
}

void ActionVolume::performTask( const unsigned& curr, MultiValue& outvals ) const {
  unsigned nref=getNumberOfAtoms()-getConstPntrToComponent(0)->getShape()[0];
  Vector wdf;
  Tensor vir;
  std::vector<Vector> refders( nref );
  double weight=calculateNumberInside( ActionAtomistic::getPosition(curr), wdf, vir, refders );

  if( not_in ) {
    weight = 1.0 - weight;
    wdf *= -1.;
    vir *=-1;
    for(unsigned i=0; i<refders.size(); ++i) {
      refders[i]*=-1;
    }
  }
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  outvals.setValue( ostrn, weight );

  if( doNotCalculateDerivatives() ) {
    return;
  }

  // Atom position
  for(unsigned i=0; i<3; ++i ) {
    outvals.addDerivative( ostrn, 3*curr+i, wdf[i] );
    outvals.updateIndex( ostrn, 3*curr+i );
  }
  // Add derivatives with respect to reference positions
  unsigned vbase = 3*(getNumberOfAtoms()-nref);
  for(unsigned i=0; i<refders.size(); ++i) {
    for(unsigned j=0; j<3; ++j ) {
      outvals.addDerivative( ostrn, vbase, refders[i][j] );
      outvals.updateIndex( ostrn, vbase );
      vbase++;
    }
  }
  // Add virial
  for(unsigned i=0; i<3; ++i) {
    for(unsigned j=0; j<3; ++j) {
      outvals.addDerivative( ostrn, vbase, vir(i,j) );
      outvals.updateIndex( ostrn, vbase );
      vbase++;
    }
  }
}

}
}
