/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "ActionWithInputMatrix.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "vesselbase/ActionWithVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void ActionWithInputMatrix::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calculates the adjacency matrix vessel we would like to analyze");
}


ActionWithInputMatrix::ActionWithInputMatrix(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao),
  mymatrix(NULL)
{
  matsums=true;
  if( keywords.exists("MATRIX") ) {
    std::vector<AtomNumber> fake_atoms;
    if( !parseMultiColvarAtomList("MATRIX",-1,fake_atoms ) ) error("unable to interpret input matrix");
    if( mybasemulticolvars.size()!=1 ) error("should be exactly one matrix input");

    // Retrieve the adjacency matrix of interest
    for(unsigned i=0; i<mybasemulticolvars[0]->getNumberOfVessels(); ++i) {
      mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( mybasemulticolvars[0]->getPntrToVessel(i) );
      if( mymatrix ) break ;
    }
    if( !mymatrix ) error( mybasemulticolvars[0]->getLabel() + " does not calculate an adjacency matrix");

    atom_lab.resize(0); unsigned nnodes; // Delete all the atom labels that have been created
    if( mymatrix->undirectedGraph() ) nnodes = (mymatrix->function)->ablocks[0].size();
    else nnodes = (mymatrix->function)->ablocks[0].size() + (mymatrix->function)->ablocks[1].size();
    for(unsigned i=0; i<nnodes; ++i) atom_lab.push_back( std::pair<unsigned,unsigned>( 1, i ) );
  }
}

unsigned ActionWithInputMatrix::getNumberOfDerivatives() {
  return (mymatrix->function)->getNumberOfDerivatives();
}

unsigned ActionWithInputMatrix::getNumberOfNodes() const {
  return (mymatrix->function)->ablocks[0].size();
}

AdjacencyMatrixVessel* ActionWithInputMatrix::getAdjacencyVessel() const {
  return mymatrix;
}

AtomNumber ActionWithInputMatrix::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
  return (mymatrix->function)->getAbsoluteIndexOfCentralAtom(i);
}

double ActionWithInputMatrix::retrieveConnectionValue( const unsigned& i, const unsigned& j, std::vector<double>& vals ) const {
  if( !mymatrix->matrixElementIsActive( i, j ) ) return 0;
  unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( i, j );

  // unsigned vi; double df;
  mymatrix->retrieveValueWithIndex( myelem, false, vals );
  return vals[0]*vals[1];       // (mymatrix->function)->transformStoredValues( vals, vi, df );
}

void ActionWithInputMatrix::getInputData( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms, std::vector<double>& orient0 ) const {
  if( (mymatrix->function)->mybasemulticolvars.size()==0  ) {
    std::vector<double> tvals( mymatrix->getNumberOfComponents() ); orient0.assign(orient0.size(),0);
    for(unsigned i=0; i<mymatrix->getNumberOfColumns(); ++i) {
      if( mymatrix->undirectedGraph() && ind==i ) continue;
      orient0[1]+=retrieveConnectionValue( ind, i, tvals );
    }
    orient0[0]=1.0; return;
  }
  (mymatrix->function)->getInputData( ind, normed, myatoms, orient0 );
}

void ActionWithInputMatrix::addConnectionDerivatives( const unsigned& i, const unsigned& j, MultiValue& myvals, MultiValue& myvout ) const {
  if( !mymatrix->matrixElementIsActive( i, j ) ) return;
  unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( i, j );
  // Get derivatives and add
  mymatrix->retrieveDerivatives( myelem, false, myvals );
  for(unsigned jd=0; jd<myvals.getNumberActive(); ++jd) {
    unsigned ider=myvals.getActiveIndex(jd);
    myvout.addDerivative( 1, ider, myvals.getDerivative( 1, ider ) );
  }
}

MultiValue& ActionWithInputMatrix::getInputDerivatives( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const {
  if( (mymatrix->function)->mybasemulticolvars.size()==0  ) {
    MultiValue& myder=mymatrix->getTemporyMultiValue(0);
    if( myder.getNumberOfValues()!=2 || myder.getNumberOfDerivatives()!=(mymatrix->function)->getNumberOfDerivatives() ) {
      myder.resize( 2, (mymatrix->function)->getNumberOfDerivatives() );
    }
    myder.clearAll();
    MultiValue myvals( (mymatrix->function)->getNumberOfQuantities(), (mymatrix->function)->getNumberOfDerivatives() );
    for(unsigned i=0; i<mymatrix->getNumberOfColumns(); ++i) {
      if( mymatrix->undirectedGraph() && ind==i ) continue;
      addConnectionDerivatives( ind, i, myvals, myder );
    }
    myder.updateDynamicList(); return myder;
  }
  return (mymatrix->function)->getInputDerivatives( ind, normed, myatoms );
}

unsigned ActionWithInputMatrix::getNumberOfNodeTypes() const {
  unsigned size = (mymatrix->function)->mybasemulticolvars.size();
  if( size==0 ) return 1;
  return size;
}

unsigned ActionWithInputMatrix::getNumberOfQuantities() const {
  if( (mymatrix->function)->mybasemulticolvars.size()==0 ) return 2;
  return (mymatrix->function)->mybasemulticolvars[0]->getNumberOfQuantities();
}

unsigned ActionWithInputMatrix::getNumberOfAtomsInGroup( const unsigned& igrp ) const {
  plumed_dbg_assert( igrp<(mymatrix->function)->mybasemulticolvars.size() );
  return (mymatrix->function)->mybasemulticolvars[igrp]->getFullNumberOfTasks();
}

multicolvar::MultiColvarBase* ActionWithInputMatrix::getBaseMultiColvar( const unsigned& igrp ) const {
  plumed_dbg_assert( igrp<(mymatrix->function)->mybasemulticolvars.size() );
  return (mymatrix->function)->mybasemulticolvars[igrp];
}

Vector ActionWithInputMatrix::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return (getAdjacencyVessel()->function)->getPositionOfAtomForLinkCells( iatom );
}

}
}
