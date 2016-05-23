/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "MatrixSummationBase.h"
#include "vesselbase/ActionWithVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void ActionWithInputMatrix::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys ); keys.remove("DATA");
  keys.add("compulsory","MATRIX","the action that calcualtes the adjacency matrix vessel we would like to analyse"); 
}


ActionWithInputMatrix::ActionWithInputMatrix(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao),
mymatrix(NULL)
{
  if( keywords.exists("MATRIX") ){
      // Find the object that calculates our adjacency matrix
      std::vector<std::string> matname(1); parse("MATRIX",matname[0]);
      bool found_matrix=interpretInputMultiColvars(matname,0.0);
      if(!found_matrix) error("unable to interpret input matrix " + matname[0] );

      // Retrieve the adjacency matrix of interest
      for(unsigned i=0;i<mybasemulticolvars[0]->getNumberOfVessels();++i){
          mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( mybasemulticolvars[0]->getPntrToVessel(i) );
          if( mymatrix ) break ;
      }

      if( !mymatrix ){
         MatrixSummationBase* mybase = dynamic_cast<MatrixSummationBase*>( mybasemulticolvars[0] );
         if( !mybase ) error( matname[0] + " does not calculate an adjacency matrix");
      }
      log.printf("  using matrix calculated by action %s \n",(mymatrix->function)->getLabel().c_str() );

      // And now finish the setup of everything in the base
      setupAtomLists();
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

bool ActionWithInputMatrix::isCurrentlyActive( const unsigned& ind ) const {
  return (mymatrix->function)->isCurrentlyActive( 0, ind );  
}

void ActionWithInputMatrix::getVectorForTask( const unsigned& ind, const bool& normed, std::vector<double>& orient0 ) const {
  if( (mymatrix->function)->colvar_label.size()==0  ){
     double df, sum=0.0; std::vector<double> tvals( mymatrix->getNumberOfComponents() ); 
     unsigned vin; unsigned ncols = mymatrix->getNumberOfColumns(); orient0.assign(orient0.size(),0);
     for(unsigned i=0;i<ncols;++i){
         if( mymatrix->isSymmetric() && ind==i ) continue;
         if( !mymatrix->matrixElementIsActive( ind, i ) ) continue; 
         unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( ind, i );
         mymatrix->retrieveValueWithIndex( myelem, false, tvals ); 
         orient0[1]+=(mymatrix->function)->transformStoredValues( tvals, vin, df);
     } 
     orient0[0]=1.0;
  } else {
     plumed_dbg_assert( isCurrentlyActive( ind ) );
     plumed_dbg_assert( ind<(mymatrix->function)->colvar_label.size() ); unsigned mmc=(mymatrix->function)->colvar_label[ind];
     plumed_dbg_assert( ((mymatrix->function)->mybasedata[mmc])->storedValueIsActive( (mymatrix->function)->convertToLocalIndex(ind,mmc) ) );
     ((mymatrix->function)->mybasedata[mmc])->retrieveValueWithIndex( (mymatrix->function)->convertToLocalIndex(ind,mmc), normed, orient0 );
  }
}

void ActionWithInputMatrix::getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder ) const {
  if( (mymatrix->function)->colvar_label.size()==0  ){
     MultiValue myvals( 2, myder.getNumberOfDerivatives() ); 
     double df, sum=0.0; std::vector<double> tvals( mymatrix->getNumberOfComponents() ); 
     unsigned vin; unsigned ncols = mymatrix->getNumberOfColumns();
     for(unsigned i=0;i<ncols;++i){
         if( mymatrix->isSymmetric() && ind==i ) continue;
         if( !mymatrix->matrixElementIsActive( ind, i ) ) continue; 
         unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( ind, i );
         mymatrix->retrieveValueWithIndex( myelem, false, tvals );
         double dum=(mymatrix->function)->transformStoredValues( tvals, vin, df);
         mymatrix->retrieveDerivatives( myelem, false, myvals ); 
         for(unsigned jd=0;jd<myvals.getNumberActive();++jd){
             unsigned ider=myvals.getActiveIndex(jd);
             myder.addDerivative( 1, ider, df*myvals.getDerivative( vin, ider ) );
         }
     }
  } else { 
     plumed_dbg_assert( isCurrentlyActive( ind ) ); 
     plumed_dbg_assert( ind<(mymatrix->function)->colvar_label.size() ); unsigned mmc=(mymatrix->function)->colvar_label[ind];
     plumed_dbg_assert( ((mymatrix->function)->mybasedata[mmc])->storedValueIsActive( (mymatrix->function)->convertToLocalIndex(ind,mmc) ) );
     if( myder.getNumberOfValues()!=(mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfQuantities() ||
         myder.getNumberOfDerivatives()!=(mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
             myder.resize( (mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfQuantities(), (mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfDerivatives() );
     }
     (mymatrix->function)->mybasedata[mmc]->retrieveDerivatives( (mymatrix->function)->convertToLocalIndex(ind,mmc), normed, myder );
  }
}

unsigned ActionWithInputMatrix::getNumberOfNodeTypes() const {
  unsigned size = (mymatrix->function)->mybasemulticolvars.size();
  if( size==0 ) return 1;
  return size; 
}

unsigned ActionWithInputMatrix::getNumberOfQuantities() const {
  if( (mymatrix->function)->colvar_label.size()==0 ) return 2;
  else return (mymatrix->function)->mybasemulticolvars[0]->getNumberOfQuantities();
}

unsigned ActionWithInputMatrix::getNumberOfAtomsInGroup( const unsigned& igrp ) const {
 plumed_dbg_assert( igrp<(mymatrix->function)->mybasemulticolvars.size() );
 return (mymatrix->function)->mybasemulticolvars[igrp]->getFullNumberOfTasks(); 
}

multicolvar::MultiColvarBase* ActionWithInputMatrix::getBaseMultiColvar( const unsigned& igrp ) const {
 plumed_dbg_assert( igrp<(mymatrix->function)->mybasemulticolvars.size() );
 return (mymatrix->function)->mybasemulticolvars[igrp];
}

}
}
