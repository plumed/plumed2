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
#include "MatrixSummationBase.h"
#include "multicolvar/AtomValuePack.h"
#include "MatrixSummationBase.h"
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void MatrixSummationBase::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calcualtes the adjacency matrix vessel we would like to analyse");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("MEAN");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN"); 
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

MatrixSummationBase::MatrixSummationBase(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // This ensures that multicolvar base does central atoms correctly
  usespecies=true;
  // Find the object that calculates our adjacency matrix
  std::string matname; parse("MATRIX",matname);
  ActionWithVessel* myvess = plumed.getActionSet().selectWithLabel<ActionWithVessel*>( matname );
  if( !myvess ) error( matname + " does not calculate an adjacency matrix");

  // Retrieve the adjacency matrix of interest
  for(unsigned i=0;i<myvess->getNumberOfVessels();++i){
      mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( myvess->getPntrToVessel(i) );
      if( mymatrix ) break ;
  }
  if( !mymatrix )  error( matname + " does not calculate an adjacency matrix");
  log.printf("  using matrix calculated by action %s \n",matname.c_str() );

}

void MatrixSummationBase::updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const {
  myatoms.updateDynamicList();
} 

bool MatrixSummationBase::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  return (mymatrix->function)->isCurrentlyActive( bno, code );
}

Vector MatrixSummationBase::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return (mymatrix->function)->getPositionOfAtomForLinkCells(iatom);
}

AtomNumber MatrixSummationBase::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
  return (mymatrix->function)->getAbsoluteIndexOfCentralAtom(i);
}

double MatrixSummationBase::retrieveConnectionValue( const unsigned& i, const unsigned& j, std::vector<double>& vals ) const {
  if( !mymatrix->matrixElementIsActive( i, j ) ) return 0;
  unsigned vi, myelem = mymatrix->getStoreIndexFromMatrixIndices( i, j );
 
  mymatrix->retrieveValueWithIndex( myelem, false, vals ); double df;
  return (mymatrix->function)->transformStoredValues( vals, vi, df );
}

void MatrixSummationBase::addConnectionDerivatives( const unsigned& i, const unsigned& j, std::vector<double>& vals, MultiValue& myvals, MultiValue& myvout ) const {
  if( !mymatrix->matrixElementIsActive( i, j ) ) return;
  unsigned vi, myelem = mymatrix->getStoreIndexFromMatrixIndices( i, j );

  mymatrix->retrieveValueWithIndex( myelem, false, vals ); double df;
  double vv = (mymatrix->function)->transformStoredValues( vals, vi, df );
  mymatrix->retrieveDerivatives( myelem, false, myvals );
  for(unsigned jd=0;jd<myvals.getNumberActive();++jd){
      unsigned ider=myvals.getActiveIndex(jd);
      myvout.addDerivative( 1, ider, df*myvals.getDerivative( vi, ider ) );
  }
}

}
}
