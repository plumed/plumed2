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
#include "MatrixSummationBase.h"
#include "multicolvar/AtomValuePack.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC MATRIXF ROWSUMS 
/*
Sum the rows of a contact matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class MatrixRowSums : public MatrixSummationBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixRowSums(const ActionOptions&);
  double compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const ; 
};

PLUMED_REGISTER_ACTION(MatrixRowSums,"ROWSUMS")

void MatrixRowSums::registerKeywords( Keywords& keys ){
  MatrixSummationBase::registerKeywords( keys );
}

MatrixRowSums::MatrixRowSums(const ActionOptions& ao):
Action(ao),
MatrixSummationBase(ao)
{
  // Setup the tasks
  unsigned nrows = mymatrix->getNumberOfRows();  
  usespecies=false; ablocks.resize(1); ablocks[0].resize( nrows );
  for(unsigned i=0;i<nrows;++i){ ablocks[0][i]=i; addTaskToList( i ); }
  // Setup the underlying multicolvar
  ActionAtomistic* matoms = dynamic_cast<ActionAtomistic*>( mymatrix->getMatrixAction() );
  plumed_assert( matoms ); setupMultiColvarBase( matoms->getAbsoluteIndexes(), true );
  addDependency( mymatrix->getMatrixAction() );
}

double MatrixRowSums::compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const {
  double sum=0.0; std::vector<double> tvals(2);
  unsigned ncols = mymatrix->getNumberOfColumns();   
  for(unsigned i=0;i<ncols;++i){
     if( mymatrix->isSymmetric() && tinded==i ) continue;
     unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( tinded, i );
     mymatrix->retrieveValue( myelem, false, tvals ); 
     sum+=tvals[1]; 
  }

  if( !doNotCalculateDerivatives() ){
      MultiValue myvals( 2, myatoms.getNumberOfDerivatives() ); 
      MultiValue& myvout=myatoms.getUnderlyingMultiValue();
      for(unsigned i=0;i<ncols;++i){
          if( mymatrix->isSymmetric() && tinded==i ) continue ;
          unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( tinded, i );
          if( !mymatrix->storedValueIsActive( myelem ) ) continue ;
          mymatrix->retrieveDerivatives( myelem, false, myvals );
          for(unsigned jd=0;jd<myvals.getNumberActive();++jd){
              unsigned ider=myvals.getActiveIndex(jd);
              myvout.addDerivative( 1, ider, myvals.getDerivative( 1, ider ) );
          }
      }
  }
  return sum;
}

}
}
