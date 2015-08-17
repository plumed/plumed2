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

//+PLUMEDOC MATRIXF COLUMNSUMS
/*
Sum the columns of a contact matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class MatrixColumnSums : public MatrixSummationBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixColumnSums(const ActionOptions&);
  double compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const ; 
  Vector getPositionOfAtomForLinkCells( const unsigned& iatom ) const ;
  bool isCurrentlyActive( const unsigned& bno, const unsigned& code );
};

PLUMED_REGISTER_ACTION(MatrixColumnSums,"COLUMNSUMS")

void MatrixColumnSums::registerKeywords( Keywords& keys ){
  MatrixSummationBase::registerKeywords( keys );
}

MatrixColumnSums::MatrixColumnSums(const ActionOptions& ao):
Action(ao),
MatrixSummationBase(ao)
{
 // Setup the tasks
  unsigned ncols = mymatrix->getNumberOfColumns();  
  usespecies=false; ablocks.resize(1); ablocks[0].resize( ncols );
  for(unsigned i=0;i<ncols;++i){ ablocks[0][i]=i; addTaskToList( i ); }
  // Setup the underlying multicolvar
  setupMultiColvarBase();
}

double MatrixColumnSums::compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const {
  double sum=0.0; std::vector<double> tvals(2);
  unsigned nrows = mymatrix->getNumberOfRows();   
  for(unsigned i=0;i<nrows;++i){
     if( mymatrix->isSymmetric() && tinded==i ) continue;
     unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( i, tinded );
     mymatrix->retrieveValue( myelem, false, tvals ); 
     sum+=tvals[1]; 
  }

  if( !doNotCalculateDerivatives() ){
      MultiValue myvals( 2, myatoms.getNumberOfDerivatives() ); 
      MultiValue& myvout=myatoms.getUnderlyingMultiValue();
      for(unsigned i=0;i<nrows;++i){
          if( mymatrix->isSymmetric() && tinded==i ) continue ;
          unsigned myelem = mymatrix->getStoreIndexFromMatrixIndices( i, tinded );
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

bool MatrixColumnSums::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  return (mymatrix->function)->myinputdata.isCurrentlyActive( bno, code );
}

Vector MatrixColumnSums::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return Vector(0.0,0.0,0.0); // (mymatrix->function)->myinputdata.getPosition(iatom);
}

}
}
