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

class MatrixRowSums : public ActionWithInputMatrix {
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixRowSums(const ActionOptions&);
  double compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const ; 
};

PLUMED_REGISTER_ACTION(MatrixRowSums,"ROWSUMS")

void MatrixRowSums::registerKeywords( Keywords& keys ){
  ActionWithInputMatrix::registerKeywords( keys );
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("MEAN");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

MatrixRowSums::MatrixRowSums(const ActionOptions& ao):
Action(ao),
ActionWithInputMatrix(ao)
{
  if( (mymatrix->getMatrixAction())->mybasemulticolvars.size()>0 ) error("matrix row sums should only be calculated when inputs are atoms");
  // Setup the tasks
  unsigned nrows = mymatrix->getNumberOfRows();  
  ablocks.resize(1); ablocks[0].resize( nrows );
  for(unsigned i=0;i<nrows;++i){ ablocks[0][i]=i; addTaskToList( i ); }
  std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
}

double MatrixRowSums::compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const {
  std::vector<double> tvals( mymatrix->getNumberOfComponents() ); 
  getInputData( tinded, false, myatoms, tvals ); double fval=tvals[1];

  if( !doNotCalculateDerivatives() ){
      tvals.assign( tvals.size(), 0 ); tvals[1]=1.0;
      mergeInputDerivatives( 1, 1, 2, tinded, tvals, getInputDerivatives( tinded, false, myatoms ), myatoms );
  }
  return fval;
}

}
}
