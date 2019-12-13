/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
Sum the rows of a adjacency matrix.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  This action allows you to calculate
the sum of the rows in this adjacency matrix and to then calculate further functions of these quantities.

\par Examples

The first instruction in the following input file tells PLUMED to compute a \f$10 \times 10\f$ matrix in which the \f$ij\f$-element
tells you whether atoms \f$i\f$ and \f$j\f$ are within 1.0 nm of each other.  The numbers in each of this rows are then added together
and the average value is computed.  As such the following input provides an alternative method for calculating the coordination numbers
of atoms 1 to 10.

\plumedfile
mat: CONTACT_MATRIX ATOMS=1-10 SWITCH={RATIONAL R_0=1.0}
rsums: ROWSUMS MATRIX=mat MEAN
PRINT ARG=rsums.* FILE=colvar
\endplumedfile

The following input demonstrates another way that an average coordination number can be computed.  This input calculates the number of atoms
with indices between 6 and 15 that are within the first coordination spheres of each of the atoms within indices between 1 and 5.  The average
coordination number is then calculated from these five coordination numbers and this quantity is output to a file.

\plumedfile
mat2: CONTACT_MATRIX ATOMSA=1-5 ATOMSB=6-15 SWITCH={RATIONAL R_0=1.0}
rsums: ROWSUMS MATRIX=mat2 MEAN
PRINT ARG=rsums.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class MatrixRowSums : public ActionWithInputMatrix {
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixRowSums(const ActionOptions&);
  double compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const override;
};

PLUMED_REGISTER_ACTION(MatrixRowSums,"ROWSUMS")

void MatrixRowSums::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrix::registerKeywords( keys );
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

MatrixRowSums::MatrixRowSums(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrix(ao)
{
  if( (mymatrix->getMatrixAction())->mybasemulticolvars.size()>0 ) warning("matrix row may be problematic when inputs are not atoms");
  // Setup the tasks
  unsigned nrows = mymatrix->getNumberOfRows();
  ablocks.resize(1); ablocks[0].resize( nrows );
  for(unsigned i=0; i<nrows; ++i) { ablocks[0][i]=i; addTaskToList( i ); }
  std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
}

double MatrixRowSums::compute( const unsigned& tinded, multicolvar::AtomValuePack& myatoms ) const {
  std::vector<double> tvals( mymatrix->getNumberOfComponents() );
  getInputData( tinded, false, myatoms, tvals ); double fval=tvals[1];

  if( !doNotCalculateDerivatives() ) {
    tvals.assign( tvals.size(), 0 ); tvals[1]=1.0;
    mergeInputDerivatives( 1, 1, 2, tinded, tvals, getInputDerivatives( tinded, false, myatoms ), myatoms );
  }
  return fval;
}

}
}
