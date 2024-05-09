/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR VORONOI
/*
Do a voronoi analysis

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class Voronoi : public ActionWithMatrix {
public:
  static void registerKeywords( Keywords& keys );
  explicit Voronoi(const ActionOptions&);
  void prepare() override ;
  unsigned getNumberOfDerivatives() override { return 0; }
  unsigned getNumberOfColumns() const override { return getConstPntrToComponent(0)->getShape()[1]; }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const override {}
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override {}
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
};

PLUMED_REGISTER_ACTION(Voronoi,"VORONOI")

void Voronoi::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys); keys.use("ARG");
  keys.setValueDescription("a matrix in which element ij is equal to one if the ij component of the input matrix is lower than all the ik elements of the matrix where k is not j and zero otherwise");
}

Voronoi::Voronoi(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao)
{
  if( getNumberOfArguments()!=1 ) error("should be one arguments to this action, a matrix");
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("argument to this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[1]>getPntrToArgument(0)->getShape()[0] ) warning("would expect number of columns in matrix to exceed number of rows");
  getPntrToArgument(0)->buildDataStore(); std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore();
}

void Voronoi::prepare() {
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(0)->getShape()[1] ) return;
  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() ); myval->setShape(shape);
}

void Voronoi::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                 const unsigned& bufstart, std::vector<double>& buffer ) const {
  Value* arg0 = getPntrToArgument(0); unsigned nv = 0; std::size_t cc=code; double minmax = arg0->get( cc*arg0->getShape()[1] );
  for(unsigned i=0; i<arg0->getShape()[1]; ++i) {
    double value = arg0->get( code*arg0->getShape()[1] + i );
    if( value<minmax ) { minmax = value; nv = i; }
  }
  buffer[bufstart + code*arg0->getShape()[1] + nv] = 1;
}

}
}
