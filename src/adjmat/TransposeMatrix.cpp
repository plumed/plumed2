/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class TransposeMatrix :
  public ActionWithArguments,
  public ActionWithValue
{
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TransposeMatrix(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return 1; }
/// Do the calculation
  void calculate();
///
  void apply();
};

PLUMED_REGISTER_ACTION(TransposeMatrix,"TRANSPOSE")

void TransposeMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.use("ARG");
}

TransposeMatrix::TransposeMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument for this action should be a matrix");
  std::vector<unsigned> shape(2); shape[0]=getPntrToArgument(0)->getShape()[1]; shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape ); getPntrToOutput(0)->alwaysStoreValues();
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
}

void TransposeMatrix::calculate() {
  // Retrieve the matrix from input
  unsigned k = 0; std::vector<unsigned> shape( getPntrToOutput(0)->getShape() );
  for(unsigned i=0; i<shape[0]; ++i) {
    for(unsigned j=0; j<shape[1]; ++j) {
      getPntrToOutput(0)->set( j*shape[1] + i, getPntrToArgument(0)->get( k ) ); k++;
    }
  }
}

void TransposeMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;

  if( getPntrToOutput(0)->forcesWereAdded() ) {
    unsigned k = 0; std::vector<unsigned> shape( getPntrToOutput(0)->getShape() );
    for(unsigned i=0; i<shape[0]; ++i) {
      for(unsigned j=0; j<shape[1]; ++j) {
        getPntrToArgument(0)->addForce( j*shape[1] + i, getPntrToOutput(0)->getForce(k) ); k++;
      }
    }
  }
}

}
}
