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
#include "OuterProduct.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "tools/LeptonCall.h"

//+PLUMEDOC COLVAR OUTER_PRODUCT
/*
Calculate the outer product matrix of two vectors

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class OuterProduct : public ActionShortcut {
public:
  static void getKeywords( Keywords& keys );
  static void registerKeywords( Keywords& keys );
  explicit OuterProduct(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(OuterProduct,"OUTER_PRODUCT")

void OuterProduct::getKeywords( Keywords& keys ) {
  keys.setDisplayName("OUTER_PRODUCT");
  keys.addInputKeyword("compulsory","ARG","vector","the labels of the two vectors from which the outer product is being computed");
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.add("compulsory","FUNC","x*y","the function of the input vectors that should be put in the elements of the outer product");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
  keys.setValueDescription("matrix","a matrix containing the outer product of the two input vectors that was obtained using the function that was input");
}

void OuterProduct::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  getKeywords( keys );
  keys.addActionNameSuffix("_MIN");
  keys.addActionNameSuffix("_MAX");
  keys.addActionNameSuffix("_FUNC");
}

OuterProduct::OuterProduct(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {

  std::string func;
  parse("FUNC",func);
  if( func=="min") {
    readInputLine( getShortcutLabel() + ": OUTER_PRODUCT_MIN FUNC=" + func + " " + convertInputLineToString() );
  } else if( func=="max" ) {
    readInputLine( getShortcutLabel() + ": OUTER_PRODUCT_MAX FUNC=" + func + " " + convertInputLineToString() );
  } else {
    readInputLine( getShortcutLabel() + ": OUTER_PRODUCT_FUNC FUNC=" + func + " " + convertInputLineToString() );
  }
}

class OutputProductMin {
public:
  static void registerKeywords( Keywords& keys );
  void setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductMin>* action );
  static void calculate( bool noderiv, const OutputProductMin& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output );
};

typedef OuterProductBase<OutputProductMin> opmin;
PLUMED_REGISTER_ACTION(opmin,"OUTER_PRODUCT_MIN")

void OutputProductMin::registerKeywords( Keywords& keys ) {
  OuterProduct::getKeywords( keys );
}

void OutputProductMin::setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductMin>* action ) {
  plumed_assert( func=="min" );
  action->log.printf("  taking minimum of two input vectors \n");
}

void OutputProductMin::calculate( bool noderiv, const OutputProductMin& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output ) {
  if( vals[0]<vals[1] ) {
    output.derivs[0][0] = 1;
    output.derivs[0][1] = 0;
    output.values[0] = vals[0];
    return;
  }
  output.derivs[0][0] = 0;
  output.derivs[0][1] = 1;
  output.values[0] = vals[1];
}

class OutputProductMax {
public:
  static void registerKeywords( Keywords& keys );
  void setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductMax>* action );
  static void calculate( bool noderiv, const OutputProductMax& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output );
};

typedef OuterProductBase<OutputProductMax> opmax;
PLUMED_REGISTER_ACTION(opmax,"OUTER_PRODUCT_MAX")

void OutputProductMax::registerKeywords( Keywords& keys ) {
  OuterProduct::getKeywords( keys );
}

void OutputProductMax::setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductMax>* action ) {
  plumed_assert( func=="max" );
  action->log.printf("  taking maximum of two input vectors \n");
}

void OutputProductMax::calculate( bool noderiv, const OutputProductMax& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output ) {
  if( vals[0]>vals[1] ) {
    output.derivs[0][0] = 1;
    output.derivs[0][1] = 0;
    output.values[0] = vals[0];
    return;
  }
  output.derivs[0][0] = 0;
  output.derivs[0][1] = 1;
  output.values[0] = vals[1];
}

class OutputProductFunc {
public:
  std::string inputf;
  LeptonCall function;
  static void registerKeywords( Keywords& keys );
  void setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductFunc>* action );
  static void calculate( bool noderiv, const OutputProductFunc& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output );
  OutputProductFunc& operator=( const OutputProductFunc& m ) {
    inputf = m.inputf;
    std::vector<std::string> var(2);
    var[0]="x";
    var[1]="y";
    function.set( inputf, var );
    return *this;
  }
};

typedef OuterProductBase<OutputProductFunc> opfunc;
PLUMED_REGISTER_ACTION(opfunc,"OUTER_PRODUCT_FUNC")

void OutputProductFunc::registerKeywords( Keywords& keys ) {
  OuterProduct::getKeywords( keys );
}

void OutputProductFunc::setup( const std::vector<std::size_t>& shape, const std::string& func, OuterProductBase<OutputProductFunc>* action ) {
  action->log.printf("  with function : %s \n", func.c_str() );
  inputf = func;
  std::vector<std::string> var(2);
  var[0]="x";
  var[1]="y";
  function.set( func, var, action );
}

void OutputProductFunc::calculate( bool noderiv, const OutputProductFunc& actdata, View<double,helpers::dynamic_extent> vals, MatrixElementOutput& output ) {
  std::vector<double> vvv(2);
  vvv[0]=vals[0];
  vvv[1] = vals[1];
  output.values[0] = actdata.function.evaluate( vvv );
  output.derivs[0][0] = actdata.function.evaluateDeriv( 0, vvv );
  output.derivs[0][1] = actdata.function.evaluateDeriv( 1, vvv );
}

}
}
