/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "FunctionSetup.h"
#include "FunctionOfScalar.h"
#include "FunctionWithSingleArgument.h"

//+PLUMEDOC FUNCTION PRODUCT
/*
Calculate the product of the input quantities

This shortcut can be used to calculate the product of a collection of scalar inputs as illustrated by the following
simple example.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
p: PRODUCT ARG=d1,d2
```

You can also calculate the product of all the elements of a vector:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
p: PRODUCT ARG=d
```

The product is computed directly, so the input values can be negative or zero. The derivative with respect to
each input is the product of all the other inputs.

This action is used in the [DETERMINANT](DETERMINANT.md) shortcut. In that
action the determinant of the input square matrix is found by finding the eigenvalues of the
input matrix using [DIAGONALIZE](DIAGONALIZE.md). The determinant is then found by evaluating
the product of these eigenvalues.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace function {

class ProductFunction {
public:
  static void registerKeywords( Keywords& keys );
  static void read( ProductFunction& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const ProductFunction& func, bool noderiv, View<const double> args, FunctionOutput& funcout );
};

typedef FunctionOfScalar<ProductFunction> ScalarProduct;
PLUMED_REGISTER_ACTION(ScalarProduct,"PRODUCT_SCALAR")
typedef FunctionWithSingleArgument<ProductFunction> OneargProduct;
PLUMED_REGISTER_ACTION(OneargProduct,"PRODUCT_ONEARG")

void ProductFunction::registerKeywords( Keywords& keys ) {
  keys.setValueDescription("scalar","the product of the input values");
}

void ProductFunction::read( ProductFunction& func, ActionWithArguments* action, FunctionOptions& options ) {
  for(unsigned i=0; i<action->getNumberOfArguments(); ++i) {
    if( action->getPntrToArgument(i)->getRank()>1 ) {
      action->error("the product of the elements of a matrix is not implemented");
    }
  }
}

void ProductFunction::calc( const ProductFunction& func, bool noderiv, const View<const double> args, FunctionOutput& funcout ) {
  const std::size_t n=args.size();
  double prod=1.0;
  for(std::size_t i=0; i<n; ++i) {
    prod *= args[i];
  }
  funcout.values[0]=prod;
  if( noderiv ) {
    return;
  }
  // d prod / d x_i is the product of all the other inputs: it is computed
  // with a forward and a backward pass, without dividing by x_i, so zeros are fine
  double left=1.0;
  for(std::size_t i=0; i<n; ++i) {
    funcout.derivs[0][i]=left;
    left *= args[i];
  }
  double right=1.0;
  for(std::size_t i=n; i>0; --i) {
    funcout.derivs[0][i-1] *= right;
    right *= args[i-1];
  }
}

class Product : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Product(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Product,"PRODUCT")

void Product::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar/vector","the scalars (or the vector) whose elements you would like to multiply together");
  keys.setValueDescription("scalar","the product of all the input values");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_ONEARG");
  keys.needsAction("CONCATENATE");
}

Product::Product( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> args;
  parseVector("ARG",args);
  if( args.empty() ) {
    error("no arguments specified");
  }
  std::vector<Value*> vals;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, vals );
  if( vals.empty() ) {
    error("found no input arguments to product");
  }
  std::string allargs=args[0];
  for(unsigned i=1; i<args.size(); ++i) {
    allargs += "," + args[i];
  }
  bool allscalars=true;
  for(const auto & v : vals) {
    if( v->getRank()>1 || (v->getRank()>0 && v->hasDerivatives()) ) {
      error("PRODUCT only accepts scalars and vectors in input");
    }
    if( v->getRank()!=0 ) {
      allscalars=false;
    }
  }
  if( allscalars ) {
    readInputLine( getShortcutLabel() + ": PRODUCT_SCALAR ARG=" + allargs );
  } else if( vals.size()==1 ) {
    readInputLine( getShortcutLabel() + ": PRODUCT_ONEARG ARG=" + allargs );
  } else {
    readInputLine( getShortcutLabel() + "_vec: CONCATENATE ARG=" + allargs );
    readInputLine( getShortcutLabel() + ": PRODUCT_ONEARG ARG=" + getShortcutLabel() + "_vec" );
  }
}

}
}
