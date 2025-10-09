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

This action can be used to calculate the [outer product](https://en.wikipedia.org/wiki/Outer_product) of two
vectors.  As a (useless) example of what can be done with this action consider the following simple input:

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8 ATOMS3=9,10
pp: OUTER_PRODUCT ARG=d1,d2
PRINT ARG=pp FILE=colvar
```

This input outputs a $2 \times 3$ matrix. If we call the 2 dimensional vector output by the first DISTANCE action
$d$ and the 3 dimensional vector output by the second DISTANCE action $h$ then the $(i,j)$ element of the matrix
output by the action with the label `pp` is given by:

$$
p_{ij} = d_i h_j
$$

These outer product matrices are useful if you are trying to calculate an adjacency matrix that says atoms are
connected if they are within a certain distance of each other and if they satisfy a certain criterion.  For example,
consider the following input:

```plumed
# Determine if atoms are within 5.3 nm of each other
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
# Calculate the coordination numbers
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
# Now use MORE_THAN to work out which atoms have a coordination number that is bigger than six
cf: MORE_THAN ARG=cc SWITCH={RATIONAL D_0=5.5 R_0=0.5}
# Now recalculate the contact matrix above as first step towards calculating adjacency matrix that measures if
# atoms are close to each other and both have a coordination number that is bigger than six
c2: CONTACT_MATRIX GROUP=1-100 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3}
# Now make a matrix in which element i,j is one if atom i and atom j both have a coordination number that is greater than 6
cfm: OUTER_PRODUCT ARG=cf,cf MASK=c2
# And multiply this by our contact matrix to determine the desired adjacency matrix
m: CUSTOM ARG=c2,cfm FUNC=x*y PERIODIC=NO
f: SUM ARG=m PERIODIC=NO
PRINT ARG=f FILE=colvar
```

This input calculates a adjacency matrix which has element $(i,j)$ equal to one if atoms $i$ and $j$ have coordination numbers
that are greater than 6 and if they are within 5.3 nm of each other.  Notice how the `MASK` keyword is used in the input to the
OUTER_PRODUCT action here to ensure that do not calculate elements of the `cfm` matrix that will be mulitplied by elements of the
matrix `c2` that are zero.  The final quantity output is equal to two times the number of pairs of atoms that are within 5.3 nm of each
and which both have coordination numbers of six.

Notice that you can specify the function of the two input vectors that is to be calculated by using the `FUNC` keyword which accepts
mathematical expressions of $x$ and $y$.  In other words, the elements of the outer product are calculated using the lepton library
that is used in the [CUSTOM](CUSTOM.md) action.  In addition, you can set `FUNC=min` or `FUNC=max` to set the elements of the outer product equal to
the minimum of the two input variables or the maximum respectively.

## Calculating angles in the first coordination sphere

We can use OUTER_PRODUCT to calculate a matrix of angles between bonds as shown below:

```plumed
# Calculate the directors for a set of vectors
d: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,6
dm: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,6
dx: CUSTOM ARG=d.x,dm FUNC=x/y PERIODIC=NO
dy: CUSTOM ARG=d.y,dm FUNC=x/y PERIODIC=NO
dz: CUSTOM ARG=d.z,dm FUNC=x/y PERIODIC=NO
# Construct a matrix that contains all the directors of the vectors calculated
v: VSTACK ARG=dx,dy,dz
# Transpose v
vT: TRANSPOSE ARG=v
# Transform the distances by a switching functions to determine pairs of atoms that are bonded
sw: LESS_THAN ARG=dm SWITCH={RATIONAL R_0=0.2}
# Calculate the matrix of dot products between the input directors
dpmat: MATRIX_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=v,vT
# Use the transformed distances to determine which triples of atoms are bonded
swmat: OUTER_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=sw,sw
# And calculate the angles
angles: CUSTOM ARG=swmat,dpmat FUNC=x*acos(y) PERIODIC=NO
# Print the matrix of angles
PRINT ARG=angles FILE=colvar
```

Notice that we have to use the `ELEMENTS_ON_DIAGONAL_ARE_ZERO` flag here to avoid numerical issues in the calculation.

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
  void setup( const std::vector<std::size_t>& shape,
              const std::string& func,
              OuterProductBase<OutputProductMin>* action );
  static void calculate( bool noderiv,
                         const OutputProductMin& actdata,
                         View<double> vals,
                         MatrixElementOutput& output );
};

typedef OuterProductBase<OutputProductMin> opmin;
PLUMED_REGISTER_ACTION(opmin,"OUTER_PRODUCT_MIN")

void OutputProductMin::registerKeywords( Keywords& keys ) {
  OuterProduct::getKeywords( keys );
}

void OutputProductMin::setup( const std::vector<std::size_t>& shape,
                              const std::string& func,
                              OuterProductBase<OutputProductMin>* action ) {
  plumed_assert( func=="min" );
  action->log.printf("  taking minimum of two input vectors \n");
}

void OutputProductMin::calculate( bool noderiv,
                                  const OutputProductMin& actdata,
                                  View<double> vals,
                                  MatrixElementOutput& output ) {
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
  void setup( const std::vector<std::size_t>& shape,
              const std::string& func,
              OuterProductBase<OutputProductMax>* action );
  static void calculate( bool noderiv,
                         const OutputProductMax& actdata,
                         View<const double> vals,
                         MatrixElementOutput& output );
};

typedef OuterProductBase<OutputProductMax> opmax;
PLUMED_REGISTER_ACTION(opmax,"OUTER_PRODUCT_MAX")

void OutputProductMax::registerKeywords( Keywords& keys ) {
  OuterProduct::getKeywords( keys );
}

void OutputProductMax::setup( const std::vector<std::size_t>& shape,
                              const std::string& func,
                              OuterProductBase<OutputProductMax>* action ) {
  plumed_assert( func=="max" );
  action->log.printf("  taking maximum of two input vectors \n");
}

void OutputProductMax::calculate( bool noderiv,
                                  const OutputProductMax& actdata,
                                  View<const double> vals,
                                  MatrixElementOutput& output ) {
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
  void setup( const std::vector<std::size_t>& shape,
              const std::string& func,
              OuterProductBase<OutputProductFunc>* action );
  static void calculate( bool noderiv,
                         const OutputProductFunc& actdata,
                         View<const double> vals,
                         MatrixElementOutput& output );
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

void OutputProductFunc::setup( const std::vector<std::size_t>& shape,
                               const std::string& func,
                               OuterProductBase<OutputProductFunc>* action ) {
  action->log.printf("  with function : %s \n", func.c_str() );
  inputf = func;
  std::vector<std::string> var(2);
  var[0]="x";
  var[1]="y";
  function.set( func, var, action );
}

void OutputProductFunc::calculate( bool noderiv,
                                   const OutputProductFunc& actdata,
                                   View<const double> vals,
                                   MatrixElementOutput& output ) {
  output.values[0] = actdata.function.evaluate( vals );
  output.derivs[0][0] = actdata.function.evaluateDeriv( 0, vals );
  output.derivs[0][1] = actdata.function.evaluateDeriv( 1, vals );
}

}
}
