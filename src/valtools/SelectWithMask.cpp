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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC PRINTANALYSIS SELECT_WITH_MASK
/*
Use a mask to select elements of an array

Output a scalar, vector or matrix that contains a subset of the elements in the input vector or matrix.
The following example shows how we can output a scalar, `v`, that contains the distance between and 3 and 4
by using the mask vector `m` to select this element from the three element vector `d`:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
m: CONSTANT VALUES=1,0,1
v: SELECT_WITH_MASK ARG=d MASK=m
```

The value, `m`, that is passed to the keyword MASK here is a vector with the same length as `d`.
Elements of `d` that whose corresponding elements in `m` are zero are copied to the output value `v`.
When elements of `m` are non-zero the corresponding elements in `d` are not transferred to the output
value - they are masked.

If you use this action with matrices you must use the keywords `ROW_MASK` and `COLUMN_MASK`. As shown in the example
inputs below, these keywords take vectors as input.  In this first example, the output matrix is $3 \times 5$ as rows
of the matrix whose corresponding elements in `m` are non-zero are not transferred:

```plumed
d: DISTANCE_MATRIX GROUP=1-5
m: CONSTANT VALUES=0,1,1,0,0
v: SELECT_WITH_MASK ARG=d ROW_MASK=m
```

For this second example the output matrix is $5 \times 3$ as columns of the matrix whose corresponding elements in `m` are non-zero
are not transferred:

```plumed
d: DISTANCE_MATRIX GROUP=1-5
m: CONSTANT VALUES=0,1,1,0,0
v: SELECT_WITH_MASK ARG=d COLUMN_MASK=m
```

For this final example the output matrix is $3 \times 3$ as we do not transfer the rows and the columns in `d` whose corresponding
elements in `m` are non-zero.

```plumed
d: DISTANCE_MATRIX GROUP=1-5
m: CONSTANT VALUES=0,1,1,0,0
v: SELECT_WITH_MASK ARG=d ROW_MASK=m COLUMN_MASK=m
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class SelectWithMask :
  public ActionWithValue,
  public ActionWithArguments {
private:
  unsigned getOutputVectorLength( const Value* mask ) const ;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit SelectWithMask(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
///
  void getMatrixColumnTitles( std::vector<std::string>& argnames ) const override ;
///
  void prepare() override ;
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(SelectWithMask,"SELECT_WITH_MASK")

void SelectWithMask::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the label for the value upon which you are going to apply the mask");
  keys.addInputKeyword("optional","ROW_MASK","vector","an array with ones in the rows of the matrix that you want to discard");
  keys.addInputKeyword("optional","COLUMN_MASK","vector","an array with ones in the columns of the matrix that you want to discard");
  keys.addInputKeyword("compulsory","MASK","vector/matrix","an array with ones in the components that you want to discard");
  keys.setValueDescription("vector/matrix","a vector/matrix of values that is obtained using a mask to select elements of interest");
  keys.remove("NUMERICAL_DERIVATIVES");
}

SelectWithMask::SelectWithMask(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument for this action");
  }
  std::vector<std::size_t> shape;
  if( getPntrToArgument(0)->getRank()==1 ) {
    std::vector<Value*> mask;
    parseArgumentList("MASK",mask);
    if( mask.size()!=1 ) {
      error("should only be one input for mask");
    }
    if( mask[0]->getNumberOfValues()!=getPntrToArgument(0)->getNumberOfValues() ) {
      error("mismatch between size of mask and input vector");
    }
    log.printf("  creating vector from elements of %s who have a corresponding element in %s that is zero\n", getPntrToArgument(0)->getName().c_str(), mask[0]->getName().c_str() );
    std::vector<Value*> args( getArguments() );
    args.push_back( mask[0] );
    requestArguments( args );
    shape.resize(1,0);
    shape[0]=getOutputVectorLength(mask[0]);
  } else if( getPntrToArgument(0)->getRank()==2 ) {
    std::vector<Value*> rmask, cmask;
    parseArgumentList("ROW_MASK",rmask);
    parseArgumentList("COLUMN_MASK",cmask);
    if( rmask.size()==0 && cmask.size()==0 ) {
      error("no mask elements have been specified");
    } else if( cmask.size()==0 ) {
      std::string con="0";
      for(unsigned i=1; i<getPntrToArgument(0)->getShape()[1]; ++i) {
        con += ",0";
      }
      plumed.readInputWords( Tools::getWords(getLabel() + "_colmask: CONSTANT VALUES=" + con), false );
      std::vector<std::string> labs(1, getLabel() + "_colmask");
      ActionWithArguments::interpretArgumentList( labs, plumed.getActionSet(), this, cmask );
    } else if( rmask.size()==0 ) {
      std::string con="0";
      for(unsigned i=1; i<getPntrToArgument(0)->getShape()[0]; ++i) {
        con += ",0";
      }
      plumed.readInputWords( Tools::getWords(getLabel() + "_rowmask: CONSTANT VALUES=" + con), false );
      std::vector<std::string> labs(1, getLabel() + "_rowmask");
      ActionWithArguments::interpretArgumentList( labs, plumed.getActionSet(), this, rmask );
    }
    shape.resize(2);
    shape[0] = getOutputVectorLength( rmask[0] );
    shape[1] = getOutputVectorLength( cmask[0] );
    std::vector<Value*> args( getArguments() );
    args.push_back( rmask[0] );
    args.push_back( cmask[0] );
    requestArguments( args );
  } else {
    error("input should be vector or matrix");
  }

  addValue( shape );
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  if( getPntrToComponent(0)->getRank()==2 ) {
    getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  }
}

unsigned SelectWithMask::getOutputVectorLength( const Value* mask ) const  {
  unsigned l=0;
  for(unsigned i=0; i<mask->getNumberOfValues(); ++i) {
    if( fabs(mask->get(i))>0 ) {
      continue;
    }
    l++;
  }
  return l;
}

void SelectWithMask::getMatrixColumnTitles( std::vector<std::string>& argnames ) const {
  std::vector<std::string> alltitles;
  (getPntrToArgument(0)->getPntrToAction())->getMatrixColumnTitles( alltitles );
  for(unsigned i=0; i<alltitles.size(); ++i) {
    if( fabs(getPntrToArgument(2)->get(i))>0 ) {
      continue;
    }
    argnames.push_back( alltitles[i] );
  }
}

void SelectWithMask::prepare() {
  Value* arg = getPntrToArgument(0);
  Value* out = getPntrToComponent(0);
  if( arg->getRank()==1 ) {
    Value* mask = getPntrToArgument(1);
    std::vector<std::size_t> shape(1);
    shape[0]=getOutputVectorLength( mask );
    if( out->getNumberOfValues()!=shape[0] ) {
      if( shape[0]==1 ) {
        shape.resize(0);
      }
      out->setShape(shape);
    }
  } else if( arg->getRank()==2 ) {
    std::vector<std::size_t> outshape(2);
    Value* rmask = getPntrToArgument(1);
    outshape[0] = getOutputVectorLength( rmask );
    Value* cmask = getPntrToArgument(2);
    outshape[1] = getOutputVectorLength( cmask );
    if( out->getShape()[0]!=outshape[0] || out->getShape()[1]!=outshape[1] ) {
      out->setShape(outshape);
      out->reshapeMatrixStore( outshape[1] );
    }
  }
}

void SelectWithMask::calculate() {
  Value* arg = getPntrToArgument(0);
  Value* out = getPntrToComponent(0);
  if( arg->getRank()==1 ) {
    Value* mask = getPntrToArgument(1);
    unsigned n=0;
    for(unsigned i=0; i<mask->getNumberOfValues(); ++i) {
      if( fabs(mask->get(i))>0 ) {
        continue;
      }
      out->set(n, arg->get(i) );
      n++;
    }
  } else if ( arg->getRank()==2 ) {
    std::vector<std::size_t> outshape( out->getShape() );
    unsigned n = 0;
    std::vector<std::size_t> inshape( arg->getShape() );
    Value* rmask = getPntrToArgument(1);
    Value* cmask = getPntrToArgument(2);
    for(unsigned i=0; i<inshape[0]; ++i) {
      if( fabs(rmask->get(i))>0 ) {
        continue;
      }
      unsigned m = 0;
      for(unsigned j=0; j<inshape[1]; ++j) {
        if( fabs(cmask->get(j))>0 ) {
          continue;
        }
        out->set( n*outshape[1] + m, arg->get(i*inshape[1] + j) );
        m++;
      }
      n++;
    }
  }
}

void SelectWithMask::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) {
    return ;
  }

  Value* arg = getPntrToArgument(0);
  Value* out = getPntrToComponent(0);
  if( arg->getRank()==1 ) {
    unsigned n=0;
    Value* mask = getPntrToArgument(1);
    for(unsigned i=0; i<mask->getNumberOfValues(); ++i) {
      if( fabs(mask->get(i))>0 ) {
        continue;
      }
      arg->addForce(i, out->getForce(n) );
      n++;
    }
  } else if( arg->getRank()==2 ) {
    unsigned n = 0;
    std::vector<std::size_t> inshape( arg->getShape() );
    std::vector<std::size_t> outshape( out->getShape() );
    Value* rmask = getPntrToArgument(1);
    Value* cmask = getPntrToArgument(2);
    for(unsigned i=0; i<inshape[0]; ++i) {
      if( fabs(rmask->get(i))>0 ) {
        continue;
      }
      unsigned m = 0;
      for(unsigned j=0; j<inshape[1]; ++j) {
        if( fabs(cmask->get(j))>0 ) {
          continue;
        }
        arg->addForce( i*inshape[1] + j, out->getForce(n*outshape[1] + m) );
        m++;
      }
      n++;
    }
  }
}



}
}
