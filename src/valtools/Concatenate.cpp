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

//+PLUMEDOC MCOLVAR CONCATENATE
/*
Join vectors or matrices together

This action can be used to join sets of scalars, vector or matrices together into a single output value.  The following
example shows how this works with vectors and scalars:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS1=3,4 ATOMS2=5,6
v: CONCATENATE ARG=d1,d2
```

The output vector here, `v`, has three components.  The first of these is the distance between atoms 1 and 2 and the second
and third are the distances between atoms 3 and 4 and 5 and 6 respectively.

Concatenating two or more vectors is similarly straightfoward. Concatenating matrices is a little more difficult as the
user must explain the axes that each matrix should be joined to. The following input should give some idea as to how
the joining of matrices is done in practise.  In this example, we define two groups of atoms and calculate the contact
matrix between all the atoms in the two groups by first calculating three contact matrices that describe the inter and intra
group contacts.

```plumed
# Define two groups of atoms
g: GROUP ATOMS=1-5
h: GROUP ATOMS=6-20

# Calculate the CONTACT_MATRIX for the atoms in group g
cg: CONTACT_MATRIX GROUP=g SWITCH={RATIONAL R_0=0.1}

# Calculate the CONTACT_MATRIX for the atoms in group h
ch: CONTACT_MATRIX GROUP=h SWITCH={RATIONAL R_0=0.2}

# Calculate the CONTACT_MATRIX between the atoms in group g and group h
cgh: CONTACT_MATRIX GROUPA=g GROUPB=h SWITCH={RATIONAL R_0=0.15}

# Now calculate the contact matrix between the atoms in group h and group h
# Notice this is just the transpose of cgh
cghT: TRANSPOSE ARG=cgh

# And concatenate the matrices together to construct the adjacency matrix between the
# adjacency matrices
m: CONCATENATE ...
 MATRIX11=cg MATRIX12=cgh
 MATRIX21=cghT MATRIX22=ch
...
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class Concatenate :
  public ActionWithValue,
  public ActionWithArguments {
private:
  bool vectors;
  std::vector<unsigned> row_starts;
  std::vector<unsigned> col_starts;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit Concatenate(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(Concatenate,"CONCATENATE")

void Concatenate::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("optional","ARG","scalar/vector","the values that should be concatenated together to form the output vector");
  keys.addInputKeyword("numbered","MATRIX","scalar/matrix","specify the matrices that you wish to join together into a single matrix");
  keys.reset_style("MATRIX","compulsory");
  keys.setValueDescription("vector/matrix","the concatenated vector/matrix that was constructed from the input values");
}

Concatenate::Concatenate(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()>0 ) {
    vectors=true;
    std::vector<std::size_t> shape(1);
    shape[0]=0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()>1 ) {
        error("cannot concatenate matrix with vectors");
      }
      shape[0] += getPntrToArgument(i)->getNumberOfValues();
    }
    log.printf("  creating vector with %ld elements \n", shape[0] );
    addValue( shape );
    bool period=getPntrToArgument(0)->isPeriodic();
    std::string min, max;
    if( period ) {
      getPntrToArgument(0)->getDomain( min, max );
    }
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( period!=getPntrToArgument(i)->isPeriodic() ) {
        error("periods of input arguments should match");
      }
      if( period ) {
        std::string min0, max0;
        getPntrToArgument(i)->getDomain( min0, max0 );
        if( min0!=min || max0!=max ) {
          error("domains of input arguments should match");
        }
      }
    }
    if( period ) {
      setPeriodic( min, max );
    } else {
      setNotPeriodic();
    }
    if( getPntrToComponent(0)->getRank()==2 ) {
      getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
    }
  } else {
    unsigned nrows=0, ncols=0;
    std::vector<Value*> arglist;
    vectors=false;
    for(unsigned i=1;; i++) {
      unsigned nt_cols=0;
      unsigned size_b4 = arglist.size();
      for(unsigned j=1;; j++) {
        if( j==10 ) {
          error("cannot combine more than 9 matrices");
        }
        std::vector<Value*> argn;
        parseArgumentList("MATRIX", i*10+j, argn);
        if( argn.size()==0 ) {
          break;
        }
        if( argn.size()>1 ) {
          error("should only be one argument to each matrix keyword");
        }
        if( argn[0]->getRank()!=0 && argn[0]->getRank()!=2 ) {
          error("input arguments for this action should be matrices");
        }
        arglist.push_back( argn[0] );
        nt_cols++;
        if( argn[0]->getRank()==0 ) {
          log.printf("  %d %d component of composed matrix is scalar labelled %s\n", i, j, argn[0]->getName().c_str() );
        } else {
          log.printf("  %d %d component of composed matrix is %ld by %ld matrix labelled %s\n", i, j, argn[0]->getShape()[0], argn[0]->getShape()[1], argn[0]->getName().c_str() );
        }
      }
      if( arglist.size()==size_b4 ) {
        break;
      }
      if( i==1 ) {
        ncols=nt_cols;
      } else if( nt_cols!=ncols ) {
        error("should be joining same number of matrices in each row");
      }
      nrows++;
    }

    std::vector<std::size_t> shape(2);
    shape[0]=0;
    unsigned k=0;
    row_starts.resize( arglist.size() );
    col_starts.resize( arglist.size() );
    for(unsigned i=0; i<nrows; ++i) {
      unsigned cstart = 0, nr = 1;
      if( arglist[k]->getRank()==2 ) {
        nr=arglist[k]->getShape()[0];
      }
      for(unsigned j=0; j<ncols; ++j) {
        if( arglist[k]->getRank()==0 ) {
          if( nr!=1 ) {
            error("mismatched matrix sizes");
          }
        } else if( nrows>1 && arglist[k]->getShape()[0]!=nr ) {
          error("mismatched matrix sizes");
        }
        row_starts[k] = shape[0];
        col_starts[k] = cstart;
        if( arglist[k]->getRank()==0 ) {
          cstart += 1;
        } else {
          cstart += arglist[k]->getShape()[1];
        }
        k++;
      }
      if( i==0 ) {
        shape[1]=cstart;
      } else if( cstart!=shape[1] ) {
        error("mismatched matrix sizes");
      }
      if( arglist[k-1]->getRank()==0 ) {
        shape[0] += 1;
      } else {
        shape[0] += arglist[k-1]->getShape()[0];
      }
    }
    // Now request the arguments to make sure we store things we need
    requestArguments(arglist);
    addValue( shape );
    setNotPeriodic();
    if( getPntrToComponent(0)->getRank()==2 ) {
      getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
    }
  }
}

void Concatenate::calculate() {
  Value* myval = getPntrToComponent(0);
  if( vectors ) {
    unsigned k=0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      Value* myarg=getPntrToArgument(i);
      unsigned nvals=myarg->getNumberOfValues();
      for(unsigned j=0; j<nvals; ++j) {
        myval->set( k, myarg->get(j) );
        k++;
      }
    }
  } else {
    // Retrieve the matrix from input
    unsigned ncols = myval->getShape()[1];
    for(unsigned k=0; k<getNumberOfArguments(); ++k) {
      Value* argn = getPntrToArgument(k);
      if( argn->getRank()==0 ) {
        myval->set( ncols*row_starts[k]+col_starts[k], argn->get() );
      } else {
        std::vector<double> vals;
        std::vector<std::pair<unsigned,unsigned> > pairs;
        bool symmetric=getPntrToArgument(k)->isSymmetric();
        unsigned nedge=0;
        getPntrToArgument(k)->retrieveEdgeList( nedge, pairs, vals );
        for(unsigned l=0; l<nedge; ++l ) {
          unsigned i=pairs[l].first, j=pairs[l].second;
          myval->set( ncols*(row_starts[k]+i)+col_starts[k]+j, vals[l] );
          if( symmetric ) {
            myval->set( ncols*(row_starts[k]+j)+col_starts[k]+i, vals[l] );
          }
        }
      }
    }
  }
}

void Concatenate::apply() {
  if( doNotCalculateDerivatives() || !getPntrToComponent(0)->forcesWereAdded() ) {
    return;
  }

  Value* val=getPntrToComponent(0);
  if( vectors ) {
    unsigned k=0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      Value* myarg=getPntrToArgument(i);
      unsigned nvals=myarg->getNumberOfValues();
      for(unsigned j=0; j<nvals; ++j) {
        myarg->addForce( j, val->getForce(k) );
        k++;
      }
    }
  } else {
    unsigned ncols=val->getShape()[1];
    for(unsigned k=0; k<getNumberOfArguments(); ++k) {
      Value* argn=getPntrToArgument(k);
      if( argn->getRank()==0 ) {
        argn->addForce( 0, val->getForce(ncols*row_starts[k]+col_starts[k]) );
      } else {
        unsigned val_ncols=val->getNumberOfColumns();
        unsigned arg_ncols=argn->getNumberOfColumns();
        for(unsigned i=0; i<argn->getShape()[0]; ++i) {
          unsigned ncol = argn->getRowLength(i);
          for(unsigned j=0; j<ncol; ++j) {
            argn->addForce( i*arg_ncols+j, val->getForce( val_ncols*(row_starts[k]+i)+col_starts[k]+argn->getRowIndex(i,j) ), false );
          }
        }
      }
    }
  }
}

}
}
