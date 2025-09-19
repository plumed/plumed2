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
#ifndef __PLUMED_matrixtools_MatrixTimesVectorBase_h
#define __PLUMED_matrixtools_MatrixTimesVectorBase_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace matrixtools {

namespace helpers {
///true if T has a public declared "isRowSum" type
template <typename T, typename=void>
constexpr bool isRowSum=false;
template <typename T>
constexpr bool isRowSum<T,std::void_t<typename T::isRowSum>> =true;
}

class MatrixTimesVectorData {
public:
  std::size_t fshift;
  Matrix<std::size_t> pairs;
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],fshift)
    pairs.toACCDevice();
  }
  void removeFromACCDevice() const {
    pairs.removeFromACCDevice();
#pragma acc exit data delete(fshift,this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

class MatrixForceIndexInput {
public:
  std::size_t rowlen;
  View<const std::size_t> indices;
//temporary template for compilation
  template <typename ParallelActionsInput>
  MatrixForceIndexInput( std::size_t task_index,
                         std::size_t ipair,
                         const MatrixTimesVectorData& actiondata,
                         const ParallelActionsInput& input ):
    rowlen(input.bookeeping[input.bookstarts[actiondata.pairs[ipair][0]]
                            + (1+input.ncols[actiondata.pairs[ipair][0]])*task_index]),
    indices(input.bookeeping + input.bookstarts[actiondata.pairs[ipair][0]]
            + (1+input.ncols[actiondata.pairs[ipair][0]])*task_index + 1,
            rowlen) {}
};

template <typename precision>
struct MatrixTimesVectorInput {
  bool noderiv;
  std::size_t rowlen;
  View<const std::size_t> indices;
  View<const precision> matrow;
  View<const precision> vector;
//temporary template for compilation
  template <typename ParallelActionsInput>
  MatrixTimesVectorInput( std::size_t task_index,
                          std::size_t ipair,
                          const MatrixTimesVectorData& actiondata,
                          const ParallelActionsInput& input,
                          precision* argdata ):
    noderiv(input.noderiv),
    rowlen(input.bookeeping[input.bookstarts[actiondata.pairs[ipair][0]]
                            + (1+input.ncols[actiondata.pairs[ipair][0]])*task_index]),
    indices(input.bookeeping + input.bookstarts[actiondata.pairs[ipair][0]]
            + (1+input.ncols[actiondata.pairs[ipair][0]])*task_index + 1,rowlen),
    matrow(argdata + input.argstarts[actiondata.pairs[ipair][0]]
           + task_index*input.ncols[actiondata.pairs[ipair][0]],rowlen),
    vector(argdata + input.argstarts[actiondata.pairs[ipair][1]], input.shapedata[1]) {
  }
};

template <typename precision>
class MatrixTimesVectorOutput {
public:
  std::size_t rowlen;
  View<precision,1> values;
  View<precision> matrow_deriv;
  View<precision> vector_deriv;
//temporary template for compilation
  template <typename ParallelActionsInput, typename ParallelActionsOutput>
  MatrixTimesVectorOutput( std::size_t task_index,
                           std::size_t ipair,
                           std::size_t nder,
                           const MatrixTimesVectorData& actiondata,
                           const ParallelActionsInput& input,
                           ParallelActionsOutput& output ):
    rowlen(input.bookeeping[input.bookstarts[actiondata.pairs[ipair][0]]
                            + (1+input.ncols[actiondata.pairs[ipair][0]])*task_index]),
    values(output.values.data()+ipair),
    matrow_deriv(output.derivatives.data()+ipair*nder,rowlen),
    vector_deriv(output.derivatives.data()+ipair*nder+rowlen,rowlen) {
  }
};

template <class CV, typename myPTM=defaultPTM>
class MatrixTimesVectorBase : public ActionWithVector {
public:
  using input_type = MatrixTimesVectorData;
  typedef cvprecision_t<CV> precision;
  typedef MatrixTimesVectorBase<CV,myPTM> mytype;
  using PTM = typename myPTM::template PTM<mytype>;
  typedef typename PTM::ParallelActionsInput ParallelActionsInput;
  typedef typename PTM::ParallelActionsOutput ParallelActionsOutput;
private:
/// The parallel task manager
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  static void registerLocalKeywords( Keywords& keys );
  explicit MatrixTimesVectorBase(const ActionOptions&);
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  unsigned getNumberOfDerivatives() override ;
  void prepare() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  int checkTaskIsActive( const unsigned& itask ) const override ;
  /// Override this so we write the graph properly
  std::string writeInGraph() const override {
    return "MATRIX_VECTOR_PRODUCT";
  }
  static void performTask( std::size_t task_index,
                           const MatrixTimesVectorData& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const MatrixTimesVectorData& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const MatrixTimesVectorData& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.setDisplayName("MATRIX_VECTOR_PRODUCT");
  registerLocalKeywords( keys );
  ActionWithValue::useCustomisableComponents(keys);
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::registerLocalKeywords( Keywords& keys ) {
  PTM::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","matrix/vector/scalar","the label for the matrix and the vector/scalar that are being multiplied.  Alternatively, you can provide labels for multiple matrices and a single vector or labels for a single matrix and multiple vectors. In these cases multiple matrix vector products will be computed.");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.setValueDescription("vector","the vector that is obtained by taking the product between the matrix and the vector that were input");
  ActionWithValue::useCustomisableComponents(keys);
}

template <class CV, typename myPTM>
std::string MatrixTimesVectorBase<CV, myPTM>::getOutputComponentDescription( const std::string& cname,
    const Keywords& keys ) const {
  if( getPntrToArgument(1)->getRank()==1 ) {
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getName().find(cname)!=std::string::npos ) {
        return "the product of the matrix " + getPntrToArgument(0)->getName() + " and the vector " + getPntrToArgument(i)->getName();
      }
    }
  }
  for(unsigned i=0; i<getNumberOfArguments()-1; ++i) {
    if( getPntrToArgument(i)->getName().find(cname)!=std::string::npos ) {
      return "the product of the matrix " + getPntrToArgument(i)->getName() + " and the vector " + getPntrToArgument(getNumberOfArguments()-1)->getName();
    }
  }
  plumed_merror( "could not understand request for component " + cname );
  return "";
}

template <class CV, typename myPTM>
MatrixTimesVectorBase<CV, myPTM>::MatrixTimesVectorBase(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  if( getNumberOfArguments()<2 ) {
    error("Not enough arguments specified");
  }
  bool vectormask=false, derivbool = true;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->hasDerivatives() ) {
      error("arguments should be vectors or matrices");
    }
    if( getPntrToArgument(i)->getRank()<=1 ) {
      ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
      if( av && av->getNumberOfMasks()>=0 ) {
        vectormask=true;
      }
    }
    if( !getPntrToArgument(i)->isDerivativeZeroWhenValueIsZero() ) {
      derivbool = false;
    }
  }
  if( !vectormask ) {
    ignoreMaskArguments();
  }

  std::vector<std::size_t> shape(1);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  if( getNumberOfArguments()==2 ) {
    addValue( shape );
    setNotPeriodic();
    if( derivbool ) {
      getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
    }
  } else {
    unsigned namestart=1, nameend=getNumberOfArguments();
    if( getPntrToArgument(1)->getRank()==2 ) {
      namestart = 0;
      nameend = getNumberOfArguments()-1;
    }

    for(unsigned i=namestart; i<nameend; ++i) {
      std::string name = getPntrToArgument(i)->getName();
      if( name.find_first_of(".")!=std::string::npos ) {
        std::size_t dot=name.find_first_of(".");
        name = name.substr(dot+1);
      }
      addComponent( name, shape );
      componentIsNotPeriodic( name );
      if( derivbool ) {
        copyOutput( getLabel() + "." + name )->setDerivativeIsZeroWhenValueIsZero();
      }
    }
  }
  // This sets up an array in the parallel task manager to hold all the indices
  // Sets up the index list in the task manager
  std::size_t nvectors, nder = getPntrToArgument(getNumberOfArguments()-1)->getNumberOfStoredValues();
  MatrixTimesVectorData input;
  input.pairs.resize( getNumberOfArguments()-1, 2 );
  if( getPntrToArgument(1)->getRank()==2 ) {
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      input.pairs[i-1][0] = i-1;
      input.pairs[i-1][1] = getNumberOfArguments()-1;
    }
    nvectors=1;
    input.fshift=0;
  } else {
    for(unsigned i=1; i<getNumberOfArguments(); ++i) {
      input.pairs[i-1][0] = 0;
      input.pairs[i-1][1] = i;
    }
    nvectors = getNumberOfArguments()-1;
    input.fshift=nder;
  }
  taskmanager.setActionInput( input );
}

template <class CV, typename myPTM>
unsigned MatrixTimesVectorBase<CV, myPTM>::getNumberOfDerivatives() {
  unsigned nderivatives=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    nderivatives += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nderivatives;
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] ) {
    return;
  }
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  myval->setShape(shape);
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::calculate() {
  std::size_t nder = getPntrToArgument(getNumberOfArguments()-1)->getNumberOfStoredValues();
  std::size_t nvectors;
  if( getPntrToArgument(1)->getRank()==2 ) {
    nvectors = 1;
  } else {
    nvectors = getNumberOfArguments()-1;
  }
  if constexpr ( helpers::isRowSum<CV>) { //getName()=="MATRIX_VECTOR_PRODUCT_ROWSUMS" ) {
    taskmanager.setupParallelTaskManager( nder, 0 );
  } else {
    taskmanager.setupParallelTaskManager( 2*nder, nvectors*nder );
  }
  taskmanager.runAllTasks();
}

template <class CV, typename myPTM>
int MatrixTimesVectorBase<CV, myPTM>::checkTaskIsActive( const unsigned& itask ) const {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    Value* myarg = getPntrToArgument(i);
    if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
      return 0;
    } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
      unsigned ncol = myarg->getRowLength(itask);
      unsigned base = itask*myarg->getNumberOfColumns();
      for(unsigned k=0; k<ncol; ++k) {
        if( fabs(myarg->get(base+k,false))>0 ) {
          return 1;
        }
      }
    } else {
      plumed_merror("should not be in action " + getName() );
    }
  }
  return 0;
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::performTask( std::size_t task_index,
    const MatrixTimesVectorData& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {
  for(unsigned i=0; i<actiondata.pairs.nrows(); ++i) {
    MatrixTimesVectorOutput<precision> doutput( task_index,
        i,
        input.nderivatives_per_scalar,
        actiondata,
        input,
        output );
    CV::performTask( MatrixTimesVectorInput<precision>( task_index,
                     i,
                     actiondata,
                     input,
                     input.inputdata ),
                     doutput );
  }
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class CV, typename myPTM>
int MatrixTimesVectorBase<CV, myPTM>::getNumberOfValuesPerTask( std::size_t task_index,
    const MatrixTimesVectorData& actiondata ) {
  return 1;
}

template <class CV, typename myPTM>
void MatrixTimesVectorBase<CV, myPTM>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const MatrixTimesVectorData& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  for(unsigned i=0; i<actiondata.pairs.nrows(); ++i) {
    std::size_t base = input.argstarts[actiondata.pairs[i][0]]
                       + task_index*input.ncols[actiondata.pairs[i][0]];
    std::size_t n = input.bookeeping[input.bookstarts[actiondata.pairs[i][0]]
                                     + (1+input.ncols[actiondata.pairs[i][0]])*task_index];
    for(unsigned j=0; j<n; ++j) {
      force_indices.indices[i][j] = base + j;
    }
    force_indices.threadsafe_derivatives_end[i] = n;
    force_indices.tot_indices[i] = CV::getAdditionalIndices( n,
                                   input.argstarts[actiondata.pairs[i][1]],
                                   MatrixForceIndexInput( task_index, i, actiondata, input ),
                                   force_indices.indices[i] );
  }
}

template <typename T>
struct MatrixTimesVectorRowSums {
  typedef void isRowSums;
  typedef T precision;
  static void performTask( const MatrixTimesVectorInput<T>& input,
                           MatrixTimesVectorOutput<T>& output ) {
    for(unsigned i=0; i<input.rowlen; ++i) {
      output.values[0] += input.matrow[i];
      if( input.noderiv ) {
        continue;
      }
      output.matrow_deriv[i] = 1;
    }
  }
  static std::size_t getAdditionalIndices( std::size_t n,
      std::size_t vecstart,
      const MatrixForceIndexInput& fin,
      View<std::size_t> indices ) {
    return n;
  }
};


template <typename T>
struct MatrixTimesVectorProper {
  typedef T precision;
  static void performTask( const MatrixTimesVectorInput<T>& input,
                           MatrixTimesVectorOutput<T>& output ) {
    for(unsigned i=0; i<input.rowlen; ++i) {
      std::size_t ind = input.indices[i];
      output.values[0] += input.matrow[i]*input.vector[ind];
      if( input.noderiv ) {
        continue;
      }
      output.matrow_deriv[i] = input.vector[ind];
      output.vector_deriv[i] = input.matrow[i];
    }
  }
  static std::size_t getAdditionalIndices( std::size_t n,
      std::size_t vecstart,
      const MatrixForceIndexInput& fin,
      View<std::size_t> indices ) {
    for(unsigned i=0; i<fin.rowlen; ++i) {
      indices[n+i] = vecstart + fin.indices[i];
    }
    return n + fin.rowlen;
  }
};
}
}
#endif
