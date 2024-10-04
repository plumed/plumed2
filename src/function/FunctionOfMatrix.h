/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_function_FunctionOfMatrix_h
#define __PLUMED_function_FunctionOfMatrix_h

#include "core/ActionWithMatrix.h"
#include "FunctionOfVector.h"
#include "Sum.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfMatrix : public ActionWithMatrix {
private:
/// Is this the first step of the calculation
  bool firststep;
/// The function that is being computed
  T myfunc;
/// The number of derivatives for this action
  unsigned nderivatives;
/// A vector that tells us if we have stored the input value
  std::vector<bool> stored_arguments;
/// Switch off updating the arguments for this action
  std::vector<bool> update_arguments;
/// The list of actiosn in this chain
  std::vector<std::string> actionsLabelsInChain;
/// Get the shape of the output matrix
  std::vector<unsigned> getValueShapeFromArguments();
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfMatrix(const ActionOptions&);
/// Get the label to write in the graph
  std::string writeInGraph() const override {
    return myfunc.getGraphInfo( getName() );
  }
/// Make sure the derivatives are turned on
  void turnOnDerivatives() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Resize the matrices
  void prepare() override ;
/// This gets the number of columns
  unsigned getNumberOfColumns() const override ;
/// This checks for tasks in the parent class
//  void buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<AtomNumber>& otasks ) override ;
/// This ensures that we create some bookeeping stuff during the first step
  void setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) override ;
/// This sets up for the task
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
/// Calculate the full matrix
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override ;
/// This updates the indices for the matrix
//  void updateCentralMatrixIndex( const unsigned& ind, const std::vector<unsigned>& indices, MultiValue& myvals ) const override ;
  void runEndOfRowJobs( const unsigned& ind, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

template <class T>
void FunctionOfMatrix<T>::registerKeywords(Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.use("ARG");
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_MATRIX");
  keys.setDisplayName( name.substr(0,und) );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc;
  tfunc.registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("the sum of all the elements in the input matrix");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("the largest element of the input matrix");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("the smallest element in the input matrix");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("the matrix obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input matrix");
  }
}

template <class T>
FunctionOfMatrix<T>::FunctionOfMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  firststep(true) {
  if( myfunc.getArgStart()>0 ) {
    error("this has not beeen implemented -- if you are interested email gareth.tribello@gmail.com");
  }
  // Get the shape of the output
  std::vector<unsigned> shape( getValueShapeFromArguments() );
  // Check if the output matrix is symmetric
  bool symmetric=true;
  unsigned argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) {
      if( !getPntrToArgument(i)->isSymmetric() ) {
        symmetric=false;
      }
    }
  }
  // Read the input and do some checks
  myfunc.read( this );
  // Setup to do this in chain if possible
  if( myfunc.doWithTasks() ) {
    done_in_chain=true;
  }
  // Check we are not calculating a sum
  if( myfunc.zeroRank() ) {
    shape.resize(0);
  }
  // Get the names of the components
  std::vector<std::string> components( keywords.getOutputComponents() );
  // Create the values to hold the output
  std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
  for(unsigned i=0; i<components.size(); ++i) {
    if( str_ind.size()>0 ) {
      std::string compstr = components[i];
      if( components[i]==".#!value" ) {
        compstr = "";
      }
      for(unsigned j=0; j<str_ind.size(); ++j) {
        if( myfunc.zeroRank() ) {
          addComponentWithDerivatives( compstr + str_ind[j], shape );
        } else {
          addComponent( compstr + str_ind[j], shape );
          getPntrToComponent(i*str_ind.size()+j)->setSymmetric( symmetric );
        }
      }
    } else if( components[i]==".#!value" && myfunc.zeroRank() ) {
      addValueWithDerivatives( shape );
    } else if( components[i]==".#!value" ) {
      addValue( shape );
      getPntrToComponent(0)->setSymmetric( symmetric );
    } else if( components[i].find_first_of("_")!=std::string::npos ) {
      if( getNumberOfArguments()-argstart==1 ) {
        addValue( shape );
        getPntrToComponent(0)->setSymmetric( symmetric );
      } else {
        for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
          addComponent( getPntrToArgument(j)->getName() + components[i], shape );
          getPntrToComponent(i*(getNumberOfArguments()-argstart)+j-argstart)->setSymmetric( symmetric );
        }
      }
    } else {
      addComponent( components[i], shape );
      getPntrToComponent(i)->setSymmetric( symmetric );
    }
  }
  // Check if this can be sped up
  if( myfunc.getDerivativeZeroIfValueIsZero() )  {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
    }
  }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // We can't do this with if we are dividing a stack by some a product v.v^T product as we need to store the vector
  // In order to do this type of calculation.  There should be a neater fix than this but I can't see it.
  bool foundneigh=false;
  const ActionWithMatrix* chainstart = NULL;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isConstant() && getNumberOfArguments()>1 ) {
      continue ;
    }
    std::string argname=(getPntrToArgument(i)->getPntrToAction())->getName();
    if( argname=="NEIGHBORS" ) {
      foundneigh=true;
      break;
    }
    ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
    if( !av ) {
      done_in_chain=false;
    }
    if( getPntrToArgument(i)->getRank()==0 ) {
      function::FunctionOfVector<function::Sum>* as = dynamic_cast<function::FunctionOfVector<function::Sum>*>( getPntrToArgument(i)->getPntrToAction() );
      if(as) {
        done_in_chain=false;
      }
    } else if( getPntrToArgument(i)->ignoreStoredValue( getLabel() ) ) {
      // This option deals with the case when you have two adjacency matrices, A_ij and B_ij, multiplied together.  This cannot be done in the chain as the rows
      // of the two adjacency matrix are run over separately.  The value A_ij is thus not available when B_ij is calculated.
      ActionWithMatrix* am = dynamic_cast<ActionWithMatrix*>( getPntrToArgument(i)->getPntrToAction() );
      plumed_assert( am );
      const ActionWithMatrix* thischain = am->getFirstMatrixInChain();
      if( !thischain->isAdjacencyMatrix() && thischain->getName()!="VSTACK" ) {
        continue;
      }
      if( !chainstart ) {
        chainstart = thischain;
      } else if( thischain!=chainstart ) {
        done_in_chain=false;
      }
    }
  }
  // If we are working with neighbors we trick PLUMED into storing ALL the components of the other arguments
  // in this way we can ensure that the function of the neighbours matrix is in a chain starting from the
  // Neighbours matrix action.
  if( foundneigh ) {
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      ActionWithValue* av=getPntrToArgument(i)->getPntrToAction();
      if( av->getName()!="NEIGHBORS" ) {
        for(int i=0; i<av->getNumberOfComponents(); ++i) {
          (av->copyOutput(i))->buildDataStore();
        }
      }
    }
  }
  // Now setup the action in the chain if we can
  nderivatives = buildArgumentStore(myfunc.getArgStart());
}

template <class T>
void FunctionOfMatrix<T>::turnOnDerivatives() {
  if( !myfunc.derivativesImplemented() ) {
    error("derivatives have not been implemended for " + getName() );
  }
  ActionWithValue::turnOnDerivatives();
  myfunc.setup(this);
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfDerivatives() {
  return nderivatives;
}

template <class T>
void FunctionOfMatrix<T>::prepare() {
  unsigned argstart = myfunc.getArgStart();
  std::vector<unsigned> shape(2);
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) {
      shape[0] = getPntrToArgument(i)->getShape()[0];
      shape[1] = getPntrToArgument(i)->getShape()[1];
      break;
    }
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->getRank()==2 && (myval->getShape()[0]!=shape[0] || myval->getShape()[1]!=shape[1]) ) {
      myval->setShape(shape);
      if( myval->valueIsStored() ) {
        myval->reshapeMatrixStore( shape[1] );
      }
    }
  }
  ActionWithVector::prepare();
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfColumns() const {
  if( getConstPntrToComponent(0)->getRank()==2 ) {
    unsigned argstart=myfunc.getArgStart();
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==2 ) {
        ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>( getPntrToArgument(i)->getPntrToAction() );
        if( am ) {
          return am->getNumberOfColumns();
        }
        return getPntrToArgument(i)->getShape()[1];
      }
    }
  }
  plumed_error();
  return 0;
}

template <class T>
void FunctionOfMatrix<T>::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    plumed_assert( getPntrToArgument(i)->getRank()==2 );
  }
  unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(0)->getShape()[1];
  if( indices.size()!=size_v+1 ) {
    indices.resize( size_v+1 );
  }
  for(unsigned i=0; i<size_v; ++i) {
    indices[i+1] = start_n + i;
  }
  myvals.setSplitIndex( size_v + 1 );
}

// template <class T>
// void FunctionOfMatrix<T>::buildTaskListFromArgumentRequests( const unsigned& ntasks, bool& reduce, std::set<AtomNumber>& otasks ) {
//   // Check if this is the first element in a chain
//   if( actionInChain() ) return;
//   // If it is computed outside a chain get the tassks the daughter chain needs
//   propegateTaskListsForValue( 0, ntasks, reduce, otasks );
// }

template <class T>
void FunctionOfMatrix<T>::setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  if( firststep ) {
    stored_arguments.resize( getNumberOfArguments() );
    update_arguments.resize( getNumberOfArguments(), true );
    std::string control = getFirstActionInChain()->getLabel();
    for(unsigned i=0; i<stored_arguments.size(); ++i) {
      stored_arguments[i] = !getPntrToArgument(i)->ignoreStoredValue( control );
      if( !stored_arguments[i] ) {
        update_arguments[i] = true;
      } else {
        update_arguments[i] = !argumentDependsOn( headstr, this, getPntrToArgument(i) );
      }
    }
    firststep=false;
  }
  ActionWithMatrix::setupStreamedComponents( headstr, nquants, nmat, maxcol, nbookeeping );
}

template <class T>
void FunctionOfMatrix<T>::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned argstart=myfunc.getArgStart();
  std::vector<double> args( getNumberOfArguments() - argstart );
  unsigned ind2 = index2;
  if( getConstPntrToComponent(0)->getRank()==2 && index2>=getConstPntrToComponent(0)->getShape()[0] ) {
    ind2 = index2 - getConstPntrToComponent(0)->getShape()[0];
  } else if( index2>=getPntrToArgument(0)->getShape()[0] ) {
    ind2 = index2 - getPntrToArgument(0)->getShape()[0];
  }
  if( actionInChain() ) {
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==0 ) {
        args[i-argstart] = getPntrToArgument(i)->get();
      } else if( !getPntrToArgument(i)->valueHasBeenSet() ) {
        args[i-argstart] = myvals.get( getPntrToArgument(i)->getPositionInStream() );
      } else {
        args[i-argstart] = getPntrToArgument(i)->get( getPntrToArgument(i)->getShape()[1]*index1 + ind2 );
      }
    }
  } else {
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==2 ) {
        args[i-argstart]=getPntrToArgument(i)->get( getPntrToArgument(i)->getShape()[1]*index1 + ind2 );
      } else {
        args[i-argstart] = getPntrToArgument(i)->get();
      }
    }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() );
  Matrix<double> derivatives( getNumberOfComponents(), getNumberOfArguments()-argstart );
  myfunc.calc( this, args, vals, derivatives );
  // And set the values
  for(unsigned i=0; i<vals.size(); ++i) {
    myvals.addValue( getConstPntrToComponent(i)->getPositionInStream(), vals[i] );
  }
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) {
    return;
  }

  if( actionInChain() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
      for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
        if( getPntrToArgument(j)->getRank()==2 ) {
          unsigned istrn = getPntrToArgument(j)->getPositionInStream();
          if( stored_arguments[j] ) {
            unsigned task_index = getPntrToArgument(i)->getShape()[1]*index1 + ind2;
            myvals.clearDerivatives(istrn);
            myvals.addDerivative( istrn, task_index, 1.0 );
            myvals.updateIndex( istrn, task_index );
          }
          for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
            unsigned kind=myvals.getActiveIndex(istrn,k);
            myvals.addDerivative( ostrn, arg_deriv_starts[j] + kind, derivatives(i,j)*myvals.getDerivative( istrn, kind ) );
          }
        }
      }
    }
    // If we are computing a matrix we need to update the indices here so that derivatives are calcualted correctly in functions of these
    if( getConstPntrToComponent(0)->getRank()==2 ) {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
        for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
          if( !update_arguments[j] || getPntrToArgument(j)->getRank()==0 ) {
            continue ;
          }
          // Ensure we only store one lot of derivative indices
          bool found=false;
          for(unsigned k=0; k<j; ++k) {
            if( arg_deriv_starts[k]==arg_deriv_starts[j] ) {
              found=true;
              break;
            }
          }
          if( found ) {
            continue;
          }
          unsigned istrn = getPntrToArgument(j)->getPositionInStream();
          for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
            unsigned kind=myvals.getActiveIndex(istrn,k);
            myvals.updateIndex( ostrn, arg_deriv_starts[j] + kind );
          }
        }
      }
    }
  } else {
    unsigned base=0;
    ind2 = index2;
    for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
      if( getPntrToArgument(j)->getRank()!=2 ) {
        continue ;
      }
      if( index2>=getPntrToArgument(j)->getShape()[0] ) {
        ind2 = index2 - getPntrToArgument(j)->getShape()[0];
      }
      break;
    }
    for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
      if( getPntrToArgument(j)->getRank()==2 ) {
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          unsigned myind = base + getPntrToArgument(j)->getShape()[1]*index1 + ind2;
          myvals.addDerivative( ostrn, myind, derivatives(i,j) );
          myvals.updateIndex( ostrn, myind );
        }
      } else {
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          myvals.addDerivative( ostrn, base, derivatives(i,j) );
          myvals.updateIndex( ostrn, base );
        }
      }
      base += getPntrToArgument(j)->getNumberOfValues();
    }
  }
}

template <class T>
void FunctionOfMatrix<T>::runEndOfRowJobs( const unsigned& ind, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) {
    return;
  }

  unsigned argstart=myfunc.getArgStart();
  if( actionInChain() && getConstPntrToComponent(0)->getRank()==2 ) {
    // This is triggered if we are outputting a matrix
    for(int vv=0; vv<getNumberOfComponents(); ++vv) {
      unsigned nmat = getConstPntrToComponent(vv)->getPositionInMatrixStash();
      std::vector<unsigned>& mat_indices( myvals.getMatrixRowDerivativeIndices( nmat ) );
      unsigned ntot_mat=0;
      if( mat_indices.size()<nderivatives ) {
        mat_indices.resize( nderivatives );
      }
      for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
        if( !update_arguments[i] || getPntrToArgument(i)->getRank()==0 ) {
          continue ;
        }
        // Ensure we only store one lot of derivative indices
        bool found=false;
        for(unsigned j=0; j<i; ++j) {
          if( arg_deriv_starts[j]==arg_deriv_starts[i] ) {
            found=true;
            break;
          }
        }
        if( found ) {
          continue;
        }

        if( stored_arguments[i] ) {
          unsigned tbase = getPntrToArgument(i)->getShape()[1]*ind;
          for(unsigned k=1; k<indices.size(); ++k) {
            unsigned ind2 = indices[k] - getConstPntrToComponent(0)->getShape()[0];
            mat_indices[ntot_mat + k - 1] = arg_deriv_starts[i] + tbase + ind2;
          }
          ntot_mat += indices.size()-1;
        } else {
          unsigned istrn = getPntrToArgument(i)->getPositionInMatrixStash();
          std::vector<unsigned>& imat_indices( myvals.getMatrixRowDerivativeIndices( istrn ) );
          for(unsigned k=0; k<myvals.getNumberOfMatrixRowDerivatives( istrn ); ++k) {
            mat_indices[ntot_mat + k] = arg_deriv_starts[i] + imat_indices[k];
          }
          ntot_mat += myvals.getNumberOfMatrixRowDerivatives( istrn );
        }
      }
      myvals.setNumberOfMatrixRowDerivatives( nmat, ntot_mat );
    }
  } else if( actionInChain() ) {
    // This is triggered if we are calculating a single scalar in the function
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      bool found=false;
      for(unsigned j=0; j<i; ++j) {
        if( arg_deriv_starts[j]==arg_deriv_starts[i] ) {
          found=true;
          break;
        }
      }
      if( found ) {
        continue;
      }
      unsigned istrn = getPntrToArgument(i)->getPositionInMatrixStash();
      std::vector<unsigned>& mat_indices( myvals.getMatrixRowDerivativeIndices( istrn ) );
      for(unsigned k=0; k<myvals.getNumberOfMatrixRowDerivatives( istrn ); ++k) {
        for(int j=0; j<getNumberOfComponents(); ++j) {
          unsigned ostrn = getConstPntrToComponent(j)->getPositionInStream();
          myvals.updateIndex( ostrn, arg_deriv_starts[i] + mat_indices[k] );
        }
      }
    }
  } else if( getConstPntrToComponent(0)->getRank()==2 ) {
    for(int vv=0; vv<getNumberOfComponents(); ++vv) {
      unsigned nmat = getConstPntrToComponent(vv)->getPositionInMatrixStash();
      std::vector<unsigned>& mat_indices( myvals.getMatrixRowDerivativeIndices( nmat ) );
      unsigned ntot_mat=0;
      if( mat_indices.size()<nderivatives ) {
        mat_indices.resize( nderivatives );
      }
      unsigned matderbase = 0;
      for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
        if( getPntrToArgument(i)->getRank()==0 ) {
          continue ;
        }
        unsigned ss = getPntrToArgument(i)->getShape()[1];
        unsigned tbase = matderbase + ss*myvals.getTaskIndex();
        for(unsigned k=0; k<ss; ++k) {
          mat_indices[ntot_mat + k] = tbase + k;
        }
        ntot_mat += ss;
        matderbase += getPntrToArgument(i)->getNumberOfValues();
      }
      myvals.setNumberOfMatrixRowDerivatives( nmat, ntot_mat );
    }
  }
}

template <class T>
std::vector<unsigned> FunctionOfMatrix<T>::getValueShapeFromArguments() {
  unsigned argstart=myfunc.getArgStart();
  std::vector<unsigned> shape(2);
  shape[0]=shape[1]=0;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    plumed_assert( getPntrToArgument(i)->getRank()==2 || getPntrToArgument(i)->getRank()==0 );
    if( getPntrToArgument(i)->getRank()==2 ) {
      if( shape[0]>0 && (getPntrToArgument(i)->getShape()[0]!=shape[0] || getPntrToArgument(i)->getShape()[1]!=shape[1]) ) {
        error("all matrices input should have the same shape");
      } else if( shape[0]==0 ) {
        shape[0]=getPntrToArgument(i)->getShape()[0];
        shape[1]=getPntrToArgument(i)->getShape()[1];
      }
      plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
    }
  }
  myfunc.setPrefactor( this, 1.0 );
  return shape;
}

}
}
#endif
