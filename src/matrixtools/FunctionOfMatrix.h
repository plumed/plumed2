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
#ifndef __PLUMED_matrixtools_FunctionOfMatrix_h
#define __PLUMED_matrixtools_FunctionOfMatrix_h

#include "adjmat/MatrixProductBase.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace matrixtools {

template <class T>
class FunctionOfMatrix : public adjmat::MatrixProductBase {
private:
/// The function that is being computed
  T myfunc;
/// The number of derivatives for this action
  unsigned nderivatives;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfMatrix(const ActionOptions&);
/// Get the shape of the output matrix
  std::vector<unsigned> getMatrixShapeForFinalTasks() override ;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const override ;
/// This gets the number of columns
  unsigned getNumberOfColumns() const override ;
/// This is not used
  double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                               const std::vector<double>& vec1, const std::vector<double>& vec2,
                               std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const { plumed_error(); }
/// Calculate the full matrix
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override ;
/// This updates the indices for the matrix
  void updateCentralMatrixIndex( const unsigned& ind, const std::vector<unsigned>& indices, MultiValue& myvals ) const override ;
};

template <class T>
void FunctionOfMatrix<T>::registerKeywords(Keywords& keys ) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc; tfunc.registerKeywords( keys );
}

template <class T>
FunctionOfMatrix<T>::FunctionOfMatrix(const ActionOptions&ao):
Action(ao),
MatrixProductBase(ao),
nderivatives(getNumberOfScalarArguments())
{
  // Get the shape of the output
  std::vector<unsigned> shape( getMatrixShapeForFinalTasks() );
  // Create the task list
  for(unsigned i=0;i<shape[0];++i) addTaskToList(i);
  // Read the input and do some checks
  myfunc.read( this );
  // Check we are not calculating a sum
  if( myfunc.getRank()==0 ) shape.resize(0);  
  // Get the names of the components
  std::vector<std::string> components( keywords.getAllOutputComponents() );
  // Create the values to hold the output
  if( components.size()==0 && myfunc.getRank()==0 ) addValueWithDerivatives( shape );
  else if( components.size()==0 ) addValue( shape );
  else { for(unsigned i=0;i<components.size();++i) addComponent( components[i], shape ); }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // Now setup the action in the chain if we can
  if( distinct_arguments.size()>0 ) nderivatives = setupActionInChain(0); 
  // Fix any functions that are time series
  if( input_timeseries ) {
      for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->makeTimeSeries();
  }
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfDerivatives() const {
  return nderivatives;
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfColumns() const {
  if( getPntrToOutput(0)->getRank()==2 ) {
     for(unsigned i=0;i<getNumberOfArguments();++i) {
         if( getPntrToArgument(0)->getRank()==2 ) return getPntrToArgument(0)->getNumberOfColumns();
     }
  }
  plumed_error(); return 0;
}

template <class T>
bool FunctionOfMatrix<T>::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  std::vector<double> args( getNumberOfArguments() );
  if( actionInChain() ) {
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()==2 ) args[i] = myvals.get( getPntrToArgument(i)->getPositionInStream() );
          else args[i] = getPntrToArgument(i)->get();
      }
  } else {
      unsigned ind2 = index2; if( index2>=getFullNumberOfTasks() ) ind2 = index2 - getFullNumberOfTasks();
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()==2 ) args[i]=getPntrToArgument(i)->get( getPntrToArgument(i)->getShape()[1]*index2 + ind2 );
          else args[i] = getPntrToArgument(i)->get();
      }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() ); Matrix<double> derivatives( getNumberOfComponents(), getNumberOfArguments() );
  myfunc.calc( args, vals, derivatives );
  // And set the values
  for(unsigned i=0;i<vals.size();++i) myvals.addValue( getPntrToOutput(i)->getPositionInStream(), vals[i] );
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) return true;

  if( actionInChain() ) {
      for(unsigned i=0;i<getNumberOfComponents();++i) {
          unsigned ostrn=getPntrToOutput(i)->getPositionInStream();
          for(unsigned j=0;j<getNumberOfArguments();++j) {
              if( getPntrToArgument(i)->getRank()==2 ) {
                  unsigned istrn = getPntrToArgument(j)->getPositionInStream();
                  for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
                      unsigned kind=myvals.getActiveIndex(istrn,k);
                      myvals.addDerivative( ostrn, arg_deriv_starts[j] + kind, derivatives(i,j)*myvals.getDerivative( istrn, kind ) );
                  }
              } else plumed_merror("not implemented this yet");
          }
      }
  } else plumed_merror("not implemented this yet");
  return true;
}

template <class T>
void FunctionOfMatrix<T>::updateCentralMatrixIndex( const unsigned& ind, const std::vector<unsigned>& indices, MultiValue& myvals ) const {
  if( actionInChain() && getPntrToOutput(0)->getRank()==2 ) {
      // This is triggered if we are outputting a matrix
      for(unsigned vv=0; vv<getNumberOfComponents(); ++vv) {
          unsigned nmat = getPntrToOutput(vv)->getPositionInMatrixStash();
          std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( nmat ) ); unsigned ntot_mat=0;
          if( mat_indices.size()<getNumberOfDerivatives() ) mat_indices.resize( getNumberOfDerivatives() );
          for(unsigned i=0; i<getNumberOfArguments(); ++i) {
            if( getPntrToArgument(i)->getRank()==0 ) continue ;
            // Ensure we only store one lot of derivative indices
            bool found=false;
            for(unsigned j=0; j<i; ++j) {
                if( arg_deriv_starts[j]==arg_deriv_starts[i] ) { found=true; break; }
            } 
            if( found ) continue;
            unsigned istrn = getPntrToArgument(i)->getPositionInMatrixStash();
            std::vector<unsigned>& imat_indices( myvals.getMatrixIndices( istrn ) );
            for(unsigned k=0; k<myvals.getNumberOfMatrixIndices( istrn ); ++k) mat_indices[ntot_mat + k] = arg_deriv_starts[i] + imat_indices[k];
            ntot_mat += myvals.getNumberOfMatrixIndices( istrn ); 
          }
          myvals.setNumberOfMatrixIndices( nmat, ntot_mat );
      }
  } else if( actionInChain() ) {
      // This is triggered if we are calculating a single scalar in the function
      for(unsigned i=0; i<getNumberOfArguments(); ++i) {
          bool found=false;
          for(unsigned j=0; j<i; ++j) {
              if( arg_deriv_starts[j]==arg_deriv_starts[i] ) { found=true; break; }
          } 
          if( found ) continue; 
          unsigned istrn = getPntrToArgument(i)->getPositionInMatrixStash();
          std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );  
          for(unsigned k=0; k<myvals.getNumberOfMatrixIndices( istrn ); ++k) {
              for(unsigned j=0; j<getNumberOfComponents(); ++j) {
                  unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
                   myvals.updateIndex( ostrn, arg_deriv_starts[i] + mat_indices[k] );
              }

          }
      }
  } else plumed_merror("I have not implemented this as I am not sure it is needed");
}

template <class T>
std::vector<unsigned> FunctionOfMatrix<T>::getMatrixShapeForFinalTasks() { 
  std::vector<unsigned> shape(2); shape[0]=shape[1]=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      plumed_assert( getPntrToArgument(i)->getRank()==2 || getPntrToArgument(i)->getRank()==0 );
      if( getPntrToArgument(i)->getRank()==2 ) {
          if( shape[0]>0 && (getPntrToArgument(i)->getShape()[0]!=shape[0] || getPntrToArgument(i)->getShape()[1]!=shape[1]) ) error("all matrices input should have the same shape");
          else if( shape[0]==0 ) { shape[0]=getPntrToArgument(i)->getShape()[0]; shape[1]=getPntrToArgument(i)->getShape()[1]; }
          plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
      }   
  }
  myfunc.setPrefactor( this ); return shape;
}

}
}
#endif
