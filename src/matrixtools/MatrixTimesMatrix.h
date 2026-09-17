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
#ifndef __PLUMED_matrixtools_MatrixTimesMatrix_h
#define __PLUMED_matrixtools_MatrixTimesMatrix_h

#include "core/ActionWithMatrix.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace matrixtools {

namespace helpers {
template <typename T, typename=void>
constexpr bool isDissimilarities=false;
template <typename T>
constexpr bool isDissimilarities<T,std::void_t<typename T::isDissimilarities>> =true;
}

template <typename T>
struct MatrixTimesMatrixInput {
  T funcinput;
  bool no_thread_gather;
  bool gatherForceOnColumns;
  bool secondMatrixIsSparse;
  std::vector<unsigned> locations;
  RequiredMatrixElements outmat;
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    funcinput.toACCDevice();
    outmat.toACCDevice();
  }
  void removeFromACCDevice() const {
    funcinput.removeFromACCDevice();
    outmat.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

class InputVectors {
public:
  std::size_t nelem;
  View<double> arg1;
  View<double> arg2;
  InputVectors( std::size_t n,  double* b ) : nelem(n), arg1(b,n), arg2(b+n,n) {}
};

template <class CV, typename myPTM=defaultPTM>
class MatrixTimesMatrix : public ActionWithMatrix {
public:
  using input_type = MatrixTimesMatrixInput<CV>;
  using mytype = MatrixTimesMatrix<CV, myPTM>;
  using PTM = typename myPTM::template PTM<mytype>;
  typedef typename PTM::ParallelActionsInput ParallelActionsInput;
  typedef typename PTM::ParallelActionsOutput ParallelActionsOutput;
  constexpr static bool isDissimilarities=helpers::isDissimilarities<CV>;
private:
  PTM taskmanager;
  void getBasicMatrixBookeeping( ArgumentsBookkeeping& argumentsMap ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesMatrix(const ActionOptions&);
  void prepare() override ;
  unsigned getNumberOfDerivatives() override;
  void calculate() override ;
  void getInputData( std::vector<double>& inputdata, ArgumentsBookkeeping& argumentsMap ) const override ;
  void getInputData( std::vector<float>& inputdata, ArgumentsBookkeeping& argumentsMap ) const override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index,
                           const input_type& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const input_type& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const input_type& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","the label of the two matrices from which the product is calculated");
  //if( keys.getDisplayName()=="MATRIX_PRODUCT" ) {
  //  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
  //}
  CV::registerKeywords( keys );
  PTM::registerKeywords( keys );
}

template <class CV, typename myPTM>
MatrixTimesMatrix<CV, myPTM>::MatrixTimesMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  taskmanager(this) {
  int nm=getNumberOfMasks();
  if( nm<0 ) {
    nm = 0;
  }
  if( getNumberOfArguments()-nm!=2 ) {
    error("should be two arguments to this action, a matrix and a vector");
  }
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
    error("first argument to this action should be a matrix");
  }
  if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) {
    error("second argument to this action should be a matrix");
  }
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) {
    error("number of columns in first matrix does not equal number of rows in second matrix");
  }
  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  if( !isDissimilarities && getPntrToArgument(0)->isDerivativeZeroWhenValueIsZero() && getPntrToArgument(1)->isDerivativeZeroWhenValueIsZero() ) {
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }

  if( nm>0 ) {
    unsigned iarg = getNumberOfArguments()-1;
    if( getPntrToArgument(iarg)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
      error("argument passed to MASK keyword should be a matrix");
    }
    if( getPntrToArgument(iarg)->getShape()[0]!=shape[0] || getPntrToArgument(iarg)->getShape()[1]!=shape[1] ) {
      error("argument passed to MASK keyword has the wrong shape");
    }
  }
  input_type actdata;
  actdata.funcinput.setup( this, getPntrToArgument(0) );
  actdata.gatherForceOnColumns = false;
  actdata.no_thread_gather = no_thread_gather;
  taskmanager.setActionInput( actdata );
}

template <class CV, typename myPTM>
unsigned MatrixTimesMatrix<CV, myPTM>::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues();
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[1] ) {
    return;
  }
  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  myval->setShape(shape);
  myval->reshapeMatrixStore( shape[1] );
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::calculate() {
  if( !getPntrToComponent(0)->isDerivativeZeroWhenValueIsZero() ) {
    if( getPntrToArgument(0)->getNumberOfColumns()<getPntrToArgument(0)->getShape()[1] ) {
      if( !doNotCalculateDerivatives() ) {
        error("cannot calculate derivatives for this action with sparse matrices");
      } else if(isDissimilarities ) {
        error("cannot calculate dissimilarities for sparse matrices");
      }
    }
    if( getPntrToArgument(1)->getNumberOfColumns()<getPntrToArgument(1)->getShape()[1] ) {
      if( !doNotCalculateDerivatives() ) {
        error("cannot calculate derivatives for this action with sparse matrices");
      } else if( isDissimilarities ) {
        error("cannot calculate dissimilarities for sparse matrices");
      }
    }
  }
  if( getNumberOfMasks()>0 || diagzero ) {
      updateBookeepingArrays( taskmanager.getActionInput().outmat );
  } else {
      // Determine the sparsity pattern of the output matrix from the input matrices
      Value* arg2 = getPntrToArgument(1);
      std::vector<unsigned> column_lengths( arg2->getShape()[1], 0 );
      for(unsigned i=0; i<arg2->getShape()[0]; ++i) {
          unsigned nr = arg2->getRowLength(i);
          for(unsigned j=0; j<nr; ++j) {
            unsigned ind = arg2->getRowIndex( i, j );
            column_lengths[ind]++;
        }
      }
      unsigned maxrow = 0;
      Value* arg1 = getPntrToArgument(0);
      for(unsigned i=0; i<arg1->getShape()[0]; ++i) {
          if( arg1->getRowLength(i)==0 ) continue;
          unsigned rowlen = 0;
          for(unsigned j=0; j<arg2->getShape()[1]; ++j) {
              if( column_lengths[j]>0 ) {
                  rowlen++;
              }
          }
          if( rowlen>maxrow ) {
              maxrow = rowlen;
          } 
      }
      Value* myval = getPntrToComponent(0);
      myval->reshapeMatrixStore( maxrow );
      for(unsigned i=0; i<arg1->getShape()[0]; ++i) {
          unsigned rstart = i*(maxrow+1);
          if( arg1->getRowLength(i)==0 ) {
              myval->setMatrixBookeepingElement( rstart, 0 );
              continue;
          }
          unsigned rowlen=0;
          for(unsigned j=0; j<arg2->getShape()[1]; ++j) {
              if( column_lengths[j]>0 ) {
                  myval->setMatrixBookeepingElement( rstart + 1 + rowlen, j );
                  rowlen++;
              }
          }
          myval->setMatrixBookeepingElement( rstart, rowlen );
      }
      copyMatrixBookeepingFromFirstComponent( taskmanager.getActionInput().outmat );
  }
  taskmanager.getActionInput().secondMatrixIsSparse = getPntrToArgument(1)->getNumberOfColumns()<getPntrToArgument(1)->getShape()[1];
  if( no_thread_gather ) {
      taskmanager.setupParallelTaskManager( getPntrToArgument(0)->getNumberOfColumns(), 0 );
      taskmanager.setWorkspaceSize( 4*getPntrToArgument(0)->getNumberOfColumns() ); 
  } else {
      taskmanager.setupParallelTaskManager( 2*getPntrToArgument(0)->getNumberOfColumns(),
                                            getPntrToArgument(1)->getNumberOfStoredValues() );
      taskmanager.setWorkspaceSize( 2*getPntrToArgument(0)->getNumberOfColumns() );
  }
  taskmanager.runAllTasks();
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::getBasicMatrixBookeeping( ArgumentsBookkeeping& argumentsMap ) const {
   argumentsMap.nargs = 2;
   argumentsMap.ranks.resize(2);
   argumentsMap.ranks[0] = argumentsMap.ranks[1] = 2;
   argumentsMap.shapestarts.resize(2);
   argumentsMap.shapestarts[0] = 0;
   argumentsMap.shapestarts[1] = 2;
   argumentsMap.ncols.resize(2);
   argumentsMap.ncols[0] = getPntrToArgument(0)->getNumberOfColumns();
   argumentsMap.shapedata.resize(4);
   for(unsigned i=0; i<2; ++i) {
       for(unsigned j=0; j<2; ++j) {
           argumentsMap.shapedata[2*i+j] = getPntrToArgument(i)->getShape()[j];
       }
   }
   argumentsMap.argstarts.resize(2);
   argumentsMap.argstarts[0] = 0;
   argumentsMap.argstarts[1] = getPntrToArgument(0)->getNumberOfStoredValues();
            
   // Find the lengths of all the columns
   Value* arg0 = getPntrToArgument(0); 
   Value* arg1 = getPntrToArgument(1);
   argumentsMap.ncols[1] = arg1->getLengthOfLongestColumn();
   argumentsMap.bookstarts.resize(2);
   argumentsMap.bookstarts[0] = 0;
   argumentsMap.bookstarts[1] = arg0->getShape()[1]*(1+argumentsMap.ncols[0]);
   argumentsMap.booksizes.resize(2);
   argumentsMap.booksizes[0] = arg0->getShape()[1]*(1+argumentsMap.ncols[0]);
   argumentsMap.booksizes[1] = arg1->getShape()[1]*(1+argumentsMap.ncols[1]);
   argumentsMap.bookeeping.resize( argumentsMap.booksizes[0] + argumentsMap.booksizes[1] );
   for(unsigned i=0; i<arg0->getShape()[0]; ++i) {
       unsigned nr = arg0->getRowLength(i);
       unsigned base = i*(argumentsMap.ncols[0]+1);
       argumentsMap.bookeeping[base] = nr;
       for(unsigned j=0; j<nr; ++j) {
           argumentsMap.bookeeping[base+1+j] = arg0->getRowIndex( i, j );
       }
   }
   for(unsigned i=0; i<arg1->getShape()[1]; ++i) {
       argumentsMap.bookeeping[argumentsMap.bookstarts[1] + i*(argumentsMap.ncols[1]+1)] = 0;
   } 
} 

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::getInputData( std::vector<double>& inputdata, ArgumentsBookkeeping& argumentsMap ) const {
   getBasicMatrixBookeeping( argumentsMap ); 
   std::size_t total_args = getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getShape()[1]*argumentsMap.ncols[1];
   if( inputdata.size()!=total_args ) {
       inputdata.resize( total_args );
   }
   getPntrToArgument(0)->assignValues(View{&inputdata[0], getPntrToArgument(0)->getNumberOfStoredValues()});
   // The transpose of the input matrix is stored here to make it easier to retrieve the columns

   Value* arg1 = getPntrToArgument(1);
   unsigned arg_ncol = arg1->getNumberOfColumns();
   for(unsigned i=0; i<arg1->getShape()[0]; ++i) {
       unsigned nr = arg1->getRowLength(i);
       for(unsigned j=0; j<nr; ++j) {
           unsigned ind = arg1->getRowIndex( i, j );
           unsigned startpos = argumentsMap.bookstarts[1] + ind*(1+argumentsMap.ncols[1]); 
           argumentsMap.bookeeping[startpos + 1 + argumentsMap.bookeeping[startpos] ] = i; 
           inputdata[ argumentsMap.argstarts[1] + ind*argumentsMap.ncols[1] + argumentsMap.bookeeping[startpos] ] = arg1->get( i*arg_ncol+j, false );
           argumentsMap.bookeeping[startpos]++;
       }
   }
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::getInputData( std::vector<float>& inputdata, ArgumentsBookkeeping& argumentsMap ) const {
   getBasicMatrixBookeeping( argumentsMap );
   std::size_t total_args = getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getShape()[1]*argumentsMap.ncols[1];
   if( inputdata.size()!=total_args ) {
       inputdata.resize( total_args );
   }
   Value* arg0 = getPntrToArgument(0);
   for(unsigned i=0; i<arg0->getNumberOfStoredValues(); ++i){
       inputdata[i] = arg0->get(i,false);
   }
   Value* arg1 = getPntrToArgument(1);
   unsigned arg_ncol = arg1->getNumberOfColumns();
   for(unsigned i=0; i<arg1->getShape()[0]; ++i) {
       unsigned nr = arg1->getRowLength(i);
       for(unsigned j=0; j<nr; ++j) {
           unsigned ind = arg1->getRowIndex( i, j );
           unsigned startpos = argumentsMap.bookstarts[1] + ind*(1+argumentsMap.ncols[1]);
           argumentsMap.bookeeping[startpos + 1 + argumentsMap.bookeeping[startpos] ] = i;
           inputdata[ argumentsMap.argstarts[1] + ind*argumentsMap.ncols[1] + argumentsMap.bookeeping[startpos] ] = arg1->get( i*arg_ncol+j, false );
           argumentsMap.bookeeping[startpos]++;
       }
   }
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::performTask( std::size_t task_index,
    const input_type& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {
  auto arg0=ArgumentBookeepingHolder::create( 0, input );
  auto arg1=ArgumentBookeepingHolder::create( 1, input );
  // This is the number of matrix elements that we compute on this loop
  std::size_t arg0_start = task_index*arg0.ncols;
  std::size_t arg0_bookstart = task_index*(1+arg0.ncols);
  std::size_t nmult = arg0.bookeeping[arg0_bookstart];
  std::size_t outrowstart = task_index*(1+actiondata.outmat.ncols);
  std::size_t output_rowlen = actiondata.outmat[outrowstart];
  InputVectors vectors( nmult, output.buffer.data() );
  if( actiondata.secondMatrixIsSparse ) {
      if( actiondata.no_thread_gather ) {
          std::size_t arg1_start = arg1.start + task_index*arg1.ncols;
          std::size_t arg1_bookstart = task_index*(1+arg1.ncols);
          for(unsigned i=0; i<output_rowlen; ++i) {
              if( actiondata.gatherForceOnColumns ) {
                  arg0_bookstart = actiondata.outmat[outrowstart+1+i]*(1+arg0.ncols);
                  arg0_start = actiondata.outmat[outrowstart+1+i]*arg0.ncols;
              } else {
                  arg1_bookstart = actiondata.outmat[outrowstart+1+i]*(1+arg1.ncols);
                  arg1_start = arg1.start + actiondata.outmat[outrowstart+1+i]*arg1.ncols;
              }           
              std::size_t nm = 0;  
              std::size_t arg1_nelements = arg1.bookeeping[arg1_bookstart];
              for(unsigned j=0; j<nmult; ++j) {
                  std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
                  for(unsigned k=0; k<arg1_nelements; ++k) {
                      if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                          nm++;
                          break;
                      }
                  }
              }
              if( nm==0 ) {
                  continue ;
              }
              vectors.nelem = nm;
              nm = 0;
              for(unsigned j=0; j<nmult; ++j) {
                  std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
                  for(unsigned k=0; k<arg1_nelements; ++k) {
                      if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                          vectors.arg1[nm] = input.inputdata[ arg0_start + j ];
                          vectors.arg2[nm] = input.inputdata[ arg1_start + k ];
                          nm++;
                          break;
                      }
                  }
              }
              MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.buffer.data() + 2*nmult );
              CV::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
              if( actiondata.gatherForceOnColumns ) {
                  for(unsigned k=0; k<nm; ++k) {
                      output.derivatives[nmult*i + k] = elem.derivs[0][nmult+k];
                  } 
              } else {
                  for(unsigned k=0; k<nm; ++k) {
                      output.derivatives[nmult*i + k] = elem.derivs[0][k];
                  }
              }
          }
      } else {
          for(unsigned i=0; i<output_rowlen; ++i) {
              std::size_t nm = 0;
              std::size_t arg1_bookstart = actiondata.outmat[outrowstart+1+i]*(1+arg1.ncols);
              std::size_t arg1_nelements = arg1.bookeeping[arg1_bookstart];
              for(unsigned j=0; j<nmult; ++j) {
                  std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
                  for(unsigned k=0; k<arg1_nelements; ++k) {
                      if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                          nm++;
                          break;
                      }
                  }
              }
              if( nm==0 ) {
                  continue ;
              }
              vectors.nelem = nm;
              nm = 0;
              std::size_t arg1_start = arg1.start + actiondata.outmat[outrowstart+1+i]*arg1.ncols; 
              for(unsigned j=0; j<nmult; ++j) {
                  std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
                  for(unsigned k=0; k<arg1_nelements; ++k) {
                      if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                          vectors.arg1[nm] = input.inputdata[ arg0_start + j ];
                          vectors.arg2[nm] = input.inputdata[ arg1_start + k ];
                          nm++;
                          break;
                      }
                  }
              }
              MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.derivatives.data() + 2*nmult*i );
              CV::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
              for(unsigned ii=vectors.nelem; ii<nmult; ++ii) {
                  elem.derivs[0][ii] = 0;
              }
          }
      }
  } else {
      if( actiondata.no_thread_gather ) { 
         std::size_t arg1_start = arg1.start + task_index*arg1.ncols;
         for(unsigned i=0; i<output_rowlen; ++i) {
            if( actiondata.gatherForceOnColumns ) {
                arg0_bookstart = actiondata.outmat[outrowstart+1+i]*(1+arg0.ncols);
                arg0_start = actiondata.outmat[outrowstart+1+i]*arg0.ncols;
            } else {
                arg1_start = arg1.start + actiondata.outmat[outrowstart+1+i]*arg1.ncols;
            }
            for(unsigned j=0; j<nmult; ++j) {
                vectors.arg1[j] = input.inputdata[arg0_start + j];
                vectors.arg2[j] = input.inputdata[arg1_start + arg0.bookeeping[arg0_bookstart+1+j]];
            }
            MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.buffer.data() + 2*nmult );
            CV::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
            if( actiondata.gatherForceOnColumns ) {
                for(unsigned k=0; k<nmult; ++k) {
                    output.derivatives[nmult*i + k] = elem.derivs[0][nmult+k];
                }       
            } else {    
                for(unsigned k=0; k<nmult; ++k) {
                    output.derivatives[nmult*i + k] = elem.derivs[0][k];
                }   
            } 
         }
      } else {
         for(unsigned i=0; i<nmult; ++i){
             vectors.arg1[i] = input.inputdata[arg0_start + i];
         }
         for(unsigned i=0; i<output_rowlen; ++i) {
            std::size_t arg1_start = arg1.start + actiondata.outmat[outrowstart+1+i]*arg1.ncols;
            for(unsigned j=0; j<nmult; ++j) {
                vectors.arg2[j] = input.inputdata[arg1_start + arg0.bookeeping[arg0_bookstart+1+j]];
            }
            MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.derivatives.data() + 2*nmult*i );
            CV::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
         }
      }
  }
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  Value* arg1 = getPntrToArgument(1);
  if( arg1->getNumberOfColumns()<arg1->getShape()[1] ){
      std::size_t maxcol = arg1->getLengthOfLongestColumn();
      std::vector<unsigned>& loc = taskmanager.getActionInput().locations;
      if( loc.size()!=(1+maxcol)*arg1->getShape()[1] ) {
          loc.resize( (1+maxcol)*arg1->getShape()[1] );
      }
      for(unsigned i=0; i<arg1->getShape()[1]; ++i) {
          loc[i*(1+maxcol)] = 0;
      }
      unsigned argncols = arg1->getNumberOfColumns();
      for(unsigned i=0; i<arg1->getShape()[0]; ++i) {
          unsigned nr = arg1->getRowLength(i);
          for(unsigned j=0; j<nr; ++j) {
              unsigned ind = arg1->getRowIndex( i, j );
              unsigned startpos = ind*(1+maxcol);
              loc[startpos+1+loc[startpos]] = i*argncols + j;
              loc[startpos]++;
          }
      }
  }
  taskmanager.applyForces( outforces );
  if( no_thread_gather ) { 
    getColumnBookeepingArrays( taskmanager.getActionInput().outmat );
    taskmanager.getActionInput().gatherForceOnColumns = gatherForceOnColumns = true;
    taskmanager.setNForceScalars( maxcolsize );
    taskmanager.applyForces( outforces, false );
    taskmanager.setNForceScalars( getPntrToComponent(0)->getNumberOfColumns() );
    taskmanager.getActionInput().gatherForceOnColumns = gatherForceOnColumns = false;
  } 
}

template <class CV, typename myPTM>
int MatrixTimesMatrix<CV, myPTM>::getNumberOfValuesPerTask( std::size_t task_index,
    const input_type& actiondata ) {
  std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
  return actiondata.outmat[fstart];
}

template <class CV, typename myPTM>
void MatrixTimesMatrix<CV, myPTM>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const input_type& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  auto arg0=ArgumentBookeepingHolder::create( 0, input );
  auto arg1=ArgumentBookeepingHolder::create( 1, input );
  std::size_t outcol = actiondata.outmat[task_index*(1+actiondata.outmat.ncols)+1+colno];
  std::size_t arg0_start = task_index*arg0.ncols;
  std::size_t arg0_bookstart = task_index*(1+arg0.ncols);
  std::size_t arg1_bookstart = outcol*(1+arg1.ncols);
  if( actiondata.gatherForceOnColumns ) {
      arg0_start = outcol*arg0.ncols;
      arg0_bookstart = outcol*(1+arg0.ncols);
      arg1_bookstart = task_index*(1+arg1.ncols);
  }
  std::size_t nmult = arg0.bookeeping[arg0_bookstart];
  if( actiondata.secondMatrixIsSparse ) {
    std::size_t nm = 0;
    std::size_t arg1_nelements = arg1.bookeeping[arg1_bookstart];
    if( !actiondata.no_thread_gather ) {
        for(unsigned j=0; j<nmult; ++j) {
            std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
            for(unsigned k=0; k<arg1_nelements; ++k) {
                if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                    force_indices.indices[0][nm] = arg0_start + j;
                    force_indices.indices[0][nmult+nm] = arg1.start + actiondata.locations[ arg1_bookstart + 1 + k ];
                    nm++;
                    break;
                }
            }
        }
        force_indices.threadsafe_derivatives_end[0] = nm;
        if( nm==0 ) {
            force_indices.tot_indices[0] = 0;
        } else {
            force_indices.tot_indices[0] = nmult + nm;
        }
    } else if( actiondata.gatherForceOnColumns ) {
       for(unsigned j=0; j<nmult; ++j) {
            std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
            for(unsigned k=0; k<arg1_nelements; ++k) {
                if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                    force_indices.indices[0][nm] = arg1.start + actiondata.locations[ arg1_bookstart + 1 + k ];
                    nm++;
                    break;
                }
            }
        }
        force_indices.threadsafe_derivatives_end[0] = nm;
        force_indices.tot_indices[0] = nm;
    } else {
        for(unsigned j=0; j<nmult; ++j) {
            std::size_t arg0_ind = arg0.bookeeping[arg0_bookstart+1+j];
            for(unsigned k=0; k<arg1_nelements; ++k) {
                if( arg1.bookeeping[arg1_bookstart+1+k]==arg0_ind ) {
                    force_indices.indices[0][nm] = arg0_start + j;
                    nm++;
                    break;
                }
            }
        }
        force_indices.threadsafe_derivatives_end[0] = nm;
        force_indices.tot_indices[0] = nm;
    }
  } else {
    if( !actiondata.no_thread_gather ) { 
        for(unsigned j=0; j<nmult; ++j) {
          force_indices.indices[0][j] = task_index*arg0.ncols + j;
          force_indices.indices[0][nmult+j] = arg1.start + arg0.bookeeping[arg0_bookstart+1+j]*arg1.shape[1] + outcol;
        }
        force_indices.threadsafe_derivatives_end[0] = nmult;
        force_indices.tot_indices[0] = 2*nmult;
    } else if( actiondata.gatherForceOnColumns ) {
        for(unsigned j=0; j<nmult; ++j) {
          force_indices.indices[0][j] = arg1.start + arg0.bookeeping[arg0_bookstart+1+j]*arg1.shape[1] + task_index;  
        }
        force_indices.threadsafe_derivatives_end[0] = nmult;
        force_indices.tot_indices[0] = nmult;
    } else {
        for(unsigned j=0; j<nmult; ++j) {
          force_indices.indices[0][j] = task_index*arg0.ncols + j;   
        }
        force_indices.threadsafe_derivatives_end[0] = nmult;
        force_indices.tot_indices[0] = nmult;
    }
  }
}

}
}
#endif
