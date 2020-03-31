/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "VectorProductMatrix.h"
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void VectorProductMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("numbered","GROUP","the vectors of arguments for which you would like to calculate the vector product matrix");
  keys.add("numbered","GROUPA","");
  keys.add("numbered","GROUPB","");
  keys.reset_style("GROUP","optional"); keys.remove("NUMERICAL_DERIVATIVES");
}

VectorProductMatrix::VectorProductMatrix(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  ncol_args(0)
{
  bool readgroup=false; unsigned inargs=0; std::vector<std::string> g_name; std::vector<Value*> args;
  for(unsigned i=1;; ++i) {
    if( !parseNumberedVector("GROUP",i,g_name) ) { break; }
    readgroup=true; std::vector<Value*> gargs; interpretArgumentList( g_name, gargs );
    unsigned jnargs=0; arg_ends.push_back( args.size() );
    for(unsigned j=0;j<gargs.size();++j) {
        args.push_back( gargs[j] );
        if( i==1 ) inargs += gargs[j]->getNumberOfValues( getLabel() );
        else jnargs += gargs[j]->getNumberOfValues( getLabel() );
    }
    if( i!=1 && inargs!=jnargs ) error("all arguments should have same shape");
  }
  if( readgroup ) {
    arg_ends.push_back( args.size() );
    log.printf("  calculating square vector product matrix \n");
    for(unsigned i=0; i<args.size(); ++i) log.printf("  %dth component of vectors for vector product is %s\n", i+1,args[i]->getOutputDescription( getLabel() ).c_str() );
  }
  if( !readgroup ) {
    std::vector<std::string> ga_name, gb_name; unsigned inargs=0;
    log.printf("  calculating rectangular vector product matrix \n");
    for(unsigned i=1;; ++i) {
      if( !parseNumberedVector("GROUPA",i,ga_name) ) { break; }
      std::vector<Value*> gargs; interpretArgumentList( ga_name, gargs );
      unsigned jnargs=0; arg_ends.push_back( args.size() );
      for(unsigned j=0;j<gargs.size();++j) {
          args.push_back( gargs[j] );
          if( i==1 ) inargs += gargs[j]->getNumberOfValues( getLabel() );
          else jnargs += gargs[j]->getNumberOfValues( getLabel() );
      }
      if( i!=1 && inargs!=jnargs ) error("all arguments to GROUPA should have same shape");
      log.printf("  %dth component of vectors in rows of vector product matrix is %s", i, gargs[0]->getOutputDescription( getLabel() ).c_str() );
      for(unsigned j=1;j<gargs.size();++j) log.printf(" %s", gargs[j]->getOutputDescription( getLabel() ).c_str() );
      log.printf("\n");
    }
    log.printf("\n"); ncol_args = arg_ends.size(); inargs=0;
    log.printf("  calculating dot matrix between with columns : \n");
    for(unsigned i=0; i<ncol_args; ++i) {
      if( !parseNumberedVector("GROUPB",i+1,gb_name) ) error("every GROUPA must have a matching GROUPB keyword");
      std::vector<Value*> gargs; interpretArgumentList( gb_name, gargs );
      unsigned jnargs=0; arg_ends.push_back( args.size() );
      for(unsigned j=0;j<gargs.size();++j) {
          args.push_back( gargs[j] );
          if( i==0 ) inargs += gargs[j]->getNumberOfValues( getLabel() );
          else jnargs += gargs[j]->getNumberOfValues( getLabel() );
      }
      if( i!=0 && inargs!=jnargs ) error("all arguments to GROUPB should have same shape");
      log.printf("  %dth component of vectors in columns of vector product matrix is %s", i+1, gargs[0]->getOutputDescription( getLabel() ).c_str() );
      for(unsigned j=1;j<gargs.size();++j) log.printf(" %s", gargs[j]->getOutputDescription( getLabel() ).c_str() );
      log.printf("\n"); 
    }
    arg_ends.push_back( args.size() );
  }
  if( args.size()==0 ) error("no arguments were read in use GROUP or GROUPA and GROUPB");
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  std::vector<unsigned> shape(2);
  unsigned tnum=0;
  for(unsigned j=arg_ends[0];j<arg_ends[1];++j) {
      for(unsigned i=0; i<args[j]->getNumberOfValues( getLabel() ); ++i){ addTaskToList(tnum); tnum++; }
  }
  // And create the matrix to hold the dot products
  unsigned cnum=0;
  for(unsigned j=arg_ends[ncol_args];j<arg_ends[ncol_args+1];++j) cnum += args[j]->getNumberOfValues( getLabel() );
  shape[0]=tnum; shape[1]=cnum; addValue( shape ); requestArguments( args, false );

  if( ncol_args>0 ) narg_derivatives = ( tnum + cnum ) * ncol_args; 
  else narg_derivatives = tnum * ( arg_ends.size() - 1 );

  // Now do some stuff for time series
  bool timeseries=getPntrToArgument(0)->isTimeSeries();
  if( timeseries ) {
      for(unsigned i=1;i<getNumberOfArguments();++i) {
          if( !getPntrToArgument(i)->isTimeSeries() ) error("all arguments should either be time series or not time series");
      }
      getPntrToOutput(0)->makeTimeSeries();
  } else {
      for(unsigned i=1;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->isTimeSeries() ) error("all arguments should either be time series or not time series");
      }
  }
}

unsigned VectorProductMatrix::getNumberOfDerivatives() const {
  if( getNumberOfAtoms()>0 ) return 3*getNumberOfAtoms() + 9 + narg_derivatives;
  return narg_derivatives;
}

bool VectorProductMatrix::canBeAfterInChain( ActionWithValue* av ) const {
  AdjacencyMatrixBase* vp = dynamic_cast<AdjacencyMatrixBase*>( av );
  if( vp ) {
      if( vp->getFullNumberOfTasks()<(vp->getNumberOfAtoms()-vp->threeblocks.size()) && ncol_args==0 ) {
          error("cannot mix GROUPA/GROUPB actions with GROUP actions");
      } else if( ncol_args>0 && vp->getFullNumberOfTasks()==(vp->getNumberOfAtoms()-vp->threeblocks.size()) ) {
          error("cannot mix GROUPA/GROUPB actions with GROUP actions");
      }
  }
  return true;
}

void VectorProductMatrix::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void VectorProductMatrix::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void VectorProductMatrix::calculateNumericalDerivatives( ActionWithValue* a ) { plumed_error(); }

void VectorProductMatrix::calculate() {
  if( actionInChain() || skipCalculate() ) return;
  runAllTasks();
}

void VectorProductMatrix::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  if( getFullNumberOfTasks()>0 ) runAllTasks();
}

void VectorProductMatrix::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  resizeForFinalTasks(); runAllTasks();
}

void VectorProductMatrix::updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;

  unsigned nargs=arg_ends.size()-1; if( ncol_args>0 ) nargs = ncol_args; 

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  unsigned invals=0; for(unsigned i=arg_ends[0];i<arg_ends[1];++i) invals += getPntrToArgument(i)->getNumberOfValues( getLabel() );

  for(unsigned i=0; i<nargs; ++i) { matrix_indices[nmat_ind] = ind + i*invals; nmat_ind++; }
  if( getNumberOfAtoms()>0 ) {
    matrix_indices[nmat_ind+0]=narg_derivatives + 3*ind+0;
    matrix_indices[nmat_ind+1]=narg_derivatives + 3*ind+1;
    matrix_indices[nmat_ind+2]=narg_derivatives + 3*ind+2;
    nmat_ind+=3; unsigned virbase = narg_derivatives + 3*getNumberOfAtoms();
    for(unsigned i=0; i<9; ++i) matrix_indices[nmat_ind+i]=virbase+i;
    nmat_ind+=9;
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
}

void VectorProductMatrix::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( actionInChain() ) {
    updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
    return ;
  }

  // Now loop over all atoms in coordination sphere
  unsigned start_n=0, size_v = 0; 
  if(  ncol_args>0 ) { start_n = getFullNumberOfTasks(); myvals.setNumberOfIndicesInFirstBlock( start_n ); }
  for(unsigned i=arg_ends[ncol_args];i<arg_ends[ncol_args+1];++i) size_v += getPntrToArgument(i)->getNumberOfValues( getLabel() ); 

  for(unsigned i=0; i<size_v; ++i) {
    // Don't do i==j
    if( ncol_args==0 && myvals.getTaskIndex()==i ) continue;
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), myvals.getTaskIndex(), current, start_n + i, myvals );
    // Now clear only elements that are not accumulated over whole row
    clearMatrixElements( myvals );
  }
  // Update the matrix index for the central atom
  updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
}

bool VectorProductMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned invals, jnvals; 
  invals=0; for(unsigned i=arg_ends[0];i<arg_ends[1];++i) invals += getPntrToArgument(i)->getNumberOfValues( getLabel() ); 

  unsigned jindex=index2, jind_start = 0, nargs=arg_ends.size()-1; if( ncol_args>0 ) nargs = ncol_args;
  if( ncol_args>0 ) {
    nargs = ncol_args; jindex = index2 - getFullNumberOfTasks(); jind_start = nargs*invals;  
    jnvals = 0; for(unsigned j=arg_ends[ncol_args];j<arg_ends[ncol_args+1];++j) jnvals += getPntrToArgument(j)->getNumberOfValues( getLabel() );
  } else jnvals = invals;

  std::vector<double> args1(nargs), args2(nargs), der1(nargs), der2(nargs);
  for(unsigned i=0; i<nargs; ++i) {
    args1[i] = retrieveRequiredArgument( i, index1 ); args2[i] = retrieveRequiredArgument( ncol_args + i, jindex );
  }

  double val = computeVectorProduct( index1, index2, args1, args2, der1, der2, myvals );

  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  myvals.setValue( ostrn, val );
  // Return after calculation of value if we do not need derivatives
  if( doNotCalculateDerivatives() ) return true;

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  plumed_dbg_assert( matrix_indices.size()>=getNumberOfDerivatives() );
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  for(unsigned i=0; i<nargs; ++i) {
    plumed_dbg_assert( index1 + i*invals<getNumberOfDerivatives() );
    myvals.addDerivative( ostrn, index1 + i*invals, der1[i] );
    myvals.updateIndex( ostrn, index1 + i*invals );
    myvals.addDerivative( ostrn, jind_start + jindex + i*jnvals, der2[i] );
    myvals.updateIndex( ostrn, jind_start + jindex + i*jnvals );
    matrix_indices[nmat_ind] = jind_start + jindex + i*jnvals;
    nmat_ind++;
  }
  if( getNumberOfAtoms()>0 ) {
    matrix_indices[nmat_ind+0]=narg_derivatives + 3*index2+0;
    matrix_indices[nmat_ind+1]=narg_derivatives + 3*index2+1;
    matrix_indices[nmat_ind+2]=narg_derivatives + 3*index2+2;
    nmat_ind+=3;
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
  return true;
}

void VectorProductMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) {
    setForcesOnAtoms( forcesToApply, mm );
    setForcesOnArguments( 0, forcesToApply, mm );
  }
}

}
}
