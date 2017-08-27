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
  bool readgroup=false; std::string g_name; std::vector<Value*> args;
  for(unsigned i=1;;++i) {
      if( !parseNumbered("GROUP",i,g_name) ){ break; }
      readgroup=true; args.push_back( convertStringToValue(g_name) );
      if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
      if( args[args.size()-1]->getShape()[0]!=args[0]->getShape()[0] ) error("all arguments should have same shape");
  }
  if( readgroup ) {
      log.printf("  calculating square vector product matrix \n");
      for(unsigned i=0;i<args.size();++i) log.printf("  %dth component of vectors for vector product is %s\n", i+1,args[i]->getName().c_str() );
  }
  if( !readgroup ){
      std::string ga_name, gb_name;
      log.printf("  calculating rectangular vector product matrix \n"); 
      for(unsigned i=1;;++i) {
          if( !parseNumbered("GROUPA",i,ga_name) ){ break; }
          args.push_back( convertStringToValue(ga_name) );
          if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
          if( args[args.size()-1]->getShape()[0]!=args[0]->getShape()[0] ) error("all arguments to GROUPA should have same shape");
          log.printf("  %dth component of vectors in rows of vector product matrix is %s \n", i, ga_name.c_str() );
      }
      log.printf("\n"); ncol_args = args.size();
      log.printf("  calculating dot matrix between with columns : \n"); 
      for(unsigned i=0;i<ncol_args;++i){ 
          if( !parseNumbered("GROUPB",i+1,gb_name) ) error("every GROUPA must have a matching GROUPB keyword");
          args.push_back( convertStringToValue(gb_name) );
          if( args[args.size()-1]->getRank()!=1 ) error("all arguments should be vectors");
          if( args[args.size()-1]->getShape()[0]!=args[ncol_args]->getShape()[0] ) error("all arguments to GROUPB should have same shape");
          log.printf("  %dth component of vectors in columns of vector product matrix is %s\n", i+1, gb_name.c_str() );
      }
  }
  if( args.size()==0 ) error("no arguments were read in use GROUP or GROUPA and GROUPB");
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  for(unsigned i=0;i<args[0]->getShape()[0];++i) addTaskToList(i);
  requestArguments( args, false ); std::vector<unsigned> shape(2); 
  shape[0]=args[0]->getShape()[0]; shape[1]=args[ncol_args]->getShape()[0];
  // And create the matrix to hold the dot products 
  addValue( shape ); 
}

Value* VectorProductMatrix::convertStringToValue( const std::string& name ) {
  std::size_t dot=name.find_first_of("."); std::vector<Value*> args;
  if( dot!=std::string::npos ) {
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name.substr(0,dot) );
      if( !action ){
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + name + str);
      }
      action->interpretDataLabel( name, this, args );
  } else {
      ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name );
      if( !action ){
          std::string str=" (hint! the actions in this ActionSet are: ";
          str+=plumed.getActionSet().getLabelList()+")";
          error("cannot find action named " + name + str);
      }
      action->interpretDataLabel( name, this, args );
  }
  plumed_assert( args.size()==1 );
  return args[0];
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

void VectorProductMatrix::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  tflags.assign(tflags.size(),1);
}

void VectorProductMatrix::calculate() {
  if( actionInChain() ) return;
  runAllTasks();
}

void VectorProductMatrix::updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;

  unsigned nat_indices = 0; 
  unsigned nargs=getNumberOfArguments(); if( ncol_args>0 ) nargs /= 2;
  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  if( getNumberOfAtoms()>0 ) {
      matrix_indices[nmat_ind+0]=3*ind+0; matrix_indices[nmat_ind+1]=3*ind+1; matrix_indices[nmat_ind+2]=3*ind+2; nmat_ind+=3; 
      unsigned virbase = 3*getNumberOfAtoms(); nat_indices = virbase+9;
      for(unsigned i=0;i<9;++i) matrix_indices[nmat_ind+i]=virbase+i;
      nmat_ind+=9; 
  }   
  unsigned invals = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0;i<nargs;++i) matrix_indices[nmat_ind+i] = nat_indices + ind + i*invals; 
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind + nargs );
}

void VectorProductMatrix::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( actionInChain() ){
     updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
     return ;
  }

  // Now loop over all atoms in coordination sphere
  unsigned start_n=0; if(  ncol_args>0 ){ start_n = getFullNumberOfTasks(); myvals.setNumberOfIndicesInFirstBlock( start_n ); }
  for(unsigned i=0;i<getPntrToArgument(ncol_args)->getShape()[0];++i){
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

void VectorProductMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned invals, jnvals; invals = jnvals = getPntrToArgument(0)->getShape()[0]; 
  unsigned jindex=index2, jind_start = 0, nargs=getNumberOfArguments(); 
  if( ncol_args>0 ) { 
      nargs /= 2; jnvals = getPntrToArgument(ncol_args)->getShape()[0];
      jindex = index2 - getFullNumberOfTasks();
      jind_start = nargs*invals;
  }
 
  std::vector<double> args1(nargs), args2(nargs), der1(nargs), der2(nargs);
  for(unsigned i=0;i<nargs;++i){ 
      args1[i] = getPntrToArgument(i)->get( index1 ); 
      args2[i] = getPntrToArgument(ncol_args+i)->get( jindex ); 
  }
  
  double val = computeVectorProduct( index1, index2, args1, args2, der1, der2, myvals );

  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  myvals.setValue( ostrn, val ); 
  // Return after calculation of value if we do not need derivatives
  if( doNotCalculateDerivatives() ) return;

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  if( matrix_indices.size()<getNumberOfDerivatives() ) matrix_indices.resize( getNumberOfDerivatives() );
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat ), nat_indices=0;
  if( getNumberOfAtoms()>0 ) {
      unsigned nat_indices=3*getNumberOfAtoms()+9;
      myvals.updateIndex( ostrn, 3*index1+0 ); myvals.updateIndex( ostrn, 3*index1+1 ); myvals.updateIndex( ostrn, 3*index1+2 );
      myvals.updateIndex( ostrn, 3*index2+0 ); myvals.updateIndex( ostrn, 3*index2+1 ); myvals.updateIndex( ostrn, 3*index2+2 );
      unsigned virbase=3*getNumberOfAtoms(); for(unsigned i=0;i<9;++i) myvals.updateIndex( ostrn, virbase + i );
      matrix_indices[nmat_ind+0]=3*index2+0; matrix_indices[nmat_ind+1]=3*index2+1; matrix_indices[nmat_ind+2]=3*index2+2; nmat_ind+=3;
  }
  for(unsigned i=0;i<nargs;++i) {
      plumed_dbg_assert( index1 + i*invals<getNumberOfDerivatives() );
      myvals.addDerivative( ostrn, nat_indices + index1 + i*invals, der1[i] ); 
      myvals.updateIndex( ostrn, nat_indices + index1 + i*invals );
      myvals.addDerivative( ostrn, nat_indices + jind_start + jindex + i*jnvals, der2[i] ); 
      myvals.updateIndex( ostrn, nat_indices + jind_start + jindex + i*jnvals );
      matrix_indices[nmat_ind+i] = nat_indices + jind_start + jindex + i*jnvals; 
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind + nargs );
}

void VectorProductMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ){ 
      setForcesOnAtoms( forcesToApply, mm ); 
      setForcesOnArguments( forcesToApply, mm ); 
  }
}

}
}
