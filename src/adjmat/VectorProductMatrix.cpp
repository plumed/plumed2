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
  ncol_args(0),
  nspAa(0),
  nspAb(0)
{
  bool readgroup=false; unsigned inargs=0; std::vector<std::string> g_name; std::vector<Value*> args;
  for(unsigned i=1;; ++i) {
    if( !parseNumberedVector("GROUP",i,g_name) ) { break; }
    if( i==1 ) nspAa=nspAb=g_name.size();
    else if( g_name.size()!=nspAa ) error("input to all GROUP keywords should have the same number of arguments");
    readgroup=true; unsigned jnargs=0;
    for(unsigned j=0;j<g_name.size();++j) { 
        args.push_back( convertStringToValue(g_name[j]) );
        if( args[args.size()-1]->getRank()>1 ) error("all arguments should be vectors");
        if( args[0]->getRank()==0 ) {
          if( args[args.size()-1]->getRank()!=0 ) error("all arguments should have same shape");
        } else if( i==1 ) inargs += args[args.size()-1]->getShape()[0];
        else jnargs += args[args.size()-1]->getShape()[0]; 
    }
    if( i!=1 && inargs!=jnargs ) error("all arguments should have same shape");
  }
  if( readgroup ) {
    log.printf("  calculating square vector product matrix \n");
    for(unsigned i=0; i<args.size(); ++i) log.printf("  %dth component of vectors for vector product is %s\n", i+1,args[i]->getName().c_str() );
  }
  if( !readgroup ) {
    std::vector<std::string> ga_name, gb_name; unsigned inargs=0;
    log.printf("  calculating rectangular vector product matrix \n");
    for(unsigned i=1;; ++i) {
      if( !parseNumberedVector("GROUPA",i,ga_name) ) { break; }
      if( i==1 ) nspAa=ga_name.size(); 
      else if( ga_name.size()!=nspAa ) error("input to all GROUPA keywords should have the same number of arguments");
      unsigned jnargs=0;
      for(unsigned j=0;j<ga_name.size();++j) {
          args.push_back( convertStringToValue(ga_name[j]) );
          if( args[args.size()-1]->getRank()>1 ) error("all arguments should be vectors");
          if( args[0]->getRank()==0 ) {
            if( args[args.size()-1]->getRank()!=0 ) error("all arguments to GROUPA should have same shape");
          } else if( i==1 ) inargs += args[args.size()-1]->getShape()[0];
          else jnargs += args[args.size()-1]->getShape()[0];
      }
      if( i!=1 && inargs!=jnargs ) error("all arguments to GROUPA should have same shape");
      log.printf("  %dth component of vectors in rows of vector product matrix is %s", i, ga_name[0].c_str() );
      for(unsigned j=1;j<ga_name.size();++j) log.printf(" %s", ga_name[j].c_str() );
      log.printf("\n");
    }
    log.printf("\n"); ncol_args = args.size(); inargs=0;
    log.printf("  calculating dot matrix between with columns : \n");
    for(unsigned i=0; i<ncol_args; ++i) {
      if( !parseNumberedVector("GROUPB",i+1,gb_name) ) error("every GROUPA must have a matching GROUPB keyword");
      if( i==0 ) nspAb=gb_name.size(); 
      else if( gb_name.size()!=nspAb ) error("input to all GROUPB keywords should have the same number of arguments");
      unsigned jnargs=0;
      for(unsigned j=0;j<gb_name.size();++j) {
          args.push_back( convertStringToValue(gb_name[j]) );
          if( args[args.size()-1]->getRank()>1 ) error("all arguments should be vectors");
          if( args[ncol_args]->getRank()==0 ) {
            if( args[args.size()-1]->getRank()!=0 ) error("all arguments to GROUPB should have same shape");
          } else if( i==0 ) inargs += args[args.size()-1]->getShape()[0];
          else jnargs += args[args.size()-1]->getShape()[0];
      }
      if( i!=0 && inargs!=jnargs ) error("all arguments to GROUPB should have same shape");
      log.printf("  %dth component of vectors in columns of vector product matrix is %s", i+1, gb_name[0].c_str() );
      for(unsigned j=1;j<gb_name.size();++j) log.printf(" %s", gb_name[j].c_str() );
      log.printf("\n"); 
    }
  }
  if( args.size()==0 ) error("no arguments were read in use GROUP or GROUPA and GROUPB");
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  std::vector<unsigned> shape(2);
  if( args[0]->getRank()==0 ) {
    addTaskToList(0); shape[0]=shape[1]=1;
  } else {
    unsigned tnum=0;
    for(unsigned j=0;j<nspAa;++j) {
        for(unsigned i=0; i<args[j]->getShape()[0]; ++i){ addTaskToList(tnum); tnum++; }
    }
    unsigned cnum=0;
    for(unsigned j=nspAa*ncol_args;j<nspAa*ncol_args+nspAb;++j) cnum += args[j]->getShape()[0];
    shape[0]=tnum; shape[1]=cnum;
  }
  requestArguments( args, false );
  // And create the matrix to hold the dot products
  addValue( shape ); // if( readgroup ) getPntrToComponent(0)->setSymmetric( true );    /// Is this right though???
  if( args[0]->getRank()==0 ) narg_derivatives = getNumberOfArguments();
  else if( ncol_args>0 ) {
    unsigned tnum=0; for(unsigned j=0;j<nspAa;++j) tnum += args[j]->getShape()[0];
    unsigned cnum=0; for(unsigned j=nspAa*ncol_args;j<nspAa*ncol_args+nspAb;++j) cnum += args[j]->getShape()[0];
    narg_derivatives = ( tnum + cnum ) * getNumberOfArguments() / ( nspAa + nspAb ); 
  } else { 
    unsigned tnum=0; for(unsigned j=0;j<nspAa;++j) tnum += args[j]->getShape()[0];
    narg_derivatives = tnum * getNumberOfArguments() / nspAa;
  }
}

Value* VectorProductMatrix::convertStringToValue( const std::string& name ) {
  std::size_t dot=name.find_first_of("."); unsigned nargs=0; std::vector<Value*> args;
  if( dot!=std::string::npos ) {
    ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name.substr(0,dot) );
    if( !action ) {
      std::string str=" (hint! the actions in this ActionSet are: ";
      str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
      error("cannot find action named " + name + str);
    }
    action->interpretDataLabel( name, this, nargs, args );
  } else {
    ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( name );
    if( !action ) {
      std::string str=" (hint! the actions in this ActionSet are: ";
      str+=plumed.getActionSet().getLabelList<ActionWithValue*>()+")";
      error("cannot find action named " + name + str);
    }
    action->interpretDataLabel( name, this, nargs, args );
  }
  plumed_assert( args.size()==1 );
  return args[0];
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
  if( actionInChain() ) return;
  runAllTasks();
}

void VectorProductMatrix::updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;

  unsigned nargs=getNumberOfArguments() / nspAa; 
  if( ncol_args>0 ) nargs = getNumberOfArguments() / (nspAa+nspAb);
  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  unsigned invals;
  if( getPntrToArgument(0)->getRank()==0 ) invals = 1;
  else {
    invals = 0; for(unsigned i=0;i<nspAa;++i) invals += getPntrToArgument(i)->getShape()[0];
  }
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
  unsigned start_n=0; if(  ncol_args>0 ) { start_n = getFullNumberOfTasks(); myvals.setNumberOfIndicesInFirstBlock( start_n ); }
  unsigned size_v = 1; if( getPntrToArgument(ncol_args)->getRank()>0 ) size_v = getPntrToArgument(ncol_args)->getShape()[0];

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
  unsigned invals, jnvals; invals = jnvals = 1;
  if( getPntrToArgument(0)->getRank()>0 ) {
      invals=0; for(unsigned i=0;i<nspAa;++i) invals +=  getPntrToArgument(i)->getShape()[0];
      jnvals=invals; 
  }

  unsigned jindex=index2, jind_start = 0, nargs=getNumberOfArguments() / nspAa;
  if( ncol_args>0 ) {
    nargs = getNumberOfArguments() / (nspAa+nspAb); jnvals = 1; 
    if( getPntrToArgument(nspAa*ncol_args)->getRank()>0 ) {
        jnvals = 0; for(unsigned j=nspAa*ncol_args;j<nspAa*ncol_args+nspAb;++j) jnvals += getPntrToArgument(j)->getShape()[0];
    }
    jindex = index2 - getFullNumberOfTasks();
    jind_start = nargs*invals;
  }

  std::vector<double> args1(nargs), args2(nargs), der1(nargs), der2(nargs);
  for(unsigned i=0; i<nargs; ++i) {
    unsigned jindex1=index1, kindex=nspAa*i; 
    for(unsigned j=0;j<nspAa;++j) {
        if( getPntrToArgument(kindex)->getRank()==0 ) {
            if( jindex1<1 ){ break; }
            jindex1 -= 1; kindex++;
        } else {
            if( jindex1<getPntrToArgument(kindex)->getShape()[0] ){ break; }
            jindex1 -= getPntrToArgument(kindex)->getShape()[0]; kindex++;
        }
    }
    args1[i] = getPntrToArgument(kindex)->get( jindex1 );
    unsigned jindex2=jindex, kindex2=nspAa*ncol_args + nspAb*i;
    for(unsigned j=0;j<nspAb;++j) {
        if( getPntrToArgument(kindex2)->getRank()==0 ) {
            if( jindex2<1 ){ break; }
            jindex2 -= 1; kindex2++;
        } else {
            if( jindex2<getPntrToArgument(kindex2)->getShape()[0] ){ break; }
            jindex2 -= getPntrToArgument(kindex2)->getShape()[0]; kindex2++;
        }
    }
    args2[i] = getPntrToArgument(kindex2)->get( jindex2 );
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
