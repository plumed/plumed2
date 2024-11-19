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
#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {

void ActionWithVector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES");
  ActionWithArguments::registerKeywords( keys );
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.reserve("optional","MASK","the label for a sparse matrix that should be used to determine which elements of the matrix should be computed");
}

ActionWithVector::ActionWithVector(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  nmask(-1),
  serial(false),
  forwardPass(false)
{
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
    if( av && av->getNumberOfMasks()>=0 ) nmask=0;
  }

  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  if( keywords.exists("MASK") ) {
    std::vector<Value*> mask; parseArgumentList("MASK",mask);
    if( mask.size()>0 ) {
      if( nmask>=0 && getNumberOfArguments()==1 ) {
        ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(0)->getPntrToAction() );
        plumed_massert( av, "input should be a vector from ActionWithVector" ); unsigned j=0, nargs = av->getNumberOfArguments();
        for(unsigned i=nargs-av->nmask; i<nargs; ++i) {
          if( av->getPntrToArgument(i)!=mask[j] ) error("the masks in subsequent actions do not match");
          j++;
        }
      }
      if( getNumberOfArguments()>0 && getPntrToArgument(0)->hasDerivatives() ) error("input for mask should be vector or matrix");
      else if( mask[0]->getRank()==2 ) {
        if( mask.size()>1 ) error("MASK should only have one argument");
        log.printf("  only computing elements of matrix that correspond to non-zero elements of matrix %s \n", mask[0]->getName().c_str() );
      } else if( mask[0]->getRank()==1 ) {
        log.printf("  only computing elements of vector that correspond to non-zero elements of vectors %s", mask[0]->getName().c_str() );
        for(unsigned i=1; i<mask.size(); ++i) {
          if( mask[i]->getRank()!=1 ) { log.printf("\n"); error("input to mask should be vector"); }
          log.printf(", %s", mask[i]->getName().c_str() );
        }
        log.printf("\n");
      }
      std::vector<Value*> allargs( getArguments() ); nmask=mask.size();
      for(unsigned i=0; i<mask.size(); ++i) allargs.push_back( mask[i] );
      requestArguments( allargs );
    }
  }
}

void ActionWithVector::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void ActionWithVector::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithVector::calculateNumericalDerivatives(ActionWithValue* av) {
  plumed_merror("cannot calculate numerical derivative for action " + getName() + " with label " + getLabel() );
}

void ActionWithVector::prepare() {
  active_tasks.resize(0);
}

int ActionWithVector::checkTaskIsActive( const unsigned& itask ) const {
  unsigned nargs = getNumberOfArguments();
  if( nargs==0 ) {
    return 1;
  } else if( nmask>0 ) {
    for(unsigned j=nargs-nmask; j<nargs; ++j) {
      Value* myarg = getPntrToArgument(j);
      if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0 ) return 1;
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        unsigned ncol = myarg->getRowLength(itask);
        unsigned base = itask*myarg->getNumberOfColumns();
        for(unsigned k=0; k<ncol; ++k) {
          if( fabs(myarg->get(base+k,false))>0 ) return 1;
        }
      } else plumed_merror("only matrices and vectors should be used as masks");
    }
  } else {
    for(unsigned i=0; i<nargs; ++i) {
      if( getName()=="OUTER_PRODUCT" && i>0 ) return -1;

      Value* myarg = getPntrToArgument(i);
      if( !myarg->isDerivativeZeroWhenValueIsZero() ) return 1;

      if( myarg->getRank()==0 ) {
        return 1;
      } else if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0 ) return 1;
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        unsigned ncol = myarg->getRowLength(itask);
        unsigned base = itask*myarg->getNumberOfColumns();
        for(unsigned k=0; k<ncol; ++k) {
          if( fabs(myarg->get(base+k,false))>0 ) return 1;
        }
      } else if( myarg->getRank()>0 ) {
        return 1;
      } else plumed_merror("should not be in action " + getName() );
    }
  }
  return -1;
}

std::vector<unsigned>& ActionWithVector::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) return active_tasks;
  unsigned ntasks=0; getNumberOfTasks( ntasks );

  std::vector<int> taskFlags( ntasks, -1 );
  for(unsigned i=0; i<ntasks; ++i) taskFlags[i] = checkTaskIsActive(i);
  unsigned nt=0;
  for(unsigned i=0; i<ntasks; ++i) {
    if( taskFlags[i]>0 ) nt++;
  }
  active_tasks.resize(nt); nt=0;
  for(unsigned i=0; i<ntasks; ++i) {
    if( taskFlags[i]>0 ) { active_tasks[nt]=i; nt++; }
  }
  return active_tasks;
}

bool ActionWithVector::doNotCalculateDerivatives() const {
  if( forwardPass ) return true;
  return ActionWithValue::doNotCalculateDerivatives();
}

void ActionWithVector::clearMatrixBookeeping() {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->getRank()==2 && myval->getNumberOfColumns()<myval->getShape()[1] ) {
      std::fill(myval->matrix_bookeeping.begin(), myval->matrix_bookeeping.end(), 0);
    }
    myval->set(0);
  }
}

void ActionWithVector::runAllTasks() {
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Clear matrix bookeeping arrays
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(this);
  if( am && stride>1 ) clearMatrixBookeeping();

  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( getListOfActiveTasks( this ) );
  unsigned nactive_tasks=partialTaskList.size();

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 ) nt=1;
  if( myvals.size()!=nt ) myvals.resize(nt);

  // Get the total number of streamed quantities that we need
  // Get size for buffer
  unsigned bufsize=0, nderivatives = 0; bool gridsInStream=false; forwardPass=true;
  for(int i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->bufstart=bufsize;
    if( getPntrToComponent(i)->hasDerivatives() || getPntrToComponent(i)->getRank()==0 ) {
      forwardPass=false; bufsize += getPntrToComponent(i)->data.size();
    }
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->hasDerivatives() ) {
      nderivatives=getConstPntrToComponent(i)->getNumberOfGridDerivatives(); gridsInStream=true;
    }
  }
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Recover the number of derivatives we require
  if( !doNotCalculateDerivatives() && !gridsInStream ) {
    unsigned nargs = getNumberOfArguments(); int nmasks=getNumberOfMasks();
    if( nargs>=nmasks && nmasks>0 ) nargs = nargs - nmasks;
    if( getNumberOfAtoms()>0 ) nderivatives += 3*getNumberOfAtoms() + 9;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) nderivatives += getPntrToArgument(i)->getNumberOfValues();
  }

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    const unsigned t=OpenMP::getThreadNum();
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    if( myvals[t].getNumberOfValues()!=getNumberOfComponents() || myvals[t].getNumberOfDerivatives()!=nderivatives ) myvals[t].resize( getNumberOfComponents(), nderivatives );
    myvals[t].clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( partialTaskList[i], myvals[t] );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) gatherAccumulators( partialTaskList[i], myvals[t], omp_buffer );
      else gatherAccumulators( partialTaskList[i], myvals[t], buffer );

      // Clear the value
      myvals[t].clearAll();
    }
    #pragma omp critical
    if( nt>1 ) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }

  // MPI Gather everything
  if( !serial ) {
    if( buffer.size()>0 ) comm.Sum( buffer );
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      Value* myval = getPntrToComponent(i);
      if( !myval->hasDeriv ) {
        comm.Sum( myval->data ); if( am && myval->getRank()==2 && myval->getNumberOfColumns()<myval->getShape()[1] ) comm.Sum( myval->matrix_bookeeping );
      }
    }
  }
  finishComputations( buffer ); forwardPass=false;
}

void ActionWithVector::getNumberOfTasks( unsigned& ntasks ) {
  if( ntasks==0 ) {
    if( getNumberOfArguments()==1 && getNumberOfComponents()==1 && getPntrToComponent(0)->getRank()==0 ) {
      if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()==2 ) ntasks = getPntrToArgument(0)->getShape()[0];
      else ntasks = getPntrToArgument(0)->getNumberOfValues();
    } else {
      plumed_assert( getNumberOfComponents()>0 && getPntrToComponent(0)->getRank()>0 );
      if( getPntrToComponent(0)->hasDerivatives() ) ntasks = getPntrToComponent(0)->getNumberOfValues();
      else ntasks = getPntrToComponent(0)->getShape()[0];
    }
  }
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==0 ) {
      if( getNumberOfArguments()!=1 ) error("mismatched numbers of tasks in streamed quantities");
      if( getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
      else if ( !getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
    } else if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
    else if( !getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
  }
}

void ActionWithVector::runTask( const unsigned& current, MultiValue& myvals ) const {
  const ActionWithMatrix* am = dynamic_cast<const ActionWithMatrix*>(this);
  myvals.setTaskIndex(current);
  myvals.vector_call=true;
  performTask( current, myvals );
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    const Value* myval = getConstPntrToComponent(i);
    if( am || myval->hasDerivatives() )
      continue;
    Value* myv = const_cast<Value*>( myval );
    if( getName()=="RMSD_VECTOR" && myv->getRank()==2 )
      continue;
    myv->set( current, myvals.get( i ) );
  }
}

void ActionWithVector::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    unsigned bufstart = getConstPntrToComponent(i)->bufstart;
    // This looks after storing of scalars that are summed from vectors/matrices
    if( getConstPntrToComponent(i)->getRank()==0 ) {
      plumed_dbg_massert( bufstart<buffer.size(), "problem in " + getLabel() );
      buffer[bufstart] += myvals.get(i);
      if( getConstPntrToComponent(i)->hasDerivatives() ) {
        for(unsigned k=0; k<myvals.getNumberActive(i); ++k) {
          unsigned kindex = myvals.getActiveIndex(i,k);
          plumed_dbg_massert( bufstart+1+kindex<buffer.size(), "problem in " + getLabel()  );
          buffer[bufstart + 1 + kindex] += myvals.getDerivative(i,kindex);
        }
      }
      // This looks after storing of vectors
    } else gatherStoredValue( i, taskCode, myvals, bufstart, buffer );
  }
}

void ActionWithVector::finishComputations( const std::vector<double>& buf ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    // This gathers vectors and grids at the end of the calculation
    unsigned bufstart = getPntrToComponent(i)->bufstart;
    if( (getPntrToComponent(i)->getRank()>0 && getPntrToComponent(i)->hasDerivatives()) ) {
      unsigned sz_v = getPntrToComponent(i)->data.size(); getPntrToComponent(i)->data.assign( getPntrToComponent(i)->data.size(), 0 );
      for(unsigned j=0; j<sz_v; ++j) {
        plumed_dbg_assert( bufstart+j<buf.size() );
        getPntrToComponent(i)->add( j, buf[bufstart+j] );
      }
      // Make sure single values are set
    } else if( getPntrToComponent(i)->getRank()==0 ) getPntrToComponent(i)->set( buf[bufstart] );
    // This gathers derivatives of scalars
    if( !doNotCalculateDerivatives() && getPntrToComponent(i)->hasDeriv && getPntrToComponent(i)->getRank()==0 ) {
      for(unsigned j=0; j<getPntrToComponent(i)->getNumberOfDerivatives(); ++j) getPntrToComponent(i)->setDerivative( j, buf[bufstart+1+j] );
    }
  }
}

bool ActionWithVector::checkChainForNonScalarForces() const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->forcesWereAdded() ) return true;
  }
  return false;
}

void ActionWithVector::getNumberOfForceDerivatives( unsigned& nforces, unsigned& nderiv ) const {
  nforces=0; unsigned nargs = getNumberOfArguments(); int nmasks = getNumberOfMasks();
  if( nargs>=nmasks && nmasks>0 ) nargs = nargs - nmasks;
  if( getNumberOfAtoms()>0 ) nforces += 3*getNumberOfAtoms() + 9;
  for(unsigned i=0; i<nargs; ++i) {
    nforces += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  nderiv = nforces;
}

bool ActionWithVector::checkForForces() {
  if( getPntrToComponent(0)->getRank()==0 ) return ActionWithValue::checkForForces();

  // Check if there are any forces
  if( !checkChainForNonScalarForces() ) return false;

  // Setup MPI parallel loop
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Get the number of tasks
  std::vector<unsigned> force_tasks; std::vector<unsigned> & partialTaskList( getListOfActiveTasks( this ) );
  for(unsigned i=0; i<partialTaskList.size(); ++i) force_tasks.push_back( partialTaskList[i] );
  unsigned nf_tasks=force_tasks.size();
  if( nf_tasks==0 ) return false;

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nf_tasks ) nt=nf_tasks/stride/10;
  if( nt==0 ) nt=1;
  if( myvals.size()!=nt ) myvals.resize(nt);
  if( omp_forces.size()!=nt ) omp_forces.resize(nt);

  // Recover the number of derivatives we require (this should be equal to the number of forces)
  unsigned nderiv, nforces; getNumberOfForceDerivatives( nforces, nderiv );
  if( forcesForApply.size()!=nforces ) forcesForApply.resize( nforces );
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );

  #pragma omp parallel num_threads(nt)
  {
    const unsigned t=OpenMP::getThreadNum();
    if( nt>1 ) {
      if( omp_forces[t].size()!=forcesForApply.size() ) omp_forces[t].resize( forcesForApply.size(), 0.0 );
      else omp_forces[t].assign( forcesForApply.size(), 0.0 );
    }
    if( myvals[t].getNumberOfValues()!=getNumberOfComponents() || myvals[t].getNumberOfDerivatives()!=nderiv ) myvals[t].resize( getNumberOfComponents(), nderiv );
    myvals[t].clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nf_tasks; i+=stride) {
      runTask( force_tasks[i], myvals[t] );

      // Now get the forces
      if( nt>1 ) gatherForces( force_tasks[i], myvals[t], omp_forces[t] );
      else gatherForces( force_tasks[i], myvals[t], forcesForApply );

      myvals[t].clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<forcesForApply.size(); ++i) forcesForApply[i]+=omp_forces[t][i];
  }
  // MPI Gather on forces
  if( !serial ) comm.Sum( forcesForApply );
  return true;
}

bool ActionWithVector::checkComponentsForForce() const {
  for(unsigned i=0; i<values.size(); ++i) {
    if( getConstPntrToComponent(i)->forcesWereAdded() ) return true;
  }
  return false;
}

bool ActionWithVector::checkForTaskForce( const unsigned& itask, const Value* myval ) const {
  plumed_dbg_assert( !(myval->getRank()==2 && !myval->hasDerivatives()) );
  return fabs(myval->getForce(itask))>epsilon;
}

void ActionWithVector::gatherForcesOnStoredValue( const unsigned& ival, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  const Value* myval = getConstPntrToComponent(ival); double fforce = myval->getForce(itask);
  for(unsigned j=0; j<myvals.getNumberActive(ival); ++j) {
    unsigned jder=myvals.getActiveIndex(ival, j); plumed_dbg_assert( jder<forces.size() );
    forces[jder] += fforce*myvals.getDerivative( ival, jder );
  }
}

void ActionWithVector::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( checkComponentsForForce() ) {
    for(unsigned k=0; k<getNumberOfComponents(); ++k) {
      const Value* myval=getConstPntrToComponent(k);
      if( myval->getRank()>0 && myval->forcesWereAdded() ) gatherForcesOnStoredValue( k, itask, myvals, forces );
    }
  }
}

void ActionWithVector::apply() {
  if( !checkForForces() ) return;
  // Find the top of the chain and add forces
  unsigned ind=0; addForcesOnArguments( 0, forcesForApply, ind ); setForcesOnAtoms( forcesForApply, ind );
}

}
