/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Function.h"
#include "core/Average.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

using namespace std;
namespace PLMD {
namespace function {

void Function::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
}

Function::Function(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true),
  nderivatives(getNumberOfScalarArguments()),
  forcesToApply(getNumberOfScalarArguments())
{
  plumed_dbg_assert( getNumberOfArguments()>0 );
  bool gridinput=false; unsigned npoints=0;
  // Method for if input to function is a function on a grid
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ){ 
          gridinput=true; npoints=getPntrToArgument(i)->getNumberOfValues( getLabel() ); 
          nderivatives = getPntrToArgument(i)->getRank(); break; 
      }
  }
  if( gridinput ) {
       unsigned nscalars=0; 
       if( arg_ends.size()==0 && getNumberOfArguments()==1 ){ arg_ends.push_back(0); arg_ends.push_back(1); }
       for(unsigned j=0;j<getNumberOfArguments();++j) {
           if( getPntrToArgument(j)->getRank()!=0 ) {
               if( getPntrToArgument(j)->getNumberOfValues( getLabel() )!=npoints || !getPntrToArgument(j)->hasDerivatives() ) error("mismatch in input arguments");
           } else { nscalars++; nderivatives++; }
       }
       if( nscalars>1 ) error("can only multiply/divide grid by one scalar at a time");
       // Now create a task list for the function
       for(unsigned j=0;j<npoints;++j) addTaskToList(j);
  } else {
       bool hasscalar=false, hasrank=false; nderivatives = getNumberOfScalarArguments(); firststep=false;
       for(unsigned i=0;i<getNumberOfArguments();++i) {
           if( getPntrToArgument(i)->getRank()==0 ) hasscalar=true;
           else hasrank=true;
       }
       if( hasscalar && hasrank ) {
           unsigned nscalars=0, nranks=0;
           for(unsigned i=0;i<getNumberOfArguments();++i) {
               if( getPntrToArgument(i)->getRank()==0 ) nscalars++;
               else {
                  nranks++; npoints=getPntrToArgument(i)->getNumberOfValues( getLabel() );
               }
           }
           if( nscalars>1 ) error("can only multiply/divide grid by one scalar at a time");
           // Now create a task list for the function
           for(unsigned j=0;j<npoints;++j) addTaskToList(j);
       } else {
           createTasksFromArguments(); 
           // Now create the stream of jobs to work through 
           if( distinct_arguments.size()>0 ){   // This is for if we have a function that needs to store - needs though GAT
               std::vector<std::string> alabels;
               for(unsigned i=0;i<getNumberOfArguments();++i){
                   bool found=false; std::string mylab = (getPntrToArgument(i)->getPntrToAction())->getLabel();
                   for(unsigned j=0;j<alabels.size();++j){
                       if( alabels[j]==mylab ){ found=true; break; }
                   }
                   if( !found ) alabels.push_back( mylab );
               }
               
               bool added=false;
               for(unsigned i=0;i<getNumberOfArguments();++i){
                   // Add this function to jobs to do in recursive loop in previous action
                   if( getPntrToArgument(i)->getRank()>0 ){
                       if( (getPntrToArgument(i)->getPntrToAction())->addActionToChain( alabels, this ) ){ added=true; break; } 
                   }
               }  
               plumed_massert(added, "could not add action " + getLabel() + " to chain of any of its arguments");

               // Now make sure we have the derivative size correct
               nderivatives=0; 
               for(unsigned i=0;i<distinct_arguments.size();++i) {
                   if( distinct_arguments[i].second==0 ) nderivatives += distinct_arguments[i].first->getNumberOfDerivatives();
                   else nderivatives += distinct_arguments[i].first->getFullNumberOfTasks();
               }
               // Set forces to apply to correct size
               forcesToApply.resize( nderivatives );
           }
       }
  }
  // This creates a group of atoms that have these weights -- not entirely foolproof and could be improved GAT
  bool checkforrank=false;
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()>0 && !getPntrToArgument(i)->hasDerivatives() ){ checkforrank=true; break; }
  }
  if( checkforrank ) {
      std::string myat_group="none";
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()>0 && !getPntrToArgument(i)->hasDerivatives() ) {
              if( plumed.getAtoms().getAllGroups().count(getPntrToArgument(i)->getPntrToAction()->getLabel()) ) {
                  myat_group = getPntrToArgument(i)->getPntrToAction()->getLabel(); break;
              } 
          }
      }
      if( myat_group!="none" ) {
          const auto m=plumed.getAtoms().getAllGroups().find(myat_group ); 
          plumed.getAtoms().insertGroup( getLabel(), m->second );
      }
  }
}

std::vector<unsigned> Function::getShape() {
  unsigned maxrank=0, rmax=0;
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()>maxrank ){ maxrank=getPntrToArgument(i)->getRank(); rmax=i; }
  } 
  std::vector<unsigned> shape;
  if( hasGridOutput() ) {
     shape.resize( maxrank );
     for(unsigned i=0;i<shape.size();++i) shape[i] = getPntrToArgument(rmax)->getShape()[i];
  } else if( maxrank==0 || !numberedkeys ){ 
     shape.resize(0);
  } else { 
     shape.resize( maxrank );
     for(unsigned i=0;i<shape.size();++i) shape[i]=getPntrToArgument(rmax)->getShape()[i]; 
  }
  return shape;
}

bool Function::hasGridOutput() const {
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) return true;
  }
  return false; 
}

void Function::addValueWithDerivatives() {
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");

  std::vector<std::string> period; 
  if( keywords.exists("PERIODIC") ){
     parseVector("PERIODIC",period);
     if( period.size()==1 ){
        if( period[0]!="NO") error("input to PERIODIC keyword does not make sense");
     } else if( period.size()!=2 ) error("input to PERIODIC keyword does not make sense"); 
  } else { period.resize(1); period[0]="NO"; }

  std::vector<unsigned> shape( getShape() ); 
  if( arg_ends.size()==0 ) {
      if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addValueWithDerivatives( shape );
      else if( hasGridOutput() ) ActionWithValue::addValueWithDerivatives( shape ); 
      else if( actionInChain() && shape.size()>0 ) ActionWithValue::addValue( shape ); 
      else if( shape.size()==0 ) ActionWithValue::addValueWithDerivatives( shape );
      else ActionWithValue::addValue( shape ); 
      if(period.size()==1 && period[0]=="NO") setNotPeriodic(); 
      else if(period.size()==2) setPeriodic(period[0],period[1]);
  } else if( arg_ends[1]-arg_ends[0]==1 ) {
      if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addValueWithDerivatives( shape );
      else if( hasGridOutput() ) ActionWithValue::addValueWithDerivatives( shape );
      else if( actionInChain() && shape.size()>0 ) ActionWithValue::addValue( shape ); 
      else if( shape.size()==0 ) ActionWithValue::addValueWithDerivatives( shape );
      else ActionWithValue::addValue( shape );
      if(period.size()==1 && period[0]=="NO") setNotPeriodic();
      else if(period.size()==2) setPeriodic(period[0],period[1]);
  } else {
      std::string num; 
      for(unsigned i=0;i<arg_ends.size()-1;++i){
          Tools::convert(i+1,num); 
          if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addComponentWithDerivatives( "arg_" + num, shape );
          else if( hasGridOutput() ) ActionWithValue::addComponentWithDerivatives( "arg_" + num, shape );
          else if( actionInChain() && shape.size()>0 ) ActionWithValue::addComponent( "arg_" + num, shape ); 
          else if( shape.size()==0 ) ActionWithValue::addComponentWithDerivatives( "arg_" + num, shape );
          else ActionWithValue::addComponent( "arg_" + num, shape );
          if(period.size()==1 && period[0]=="NO") componentIsNotPeriodic( "arg_" + num ); 
          else if(period.size()==2) componentIsPeriodic("arg_" + num, period[0], period[1]);
      } 
  }
}

void Function::addComponentWithDerivatives( const std::string& name ) {
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");

  std::vector<unsigned> shape( getShape() );
  if( arg_ends.size()==0 ){ 
      if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name,shape);
      else if( hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name,shape );
      else if( actionInChain() && shape.size()>0 ) ActionWithValue::addComponent(name,shape); 
      else if( shape.size()==0 ) ActionWithValue::addComponentWithDerivatives(name,shape); 
      else ActionWithValue::addComponent(name,shape);
  } else if( arg_ends[1]-arg_ends[0]==1 ) {
      if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name,shape);
      else if( hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name,shape );
      else if( actionInChain() && shape.size()>0 ) ActionWithValue::addComponent(name,shape);
      else if( shape.size()==0 ) ActionWithValue::addComponentWithDerivatives(name,shape);
      else ActionWithValue::addComponent(name,shape);
  } else { 
      std::string num; 
      for(unsigned i=0;i<arg_ends.size()-1;++i){ 
          Tools::convert(i+1,num);
          if( actionInChain() && shape.size()>0 && hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name + "_arg_" + num, shape);
          else if( hasGridOutput() ) ActionWithValue::addComponentWithDerivatives(name + "_arg_" + num, shape);
          else if( actionInChain() && shape.size()>0 ) ActionWithValue::addComponent( name + "_arg_" + num, shape ); 
          else if( shape.size()==0 ) ActionWithValue::addComponentWithDerivatives(name + "_arg_" + num, shape);
          else ActionWithValue::addComponent( name + "_arg_" + num, shape ); 
      } 
  }
}

void Function::evaluateAllFunctions() {
  if( firststep ) {
      std::vector<unsigned> shape( getShape() );
      unsigned ival = getPntrToOutput(0)->getNumberOfValues( getLabel() );
      getPntrToOutput(0)->setShape( shape ); firststep=false;
      if( ival<getPntrToOutput(0)->getNumberOfValues( getLabel() ) ) {
          for(unsigned j=ival;j<getPntrToOutput(0)->getNumberOfValues( getLabel() );++j) addTaskToList(j);
      }
  }
  runAllTasks();
}

void Function::calculate(){
  // Everything is done elsewhere
  if( hasAverageAsArgument() || actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  evaluateAllFunctions(); 
}

void Function::update() {
  if( !hasAverageAsArgument() ) return;
  plumed_dbg_assert( !actionInChain() && getFullNumberOfTasks()>0 ); 
  evaluateAllFunctions();
}

void Function::runFinalJobs() {
  if( !hasAverageAsArgument() ) return;
  plumed_dbg_assert( !actionInChain() && getFullNumberOfTasks()>0 );
  evaluateAllFunctions();
}

void Function::buildCurrentTaskList( std::vector<unsigned>& tflags ) {
  if( !actionInChain() ) tflags.assign(tflags.size(),1);
}

void Function::getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min, 
                                     std::vector<std::string>& max, std::vector<unsigned>& nbin, 
                                     std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const { 
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( argn, min, max, nbin, spacing, pbc, dumpcube );
}

void Function::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const { 
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointIndicesAndCoordinates( ind, indices, coords );
}

void Function::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  (getPntrToArgument(0)->getPntrToAction())->getGridPointAsCoordinate( ind, false, coords );
  if( coords.size()==(getPntrToOutput(0)->getRank()+1) ) coords[getPntrToOutput(0)->getRank()] = getPntrToOutput(0)->get(ind);
  else if( setlength ) {
    double val=getPntrToOutput(0)->get(ind);
    for(unsigned i=0;i<coords.size();++i) coords[i] = val*coords[i];
  } 
}

void Function::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Get the values of all the arguments
  bool matout=false, matinp=false;
  if( actionInChain() ) {
      matinp=getPntrToArgument(0)->getRank()==2 && !getPntrToArgument(0)->hasDerivatives();
#ifdef DNDEBUG
      if( matinp ){
          for(unsigned i=1;i<getNumberOfArguments();++i) plumed_dbg_assert( getPntrToArgument(i)->getRank()==2 && !getPntrToArgument(0)->hasDerivatives() );
      }
#endif
      if( matinp ) {
          matout=getPntrToOutput(0)->getRank()==2;
#ifdef DNDEBUG
          if( matout ){
              for(unsigned i=1;i<getNumberOfComponents();++i) plumed_dbg_assert( getPntrToOutput(i)->getRank()==2 );
          }
#endif
      }
  }
  // Calculate whatever we are calculating
  if( (matinp && !myvals.inVectorCall()) || !matinp ){
       std::vector<double> args( getNumberOfArgumentsPerTask() ); retrieveArguments( myvals, args );
       calculateFunction( args, myvals );
       // Make sure grid derivatives are updated 
       if( getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() ) {
           unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
           for(unsigned i=0;i<nderivatives;++i) myvals.updateIndex( ostrn, i ); 
           return;
       } 
  }
  // And update the dynamic list
  if( doNotCalculateDerivatives() ) return ;
  if( actionInChain() ) {
      if( (matinp && matout && !myvals.inVectorCall()) || !matinp ) {
           unsigned der_start=0;
           for(unsigned i=0;i<distinct_arguments.size();++i){
               unsigned istrn, jvalind;
               for(unsigned j=0;j<getNumberOfArguments();++j) {
                   if( getPntrToArgument(j)->getPntrToAction()==distinct_arguments[i].first ){ 
                       istrn = getArgumentPositionInStream(j,myvals); 
                       jvalind = j; break; 
                   }
               }
               for(unsigned k=0;k<myvals.getNumberActive(istrn);++k){
                   unsigned kind = myvals.getActiveIndex(istrn,k);
                   for(unsigned j=0;j<getNumberOfComponents();++j){
                       unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
                       myvals.updateIndex( ostrn, der_start + kind );
                   }
               }
               if( distinct_arguments[i].second==0 ) der_start += distinct_arguments[i].first->getNumberOfDerivatives();
               else der_start += getPntrToArgument(jvalind)->getNumberOfValues(getLabel());
           }
      } else if( (matinp && matout && myvals.inVectorCall()) ) {
           unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
           std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( nmat ) ); unsigned der_start=0, ntot_mat=0;
           if( mat_indices.size()<getNumberOfDerivatives() ) mat_indices.resize( getNumberOfDerivatives() );
           for(unsigned i=0;i<distinct_arguments.size();++i){
               unsigned istrn = (distinct_arguments[i].first->copyOutput(0))->getPositionInMatrixStash();
               std::vector<unsigned>& imat_indices( myvals.getMatrixIndices( istrn ) );
               for(unsigned k=0;k<myvals.getNumberOfMatrixIndices( istrn );++k) mat_indices[ntot_mat + k] = der_start + imat_indices[k]; 
               ntot_mat += myvals.getNumberOfMatrixIndices( istrn ); der_start += distinct_arguments[i].first->getNumberOfDerivatives();
           }
           myvals.setNumberOfMatrixIndices( nmat, ntot_mat ); 
      } else if( myvals.inVectorCall() ) { 
           for(unsigned i=0;i<distinct_arguments.size();++i){
               unsigned der_start = 0;
               unsigned istrn = (distinct_arguments[i].first->copyOutput(0))->getPositionInMatrixStash();
               std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
               for(unsigned k=0;k<myvals.getNumberOfMatrixIndices( istrn );++k){
                   for(unsigned j=0;j<getNumberOfComponents();++j){
                       unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
                       myvals.updateIndex( ostrn, der_start + mat_indices[k] );
                   }
               }
               der_start += distinct_arguments[i].first->getNumberOfDerivatives();
           }
      }
  } else {
      if( arg_ends.size()>0 ) {
          unsigned base=0;
          for(unsigned i=0;i<arg_ends.size()-1;++i) {
              for(unsigned j=0;j<getNumberOfComponents();++j){
                  unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
                  if( arg_ends[i+1]==(arg_ends[i]+1) && getPntrToArgument(arg_ends[i])->getRank()==0 ) {
                      myvals.updateIndex( ostrn, base );  
                  } else { 
                      myvals.updateIndex( ostrn, base + myvals.getTaskIndex() );
                  }
              }
              for(unsigned k=arg_ends[i];k<arg_ends[i+1];++k) base += getPntrToArgument(k)->getSize();
          } 
      } else {
          for(unsigned j=0;j<getNumberOfComponents();++j){
              unsigned ostrn = getPntrToOutput(j)->getPositionInStream();
              for(unsigned i=0;i<nderivatives;++i) myvals.updateIndex( ostrn, i ); 
          }
      }
  }
}

void Function::gatherGridAccumulators( const unsigned& code, const MultiValue& myvals,
                                       const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( getNumberOfComponents()==1 && getPntrToOutput(0)->getRank()>0 && getPntrToOutput(0)->hasDerivatives() );
  unsigned nder = getPntrToOutput(0)->getRank(), ostr = getPntrToOutput(0)->getPositionInStream();
  unsigned kp = bufstart + code*(1+nderivatives); buffer[kp] += myvals.get( ostr );
  for(unsigned i=0;i<nderivatives;++i) buffer[kp + 1 + i] += myvals.getDerivative( ostr, i ); 
}

void Function::apply()
{
  // Everything is done elsewhere
  if( doNotCalculateDerivatives() ) return;
  // And add forces
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( forcesToApply, ss ); 

//   const unsigned noa=getNumberOfArguments();
//   const unsigned ncp=getNumberOfComponents();
//   const unsigned cgs=comm.Get_size();
// 
//   vector<double> f(noa,0.0);
// 
//   unsigned stride=1;
//   unsigned rank=0;
//   if(ncp>4*cgs) {
//     stride=comm.Get_size();
//     rank=comm.Get_rank();
//   }
// 
//   unsigned at_least_one_forced=0;
//   #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(f)
//   {
//     vector<double> omp_f(noa,0.0);
//     vector<double> forces(noa);
//     #pragma omp for reduction( + : at_least_one_forced)
//     for(unsigned i=rank; i<ncp; i+=stride) {
//       if(getPntrToComponent(i)->applyForce(forces)) {
//         at_least_one_forced+=1;
//         for(unsigned j=0; j<noa; j++) omp_f[j]+=forces[j];
//       }
//     }
//     #pragma omp critical
//     for(unsigned j=0; j<noa; j++) f[j]+=omp_f[j];
//   }
// 
//   if(noa>0&&ncp>4*cgs) { comm.Sum(&f[0],noa); comm.Sum(at_least_one_forced); }

//  if(at_least_one_forced>0) for(unsigned i=0; i<noa; ++i) getPntrToArgument(i)->addForce(f[i]);
}

}
}
