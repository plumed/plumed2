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
  ActionWithArguments(ao)
{
  createTasksFromArguments(); nderivatives = getNumberOfArguments();
  // Now create the stream of jobs to work through 
  if( done_over_stream ){   // getFullNumberOfTasks()>0   // This is for if we have a function that needs to store - needs though GAT
      nderivatives = 0; std::vector<ActionWithValue*> tvals; 
      for(unsigned i=0;i<getNumberOfArguments();++i){
          // Add this function to jobs to do in recursive loop in previous action
          if( getPntrToArgument(i)->getRank()>0 ) (getPntrToArgument(i)->getPntrToAction())->addActionToChain( this ); 
          // Check for number of derivatives 
          bool found=false; std::string mylabstr = (getPntrToArgument(i)->getPntrToAction())->getLabel();
          for(unsigned j=0;j<tvals.size();++j){
             if( mylabstr==tvals[j]->getLabel() ){ found=true; break; }
          }
          if( !found ){ 
              tvals.push_back( getPntrToArgument(i)->getPntrToAction() );
              nderivatives += (getPntrToArgument(i)->getPntrToAction())->getNumberOfDerivatives(); 
          }
      }  
  }
}

std::vector<unsigned> Function::getShape() {
  bool rank0=false; unsigned size=0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
      if( getPntrToArgument(i)->getRank()==0 ) rank0=true;
      else if( size>0 && getPntrToArgument(i)->getRank()>1 ) error("cannot have more than one 2 or more rank inputs to function");
      else size += getPntrToArgument(i)->getSize();
  }
  std::vector<unsigned> shape;
  if( rank0 ){ 
     shape.resize(0);
  } else if( getPntrToArgument(0)->getRank()>1 ){ 
     plumed_assert( getNumberOfArguments()==1 ); shape.resize( getPntrToArgument(0)->getRank() );
     for(unsigned i=0;i<shape.size();++i) shape[i]=getPntrToArgument(0)->getShape()[i]; 
  } else {
     shape.resize(1); shape[0]=size;
  }
  return shape;
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
  if( arg_ends[1]-arg_ends[0]==1 ){ 
      if( done_over_stream ) ActionWithValue::addValue( shape ); 
      else ActionWithValue::addValueWithDerivatives( shape ); 
      if(period.size()==1 && period[0]=="NO") setNotPeriodic(); 
      else if(period.size()==2) setPeriodic(period[0],period[1]);
  } else {
      std::string num; 
      for(unsigned i=0;i<arg_ends.size()-1;++i){
          Tools::convert(i+1,num); 
          if( done_over_stream ) ActionWithValue::addComponent( "arg_" + num, shape ); 
          else ActionWithValue::addComponentWithDerivatives( "arg_" + num, shape );
          if(period.size()==1 && period[0]=="NO") componentIsNotPeriodic( "arg_" + num ); 
          else if(period.size()==2) componentIsPeriodic("arg_" + num, period[0], period[1]);
      } 
  }
}

void Function::addComponentWithDerivatives( const std::string& name ) {
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");

  std::vector<unsigned> shape( getShape() );
  if( arg_ends[1]-arg_ends[0]==1 ){ 
      shape.resize(0); 
      if( done_over_stream ) ActionWithValue::addComponent(name,shape); 
      else ActionWithValue::addComponentWithDerivatives(name,shape); 
  } else { 
      std::string num; shape.resize(0); 
      for(unsigned i=0;i<arg_ends.size()-1;++i){ 
          Tools::convert(i+1,num);
          if( done_over_stream ) ActionWithValue::addComponent( name + "_arg_" + num, shape ); 
          else ActionWithValue::addComponentWithDerivatives( name + "_arg_" + num, shape ); 
      } 
  }
}

void Function::calculate(){
  // Everything is done elsewhere
  if( done_over_stream ) return;
  // This is done if we are calculating a function of multiple cvs
  plumed_dbg_assert( getFullNumberOfTasks()>0 ); runAllTasks(); 
}

void Function::buildCurrentTaskList( std::vector<unsigned>& tflags ) const {
  if( !done_over_stream ) tflags.assign(tflags.size(),1);
}

void Function::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Get the values of all the arguments
  std::vector<double> args( arg_ends.size()-1 ); retrieveArguments( myvals, args );
  // Calculate whatever we are calculating
  calculateFunction( args, myvals ); 
  // And update the dynamic list
  if( !doNotCalculateDerivatives() && !done_over_stream ) myvals.updateDynamicList();
}

void Function::apply()
{
  // Everything is done elsewhere
  if( done_over_stream ) return;
  const unsigned noa=getNumberOfArguments();
  const unsigned ncp=getNumberOfComponents();
  const unsigned cgs=comm.Get_size();

  vector<double> f(noa,0.0);

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>4*cgs) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned at_least_one_forced=0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(f)
  {
    vector<double> omp_f(noa,0.0);
    vector<double> forces(noa);
    #pragma omp for reduction( + : at_least_one_forced)
    for(unsigned i=rank; i<ncp; i+=stride) {
      if(getPntrToComponent(i)->applyForce(forces)) {
        at_least_one_forced+=1;
        for(unsigned j=0; j<noa; j++) omp_f[j]+=forces[j];
      }
    }
    #pragma omp critical
    for(unsigned j=0; j<noa; j++) f[j]+=omp_f[j];
  }

  if(noa>0&&ncp>4*cgs) { comm.Sum(&f[0],noa); comm.Sum(at_least_one_forced); }

  if(at_least_one_forced>0) for(unsigned i=0; i<noa; ++i) getPntrToArgument(i)->addForce(f[i]);
}

}
}
