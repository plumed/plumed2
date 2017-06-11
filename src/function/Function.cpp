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
  createTasksFromArguments();
}

void Function::addValueWithDerivatives() {
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");

  std::vector<unsigned> shape(0);
  if( arg_ends[1]-arg_ends[0]==1 ){ shape.resize(0); ActionWithValue::addValueWithDerivatives( shape ); }
  else { shape.resize(1); shape[0]=getFullNumberOfTasks(); ActionWithValue::addValue( shape );} 
        
  if( keywords.exists("PERIODIC") ) {
      std::vector<std::string> period; parseVector("PERIODIC",period);
      if(period.size()==1 && period[0]=="NO") setNotPeriodic();
      else if(period.size()==2) setPeriodic(period[0],period[1]);
      else error("missing PERIODIC keyword");
  }
}

void Function::addComponentWithDerivatives( const std::string& name ) {
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");

  std::vector<unsigned> shape(0);
  if( arg_ends[1]-arg_ends[0]==1 ){ shape.resize(0); }
  else { shape.resize(1); shape[0]=getFullNumberOfTasks(); } 
  ActionWithValue::addComponentWithDerivatives(name,shape);
}

void Function::calculate(){
  // Everything is done elsewhere
  if( done_over_stream ) return;
  // Set all values equal to zero prior to computing
  // for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->set(0);
  // This is done if we are calculating a function of multiple cvs
  plumed_dbg_assert( getFullNumberOfTasks()>0 ); runAllTasks(); 
}

void Function::activateTasks( std::vector<unsigned>& tflags ) const {
  tflags.assign(tflags.size(),1);
}

void Function::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Get the values of all the arguments
  std::vector<double> args( arg_ends.size()-1 ); retrieveArguments( myvals, args );
  // Calculate whatever we are calculating
  calculateFunction( args, myvals ); 
  // And update the dynamic list
  if( !doNotCalculateDerivatives() ) myvals.updateDynamicList();
}

// void Function::accumulate( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
//   for(unsigned i=0;i<getNumberOfComponents();++i){
//       Value* val=getPntrToOutput(i); if( val->isDataStream() ) continue ;
//       unsigned sind = val->getPositionInStream(), bufstart = val->getBufferPosition(); buffer[bufstart] += myvals.get(sind);
//       for(unsigned k=0;k<myvals.getNumberActive();++k){ 
//           unsigned kindex = myvals.getActiveIndex(k); buffer[bufstart + 1 + kindex] += myvals.getDerivative(sind,kindex); 
//       }
//   }
// }

// void Function::finish( const std::vector<double>& buffer ){
//   for(unsigned i=0;i<getNumberOfComponents();++i){
//       Value* val=getPntrToComponent(i); if( val->isDataStream() ) continue ;
//       unsigned bufstart = getPntrToComponent(i)->getBufferPosition(); val->set( buffer[bufstart] );
//       for(unsigned j=0;j<val->getNumberOfDerivatives();++j) val->setDerivative( j, buffer[bufstart+1+j] ); 
//   }
// }

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
