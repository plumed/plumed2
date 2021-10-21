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
#ifndef __PLUMED_function_FunctionOfScalar_h
#define __PLUMED_function_FunctionOfScalar_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"

namespace PLMD {
namespace function {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new CV function, within it there is
\ref AddingAFunction "information" as to how to go about implementing a new function.
*/

template <class T>
class FunctionOfScalar :
  public ActionWithValue,
  public ActionWithArguments
{
private:
/// The function that is being computed
  T myfunc;
protected:
  void setDerivative(int,double);
  void setDerivative(Value*,int,double);
public:
  explicit FunctionOfScalar(const ActionOptions&);
  virtual ~FunctionOfScalar() {}
/// Get the label to write in the graph
  std::string writeInGraph() const override { return myfunc.getGraphInfo( getName() ); }
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
  void apply() override;
  static void registerKeywords(Keywords&);
  unsigned getNumberOfDerivatives() const override;
  void turnOnDerivatives() override;
};

template <class T>
void FunctionOfScalar<T>::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc; tfunc.registerKeywords( keys );
}

template <class T>
FunctionOfScalar<T>::FunctionOfScalar(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  myfunc.read( this ); 
  // Get the names of the components
  std::vector<std::string> components( keywords.getAllOutputComponents() );
  // Create the values to hold the output
  std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
  if( components.size()==0 && str_ind.size()==0 ) addValueWithDerivatives();
  else if ( components.size()==0 ) {
    for(unsigned j=0;j<str_ind.size();++j) addComponentWithDerivatives( str_ind[j] ); 
  } else { 
    std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
    for(unsigned i=0;i<components.size();++i) {
        if( str_ind.size()>0 ) {
            for(unsigned j=0;j<str_ind.size();++j) addComponentWithDerivatives( components[i] + str_ind[j] );
        } else if( components[i].find_first_of("_")!=std::string::npos ) {
            if( getNumberOfArguments()==1 ) addValueWithDerivatives(); 
            else { for(unsigned j=0; j<getNumberOfArguments(); ++j) addComponentWithDerivatives( getPntrToArgument(j)->getName() + components[i] ); }
        } else addComponentWithDerivatives( components[i] );
    } 
  }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this ); myfunc.setPrefactor( this, 1.0 );
}

template <class T>
void FunctionOfScalar<T>::setDerivative(Value*v,int i,double d) {
  v->addDerivative(i,d);
}

template <class T>
void FunctionOfScalar<T>::setDerivative(int i,double d) {
  setDerivative(getPntrToValue(),i,d);
}

template <class T>
void FunctionOfScalar<T>::turnOnDerivatives() {
  if( !myfunc.derivativesImplemented() ) error("derivatives have not been implemended for " + getName() );
  ActionWithValue::turnOnDerivatives(); 
}

template <class T>
unsigned FunctionOfScalar<T>::getNumberOfDerivatives() const {
  return getNumberOfArguments();
}

template <class T>
void FunctionOfScalar<T>::calculate() {
  std::vector<double> args( getNumberOfArguments() ); for(unsigned i=0;i<args.size();++i) args[i]=getPntrToArgument(i)->get();
  std::vector<double> vals( getNumberOfComponents() ); Matrix<double> derivatives( getNumberOfComponents(), getNumberOfArguments() );
  myfunc.calc( this, args, vals, derivatives );
  for(unsigned i=0;i<vals.size();++i) getPntrToOutput(i)->set(vals[i]);
  if( doNotCalculateDerivatives() ) return;

  for(unsigned i=0;i<vals.size();++i) { 
      Value* val = getPntrToComponent(i);
      for(unsigned j=0;j<args.size();++j) setDerivative( val, j, derivatives(i,j) ); 
  }
}

template <class T>
void FunctionOfScalar<T>::update() {
  if( skipUpdate() ) return;
  calculate();
}

template <class T>
void FunctionOfScalar<T>::runFinalJobs() {
  if( skipUpdate() ) return;
  calculate();
}

template <class T>
void FunctionOfScalar<T>::apply() {
  const unsigned noa=getNumberOfArguments();
  const unsigned ncp=getNumberOfComponents();
  const unsigned cgs=comm.Get_size();

  std::vector<double> f(noa,0.0);

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>4*cgs) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned at_least_one_forced=0;
  #pragma omp parallel num_threads(OpenMP::getNumThreads()) shared(f)
  {
    std::vector<double> omp_f(noa,0.0);
    std::vector<double> forces(noa);
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

  if(at_least_one_forced>0) for(unsigned i=0; i<noa; ++i) getPntrToArgument(i)->addForce(0, f[i]);
}

}
}
#endif
