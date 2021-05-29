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
#ifndef __PLUMED_core_ActionWithArguments_h
#define __PLUMED_core_ActionWithArguments_h

#include "Action.h"
#include "Value.h"
#include "tools/MultiValue.h"
#include <vector>

namespace PLMD {

/**
\ingroup MULTIINHERIT
This is used to create PLMD::Action objects that take the output from some other Action as input.
This is used in PLMD::Function and PLMD::Bias
 PLMD::Action objects that inherit from PLMD::ActionWithArguments take
 values and components calculated in other PLMD::Action objects and
 use this information to calculate some new function.  If you have
 only one list of arguments you should use the reserved keyword <b> ARG </b>
 when you use parseArgumentList.
*/

class AverageBase;
class ReweightBase;

class ActionWithArguments : public virtual Action {
  friend class ActionWithValue;
  friend class CollectFrames; 
private:
  bool allrankzero;
  std::vector<bool> usingAllArgs;
  std::vector<Value*> arguments;
  bool lockRequestArguments;
  AverageBase* theAverageInArguments;
  ReweightBase* theReweightBase;
  const ActionWithValue* thisAsActionWithValue;
  ActionWithValue* getFirstNonStream();
  Value* getArgumentForScalar(const unsigned n) const ;
protected:
  bool numberedkeys, done_over_stream;
  std::vector< std::pair<ActionWithValue*,unsigned> > distinct_arguments;
  std::vector<unsigned> arg_ends, arg_deriv_starts;
/// This changes the arg keyword in the pdb file
  void expandArgKeywordInPDB( const PDB& pdb );
/// Create a list of tasks from the argument streams
  void createTasksFromArguments();
/// Get the total number of input arguments
  unsigned getNumberOfScalarArguments() const ;
/// This is used to create a chain of actions that can be used to calculate a function/multibias
  unsigned setupActionInChain( const unsigned& argstart ) ;
/// Should we skip running calculate for this action
  bool skipCalculate() const ;
/// Should we skip running update for this action
  bool skipUpdate() const ;
/// Resize all the values for the final task
  void resizeForFinalTasks();
public:
/// Get the scalar product between the gradients of two variables
  double getProjection(unsigned i,unsigned j)const;
/// Registers the list of keywords
  static void registerKeywords( Keywords& keys );
/// Returns the value of an argument
  double getArgumentScalar( const unsigned n ) const;
/// Return a pointer to specific argument
  Value* getPntrToArgument( const unsigned n ) const ;
/// Set the force on the nth scalar argumetn
  void setForceOnScalarArgument(const unsigned n, const double& ff);
/// Returns the number of arguments
  virtual unsigned getNumberOfArguments() const ;
/// Get the number of arguments in each task
  unsigned getNumberOfArgumentsPerTask() const ;
/// Takes the difference taking into account pbc for arg i
  double difference(int, double, double) const;
/// Takes one value and brings it back into the pbc of argument i
  double bringBackInPbc(int i,double d1)const;
/// Parse a list of arguments
  void parseArgumentList(const std::string&key,std::vector<Value*>&args);
/// Parse a numbered list of arguments
  bool parseArgumentList(const std::string&key,int i,std::vector<Value*>&args);
/// Setup the dependencies
  void requestArguments(const std::vector<Value*> &arg, const bool& allow_streams, const unsigned& argstart=0 );
/// Set the forces on the arguments
  void setForcesOnArguments( const unsigned& argstart, const std::vector<double>& forces, unsigned& start );
/// This sets the forces on if the action is in a chain
  static void setForcesOnActionChain( const std::vector<double>& forces, unsigned& start, ActionWithValue* av );
/// This gets the component jcomp of argument iarg
  double retrieveRequiredArgument( const unsigned& iarg, const unsigned& jcomp ) const ;
public:
  explicit ActionWithArguments(const ActionOptions&);
  virtual ~ActionWithArguments() {}
/// Calculate the numerical derivatives
/// N.B. only pass an ActionWithValue to this routine if you know exactly what you
/// are doing.  The default will be correct for the vast majority of cases
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;
  void lockRequests() override;
  void unlockRequests() override;
/// Returns an array of pointers to the arguments
  virtual const std::vector<Value*>    & getArguments() const ;
/// Convert a list of argument names into a list of pointers to the values
  void interpretArgumentList(const std::vector<std::string>& c, std::vector<Value*>&arg);
///
  virtual void getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
/// Retrieve the argument values
  void retrieveArguments( const MultiValue& myvals, std::vector<double>& args, const unsigned& argstart ) const ;
/// This tells us which arguments must be treated as distinct in functions
  virtual bool mustBeTreatedAsDistinctArguments() const ;
/// Get the number of quantities that must be stored in input
  void getNumberOfStashedInputArguments( unsigned& nquants ) const ;
/// This gets the position of the argument in the stream
  unsigned getArgumentPositionInStream( const unsigned& jder, MultiValue& myvals ) const ;
};


inline
Value* ActionWithArguments::getPntrToArgument( const unsigned n ) const {
  return arguments[n];
}

inline
unsigned ActionWithArguments::getNumberOfScalarArguments() const {
  unsigned nscalars=0;
  for(unsigned i=0; i<arguments.size(); ++i) nscalars += arguments[i]->getNumberOfValues( getLabel() );
  return nscalars;
}

inline 
Value* ActionWithArguments::getArgumentForScalar(const unsigned n) const {
  unsigned nt = 0, nn = 0, j=0;
  for(unsigned i=0; i<arguments.size(); ++i) {
    nt += arguments[i]->getNumberOfValues( getLabel() );
    if( n<nt ) { j=i; break ; }
    nn += arguments[i]->getNumberOfValues( getLabel() );
  }
  return arguments[j];
}

inline
double ActionWithArguments::getArgumentScalar(const unsigned n) const {
  unsigned nt = 0, nn = 0, j=0;
  for(unsigned i=0; i<arguments.size(); ++i) {
    nt += arguments[i]->getNumberOfValues( getLabel() );
    if( n<nt ) { j=i; break ; }
    nn += arguments[i]->getNumberOfValues( getLabel() );
  }
  return arguments[j]->get( n - nn );
}

inline
unsigned ActionWithArguments::getNumberOfArguments()const {
  return arguments.size();
}

inline
double ActionWithArguments::difference(int i,double d1,double d2)const {
  return getArgumentForScalar(i)->difference(d1,d2);
}

inline
double ActionWithArguments::bringBackInPbc(int i,double d1)const {
  return getArgumentForScalar(i)->bringBackInPbc(d1);
}

inline
void ActionWithArguments::lockRequests() {
  lockRequestArguments=true;
}

inline
void ActionWithArguments::unlockRequests() {
  lockRequestArguments=false;
}

inline
const std::vector<Value*> & ActionWithArguments::getArguments() const {
  return arguments;
}

}

#endif
