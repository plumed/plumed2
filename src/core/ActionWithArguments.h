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

class ActionWithArguments:
  public virtual Action
{
  std::vector<Value*> arguments;
  bool lockRequestArguments;
protected:
/// This changes the arg keyword in the pdb file
  void expandArgKeywordInPDB( PDB& pdb );
public:
/// Get the scalar product between the gradients of two variables
  double getProjection(unsigned i,unsigned j)const;
/// Registers the list of keywords
  static void registerKeywords( Keywords& keys );
/// Returns the value of an argument
  double getArgument( const unsigned n ) const;
/// Return a pointer to specific argument
  Value* getPntrToArgument( const unsigned n );
/// Returns the number of arguments
  virtual unsigned getNumberOfArguments() const ;
/// Takes the difference taking into account pbc for arg i
  double difference(int, double, double) const;
/// Takes one value and brings it back into the pbc of argument i
  double bringBackInPbc(int i,double d1)const;
/// Parse a list of arguments
  void parseArgumentList(const std::string&key,std::vector<Value*>&args);
/// Parse a numbered list of arguments
  bool parseArgumentList(const std::string&key,int i,std::vector<Value*>&args);
/// Setup the dependencies
  void requestArguments(const std::vector<Value*> &arg);
/// Add forces to arguments (used in apply)
  void addForcesOnArguments( const std::vector<double>& forces );
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
};


inline
Value* ActionWithArguments::getPntrToArgument( const unsigned n ) {
  return arguments[n];
}

inline
double ActionWithArguments::getArgument(const unsigned n) const {
  return arguments[n]->get();
}

inline
unsigned ActionWithArguments::getNumberOfArguments()const {
  return arguments.size();
}

inline
double ActionWithArguments::difference(int i,double d1,double d2)const {
  return arguments[i]->difference(d1,d2);
}

inline
double ActionWithArguments::bringBackInPbc(int i,double d1)const {
  return arguments[i]->bringBackInPbc(d1);
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
