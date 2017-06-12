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
#ifndef __PLUMED_function_Function_h
#define __PLUMED_function_Function_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace function {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new CV function, within it there is
\ref AddingAFunction "information" as to how to go about implementing a new function.
*/

class Function:
  public ActionWithValue,
  public ActionWithArguments
{
private:
  unsigned nderivatives;
  std::vector<unsigned> getShape();
protected:
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name );
  void setValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const ;
public:
  static void registerKeywords(Keywords&);
  explicit Function(const ActionOptions&);
  virtual ~Function() {}
  void calculate();
  void buildCurrentTaskList( std::vector<unsigned>& tflags ) const ;
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  virtual void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const = 0;
  void apply();
  unsigned getNumberOfDerivatives();
};

inline
unsigned Function::getNumberOfDerivatives() {
  return nderivatives;
}

inline
void Function::setValue( const unsigned& ival, const double& val, MultiValue& myvals ) const {
  myvals.setValue( getPntrToOutput(ival)->getPositionInStream(), val ); 
}

inline
void Function::addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return ;

  if( done_over_stream ){ plumed_error(); return; }
  if( getPntrToArgument(0)->getRank()>0 ){ plumed_error(); return; }
  myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), arg_ends[jder] + myvals.getTaskIndex(), der ); 
}

}
}

#endif

