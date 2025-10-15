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
  public ActionWithArguments {
protected:
  void setDerivative(int,double);
  void setDerivative(Value*,int,double);
//overriding explicitly the two functions avoids the [-Woverloaded-virtual] warning
/// the shape argument will be ignored
  void addValueWithDerivatives( const std::vector<std::size_t>& =std::vector<std::size_t>() ) override;
/// the shape will be ignored
  void addComponentWithDerivatives( const std::string& valname, const std::vector<std::size_t>& =std::vector<std::size_t>() )override;
public:
  explicit Function(const ActionOptions&);
  virtual ~Function() {}
  double getArgument( const unsigned& iarg );
  void apply() override;
  static void registerKeywords(Keywords&);
  unsigned getNumberOfDerivatives() override;
};

inline
void Function::setDerivative(Value*v,int i,double d) {
  v->addDerivative(i,d);
}

inline
void Function::setDerivative(int i,double d) {
  setDerivative(getPntrToValue(),i,d);
}

inline
unsigned Function::getNumberOfDerivatives() {
  unsigned narg=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    narg += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return narg;
}

inline
double Function::getArgument( const unsigned& iarg ) {
  return getPntrToArgument(iarg)->get();
}

}
}

#endif

