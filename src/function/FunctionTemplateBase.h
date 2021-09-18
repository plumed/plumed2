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
#ifndef __PLUMED_function_FunctionTemplateBase_h
#define __PLUMED_function_FunctionTemplateBase_h

#include "FunctionBase.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace function {

class FunctionTemplateBase { 
protected:
/// Parse a keyword from the input as a value
  template<class T>
  void parse( Action* action, const std::string&key, T&t );
/// Parse a keyword from the input as a vector
  template<class T>
  void parseVector( Action* action, const std::string&key,std::vector<T>&t);
/// Parse a keyword from the input as a flag
  void parseFlag( Action* action, const std::string&key, bool&t );
public:
/// Override this function if you have not implemented the derivatives
  virtual bool derivativesImplemented() { return true; }
  virtual void registerKeywords( Keywords& keys ) = 0;
  virtual void read( ActionWithArguments* action ) = 0;
  virtual unsigned getRank() = 0;
  virtual void setPeriodicityForOutputs( ActionWithValue* action ) = 0;
  virtual void setPrefactor( ActionWithArguments* action ) {}
  virtual void calc( const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const = 0;
};

template<class T>
void FunctionTemplateBase::parse( Action* action, const std::string&key, T&t ) {
  action->parse(key,t);
}

template<class T>
void FunctionTemplateBase::parseVector( Action* action, const std::string&key,std::vector<T>&t) {
  action->parseVector(key,t);
}

inline
void FunctionTemplateBase::parseFlag( Action* action, const std::string&key, bool&t ) {
  action->parseFlag(key,t);
}

}
}
#endif
