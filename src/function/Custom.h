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
#ifndef __PLUMED_function_Custom_h
#define __PLUMED_function_Custom_h

#include "FunctionSetup.h"
#include "tools/LeptonCall.h"

namespace PLMD {
namespace function {

class Custom {
public:
  std::string func;
  std::vector<std::string> var;
  bool zerowhenallzero;
  LeptonCall function;
  std::vector<unsigned> check_multiplication_vars;
  static void registerKeywords( Keywords& keys );
  static std::string getFunctionString( const std::string& func );
  static void read( Custom& f, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Custom& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout );
  Custom& operator=( const Custom& m ) {
    func = m.func;
    var = m.var;
    function.set( func, var, NULL );
    zerowhenallzero = m.zerowhenallzero;
    check_multiplication_vars = m.check_multiplication_vars;
    return *this;
  }
};

}
}
#endif
