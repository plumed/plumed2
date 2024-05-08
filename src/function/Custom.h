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

#include "FunctionTemplateBase.h"
#include "tools/LeptonCall.h"

namespace PLMD {
namespace function {

class Custom : public FunctionTemplateBase {
  std::string func;
  bool zerowhenallzero;
  LeptonCall function;
/// Check if only multiplication is done in function.  If only multiplication is done we can do some tricks
/// to speed things up
  std::vector<unsigned> check_multiplication_vars;
public:
  void registerKeywords( Keywords& keys ) override;
  std::string getGraphInfo( const std::string& lab ) const override;
  void read( ActionWithArguments* action ) override;
  bool getDerivativeZeroIfValueIsZero() const override;
  std::vector<Value*> getArgumentsToCheck( const std::vector<Value*>& args ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

}
}
#endif
