/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_tools_LeptonCall_h
#define __PLUMED_tools_LeptonCall_h

#include "core/Action.h"
#include "lepton/Lepton.h"
#include "View.h"

namespace PLMD {

/// \ingroup TOOLBOX
class LeptonCall {
private:
  unsigned nargs;
  bool allow_extra_args;
/// Lepton expression.
/// \warning Since lepton::CompiledExpression is mutable, a vector is necessary for multithreading!
  std::vector<lepton::CompiledExpression> expression;
/// Lepton expression for derivative
/// \warning Since lepton::CompiledExpression is mutable, a vector is necessary for multithreading!
  std::vector<std::vector<lepton::CompiledExpression> > expression_deriv;
  std::vector<double*> lepton_ref;
  std::vector<double*> lepton_ref_deriv;
public:
  void set(const std::string & func,
           const std::vector<std::string>& var,
           Action* action=nullptr,
           bool extraArgs=false );
  unsigned getNumberOfArguments() const ;
  double evaluate( const std::vector<double>& args ) const ;
  double evaluate( View<const double> args ) const ;
  double evaluateDeriv( unsigned ider, View<const double> args ) const ;
  double evaluateDeriv( unsigned ider, const std::vector<double>& args ) const ;
};

inline
unsigned LeptonCall::getNumberOfArguments() const {
  return nargs;
}

}

#endif

