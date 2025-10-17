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
#include "LeptonCall.h"
#include "OpenMP.h"

namespace PLMD {

void LeptonCall::set(const std::string & func,
                     const std::vector<std::string>& var,
                     Action* action,
                     const bool extraArgs ) {
  allow_extra_args=extraArgs;
  unsigned nth=OpenMP::getNumThreads();
  expression.resize(nth);
  expression_deriv.resize(var.size());
  // Resize the expression for the derivatives
  for(unsigned i=0; i<expression_deriv.size(); ++i) {
    expression_deriv[i].resize(OpenMP::getNumThreads());
  }
  nargs=var.size();

  lepton_ref.resize(nth*nargs,nullptr);
  lepton::ParsedExpression pe=lepton::Parser::parse(func).optimize(lepton::Constants());
  unsigned nt=0;
  if( action ) {
    action->log<<"  function as parsed by lepton: "<<pe<<"\n";
  }
  for(auto & e : expression) {
    e=pe.createCompiledExpression();
    for(unsigned j=0; j<var.size(); ++j) {
      try {
        lepton_ref[nt*var.size()+j]=&const_cast<lepton::CompiledExpression*>(&expression[nt])->getVariableReference(var[j]);
      } catch(const PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
      }
    }
    nt++;
  }
  for(auto & p : expression[0].getVariables()) {
    if(std::find(var.begin(),var.end(),p)==var.end()) {
      if( action ) {
        action->error("variable " + p + " is not defined");
      } else {
        plumed_merror("variable " + p + " is not defined in lepton function");
      }
    }
  }
  if( action ) {
    action->log<<"  derivatives as computed by lepton:\n";
  }
  lepton_ref_deriv.resize(nth*nargs*nargs,nullptr);
  for(unsigned i=0; i<var.size(); i++) {
    lepton::ParsedExpression pevar=lepton::Parser::parse(func).differentiate(var[i]).optimize(lepton::Constants());
    nt=0;
    if( action ) {
      action->log<<"    "<<pevar<<"\n";
    }
    for(auto & e : expression_deriv[i]) {
      e=pevar.createCompiledExpression();
      for(unsigned j=0; j<var.size(); ++j) {
        try {
          lepton_ref_deriv[i*OpenMP::getNumThreads()*var.size() + nt*var.size()+j]=&const_cast<lepton::CompiledExpression*>(&expression_deriv[i][nt])->getVariableReference(var[j]);
        } catch(const PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
        }
      }
      nt++;
    }
  }
}

double LeptonCall::evaluate( const std::vector<double>& args ) const {
  return evaluate( View<const double>( args.data(), args.size() ) );
}

double LeptonCall::evaluate(const View<const double> args ) const {
  plumed_dbg_assert( allow_extra_args || args.size()==nargs );
  const unsigned t=OpenMP::getThreadNum(), tbas=t*nargs;
  for(unsigned i=0; i<nargs; ++i) {
    if( lepton_ref[tbas+i] ) {
      *lepton_ref[tbas+i] = args[i];
    }
  }
  return expression[t].evaluate();
}

double LeptonCall::evaluateDeriv( const unsigned ider, const std::vector<double>& args ) const {
  return evaluateDeriv( ider, View<const double>( args.data(), args.size() ) );
}

double LeptonCall::evaluateDeriv( const unsigned ider, const View<const double> args ) const {
  plumed_dbg_assert( allow_extra_args || args.size()==nargs );
  plumed_dbg_assert( ider<nargs );
  const unsigned t=OpenMP::getThreadNum(), dbas = ider*OpenMP::getNumThreads()*nargs + t*nargs;
  for(unsigned j=0; j<nargs; j++) {
    if(lepton_ref_deriv[dbas+j] ) {
      *lepton_ref_deriv[dbas+j] = args[j];
    }
  }
  return expression_deriv[ider][t].evaluate();
}

}
