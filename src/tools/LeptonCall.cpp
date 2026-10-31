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

LeptonCall::LeptonCall(const LeptonCall& other):
  nargs(other.nargs),
  allow_extra_args(other.allow_extra_args),
  varnames(other.varnames),
  expression(other.expression),
  expression_deriv(other.expression_deriv) {
  linkVariables();
}

LeptonCall& LeptonCall::operator=(const LeptonCall& other) {
  if( this!=&other ) {
    nargs=other.nargs;
    allow_extra_args=other.allow_extra_args;
    varnames=other.varnames;
    expression=other.expression;
    expression_deriv=other.expression_deriv;
    linkVariables();
  }
  return *this;
}

void LeptonCall::linkVariables() {
  const unsigned nth=expression.size();
  lepton_ref.assign(nth*nargs,nullptr);
  for(unsigned t=0; t<nth; ++t) {
    for(unsigned j=0; j<nargs; ++j) {
      try {
        lepton_ref[t*nargs+j]=&expression[t].getVariableReference(varnames[j]);
      } catch(const PLMD::lepton::Exception& exc) {
// this is necessary since in some cases lepton things a variable is not present even though it is present
// e.g. func=0*x
      }
    }
  }
  lepton_ref_deriv.assign(nth*nargs*nargs,nullptr);
  for(unsigned i=0; i<expression_deriv.size(); i++) {
    plumed_assert( expression_deriv[i].size()==nth );
    for(unsigned t=0; t<nth; ++t) {
      for(unsigned j=0; j<nargs; ++j) {
        try {
          lepton_ref_deriv[i*nth*nargs + t*nargs+j]=&expression_deriv[i][t].getVariableReference(varnames[j]);
        } catch(const PLMD::lepton::Exception& exc) {
        }
      }
    }
  }
}

void LeptonCall::set(const std::string & func,
                     const std::vector<std::string>& var,
                     Action* action,
                     const bool extraArgs ) {
  allow_extra_args=extraArgs;
  varnames=var;
  nargs=var.size();
  const unsigned nth=OpenMP::getNumThreads();
  lepton::ParsedExpression pe=lepton::Parser::parse(func).optimize(lepton::Constants());
  if( action ) {
    action->log<<"  function as parsed by lepton: "<<pe<<"\n";
  }
  expression.assign(nth,pe.createCompiledExpression());
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
  expression_deriv.resize(var.size());
  for(unsigned i=0; i<var.size(); i++) {
    lepton::ParsedExpression pevar=lepton::Parser::parse(func).differentiate(var[i]).optimize(lepton::Constants());
    if( action ) {
      action->log<<"    "<<pevar<<"\n";
    }
    expression_deriv[i].assign(nth,pevar.createCompiledExpression());
  }
  linkVariables();
}

double LeptonCall::evaluate( const std::vector<double>& args ) const {
  return evaluate( View<const double>( args.data(), args.size() ) );
}

double LeptonCall::evaluate(const View<const double> args ) const {
  plumed_dbg_assert( allow_extra_args || args.size()==nargs );
  const unsigned t=OpenMP::getThreadNum(), tbas=t*nargs;
  plumed_dbg_assert( t<expression.size() );
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
  const unsigned t=OpenMP::getThreadNum(), dbas = ider*expression.size()*nargs + t*nargs;
  for(unsigned j=0; j<nargs; j++) {
    if(lepton_ref_deriv[dbas+j] ) {
      *lepton_ref_deriv[dbas+j] = args[j];
    }
  }
  return expression_deriv[ider][t].evaluate();
}

}
