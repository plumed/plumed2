/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

class ComposeMatrix :
  public ActionWithValue,
  public ActionWithArguments 
{
public:
  static void registerKeywords(Keywords&);
  explicit ComposeMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  unsigned getNumberOfColumns() const override;
  bool canChainFromThisAction() const { return false; }
  void calculate();
  void apply();
};

PLUMED_REGISTER_ACTION(ComposeMatrix,"COMPOSE_MATRIX")

void ComposeMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); 
  ActionWithArguments::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES"); keys.use("ARG");
}

ComposeMatrix::ComposeMatrix(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{ 
  std::vector<unsigned> shape(2); shape[0]=arg_ends.size()-1; shape[1]=0;
  for(unsigned i=arg_ends[0];i<arg_ends[1];++i) shape[1] += getPntrToArgument(i)->getNumberOfValues(getLabel());
  for(unsigned j=1;j<shape[0];++j) {
      unsigned nv=0; for(unsigned i=arg_ends[j];i<arg_ends[j+1];++i) nv += getPntrToArgument(i)->getNumberOfValues(getLabel());
      if( nv!=shape[1] ) error("mismatch in numbers of arguments for rows of matrix");
  }
  for(unsigned i=0;i<getNumberOfArguments();++i) getPntrToArgument(i)->buildDataStore( getLabel() );
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();
}

unsigned ComposeMatrix::getNumberOfDerivatives() const {
  return 1;
}

unsigned ComposeMatrix::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[1];
}

void ComposeMatrix::calculate() {
  unsigned nn=0; Value* val=getPntrToOutput(0);
  for(unsigned j=0;j<arg_ends.size()-1;++j) {
      for(unsigned i=arg_ends[j];i<arg_ends[j+1];++i) {
          unsigned nvals = getPntrToArgument(i)->getNumberOfValues(getLabel());
          for(unsigned k=0;k<nvals;++k) { val->set( nn, getPntrToArgument(i)->get(k) ); nn++; }
      }
  }
}

void ComposeMatrix::apply() {
  unsigned nn=0; Value* val=getPntrToOutput(0);
  for(unsigned j=0;j<arg_ends.size()-1;++j) {
      for(unsigned i=arg_ends[j];i<arg_ends[j+1];++i) {
          unsigned nvals = getPntrToArgument(i)->getNumberOfValues(getLabel());
          for(unsigned k=0;k<nvals;++k) { getPntrToArgument(i)->addForce( k, val->getForce(nn) ); nn++; }
      }
  }
}

}
}

