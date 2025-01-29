/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC GRIDCALC ACCUMULATE
/*
Sum the elements of this value over the course of the trajectory

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Accumulate :
  public ActionWithValue,
  public ActionWithArguments,
  public ActionPilot {
private:
  bool clearnextstep;
  unsigned clearstride;
public:
  static void registerKeywords( Keywords& keys );
  Accumulate( const ActionOptions& );
  unsigned getNumberOfDerivatives();
  bool calculateOnUpdate() override {
    return false;
  }
  bool calculateConstantValues( const bool& have_atoms ) override {
    return false;
  }
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(Accumulate,"ACCUMULATE")

void Accumulate::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.addInputKeyword("compulsory","ARG","scalar/grid","the label of the argument that is being added to on each timestep");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.setValueDescription("scalar/grid","a sum calculated from the time series of the input quantity");
}

Accumulate::Accumulate( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  clearnextstep(true) {
  if( getNumberOfArguments()!=1 ) {
    error("there should only be one argument to this action");
  }
  if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()!=0 ) {
    error("input to the accumulate action should be a scalar or a grid");
  }

  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) {
      error("CLEAR parameter must be a multiple of STRIDE");
    }
    log.printf("  clearing average every %u steps \n",clearstride);
  }
  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
  addValueWithDerivatives( shape );
  setNotPeriodic();
  if( getPntrToArgument(0)->isPeriodic() ) {
    error("you cannot accumulate a periodic quantity");
  }
}

unsigned Accumulate::getNumberOfDerivatives() {
  if( getPntrToArgument(0)->getRank()>0 ) {
    return getPntrToArgument(0)->getNumberOfGridDerivatives();
  }
  return getPntrToArgument(0)->getNumberOfDerivatives();
}

void Accumulate::update() {
  if( clearnextstep ) {
    if( getPntrToComponent(0)->getNumberOfValues()!=getPntrToArgument(0)->getNumberOfValues() ) {
      getPntrToComponent(0)->setShape( getPntrToArgument(0)->getShape() );
    }
    clearnextstep=false;
    getPntrToComponent(0)->set(0,0.0);
    getPntrToComponent(0)->clearDerivatives(true);
  }
  if( getStep()==0 ) {
    return;
  }

  Value* myarg=getPntrToArgument(0);
  Value* myout = getPntrToComponent(0);
  if( getPntrToArgument(0)->getRank()>0 ) {
    unsigned nvals = myarg->getNumberOfValues(), nder = myarg->getNumberOfGridDerivatives();
    for(unsigned i=0; i<nvals; ++i) {
      myout->set( i, myout->get(i) + myarg->get(i) );
      for(unsigned j=0; j<nder; ++j) {
        myout->addGridDerivatives( i, j, myarg->getGridDerivative( i, j ) );
      }
    }
  } else {
    getPntrToComponent(0)->add( getPntrToArgument(0)->get() );
  }

  // Clear if required
  if( clearstride>0 && getStep()%clearstride==0 ) {
    clearnextstep=true;
  }
}

}
}
