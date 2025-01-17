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

//+PLUMEDOC ANALYSIS COLLECT
/*
Collect data from the trajectory for later analysis

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Collect :
  public ActionWithValue,
  public ActionWithArguments,
  public ActionPilot {
private:
  bool usefirstconf;
  unsigned clearstride;
public:
  static void registerKeywords( Keywords& keys );
  Collect( const ActionOptions& );
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

PLUMED_REGISTER_ACTION(Collect,"COLLECT")

void Collect::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.use("ARG");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("compulsory","TYPE","auto","required if you are collecting an object with rank>0. Should be vector/matrix and determines how data is stored.  If rank==0 then data has to be stored as a vector");
  keys.setValueDescription("the time series for the input quantity");
}

Collect::Collect( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  usefirstconf(false) {
  if( getNumberOfArguments()!=1 ) {
    error("there should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
    error("input to the collect argument cannot be a grid");
  }

  std::string type;
  parse("TYPE",type);
  if( getPntrToArgument(0)->getNumberOfValues()==1 && (type=="auto" || type=="vector") ) {
    type="vector";
  } else if( getPntrToArgument(0)->getNumberOfValues()==1 && type=="matrix" ) {
    error("invalid type specified. Cannot construct a matrix by collecting scalars");
  } else if(  getPntrToArgument(0)->getNumberOfValues()!=1 && type=="auto" ) {
    error("missing TYPE keyword.  TYPE should specify whether data is to be stored as a vector or a matrix");
  } else if( type!="vector" && type!="matrix" ) {
    error("invalid TYPE specified. Should be matrix/scalar found " + type);
  }

  if( type=="vector" ) {
    log.printf("  adding %d elements to stored vector each time we collect\n", getPntrToArgument(0)->getNumberOfValues() );
  } else {
    log.printf("  constructing matrix with rows of length %d from input data\n", getPntrToArgument(0)->getNumberOfValues() );
  }

  parse("CLEAR",clearstride);
  unsigned nvals=0;
  if( clearstride==getStride() ) {
    nvals=1;
    usefirstconf=(getStride()==0);
  } else if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) {
      error("CLEAR parameter must be a multiple of STRIDE");
    }
    log.printf("  clearing collected data every %u steps \n",clearstride);
    nvals=(clearstride/getStride());
  }

  std::vector<unsigned> shape(1);
  shape[0]=nvals;
  getPntrToArgument(0)->buildDataStore();
  if( type=="matrix" ) {
    shape.resize(2);
    shape[1] = getPntrToArgument(0)->getNumberOfValues();
  }
  if( type=="vector" ) {
    shape[0] = nvals*getPntrToArgument(0)->getNumberOfValues();
  }
  addValue( shape );
  if( shape.size()==2 ) {
    getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  }
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
}

unsigned Collect::getNumberOfDerivatives() {
  return 0;
}

void Collect::update() {
  if( getStep()==0 || (!onStep() && !usefirstconf) ) {
    return ;
  }
  usefirstconf=false;

  Value* myin=getPntrToArgument(0);
  Value* myout=getPntrToComponent(0);
  unsigned nargs=myin->getNumberOfValues();
  if( clearstride==getStride() ) {
    for(unsigned i=0; i<nargs; ++i) {
      myout->set( i, myin->get(i) );
    }
  } else if( clearstride>0 ) {
    unsigned step = getStep() - clearstride*std::floor( getStep() / clearstride );
    if( getStep()%clearstride==0 ) {
      step = step + clearstride;
    }
    unsigned base = (step/getStride()-1)*nargs;
    for(unsigned i=0; i<nargs; ++i) {
      myout->set( base+i, myin->get(i) );
    }
  } else {
    for(unsigned i=0; i<nargs; ++i) {
      myout->push_back( myin->get(i) );
    }
    if( myout->getRank()==2 ) {
      myout->reshapeMatrixStore( nargs );
    }
  }
}

}
}
