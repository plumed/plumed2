/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "ActionToGetData.h"
#include "ActionRegister.h"
#include "PlumedMain.h"

//+PLUMEDOC ANALYSIS GET
/*
Get data from PLUMED for another code

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {

PLUMED_REGISTER_ACTION(ActionToGetData,"GET")

void ActionToGetData::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be stored");
  keys.add("compulsory","TYPE","value","what do you want to collect for the value can be derivative/force");
  keys.use("ARG");
  keys.setValueDescription("a copy of the data in the value specified by the ARG keyword");
}

ActionToGetData::ActionToGetData(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  mydata(DataPassingObject::create(plumed.getRealPrecision())) {
  std::string type;
  parse("TYPE",type);
  if( type=="value" ) {
    gtype=val;
  } else if( type=="derivatives" ) {
    gtype=deriv;
  } else if( type=="forces" ) {
    gtype=force;
  } else {
    plumed_merror("cannot get " + type + " for value TYPE should be value/derivative/force");
  }

  if( gtype!=val ) {
    error("not implemented functionality to pass derviatives or forces to python.  Email gareth.tribello@gmail.com if you want this.");
  }

  if( getNumberOfArguments()!=1 ) {
    error("python interface works best when you ask for one argument at a time");
  }
  if( getPntrToArgument(0)->getNumberOfValues()==0 ) {
    error("cannot get data as shape of value " + getPntrToArgument(0)->getName() + " has not been set");
  }
  getPntrToArgument(0)->buildDataStore();
  data.resize( getPntrToArgument(0)->getNumberOfValues() );
}

void ActionToGetData::get_rank( const TypesafePtr & dims ) {
  if( getPntrToArgument(0)->getRank()==0 ) {
    dims.set(long(1));
    return;
  }
  dims.set(long(getPntrToArgument(0)->getRank()));
}

void ActionToGetData::get_shape( const TypesafePtr & dims ) {
  if( getPntrToArgument(0)->getRank()==0 ) {
    dims.set(long(1));
    return;
  }
  auto dims_=dims.get<long*>( { getPntrToArgument(0)->getRank() } );
  for(unsigned j=0; j<getPntrToArgument(0)->getRank(); ++j) {
    dims_[j] = getPntrToArgument(0)->getShape()[j];
  }
}

void ActionToGetData::set_memory( const TypesafePtr & val ) {
  mydata->setValuePointer(val,getPntrToArgument(0)->getShape(),false);
}

void ActionToGetData::calculate() {
  plumed_assert( gtype==val );
  mydata->setData( getPntrToArgument(0) );
}

}
