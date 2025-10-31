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

The GET command takes in the label of a Value and transfers the contents of the Value to
to a void pointer that is accessible in the code that called PLUMED. As the calling code does not know the shape of the Value in advance,
we provide the functionality to get the data rank and shape. The following Python snippet illustrates how this works in practice:

```python
import plumed

# Create a PLUMED object
p = plumed.Plumed()
# Setup PLUMED
num_atoms = 10
p.cmd("setNatoms",num_atoms)
p.cmd("setLogFile","test.log")
p.cmd("init")
# Tell PLUMED to calculate the distance between two atoms
p.cmd("readInputLine", "d1: DISTANCE ATOMS=1,2")
# Get the rank of the PLMD::Value that holds the distance
# This command sets up the GET object
rank = np.zeros( 1, dtype=np.int_ )
p.cmd("getDataRank d1", rank )
# Now get the shape of the PLMD::Value d1 that we are asking for in the GET object
shape = np.zeros( rank, dtype=np.int_ )
p.cmd("getDataShape d1", shape )
# And now set the void pointer that the data in PLMD::Value d1 should be
# transferred to so it can be accessed in our python script when asking PLMD to do a calculation
d1 = np.zeros( shape )
p.cmd("setMemoryForData d1", data )

# if we now transfer some atomic positions to plumed and call calc the variable d1 is set equal to the distance between atom 1 and atom 2.
```

Notice that you can have as many GET actions as you need. The data is transferred from the PLMD::Value to the void pointer when the `calculate` method of GET is called.
Transferring variables is thus seamlessly integrated into the PLUMED calculation cycle.

You would only use the GET command if you were calling PLUMED from python or an MD code. The equivalent commands that you would use for this action in a conventional PLUMED input file as follows.

```plumed
d: DISTANCE ATOMS=1,2
GET ARG=d STRIDE=1 TYPE=value
```

!!! warning "TYPE not fully implemented"

    At the moment you can only use GET to extract the values. The TYPE keyword was added so that we could support passing of the forces
    on a value or the derivatives of the value.  Currently, we do not support these options. If you are interested in using this feature
    please let us know as we could implement this functionality here relatively easily.

*/
//+ENDPLUMEDOC

namespace PLMD {

PLUMED_REGISTER_ACTION(ActionToGetData,"GET")

void ActionToGetData::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("optional","ARG","scalar/vector/matrix/grid","the label of the value that you would like to GET");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be stored");
  keys.add("compulsory","TYPE","value","what do you want to collect for the value can be derivative/force");
  keys.setValueDescription("scalar/vector/matrix/grid","a copy of the data in the value specified by the ARG keyword");
}

ActionToGetData::ActionToGetData(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  mydata(DataPassingObject::create(plumed.getRealPrecision())) {
  std::string type;
  parse("TYPE",type);
  if( type=="value" ) {
    gtype=dataType::val;
  } else if( type=="derivatives" ) {
    gtype=dataType::deriv;
  } else if( type=="forces" ) {
    gtype=dataType::force;
  } else {
    plumed_merror("cannot get " + type + " for value TYPE should be value/derivative/force");
  }

  if( gtype!=dataType::val ) {
    error("not implemented functionality to pass derviatives or forces to python.  Email gareth.tribello@gmail.com if you want this.");
  }

  if( getNumberOfArguments()!=1 ) {
    error("python interface works best when you ask for one argument at a time");
  }
  if( getPntrToArgument(0)->getNumberOfValues()==0 ) {
    error("cannot get data as shape of value " + getPntrToArgument(0)->getName() + " has not been set");
  }
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
  plumed_assert( gtype==dataType::val );
  mydata->setData( getPntrToArgument(0) );
}

}
