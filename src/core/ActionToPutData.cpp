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
#include "ActionToPutData.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"

//+PLUMEDOC ANALYSIS PUT
/*
Pass data into PLUMED

The PUT command transfers the data from a void pointer passed to PLUMED from the calling code to a PLMD::Value. The calling
object knows the shapes of the variables it is passing, so if you want to pass a 3x3 matrix from the MD code to PLUMED, you create the space to do so as follows:

```c++
plumed.cmd("readInputLine n: PUT SHAPE=3,3 UNIT=length PERIODIC=NO");
```

This command then creates a PLMD::Value called `n` that you can refer to later in your PLUMED input file. To transfer data from the void pointer called val into the PLMD::Value
called `n`, you would then use the following command:

```c++
plumed.cmd("setValue n", val);
```

Notice also that if you expect PLUMED to try to apply forces on `n`, you can pass a void pointer called `force` to get the forces that PLUMED has applied on the elements of n as follows:

```c++
plumed.cmd("setValueForces n", force);
```

Within the PLMD::Value `n`, the forces that PLUMED wishes to apply on the components of the input object are stored in the std::vector called `inputForce`. Furthermore, whenever a PLMD::Value
is created from a PUT action `storedata` is set to true. PUT also has a CONSTANT flag that allows you to transfer variables such as the value of the timestep that is set only once during the
simulation (i.e. during startup).

Data is transferred from the input void pointers to the PLMD value when the `share` and `wait` methods are called. Vectors, e.g. positions, that are split between the domains
are transferred when the share and wait methods of the [DOMAIN_DECOMPOSITION](DOMAIN_DECOMPOSITION.md) action are called.

You would only use the PUT command if you were calling PLUMED from python or an MD code. The equivalent commands in a convetional PLUMED input file would look like this.

```plumed
# This is how you create a value to hold the energy the MD code passes energy in plumed
eng: PUT UNIT=energy PERIODIC=NO
# This is how you create an vector of the 100 x positions to plumed
# Notice how we use ROLE here in order to tell PLUMED that these are the x coordinates.
# Further note that we need to use the FROM_DOMAINS flag if the MD code we are using uses
# domain decomposition.
xpos: PUT SHAPE=100 UNIT=length PERIODIC=NO ROLE=x FROM_DOMAINS
# This is how you create a scalar to hold the timestep
# The constant flag indicates that the value of the timestep doesn't change during the simulation
tstep: PUT CONSTANT UNIT=time PERIODIC=NO
# This is how you create a value to hold a 10 x 10 matrix in plumed whose elements are unitless
matrix: PUT SHAPE=10,10 UNIT=number PERIODIC=NO
# Lastly, if you want to pass a value that has a periodic domain you can do so as follows
# By adding the MUTABLE flag here we pass the data to PLUMED in a way that ensures that PLUMED can modify
# the value that was passed from the MD code and thus pass back a different value to the underlying MD code.
tor: PUT UNIT=number PERIODIC=-pi,pi MUTABLE
```

*/
//+ENDPLUMEDOC

namespace PLMD {

PLUMED_REGISTER_ACTION(ActionToPutData,"PUT")

void ActionToPutData::registerKeywords(Keywords& keys) {
  ActionForInterface::registerKeywords( keys );
  keys.add("compulsory","SHAPE","0","the shape of the value that is being passed to PLUMED");
  keys.add("compulsory","UNIT","the unit of the quantity that is being passed to PLUMED through this value.  Can be either number, energy, time, length, mass or charge");
  keys.add("compulsory","FORCE_UNIT","default","the units to use for the force");
  keys.add("compulsory","PERIODIC","if the value being passed to plumed is periodic then you should specify the periodicity of the function.  If the value "
           "is not periodic you must state this using PERIODIC=NO.  Positions are passed with PERIODIC=NO even though special methods are used "
           "to deal with pbc");
  keys.addFlag("CONSTANT",false,"does this quantity not depend on time");
  keys.addFlag("FROM_DOMAINS",false,"is this quantity passed through the domain decomposition object");
  keys.addFlag("MUTABLE",false,"can plumed change the value of the pointer that is passed from the MD code");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.setValueDescription("scalar/vector/matrix/grid","the data that was passed from the MD code");
}

ActionToPutData::ActionToPutData(const ActionOptions&ao):
  Action(ao),
  ActionForInterface(ao),
  noforce(false),
  fixed(false),
  from_domains(false),
  resetable(false),
  dataCanBeSet(true),
  unit(n),
  mydata(DataPassingObject::create(plumed.getRealPrecision())) {
  if( getName()!="ENERGY" && getName()!="PBC" ) {
    std::vector<std::size_t> shape;
    parseVector("SHAPE",shape);
    if( shape.size()==1 && shape[0]==0 ) {
      shape.resize(0);
      addValue( shape );
    } else {
      addValue( shape );
    }

    std::string unitstr, funitstr;
    parse("UNIT",unitstr);
    parse("FORCE_UNIT",funitstr);
    setUnit( unitstr, funitstr );

    // Now sort out period
    std::vector<std::string> period;
    parseVector("PERIODIC",period);
    if( period.size()==1 ) {
      if( period[0]!="NO") {
        error("input to PERIODIC keyword does not make sense");
      }
      setNotPeriodic();
    } else if( period.size()==2 ) {
      setPeriodic( period[0], period[1] );
    } else {
      error("input to PERIODIC keyword does not make sense");
    }

    parseFlag("CONSTANT",fixed);
    if( fixed ) {
      noforce=true;
      copyOutput(0)->setConstant();
    }
    parseFlag("FROM_DOMAINS",from_domains);
    parseFlag("MUTABLE",resetable);
  }
}

void ActionToPutData::setUnit( const std::string& unitstr, const std::string& funitstr ) {
  if( unitstr=="number" ) {
    unit=n;
  } else if( unitstr=="energy" ) {
    unit=e;
  } else if( unitstr=="length" ) {
    unit=l;
  } else if( unitstr=="mass" ) {
    unit=m;
  } else if( unitstr=="charge" ) {
    unit=q;
  } else if( unitstr=="time" ) {
    unit=t;
  } else {
    error( unitstr + " is not a valid input unit");
  }
  // Set the force units
  if( funitstr=="default" ) {
    funit=d;
  } else if( funitstr=="energy" ) {
    funit=eng;
  } else {
    error( funitstr + " is not a valid input force unit");
  }
}

std::string ActionToPutData::getUnitName() const {
  if( unit==e ) {
    return "energy";
  }
  if( unit==l ) {
    return "length";
  }
  if( unit==m ) {
    return "mass";
  }
  if( unit==q ) {
    return "charge";
  }
  if( unit==t ) {
    return "time";
  }
  plumed_error();
}

void ActionToPutData::setStart( const std::string& actname, const unsigned& sss) {
  plumed_assert( actname==getLabel() );
  mydata->setStart(sss);
}

void ActionToPutData::setStride( const std::string& actname, const unsigned& sss ) {
  plumed_assert( actname==getLabel() );
  mydata->setStride(sss);
}

void ActionToPutData::updateUnits( DataPassingTools* passtools ) {
  // Don't need to do anythign if this is just a number
  if( unit==n ) {
    return ;
  }

  double vunits=passtools->getUnitConversion( getUnitName() );
  mydata->setUnit(vunits);
  if( fixed && wasset ) {
    mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
  }
  if( funit==eng ) {
    mydata->setForceUnit( 1/passtools->getUnitConversion("energy"));
  } else if( funit==d ) {
    mydata->setForceUnit(1/passtools->getUnitConversion("energy")*vunits);
  }
}

bool ActionToPutData::setValuePointer( const std::string& actname, const TypesafePtr & val ) {
  if( actname!=getLabel() ) {
    return false;
  }
  wasset=true;
  plumed_massert( dataCanBeSet, "set " + getLabel() + " cannot be set at this time");
  if( !from_domains ) {
    if( !resetable && getPntrToComponent(0)->getRank()==0 ) {
      mydata->saveValueAsDouble( val );
      if( fixed ) {
        mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
      }
    } else {
      mydata->setValuePointer(val,getPntrToComponent(0)->getShape(), !resetable);
    }
  } else {
    mydata->setValuePointer(val,std::vector<std::size_t>(), !resetable);
  }
  return true;
}

bool ActionToPutData::setForcePointer( const std::string& actname, const TypesafePtr & val ) {
  if( actname!=getLabel() ) {
    return false;
  }
  plumed_massert( dataCanBeSet, "force on " + getLabel() + " cannot be set at this time");
  if( !from_domains ) {
    mydata->setForcePointer(val,getPntrToComponent(0)->getShape());
  } else {
    mydata->setForcePointer(val,std::vector<std::size_t>());
  }
  return true;
}

void ActionToPutData::getLocalValues( std::vector<double>& vals ) const {
  mydata->share_data( vals );
}

void ActionToPutData::wait() {
  dataCanBeSet=false;
  if( fixed || !wasset ) {
    return;
  }
  plumed_assert( wasset );
  mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
}

void ActionToPutData::apply() {
  if( getPntrToValue()->forcesWereAdded() && !noforce ) {
    if( getName()=="ENERGY" || getDependencies().size()==0 ) {
      mydata->add_force( getPntrToValue() );
    }
  }
}

unsigned ActionToPutData::getNumberOfForcesToRescale() const {
  if( getName()!="ENERGY" || getDependencies().size()>0 ) {
    return copyOutput(0)->getNumberOfValues();
  }
  plumed_assert( getDependencies().size()==1 );
  plumed_assert(getDependencies()[0]); // needed for following calls, see #1046
  ActionForInterface* ai = getDependencies()[0]->castToActionForInterface();
  return ai->getNumberOfForcesToRescale();
}

void ActionToPutData::rescaleForces( const double& alpha ) {
  if( noforce ) {
    return;
  }
  wasscaled=true;
  mydata->rescale_force( getNumberOfForcesToRescale(), alpha, getPntrToValue() );

}

void ActionToPutData::writeBinary(std::ostream&o) {
  if(!fixed) {
    getPntrToValue()->writeBinary(o);
  }
}

void ActionToPutData::readBinary(std::istream&i) {
  if(!fixed) {
    getPntrToValue()->readBinary(i);
  }
}

}
