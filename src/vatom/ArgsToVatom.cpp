/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/PbcAction.h"
#include "tools/Pbc.h"

//+PLUMEDOC VATOM ARGS2VATOM
/*
Create a virtual atom from the input scalars

This action takes five scalars that are computed by other actions in input and uses them to set the
x, y and z positions and the mass and charge of a virtual atom.  This action is used within the
[CENTER](CENTER.md) shortcut to compute a center of mass.  An example input that shows how you
can use this command to calculate the center of mass of atoms 1-10 is as follows:

```plumed
# Calculate the total mass of the atoms
m: MASS ATOMS=1-10
mass: SUM ARG=m PERIODIC=NO
# Calculate the totla charge of the atoms
q: CHARGE ATOMS=1-10
charge: SUM ARG=q PERIODIC=NO
# Now get the positions of the atoms
pos: POSITION WHOLEMOLECULES ATOMS=1-10
# Multiply each vector of positions by the masses
xwvec: CUSTOM ARG=m,pos.x FUNC=x*y PERIODIC=NO
ywvec: CUSTOM ARG=m,pos.y FUNC=x*y PERIODIC=NO
zwvec: CUSTOM ARG=m,pos.z FUNC=x*y PERIODIC=NO
# Sum the numerators in the expression for the center of mass
xnum: SUM ARG=xwvec PERIODIC=NO
ynum: SUM ARG=ywvec PERIODIC=NO
znum: SUM ARG=zwvec PERIODIC=NO
# And compute the x, y and z positions of the center of mass
x: CUSTOM ARG=xnum,mass FUNC=x/y PERIODIC=NO
y: CUSTOM ARG=ynum,mass FUNC=x/y PERIODIC=NO
z: CUSTOM ARG=znum,mass FUNC=x/y PERIODIC=NO
# And now create the virtual atom
p: ARGS2VATOM XPOS=x YPOS=y ZPOS=z MASS=mass CHARGE=charge
```

In the following input by contrast we use the PHASES method that is discussed in the documentation for [CENTER](CENTER.md) to calculate the position of the center of mass:

```plumed
# Calculate the total mass of the atoms
m: MASS ATOMS=1-10
mass: SUM ARG=m PERIODIC=NO
# Calculate the totla charge of the atoms
q: CHARGE ATOMS=1-10
charge: SUM ARG=q PERIODIC=NO
# Now get the positions of the atoms
pos: POSITION SCALED_COMPONENTS ATOMS=1-10
# Multiply the sins and cosines of the scaled positions by the weights
sina: CUSTOM ARG=m,pos.a FUNC=x*sin(2*pi*y) PERIODIC=NO
cosa: CUSTOM ARG=m,pos.a FUNC=x*cos(2*pi*y) PERIODIC=NO
sinb: CUSTOM ARG=m,pos.b FUNC=x*sin(2*pi*y) PERIODIC=NO
cosb: CUSTOM ARG=m,pos.b FUNC=x*cos(2*pi*y) PERIODIC=NO
sinc: CUSTOM ARG=m,pos.c FUNC=x*sin(2*pi*y) PERIODIC=NO
cosc: CUSTOM ARG=m,pos.c FUNC=x*cos(2*pi*y) PERIODIC=NO
# And accumulate the sums
sinsuma: SUM ARG=sina PERIODIC=NO
cossuma: SUM ARG=cosa PERIODIC=NO
sinsumb: SUM ARG=sinb PERIODIC=NO
cossumb: SUM ARG=cosb PERIODIC=NO
sinsumc: SUM ARG=sinc PERIODIC=NO
cossumc: SUM ARG=cosc PERIODIC=NO
# And get the position of the center in fractional coordinates
a: CUSTOM ARG=sinsuma,cossuma FUNC=atan2(x,y)/(2*pi) PERIODIC=NO
b: CUSTOM ARG=sinsumb,cossumb FUNC=atan2(x,y)/(2*pi) PERIODIC=NO
c: CUSTOM ARG=sinsumc,cossumc FUNC=atan2(x,y)/(2*pi) PERIODIC=NO
# And now create the virtual atom
p: ARGS2VATOM XPOS=a YPOS=b ZPOS=c MASS=mass CHARGE=charge FRACTIONAL
```

These inputs provide very slow ways of computing a center of mass - PLUMED contains faster implementations that do calculations that are equivalent to both of these inputs.
This type of input is nevertheless useful if you are using arbitary weights when computing the sums in the numerator
and denominator of the expression for the center as is detailed in the documentation for the [CENTER](CENTER.md) command.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace vatom {

class ArgsToVatom :
  public ActionWithValue,
  public ActionWithArguments {
private:
  bool fractional;
  PbcAction* pbc_action;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ArgsToVatom(const ActionOptions&);
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() override {
    return getNumberOfArguments();
  }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(ArgsToVatom,"ARGS2VATOM")

void ArgsToVatom::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","XPOS","scalar","the value to use for the x position of the atom");
  keys.addInputKeyword("compulsory","YPOS","scalar","the value to use for the y position of the atom");
  keys.addInputKeyword("compulsory","ZPOS","scalar","the value to use for the z position of the atom");
  keys.addInputKeyword("compulsory","MASS","scalar","the value to use for the mass of the atom");
  keys.addInputKeyword("compulsory","CHARGE","scalar","the value to use for the charge of the atom");
  keys.addInputKeyword("hidden","XBKP","scalar","x position to use in case PBC not set when using PHASES");
  keys.addInputKeyword("hidden","YBKP","scalar","y position to use in case PBC not set when using PHASES");
  keys.addInputKeyword("hidden","ZBKP","scalar","z position to use in case PBC not set when using PHASES");
  keys.addFlag("FRACTIONAL",false,"the input arguments are calculated in fractional coordinates so you need to multiply by the cell");
  keys.addOutputComponent("x","default","scalar","the x coordinate of the virtual atom");
  keys.addOutputComponent("y","default","scalar","the y coordinate of the virtual atom");
  keys.addOutputComponent("z","default","scalar","the z coordinate of the virtual atom");
  keys.addOutputComponent("mass","default","scalar","the mass of the virtual atom");
  keys.addOutputComponent("charge","default","scalar","the charge of the virtual atom");
}

ArgsToVatom::ArgsToVatom(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  parseFlag("FRACTIONAL",fractional);
  std::vector<Value*> xpos;
  parseArgumentList("XPOS",xpos);
  if( xpos.size()!=1 && xpos[0]->getRank()!=0 ) {
    error("invalid input argument for XPOS");
  }
  std::vector<Value*> ypos;
  parseArgumentList("YPOS",ypos);
  if( ypos.size()!=1 && ypos[0]->getRank()!=0 ) {
    error("invalid input argument for YPOS");
  }
  std::vector<Value*> zpos;
  parseArgumentList("ZPOS",zpos);
  if( zpos.size()!=1 && zpos[0]->getRank()!=0 ) {
    error("invalid input argument for ZPOS");
  }
  std::vector<Value*> mass;
  parseArgumentList("MASS",mass);
  if( mass.size()!=1 && mass[0]->getRank()!=0 ) {
    error("invalid input argument for MASS");
  }
  std::vector<Value*> charge;
  parseArgumentList("CHARGE",charge);
  if( charge.size()!=1 && charge[0]->getRank()!=0 ) {
    error("invalid input argument for CHARGE");
  }
  // Make sure we have requested everything that we need in xpos
  xpos.push_back(ypos[0]);
  xpos.push_back(zpos[0]);
  xpos.push_back(mass[0]);
  xpos.push_back(charge[0]);
  if( fractional ) {
    log.printf("  creating atom from fractional pos a=%s, b=%s and c=%s \n", xpos[0]->getName().c_str(), ypos[0]->getName().c_str(), zpos[0]->getName().c_str() );
    std::vector<Value*> xbkp;
    parseArgumentList("XBKP",xbkp);
    if( xbkp.size()>0 ) {
      if( xbkp.size()!=1 && xbkp[0]->getRank()!=0 ) {
        error("invalid input argument for XBKP");
      }
      std::vector<Value*> ybkp;
      parseArgumentList("YBKP",ybkp);
      if( ybkp.size()!=1 && ybkp[0]->getRank()!=0 ) {
        error("invalid input argument for YBKP");
      }
      std::vector<Value*> zbkp;
      parseArgumentList("ZBKP",zbkp);
      if( zbkp.size()!=1 && zpos[0]->getRank()!=0 ) {
        error("invalid input argument for ZBKP");
      }
      // Store backup for NOPBC
      xpos.push_back(xbkp[0]);
      xpos.push_back(ybkp[0]);
      xpos.push_back(zbkp[0]);
      log.printf("  using x=%s, y=%s and z=%s if PBC not set \n", xbkp[0]->getName().c_str(), ybkp[0]->getName().c_str(), zbkp[0]->getName().c_str() );
    }
  } else {
    log.printf("  creating atom at x=%s, y=%s and z=%s \n", xpos[0]->getName().c_str(), ypos[0]->getName().c_str(), zpos[0]->getName().c_str() );
  }
  log.printf("  mass of atom is %s and charge is %s \n", mass[0]->getName().c_str(), charge[0]->getName().c_str() );
  // Request the arguments
  requestArguments(xpos);
  // Create the components to hold the atom
  addComponentWithDerivatives("x");
  componentIsNotPeriodic("x");
  addComponentWithDerivatives("y");
  componentIsNotPeriodic("y");
  addComponentWithDerivatives("z");
  componentIsNotPeriodic("z");
  addComponent("mass");
  componentIsNotPeriodic("mass");
  if( mass[0]->isConstant() ) {
    getPntrToComponent(3)->setConstant();
  }
  addComponent("charge");
  componentIsNotPeriodic("charge");
  if( charge[0]->isConstant() ) {
    getPntrToComponent(4)->setConstant();
  }
  pbc_action = plumed.getActionSet().selectWithLabel<PbcAction*>("Box");
  for(unsigned i=0; i<3; ++i) {
    getPntrToComponent(i)->resizeDerivatives( getNumberOfArguments() );
  }
}

void ArgsToVatom::calculate() {
  if( fractional ) {
    if( pbc_action->getPbc().isSet() ) {
      // Get the position in fractional coordinates
      Vector fpos;
      for(unsigned i=0; i<3; ++i) {
        fpos[i] = getPntrToArgument(i)->get();
      }
      // Convert fractioanl coordinates to cartesian coordinates
      Tensor box=pbc_action->getPbc().getBox();
      Vector cpos=matmul(fpos,box);
      // Set the final position and derivatives
      for(unsigned i=0; i<3; ++i) {
        Value* vv=getPntrToComponent(i);
        vv->set( cpos[i] );
        for(unsigned j=0; j<3; ++j) {
          vv->addDerivative( j, box[j][i] );
        }
      }
    } else {
      if( getNumberOfArguments()<8 ) {
        error("cannot use PHASES option if box is not set");
      }
      // Set the values
      for(unsigned i=0; i<3; ++i) {
        getPntrToComponent(i)->set( getPntrToArgument(5+i)->get() );
      }
      // And the derivatives
      for(unsigned i=0; i<3; ++i) {
        getPntrToComponent(i)->addDerivative( 5+i, 1.0 );
      }
    }
    // Set the mass and charge
    for(unsigned i=3; i<5; ++i) {
      getPntrToComponent(i)->set( getPntrToArgument(i)->get() );
    }
  } else {
    // Set the values
    for(unsigned i=0; i<5; ++i) {
      getPntrToComponent(i)->set( getPntrToArgument(i)->get() );
    }
    // And the derivatives
    for(unsigned i=0; i<3; ++i) {
      getPntrToComponent(i)->addDerivative( i, 1.0 );
    }
  }
}

void ArgsToVatom::apply() {
  if( !checkForForces() ) {
    return ;
  }

  unsigned start=0;
  addForcesOnArguments( 0, getForcesToApply(), start );
}

}
}
