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

\par Examples

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
  unsigned getNumberOfDerivatives() override { return getNumberOfArguments(); }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
};

PLUMED_REGISTER_ACTION(ArgsToVatom,"ARGS2VATOM")

void ArgsToVatom::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys );
  keys.add("compulsory","XPOS","the x position of the atom");
  keys.add("compulsory","YPOS","the y position of the atom");
  keys.add("compulsory","ZPOS","the z position of the atom");
  keys.add("compulsory","MASS","the mass of the atom");
  keys.add("compulsory","CHARGE","the charge of the atom");
  keys.add("hidden","XBKP","x position to use in case PBC not set when using PHASES");
  keys.add("hidden","YBKP","y position to use in case PBC not set when using PHASES");
  keys.add("hidden","ZBKP","z position to use in case PBC not set when using PHASES");
  keys.addFlag("FRACTIONAL",false,"the input arguments are calculated in fractional coordinates so you need to multiply by the cell");
  keys.addOutputComponent("x","default","the x coordinate of the virtual atom");
  keys.addOutputComponent("y","default","the y coordinate of the virtual atom");
  keys.addOutputComponent("z","default","the z coordinate of the virtual atom");
  keys.addOutputComponent("mass","default","the mass of the virtual atom");
  keys.addOutputComponent("charge","default","the charge of the virtual atom");
}

ArgsToVatom::ArgsToVatom(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  parseFlag("FRACTIONAL",fractional);
  std::vector<Value*> xpos; parseArgumentList("XPOS",xpos);
  if( xpos.size()!=1 && xpos[0]->getRank()!=0 ) error("invalid input argument for XPOS");
  std::vector<Value*> ypos; parseArgumentList("YPOS",ypos);
  if( ypos.size()!=1 && ypos[0]->getRank()!=0 ) error("invalid input argument for YPOS");
  std::vector<Value*> zpos; parseArgumentList("ZPOS",zpos);
  if( zpos.size()!=1 && zpos[0]->getRank()!=0 ) error("invalid input argument for ZPOS");
  std::vector<Value*> mass; parseArgumentList("MASS",mass);
  if( mass.size()!=1 && mass[0]->getRank()!=0 ) error("invalid input argument for MASS");
  std::vector<Value*> charge; parseArgumentList("CHARGE",charge);
  if( charge.size()!=1 && charge[0]->getRank()!=0 ) error("invalid input argument for CHARGE");
  // Make sure we have requested everything that we need in xpos
  xpos.push_back(ypos[0]); xpos.push_back(zpos[0]); xpos.push_back(mass[0]); xpos.push_back(charge[0]);
  if( fractional ) {
    log.printf("  creating atom from fractional pos a=%s, b=%s and c=%s \n", xpos[0]->getName().c_str(), ypos[0]->getName().c_str(), zpos[0]->getName().c_str() );
    std::vector<Value*> xbkp; parseArgumentList("XBKP",xbkp);
    if( xbkp.size()>0 ) {
      if( xbkp.size()!=1 && xbkp[0]->getRank()!=0 ) error("invalid input argument for XBKP");
      std::vector<Value*> ybkp; parseArgumentList("YBKP",ybkp);
      if( ybkp.size()!=1 && ybkp[0]->getRank()!=0 ) error("invalid input argument for YBKP");
      std::vector<Value*> zbkp; parseArgumentList("ZBKP",zbkp);
      if( zbkp.size()!=1 && zpos[0]->getRank()!=0 ) error("invalid input argument for ZBKP");
      // Store backup for NOPBC
      xpos.push_back(xbkp[0]); xpos.push_back(ybkp[0]); xpos.push_back(zbkp[0]);
      log.printf("  using x=%s, y=%s and z=%s if PBC not set \n", xbkp[0]->getName().c_str(), ybkp[0]->getName().c_str(), zbkp[0]->getName().c_str() );
    }
  } else log.printf("  creating atom at x=%s, y=%s and z=%s \n", xpos[0]->getName().c_str(), ypos[0]->getName().c_str(), zpos[0]->getName().c_str() );
  log.printf("  mass of atom is %s and charge is %s \n", mass[0]->getName().c_str(), charge[0]->getName().c_str() );
  // Request the arguments
  requestArguments(xpos);
  // Create the components to hold the atom
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  addComponent("mass"); componentIsNotPeriodic("mass"); if( mass[0]->isConstant() ) getPntrToComponent(3)->setConstant();
  addComponent("charge"); componentIsNotPeriodic("charge"); if( charge[0]->isConstant() ) getPntrToComponent(4)->setConstant();
  pbc_action = plumed.getActionSet().selectWithLabel<PbcAction*>("Box");
  for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->resizeDerivatives( getNumberOfArguments() );
}

void ArgsToVatom::calculate() {
  if( fractional ) {
    if( pbc_action->getPbc().isSet() ) {
      // Get the position in fractional coordinates
      Vector fpos; for(unsigned i=0; i<3; ++i) fpos[i] = getPntrToArgument(i)->get();
      // Convert fractioanl coordinates to cartesian coordinates
      Tensor box=pbc_action->getPbc().getBox(); Vector cpos=matmul(fpos,box);
      // Set the final position and derivatives
      for(unsigned i=0; i<3; ++i) {
        Value* vv=getPntrToComponent(i); vv->set( cpos[i] );
        for(unsigned j=0; j<3; ++j) vv->addDerivative( j, box[j][i] );
      }
    } else {
      if( getNumberOfArguments()<8 ) error("cannot use PHASES option if box is not set");
      // Set the values
      for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->set( getPntrToArgument(5+i)->get() );
      // And the derivatives
      for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->addDerivative( 5+i, 1.0 );
    }
    // Set the mass and charge
    for(unsigned i=3; i<5; ++i) getPntrToComponent(i)->set( getPntrToArgument(i)->get() );
  } else {
    // Set the values
    for(unsigned i=0; i<5; ++i) getPntrToComponent(i)->set( getPntrToArgument(i)->get() );
    // And the derivatives
    for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->addDerivative( i, 1.0 );
  }
}

void ArgsToVatom::apply() {
  if( !checkForForces() ) return ;

  unsigned start=0; addForcesOnArguments( 0, getForcesToApply(), start, getLabel() );
}

}
}
