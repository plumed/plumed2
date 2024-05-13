/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"

//+PLUMEDOC MCOLVAR DISPLACEMENT
/*
Calculate the displacement vector between the pair of input vectors

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class Displacement : public ActionShortcut {
public:
  static std::string fixArgumentDot( const std::string& argin );
  static void registerKeywords( Keywords& keys );
  Value* getValueWithLabel( const std::string& name ) const ;
  explicit Displacement(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Displacement,"DISPLACEMENT")

void Displacement::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The point that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.setValueDescription("the differences between the input arguments");
  keys.needsAction("DIFFERENCE"); keys.needsAction("TRANSPOSE"); keys.needsAction("VSTACK");
}

std::string Displacement::fixArgumentDot( const std::string& argin ) {
  std::string argout = argin; std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  return argout;
}

Displacement::Displacement( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Read in argument names
  std::vector<std::string> arg1f, arg2f; parseVector("ARG1",arg1f); parseVector("ARG2",arg2f);
  // Check if one of the input arguments is a reference cluster
  if( arg1f.size()!=arg2f.size() ) error("number of arguments specified to ARG1 should be same as number for ARG2");

  Value* val1=getValueWithLabel( arg1f[0] );
  if( arg1f.size()==1 && val1->getRank()!=0 ) {
    Value* val2=getValueWithLabel( arg2f[0] );
    if( val1->getNumberOfValues()==val2->getNumberOfValues() ) {
      readInputLine( getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff: DIFFERENCE ARG=" + arg1f[0] + "," + arg2f[0] );
      readInputLine( getShortcutLabel() + ": TRANSPOSE ARG=" + getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff");
    } else readInputLine( getShortcutLabel() + ": DIFFERENCE ARG=" + arg1f[0] + "," + arg2f[0] );
  } else {
    for(unsigned i=0; i<arg1f.size(); ++i) readInputLine( getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff: DIFFERENCE ARG=" + arg1f[i] + "," + arg2f[i] );
    std::string argdat = "ARG=" + getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff"; for(unsigned i=1; i<arg1f.size(); ++i) argdat += "," +  getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff";
    readInputLine( getShortcutLabel() + ": VSTACK " + argdat );
  }
}

Value* Displacement::getValueWithLabel( const std::string& name ) const {
  std::size_t dot=name.find("."); std::string sname = name.substr(0,dot);
  ActionWithValue* vv=plumed.getActionSet().selectWithLabel<ActionWithValue*>( sname );
  if( !vv ) error("cannot find value with name " + name );
  if( dot==std::string::npos ) return vv->copyOutput(0);
  if( !vv->exists(name) ) error("cannot find value with name " + name );
  return vv->copyOutput( name );
}

}
}
