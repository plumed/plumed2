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
#include "setup/SetupReferenceBase.h"

namespace PLMD {
namespace function {

class Displacement : public ActionShortcut {
public:
  static std::string fixArgumentDot( const std::string& argin );
  static void registerKeywords( Keywords& keys );
  explicit Displacement(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Displacement,"DISPLACEMENT")

void Displacement::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The point that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
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
  std::vector<std::string> arg1, arg2, arg1f, arg2f; parseVector("ARG1",arg1); parseVector("ARG2",arg2);
  // Check if one of the input arguments is a reference cluster
  bool arg1iscenter=false;
  if( arg1.size()==1 && arg1[0].find("center")!=std::string::npos ) {
      std::size_t dot=arg1[0].find("."); setup::SetupReferenceBase* ss=plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( arg1[0].substr(0,dot) );
      if( ss ) {
          arg1iscenter=true; 
          if( ss->getNumberOfArguments()!=arg2.size() ) error("mismatch between number of arguments in input and number of arguments in reference");
          for(unsigned i=0;i<ss->getNumberOfArguments();++i) {
              if( ss->getPntrToArgument(i)->getName()!=arg2[i] ) error("mismatch between arguments in input and refrence cluster");
              std::string num; Tools::convert( i+1, num ); arg1f.push_back( getShortcutLabel() + "_" + fixArgumentDot(arg2[i]) + "_ref" );
              readInputLine( arg1f[i] + ": SELECT_COMPONENTS ARG=" + arg1[0] + " COMPONENTS=" + num );
          }
      }
  } else { for(unsigned i=0;i<arg1.size();++i) arg1f.push_back( arg1[i] ); }
 
  if( arg2.size()==1 && arg2[0].find("center")!=std::string::npos ) {
      std::size_t dot=arg2[0].find("."); setup::SetupReferenceBase* ss=plumed.getActionSet().selectWithLabel<setup::SetupReferenceBase*>( arg2[0].substr(0,dot) );
      if( ss ) {
          if( arg1iscenter ) error("both arguments are set constact during setup");
          if( ss->getNumberOfArguments()!=arg1.size() ) error("mismatch between number of arguments in input and number of arguments in reference");
          for(unsigned i=0;i<ss->getNumberOfArguments();++i) {
              if( ss->getPntrToArgument(i)->getName()!=arg1[i] ) error("mismatch between arguments in input and refrence cluster");
              std::string num; Tools::convert( i+1, num ); arg2f.push_back( getShortcutLabel() + "_" + fixArgumentDot(arg1[i]) + "_ref" );
              readInputLine( arg2f[i] + ": SELECT_COMPONENTS ARG=" + arg2[0] + " COMPONENTS=" + num );
          }
      }
  } else { for(unsigned i=0;i<arg2.size();++i) arg2f.push_back( arg2[i] ); }

  if( arg1f.size()==1 ) readInputLine( getShortcutLabel() + ": DIFFERENCE ARG1=" + arg1f[0] + " ARG2=" + arg2f[0] );
  else {
      for(unsigned i=0;i<arg1f.size();++i) readInputLine( getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff: DIFFERENCE ARG1=" + arg1f[i] + " ARG2=" + arg2f[i] );
      std::string argdat = "ARG=" + getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff"; for(unsigned i=1;i<arg1f.size();++i) argdat += "," +  getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff";
      readInputLine( getShortcutLabel() + ": VSTACK " + argdat );
  }
}

}
}
