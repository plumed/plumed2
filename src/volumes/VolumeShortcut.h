/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_volumes_VolumeShortcut_h
#define __PLUMED_volumes_VolumeShortcut_h

#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Group.h"

namespace PLMD {
namespace volumes {

template <const char* v>
class VolumeShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  VolumeShortcut(const ActionOptions&ao);
};

template <const char* v>
void VolumeShortcut<v>::registerKeywords( Keywords& keys ) {
  actionRegister().getKeywords( std::string(v) + "_CALC", keys );
  keys.add("hidden","IS_SHORTCUT","hidden keyword to tell if actions are shortcuts so that example generator can provide expansions of shortcuts");
  keys.addDeprecatedKeyword("DATA","");
  keys.addDeprecatedKeyword("LESS_THAN","");
  keys.addOutputComponent("lessthan","LESS_THAN","scalar","the number of cvs in the region of interest that are less than a certain threshold");
  keys.addDeprecatedKeyword("MORE_THAN","");
  keys.addOutputComponent("morethan","MORE_THAN","scalar","the number of cvs in the region of interest that are more than a certain threshold");
  keys.addDeprecatedKeyword("BETWEEN","");
  keys.addOutputComponent("between","BETWEEN","scalar","the number of cvs in the region of interest that are within a certain range");
  keys.addDeprecatedFlag("SUM","");
  keys.addOutputComponent("sum","SUM","scalar","the sum of all the colvars weighted by the function that determines if we are in the region");
  keys.addDeprecatedFlag("MEAN","");
  keys.addOutputComponent("mean","MEAN","scalar","the average values of the colvar in the region of interest");
  keys.addActionNameSuffix("_CALC");
  keys.needsAction("LESS_THAN");
  keys.needsAction("MORE_THAN");
  keys.needsAction("GROUP");
  keys.needsAction("BETWEEN");
  keys.needsAction("SUM");
  keys.needsAction("MEAN");
  keys.needsAction("CUSTOM");
  keys.setValueDescription("scalar","sum of values of input CVs in regin of interest");
}

template <const char* v>
VolumeShortcut<v>::VolumeShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string voltype(v), mc_lab;
  parse("DATA",mc_lab);
  bool dosum;
  parseFlag("SUM",dosum);
  if( mc_lab.length()>0 ) {
    Group* mygrp = plumed.getActionSet().template selectWithLabel<Group*>(mc_lab);
    Group* mygrp2 = plumed.getActionSet().template selectWithLabel<Group*>(mc_lab + "_grp");
    if( mygrp || mygrp2 ) {
      readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + mc_lab );
    }
    bool domean;
    parseFlag("MEAN",domean);
    std::string lt_input, mt_input, bt_input;
    parse("LESS_THAN",lt_input);
    parse("MORE_THAN",mt_input);
    parse("BETWEEN",bt_input);
    std::string atomsd;
    parse("ATOMS",atomsd);
    if( atomsd.length()==0 ) {
      atomsd=mc_lab;
    }
    // Create the apprpriate volume object
    readInputLine( getShortcutLabel() + ": " + voltype + "_CALC " + convertInputLineToString() + " ATOMS=" + atomsd );
    // Now create input for sums
    if( dosum || domean ) {
      readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + mc_lab + "," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
      std::string tlab = getShortcutLabel() + "_numer";
      if( dosum ) {
        tlab = getShortcutLabel() + "_sum";
      }
      readInputLine( tlab + ": SUM ARG=" + getShortcutLabel() + "_prod PERIODIC=NO");
    }
    if( domean ) {
      // Calculate denominator
      readInputLine( getShortcutLabel() + "_norm: SUM ARG=" + getShortcutLabel() + " PERIODIC=NO");
      // And calculate final quantity which is mean of these two actions
      std::string arg1_lab = getShortcutLabel() + "_numer";
      if( dosum ) {
        arg1_lab = getShortcutLabel()  + "_sum";
      }
      readInputLine( getShortcutLabel() + "_mean: CUSTOM ARG=" + arg1_lab + "," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
    }
    if( lt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_lt: LESS_THAN ARG=" + mc_lab + " SWITCH={" + lt_input +"}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_lt: CUSTOM ARG=" + mc_lab + "_" + getShortcutLabel() + "_lt," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
      // And the final sum
      readInputLine( getShortcutLabel() + "_lessthan: SUM ARG=" + getShortcutLabel() + "_lt PERIODIC=NO");
    }
    if( mt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_mt: MORE_THAN ARG=" + mc_lab + " SWITCH={" + mt_input  + "}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_mt: CUSTOM ARG=" + mc_lab + "_" + getShortcutLabel() + "_mt," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
      // And the final sum
      readInputLine( getShortcutLabel() + "_morethan: SUM ARG=" + getShortcutLabel() + "_mt PERIODIC=NO");
    }
    if( bt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_bt: BETWEEN ARG=" + mc_lab + " SWITCH={" + bt_input +"}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_bt: CUSTOM ARG=" + mc_lab + "_" + getShortcutLabel() + "_bt," + getShortcutLabel() + " FUNC=x*y PERIODIC=NO");
      // And the final sum
      readInputLine( getShortcutLabel() + "_between: SUM ARG=" + getShortcutLabel() + "_bt PERIODIC=NO");
    }
  } else if( dosum ) {
    readInputLine( getShortcutLabel() + "_vols: " + voltype + "_CALC " + convertInputLineToString() );
    readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_vols PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + ": " + voltype + "_CALC " + convertInputLineToString() );
  }
}

}
}
#endif
