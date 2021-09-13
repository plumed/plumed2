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
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","DATA","the label of an action that calculates multicolvars.  Weighted sums based on the location of the colvars calculated by this action will be calcualted");
  keys.add("optional","ATOMS","the atoms to use for the positions of the multicolvars");
  keys.add("optional","LESS_THAN","calcualte the number of colvars that are inside the region of interest and that are less than a certain threshold");
  keys.addOutputComponent("lessthan","LESS_THAN","the number of cvs in the region of interest that are less than a certain threshold");
  keys.add("optional","MORE_THAN","calcualte the number of colvars that are inside the region of interest and that are greater that a certain threshold");
  keys.addOutputComponent("morethan","MORE_THAN","the number of cvs in the region of interest that are more than a certain threshold");
  keys.add("optional","BETWEEN","calculate the number of colvars that are inside the region of interest and that have a CV value that is between a particular set of bounds");
  keys.addOutputComponent("between","BETWEEN","the number of cvs in the region of interest that are within a certain range");
  keys.addFlag("SUM",false,"calculate the sum of all the quantities.");
  keys.addOutputComponent("sum","SUM","the sum of all the colvars weighted by the function that determines if we are in the region");
  keys.addFlag("MEAN",false,"calculate the average value of the colvar inside the region of interest");
  keys.addOutputComponent("mean","MEAN","the average values of the colvar in the region of interest");
}

template <const char* v>
VolumeShortcut<v>::VolumeShortcut(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao) 
{
  std::string voltype(v), mc_lab; parse("DATA",mc_lab); bool dosum; parseFlag("SUM",dosum);
  if( mc_lab.length()>0 ) {
    bool domean; parseFlag("MEAN",domean); std::string lt_input, mt_input, bt_input; 
    parse("LESS_THAN",lt_input); parse("MORE_THAN",mt_input); parse("BETWEEN",bt_input); 
    std::string atomsd; parse("ATOMS",atomsd); if( atomsd.length()==0 ) atomsd=mc_lab;
    // Create the apprpriate volume object
    readInputLine( getShortcutLabel() + ": " + voltype + "_CALC " + convertInputLineToString() + " ATOMS=" + atomsd );
    // Now create input for sums
    if( dosum || domean ) {
      readInputLine( getShortcutLabel() + "_prod: MATHEVAL ARG1=" + mc_lab + " ARG2=" + getShortcutLabel() + " FUNC=x*y PERIODIC=NO"); 
      std::string tlab = getShortcutLabel() + "_numer"; if( dosum ) tlab = getShortcutLabel() + "_sum:";
      readInputLine( tlab + ": SUM ARG=" + getShortcutLabel() + "_prod PERIODIC=NO"); 
    }
    if( domean ) {
      // Calculate denominator
      readInputLine( getShortcutLabel() + "_norm: SUM ARG=" + getShortcutLabel() + " PERIODIC=NO"); 
      // And calculate final quantity which is mean of these two actions
      std::string arg1_lab = getShortcutLabel() + "_numer"; if( dosum ) arg1_lab = getShortcutLabel()  + "_sum"; 
      readInputLine( getShortcutLabel() + "_mean: MATHEVAL ARG1=" + arg1_lab + " ARG2=" + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO"); 
    }
    if( lt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_lt: LESS_THAN ARG1=" + mc_lab + " SWITCH={" + lt_input +"}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_lt: MATHEVAL ARG1=" + mc_lab + "_" + getShortcutLabel() + "_lt ARG2=" + getShortcutLabel() + " FUNC=x*y PERIODIC=NO"); 
      // And the final sum
      readInputLine( getShortcutLabel() + "_lessthan: SUM ARG=" + getShortcutLabel() + "_lt PERIODIC=NO"); 
    }
    if( mt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_mt: MORE_THAN ARG1=" + mc_lab + " SWITCH={" + mt_input  + "}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_mt: MATHEVAL ARG1=" + mc_lab + "_" + getShortcutLabel() + "_mt ARG2=" + getShortcutLabel() + " FUNC=x*y PERIODIC=NO"); 
      // And the final sum
      readInputLine( getShortcutLabel() + "_morethan: SUM ARG=" + getShortcutLabel() + "_mt PERIODIC=NO"); 
    }
    if( bt_input.length()>0 ) {
      // Calculate number less than
      readInputLine( mc_lab + "_" + getShortcutLabel() + "_bt: LESS_THAN ARG1=" + mc_lab + " SWITCH={" + bt_input +"}");
      // And the matheval bit
      readInputLine( getShortcutLabel() + "_bt: MATHEVAL ARG1=" + mc_lab + "_" + getShortcutLabel() + "_bt ARG2=" + getShortcutLabel() + " FUNC=x*y PERIODIC=NO"); 
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
