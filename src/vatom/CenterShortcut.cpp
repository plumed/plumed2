/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/Group.h"
#include "core/ActionWithValue.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace vatom {

class CenterShortcut : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit CenterShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CenterShortcut,"CENTER")

void CenterShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","ATOMS","the group of atoms that appear in the definition of this center");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","WEIGHTS","what weights should be used when calculating the center.  If this keyword is not present the geometric center is computed. "
           "If WEIGHTS=@Masses is used the center of mass is computed.  If WEIGHTS=@charges the center of charge is computed.  If "
           "the label of an action is provided PLUMED assumes that that action calculates a list of symmetry functions that can be used "
           "as weights. Lastly, an explicit list of numbers to use as weights can be provided");
  keys.addFlag("PHASES",false,"use trigonometric phases when computing position of center");
  keys.addFlag("SAFE_PHASES",false,"use trignomentric phases when computing position of center but also compute the center in ths usual way and use this when the pbc are not set. "
               "There are two reasons for using this option (1) you are doing something that you know is really weird or (2) you are an idiot");
  keys.addFlag("MASS",false,"calculate the center of mass");
  keys.addActionNameSuffix("_FAST");
  keys.setValueDescription("atom","the position of the center of mass");
  keys.needsAction("MASS");
  keys.needsAction("SUM");
  keys.needsAction("CHARGE");
  keys.needsAction("CONSTANT");
  keys.needsAction("CUSTOM");
  keys.needsAction("POSITION");
  keys.needsAction("ARGS2VATOM");
}

CenterShortcut::CenterShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in what we are doing with the weights
  bool usemass;
  parseFlag("MASS",usemass);
  bool safe_phases;
  parseFlag("SAFE_PHASES",safe_phases);
  std::vector<std::string> str_weights;
  parseVector("WEIGHTS",str_weights);
  if( usemass || str_weights.size()==0 || str_weights.size()>1 || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
    std::string useph = "";
    if( safe_phases ) {
      useph = " PHASES";
    }
    if( usemass && str_weights.size()!=0 ) {
      error("WEIGHTS and MASS keywords cannot not be used simultaneously");
    }
    std::string wt_str;
    if( str_weights.size()>0 ) {
      wt_str="WEIGHTS=" + str_weights[0];
      for(unsigned i=1; i<str_weights.size(); ++i) {
        wt_str += "," + str_weights[i];
      }
    }
    if( usemass || (str_weights.size()==1 && str_weights[0]=="@Masses") ) {
      wt_str = "MASS";
    }
    readInputLine( getShortcutLabel() + ": CENTER_FAST " + wt_str + " " + convertInputLineToString() + useph );
    return;
  }
  // Read in the atoms
  std::string atlist;
  parse("ATOMS",atlist);
  // Calculate the mass of the vatom
  readInputLine( getShortcutLabel() + "_m: MASS ATOMS=" + atlist );
  readInputLine( getShortcutLabel() + "_mass: SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_m" );
  // Calculate the charge of the vatom
  readInputLine( getShortcutLabel() + "_q: CHARGE ATOMS=" + atlist );
  readInputLine( getShortcutLabel() + "_charge: SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_q" );
  // Retrieve the number of atoms
  ActionWithValue* am = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_m" );
  unsigned nat=am->copyOutput(0)->getNumberOfValues();
  // Get the weights to use for each atom
  std::string wlab = getShortcutLabel() + "_w";
  if( str_weights.size()>0 ) {
    if( str_weights.size()==1 ) {
      if( str_weights[0]=="@Charges" ) {
        wlab = getShortcutLabel() + "_q";
      } else {
        wlab=str_weights[0];
      }
    } else if( str_weights.size()==nat ) {
      std::string vals=str_weights[0];
      for(unsigned i=1; i<str_weights.size(); ++i) {
        vals += "," + str_weights[i];
      }
      readInputLine( getShortcutLabel() + "_w: CONSTANT VALUES=" + vals );
    } else {
      error("invalid input for WEIGHTS keyword " + str_weights[0] );
    }
  } else {
    std::string ones="1";
    for(unsigned i=1; i<nat; ++i) {
      ones += ",1";
    }
    readInputLine( getShortcutLabel() + "_w: CONSTANT VALUES=" + ones );
  }
  // Read in the instructions on how to compute the center of mass
  bool phases, nopbc;
  parseFlag("NOPBC",nopbc);
  if( safe_phases ) {
    phases=true;
  } else {
    parseFlag("PHASES",phases);
  }
  // This computes a center in the conventional way
  if( !phases || safe_phases ) {
    // Calculate the sum of the weights
    readInputLine( getShortcutLabel() + "_wnorm: SUM PERIODIC=NO ARG=" + wlab );
    // Compute the normalised weights
    readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_wnorm FUNC=x/y PERIODIC=NO");
    // Get the positions into a multicolvar
    if( phases || nopbc ) {
      readInputLine( getShortcutLabel() + "_pos: POSITION NOPBC ATOMS=" + atlist );
    } else {
      readInputLine( getShortcutLabel() + "_pos: POSITION WHOLEMOLECULES ATOMS=" + atlist );
    }
    // Multiply each vector of positions by the weight
    readInputLine( getShortcutLabel() + "_xwvec: CUSTOM ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_pos.x FUNC=x*y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_ywvec: CUSTOM ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_pos.y FUNC=x*y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zwvec: CUSTOM ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_pos.z FUNC=x*y PERIODIC=NO");
    // And sum the weighted vectors
    readInputLine( getShortcutLabel() + "_x: SUM ARG=" + getShortcutLabel() + "_xwvec PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_y: SUM ARG=" + getShortcutLabel() + "_ywvec PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_z: SUM ARG=" + getShortcutLabel() + "_zwvec PERIODIC=NO");
  }
  // This computes a center using the trigonometric phases
  if( phases ) {
    // Get the positions into a multicolvar
    readInputLine( getShortcutLabel() + "_fpos: POSITION SCALED_COMPONENTS ATOMS=" + atlist );
    // Calculate the sines and cosines of the positions and multiply by the weights
    readInputLine( getShortcutLabel() + "_sina: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.a FUNC=x*sin(2*pi*y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cosa: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.a FUNC=x*cos(2*pi*y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sinb: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.b FUNC=x*sin(2*pi*y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cosb: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.b FUNC=x*cos(2*pi*y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sinc: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.c FUNC=x*sin(2*pi*y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cosc: CUSTOM ARG=" + wlab + "," + getShortcutLabel() + "_fpos.c FUNC=x*cos(2*pi*y) PERIODIC=NO");
    // Sum the sines and cosines
    readInputLine( getShortcutLabel() + "_sinsuma: SUM ARG=" + getShortcutLabel() + "_sina PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cossuma: SUM ARG=" + getShortcutLabel() + "_cosa PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sinsumb: SUM ARG=" + getShortcutLabel() + "_sinb PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cossumb: SUM ARG=" + getShortcutLabel() + "_cosb PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sinsumc: SUM ARG=" + getShortcutLabel() + "_sinc PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cossumc: SUM ARG=" + getShortcutLabel() + "_cosc PERIODIC=NO");
    // And get the final position in fractional coordinates
    readInputLine( getShortcutLabel() + "_a: CUSTOM ARG=" + getShortcutLabel() + "_sinsuma," + getShortcutLabel() + "_cossuma FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_b: CUSTOM ARG=" + getShortcutLabel() + "_sinsumb," + getShortcutLabel() + "_cossumb FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_c: CUSTOM ARG=" + getShortcutLabel() + "_sinsumc," + getShortcutLabel() + "_cossumc FUNC=atan2(x,y)/(2*pi) PERIODIC=NO");
    // And create the virtual atom
    if( safe_phases ) {
      readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_a YPOS=" + getShortcutLabel() + "_b ZPOS=" + getShortcutLabel() + "_c "
                     + " XBKP=" + getShortcutLabel() + "_x YBKP=" + getShortcutLabel() + "_y ZBKP=" + getShortcutLabel() + "_z "
                     + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge FRACTIONAL");
    } else {
      readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_a YPOS=" + getShortcutLabel() + "_b ZPOS=" + getShortcutLabel() + "_c "
                     + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge FRACTIONAL");
    }
  } else {
    // And create the virtual atom
    readInputLine( getShortcutLabel() + ": ARGS2VATOM XPOS=" + getShortcutLabel() + "_x YPOS=" + getShortcutLabel() + "_y ZPOS=" + getShortcutLabel() + "_z "
                   + " MASS=" + getShortcutLabel() + "_mass CHARGE=" + getShortcutLabel() + "_charge ");
  }
}

}
}
