/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "KDE.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace gridtools {

class MultiColvarDensity : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit MultiColvarDensity(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"MULTICOLVARDENS")

void MultiColvarDensity::registerKeywords( Keywords& keys ) {
  KDE::registerKeywords( keys );
  HistogramBase::histogramKeywords( keys );
  keys.add("compulsory","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","DIR","the direction in which to calculate the density profile");
  keys.add("optional","DATA","the multicolvar which you would like to calculate the density profile for");
  keys.add("optional","ATOMS","if you are calculating a atomic density you use this keyword to specify the atoms that are involved");
  keys.use("UPDATE_FROM"); keys.use("UPDATE_UNTIL");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
Action(ao),
ActionShortcut(ao)
{
    // Read in the position of the origin 
    std::string origin_str; parse("ORIGIN",origin_str);
    // Read in the quantity we are calculating the density for
    std::string atoms_str, data_str; parse("ATOMS",atoms_str); parse("DATA",data_str);
    if( atoms_str.length()==0 && data_str.length()==0 ) error("quantity to calculate the density for was not specified used DATA/ATOMS");
    // Get the information on the direction for the density
    std::string dir, direction_string; parse("DIR",dir);
    if( dir=="x" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.x";
    else if( dir=="y" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.y";
    else if( dir=="z" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.z";
    else if( dir=="xy" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.x ARG2=" + getShortcutLabel() + "_dist.y";
    else if( dir=="xz" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.x ARG2=" + getShortcutLabel() + "_dist.z";
    else if( dir=="yz" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.y ARG2=" + getShortcutLabel() + "_dist.z";
    else if( dir=="xyz" ) direction_string = "ARG1=" + getShortcutLabel() + "_dist.x ARG2=" + getShortcutLabel() + "_dist.y ARG3=" + getShortcutLabel() + "_dist.z";
    else error( dir + " is invalid dir specification use x/y/z/xy/xz/yz/xyz");

    // Parse the keymap for this averaging stuff
    std::map<std::string,std::string> keymap; HistogramBase::readHistogramKeywords( keymap, this );
    // Create distance action
    bool hasheights; std::string dist_words = getShortcutLabel() + "_dist: DISTANCES COMPONENTS ORIGIN=" + origin_str; 
    if( atoms_str.length()>0 ) { hasheights=false; dist_words += " ATOMS=" + atoms_str; }
    else { hasheights=true; dist_words += " ATOMS=" + data_str; }
    // plumed_massert( keys.count("ORIGIN"), "you must specify the position of the origin" );
    readInputLine( dist_words );

    std::string inputLine = convertInputLineToString();
    // Make the kde object for the numerator if needed
    if( hasheights ) {
      readInputLine( getShortcutLabel() + "_inumer: KDE UNORMALIZED HEIGHTS=" + data_str + " " + direction_string + " " + inputLine );
      HistogramBase::createAveragingObject( getShortcutLabel() + "_inumer", getShortcutLabel() + "_numer", keymap, this );
    }
    // Make the kde object
    readInputLine( getShortcutLabel() + "_kde: KDE " + inputLine  + " " + direction_string );
    // Make the division object if it is required
    if( hasheights ) {
      HistogramBase::createAveragingObject( getShortcutLabel() + "_kde", getShortcutLabel() + "_denom", keymap, this);
      readInputLine( getShortcutLabel() + ": MATHEVAL ARG1=" + getShortcutLabel() + "_numer ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    } else {
      HistogramBase::createAveragingObject( getShortcutLabel() + "_kde", getShortcutLabel(), keymap, this );
    }
}

}
}
