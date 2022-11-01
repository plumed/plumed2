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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR ENVIRONMENT
/*

*/
//+ENDPLUMEDOC


class EnvironmentalSimilarity : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit EnvironmentalSimilarity(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(EnvironmentalSimilarity,"ENVIRONMENT")

void EnvironmentalSimilarity::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the input adjacency matrix that you would like to use for this calculation");
  keys.add("compulsory","SIGMA","the width to use for the gaussian kernels");
  keys.add("compulsory","REFERENCE","a list of vectors that describe the environment");
}

EnvironmentalSimilarity::EnvironmentalSimilarity(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string matlab; parse("MATRIX",matlab); std::string reffile; parse("REFERENCE",reffile); PDB pdb; std::string argstr; 
  pdb.read(reffile,plumed.usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength());
  unsigned natoms=pdb.getPositions().size(); std::string str_natoms; Tools::convert( natoms, str_natoms ); 
  double sig; parse("SIGMA",sig); std::string sig2; Tools::convert( sig*sig, sig2 );
  for(unsigned i=0;i<natoms;++i) {
     std::string num, xpos, ypos, zpos; Tools::convert( i+1, num );
     Tools::convert( pdb.getPositions()[i][0], xpos ); Tools::convert( pdb.getPositions()[i][1], ypos ); Tools::convert( pdb.getPositions()[i][2], zpos );
     readInputLine( getShortcutLabel() + "_exp" + num + ": MATHEVAL ARG1=" + matlab + ".x ARG2=" + matlab + ".y ARG3=" + matlab + ".z ARG4=" + matlab + ".w VAR=x,y,z,w" + 
                                                        " FUNC=w*exp(-((x-" + xpos + ")^2+(y-" + ypos + ")^2+(z-" + zpos + ")^2)/4*" + sig2 + ")/" + str_natoms + " PERIODIC=NO");
     argstr += " ARG" + num + "=" + getShortcutLabel() + "_exp" + num;
  }
  readInputLine( getShortcutLabel() + "_tot: COMBINE PERIODIC=NO" + argstr );
  readInputLine( getShortcutLabel() + ": COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_tot");
}

}
}
