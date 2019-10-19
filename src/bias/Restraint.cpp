/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "ActionRegister.h"
#include "core/ActionShortcut.h"


using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS RESTRAINT
/*
Adds harmonic and/or linear restraints on one or more variables.

Either or both
of SLOPE and KAPPA must be present to specify the linear and harmonic force constants
respectively.  The resulting potential is given by:
\f[
  \sum_i \frac{k_i}{2} (x_i-a_i)^2 + m_i*(x_i-a_i)
\f].

The number of components for any vector of force constants must be equal to the number
of arguments to the action.

Additional material and examples can be also found in the tutorial \ref belfast-4

\par Examples

The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and to print the energy of the restraint
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
RESTRAINT ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 LABEL=restraint
PRINT ARG=restraint.bias
\endplumedfile

*/
//+ENDPLUMEDOC

class Restraint : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Restraint(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Restraint,"RESTRAINT")

void Restraint::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ARG","the arguments on which the bias is acting");
  keys.add("compulsory","SLOPE","0.0","specifies that the restraint is linear and what the values of the force constants on each of the variables are");
  keys.add("compulsory","KAPPA","0.0","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
  keys.add("compulsory","AT","the position of the restraint");
  keys.add("hidden","STRIDE","1","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}

Restraint::Restraint(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string stride; parse("STRIDE",stride);
  std::vector<std::string> at; parseVector("AT",at);
  std::vector<std::string> slope(at.size()); parseVector("SLOPE",slope);
  std::vector<double> kappa(at.size()); parseVector("KAPPA",kappa);
  // This should just be the arguments
  std::string argstr=convertInputLineToString();
  // This makes us the list of centers
  std::string atstr=" PARAMETERS=" + at[0]; for(unsigned i=1;i<at.size();++i) atstr += "," + at[i];

  // Sort the calculation of the bias
  std::string kstr; Tools::convert( 0.5*kappa[0], kstr ); 
  std::string powstr=" POWERS=2", coefs=" COEFFICIENTS=" + kstr, slopes=" COEFFICIENTS=" + slope[0];
  for(unsigned i=1;i<kappa.size();++i) { powstr += ",2"; Tools::convert( 0.5*kappa[i], kstr ); coefs += "," + kstr; slopes += "," + slope[i]; }
  readInputLine( getShortcutLabel() + "_kap: COMBINE PERIODIC=NO " + argstr + powstr + atstr + coefs );
  readInputLine( getShortcutLabel() + "_slope: COMBINE PERIODIC=NO " + argstr + atstr + slopes );
  readInputLine( getShortcutLabel() + "_bias: COMBINE ARG1=" + getShortcutLabel() + "_kap ARG2=" + getShortcutLabel() + "_slope PERIODIC=NO");
  readInputLine( getShortcutLabel() + ": BIASVALUE ARG=" + getShortcutLabel() + "_bias STRIDE=" + stride );
  // And now the force squared
  Tools::convert( kappa[0]*kappa[0], kstr ); coefs=" COEFFICIENTS=" + kstr; double sltmp, slsum; Tools::convert( slope[0], sltmp ); slsum=sltmp*sltmp;
  for(unsigned i=1;i<kappa.size();++i) { Tools::convert( kappa[i]*kappa[i], kstr ); coefs += "," + kstr; Tools::convert( slope[i], sltmp ); slsum+=sltmp*sltmp; }
  readInputLine( getShortcutLabel() + "_der2: COMBINE PERIODIC=NO " + argstr + powstr + coefs + atstr ); 
  Tools::convert( slsum, kstr ); readInputLine( getShortcutLabel() + "_force2: MATHEVAL PERIODIC=NO ARG=" + getShortcutLabel() + "_der2 FUNC=x+" + kstr );
}

}


}
