/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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

Additional material and examples can be also found in the tutorial \ref lugano-2

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
  std::vector<std::string> kappa(at.size()); parseVector("KAPPA",kappa);
  // Read in the args
  std::vector<std::string> args; parseVector("ARG",args); std::string biasinp, forceinp;
  if( args.size()==0 ) {
      args.resize(kappa.size() );
      for(unsigned i=0;i<kappa.size();++i) {
          if( !parseNumbered("ARG",i+1,args[i]) ) error("failed to find sufficient numbered ARG keywords");
      }
  }
  
  std::string biasargs, forceargs; bool non_constant_force=false;
  for(unsigned i=0;i<args.size();++i) {
      std::string argn=args[i]; std::size_t dot=argn.find_first_of("."); if(dot!=std::string::npos) argn = argn.substr(0,dot) + "_" + argn.substr(dot+1);
      readInputLine( getShortcutLabel() + "_cv_" + argn + ": COMBINE PERIODIC=NO ARG1=" + args[i] + " PARAMETERS=" + at[i] );
      if( kappa[i]!="0.0" ) {
          non_constant_force = true;
          readInputLine( getShortcutLabel() + "_harm_" + argn + ": MATHEVAL NO_WILDCARD PERIODIC=NO FUNC=0.5*" + kappa[i] + "*x^2 ARG1=" + getShortcutLabel() + "_cv_" + argn );
          readInputLine( getShortcutLabel() + "_kap_" + argn + ": SUM NO_WILDCARD PERIODIC=NO ARG=" + getShortcutLabel() + "_harm_" + argn ); 
          readInputLine( getShortcutLabel() + "_f2_" + argn + ": MATHEVAL NO_WILDCARD PERIODIC=NO FUNC=" + kappa[i] + "*2*x+" + slope[i] + "*" + slope[i] + " ARG1=" + getShortcutLabel() + "_harm_" + argn ); 
          if( i==0 ) biasargs = "ARG=" + getShortcutLabel() + "_kap_" + argn; else biasargs += "," + getShortcutLabel() + "_kap_" + argn;
          if( i==0 ) forceargs = "ARG=" + getShortcutLabel() + "_f2_" + argn; else forceargs += "," + getShortcutLabel() + "_f2_" + argn;    
     }  
      if( slope[i]!="0.0" ) {
          readInputLine( getShortcutLabel() + "_linear_" + argn + ": MATHEVAL PERIODIC=NO FUNC=" + slope[i] + "*x ARG1=" + getShortcutLabel() + "_cv_" + argn );
          readInputLine( getShortcutLabel() + "_slope_" + argn + ": SUM NO_WILDCARD PERIODIC=NO ARG=" + getShortcutLabel() + "_linear_" + argn );
          if( biasargs.length()==0 ) biasargs = "ARG=" + getShortcutLabel() + "_slope_" + argn; else biasargs += "," + getShortcutLabel() + "_slope_" + argn;
      }
  }
  // This is the bias
  readInputLine( getShortcutLabel() + "_bias: COMBINE PERIODIC=NO " + biasargs );
  readInputLine( getShortcutLabel() + ": BIASVALUE ARG=" + getShortcutLabel() + "_bias STRIDE=" + stride );
  if( non_constant_force ) readInputLine( getShortcutLabel() + "_force2: COMBINE PERIODIC=NO " + forceargs );
}

}


}
