/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "ActionRegister.h"


using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS UPPER_WALLS
/*
Defines a wall for the value of one or more collective variables,
 which limits the region of the phase space accessible during the simulation.

The restraining potential starts acting on the system when the value of the CV is greater
(in the case of UPPER_WALLS) or lower (in the case of LOWER_WALLS) than a certain limit \f$a_i\f$ (AT)
minus an offset \f$o_i\f$ (OFFSET).
The expression for the bias due to the wall is given by:

\f$
  \sum_i {k_i}((x_i-a_i+o_i)/s_i)^e_i
\f$

\f$k_i\f$ (KAPPA) is an energy constant in internal unit of the code, \f$s_i\f$ (EPS) a rescaling factor and
\f$e_i\f$ (EXP) the exponent determining the power law. By default: EXP = 2, EPS = 1.0, OFFSET = 0.


\par Examples

The following input tells plumed to add both a lower and an upper walls on the distance
between atoms 3 and 5 and the distance between atoms 2 and 4. The lower and upper limits
are defined at different values. The strength of the walls is the same for the four cases.
It also tells plumed to print the energy of the walls.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
UPPER_WALLS ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET=0,0 LABEL=uwall
LOWER_WALLS ARG=d1,d2 AT=0.0,1.0 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET=0,0 LABEL=lwall
PRINT ARG=uwall.bias,lwall.bias
\endplumedfile

*/
//+ENDPLUMEDOC

class UWalls : public ActionShortcut {
public:
  explicit UWalls(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(UWalls,"UPPER_WALLS")
PLUMED_REGISTER_ACTION(UWalls,"LOWER_WALLS")

void UWalls::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.add("numbered","ARG","the arguments on which the bias is acting");
  keys.add("compulsory","AT","the positions of the wall. The a_i in the expression for a wall.");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
}

UWalls::UWalls(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Note : the sizes of these vectors are checked automatically by parseVector
  std::vector<std::string> kappa; parseVector("KAPPA",kappa);
  std::vector<std::string> offset(kappa.size()); parseVector("OFFSET",offset);
  std::vector<std::string> eps(kappa.size()); parseVector("EPS",eps);
  std::vector<std::string> exp(kappa.size()); parseVector("EXP",exp);
  std::vector<std::string> at(kappa.size()); parseVector("AT",at);
  // Read in the args
  std::vector<std::string> args; parseVector("ARG",args); std::string biasinp, forceinp;
  if( args.size()==0 ) { 
      args.resize(kappa.size() ); 
      for(unsigned i=0;i<kappa.size();++i) {
          if( !parseNumbered("ARG",i+1,args[i]) ) error("failed to find sufficient numbered ARG keywords");
      }
  }
  
  for(unsigned i=0;i<args.size();++i) {
      std::string argn=args[i]; std::size_t dot=argn.find_first_of("."); if(dot!=std::string::npos) argn = argn.substr(0,dot) + "_" + argn.substr(dot+1);  
      readInputLine( getShortcutLabel() + "_cv_" + argn + ": COMBINE PERIODIC=NO ARG1=" + args[i] + " PARAMETERS=" + at[i] );
      if( getName()=="UPPER_WALLS" ) {
          readInputLine( getShortcutLabel() + "_scale_" + argn + ": MATHEVAL PERIODIC=NO FUNC=(x+" + offset[i] +")/" + eps[i] + " ARG1=" + getShortcutLabel() + "_cv_" + argn );
          readInputLine( getShortcutLabel() + "_pow_" + argn + ": MATHEVAL PERIODIC=NO FUNC=step(x)*x^" + exp[i] + " ARG1=" + getShortcutLabel() + "_scale_" + argn );
      } else if( getName()=="LOWER_WALLS" ) {
          readInputLine( getShortcutLabel() + "_scale_" + argn + ": MATHEVAL PERIODIC=NO FUNC=(x-" + offset[i] +")/" + eps[i] + " ARG1=" + getShortcutLabel() + "_cv_" + argn );
          readInputLine( getShortcutLabel() + "_pow_" + argn + ": MATHEVAL PERIODIC=NO FUNC=step(-x)*x^" + exp[i] + " ARG1=" + getShortcutLabel() + "_scale_" + argn );
      }
      readInputLine( getShortcutLabel() + "_wall_" + argn + ": MATHEVAL PERIODIC=NO FUNC=" + kappa[i] +"*x" + " ARG1=" + getShortcutLabel() + "_pow_" + argn );
      readInputLine( getShortcutLabel() + "_force_" + argn + ": MATHEVAL PERIODIC=NO FUNC=" + kappa[i] + "*" + exp[i] + "*x/(y*" + eps[i] + ") " +
                     "ARG1=" + getShortcutLabel() + "_pow_" + argn + " ARG2=" + getShortcutLabel() + "_scale_" + argn ); 
      readInputLine( getShortcutLabel() + "_force2_" + argn + ": MATHEVAL PERIODIC=NO FUNC=x*x ARG1=" + getShortcutLabel() + "_force_" + argn );
      if(i==0) {
         biasinp = " ARG=" + getShortcutLabel() + "_wall_" + argn; 
         forceinp = " ARG=" + getShortcutLabel() + "_force2_" + argn;
      } else {
         biasinp += "," + getShortcutLabel() + "_wall_" + argn;
         forceinp += "," + getShortcutLabel() + "_force2_" + argn;
      }
  }
  readInputLine( getShortcutLabel() + "_bias: COMBINE PERIODIC=NO " + biasinp );
  readInputLine( "BIASVALUE ARG=" + getShortcutLabel() + "_bias" );
  readInputLine( getShortcutLabel() + "_force2: COMBINE PERIODIC=NO " + forceinp  );
}

}
}
