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
#include "core/ActionShortcut.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace bias {

class Walls : public ActionShortcut {
public:
  explicit Walls(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Walls,"UPPER_WALLS")
PLUMED_REGISTER_ACTION(Walls,"LOWER_WALLS")

void Walls::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.addInputKeyword("numbered","ARG","scalar/vector","the arguments on which the bias is acting");
  keys.add("compulsory","AT","the positions of the wall. The a_i in the expression for a wall.");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.addOutputComponent("bias","default","scalar","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","scalar","the instantaneous value of the squared force due to this bias potential");
  keys.addActionNameSuffix("_SCALAR");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
  keys.needsAction("COMBINE");
  keys.needsAction("BIASVALUE");
}

Walls::Walls(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read the arguments
  std::vector<std::string> args;
  parseVector("ARG",args);
  if( args.size()==0 ) {
    error("found no input arguments");
  }
  std::string allargs=args[0];
  for(unsigned i=1; i<args.size(); ++i) {
    allargs += "," + args[i];
  }
  std::vector<Value*> vals;
  ActionWithArguments::interpretArgumentList( args, plumed.getActionSet(), this, vals );
  if( vals.size()==0 ) {
    error("found no input arguments");
  }

  // Find the rank
  unsigned rank=vals[0]->getRank();
  for(unsigned i=0; i<vals.size(); ++i) {
    if( vals[i]->getRank()>0 && vals[i]->hasDerivatives() ) {
      error("argument should not be function on grid");
    }
    if( vals[i]->getRank()!=rank ) {
      error("all arguments should have same rank");
    }
  }
  if( rank==0 ) {
    if( getName()=="UPPER_WALLS") {
      readInputLine( getShortcutLabel() + ": UPPER_WALLS_SCALAR ARG=" + allargs + " " + convertInputLineToString() );
    } else if( getName()=="LOWER_WALLS") {
      readInputLine( getShortcutLabel() + ": LOWER_WALLS_SCALAR ARG=" + allargs + " " + convertInputLineToString() );
    } else {
      plumed_merror( getName() + " is not valid");
    }
    return;
  }

  // Note : the sizes of these vectors are checked automatically by parseVector
  std::vector<std::string> kappa(args.size());
  parseVector("KAPPA",kappa);
  std::vector<std::string> offset(kappa.size());
  parseVector("OFFSET",offset);
  std::vector<std::string> eps(kappa.size());
  parseVector("EPS",eps);
  std::vector<std::string> exp(kappa.size());
  parseVector("EXP",exp);
  std::vector<std::string> at(kappa.size());
  parseVector("AT",at);

  std::string biasinp, forceinp;
  for(unsigned i=0; i<args.size(); ++i) {
    std::string argn;
    std::size_t dot=args[i].find_first_of(".");
    if(dot!=std::string::npos) {
      argn = args[i].substr(0,dot) + "_" + args[i].substr(dot+1);
    } else {
      argn = args[i];
    }
    readInputLine( getShortcutLabel() + "_cv_" + argn + ": COMBINE PERIODIC=NO ARG=" + args[i] + " PARAMETERS=" + at[i] );
    if( getName()=="UPPER_WALLS" ) {
      readInputLine( getShortcutLabel() + "_scale_" + argn + ": CUSTOM PERIODIC=NO FUNC=(x+" + offset[i] +")/" + eps[i] + " ARG=" + getShortcutLabel() + "_cv_" + argn );
      readInputLine( getShortcutLabel() + "_pow_" + argn + ": CUSTOM PERIODIC=NO FUNC=step(x)*x^" + exp[i] + " ARG=" + getShortcutLabel() + "_scale_" + argn );
    } else if( getName()=="LOWER_WALLS" ) {
      readInputLine( getShortcutLabel() + "_scale_" + argn + ": CUSTOM PERIODIC=NO FUNC=(x-" + offset[i] +")/" + eps[i] + " ARG=" + getShortcutLabel() + "_cv_" + argn );
      readInputLine( getShortcutLabel() + "_pow_" + argn + ": CUSTOM PERIODIC=NO FUNC=step(-x)*(-x)^" + exp[i] + " ARG=" + getShortcutLabel() + "_scale_" + argn );
    }
    readInputLine( getShortcutLabel() + "_v_wall_" + argn + ": CUSTOM PERIODIC=NO FUNC=" + kappa[i] +"*x" + " ARG=" + getShortcutLabel() + "_pow_" + argn );
    readInputLine( getShortcutLabel() + "_wall_" + argn + ": SUM ARG=" + getShortcutLabel() + "_v_wall_" + argn + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_force_" + argn + ": CUSTOM PERIODIC=NO FUNC=" + kappa[i] + "*" + exp[i] + "*x/(y*" + eps[i] + ") " +
                   "ARG=" + getShortcutLabel() + "_pow_" + argn + "," + getShortcutLabel() + "_scale_" + argn );
    readInputLine( getShortcutLabel() + "_v_force2_" + argn + ": CUSTOM PERIODIC=NO FUNC=x*x ARG=" + getShortcutLabel() + "_force_" + argn );
    readInputLine( getShortcutLabel() + "_force2_" + argn + ": SUM ARG=" + getShortcutLabel() + "_v_force2_" + argn + " PERIODIC=NO");
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
