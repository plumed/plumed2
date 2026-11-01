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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace bias {

class RestraintShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit RestraintShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(RestraintShortcut,"RESTRAINT")

void RestraintShortcut::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("numbered","ARG","scalar/vector","the values the harmonic restraint acts upon");
  keys.add("compulsory","SLOPE","0.0","specifies that the restraint is linear and what the values of the force constants on each of the variables are");
  keys.add("compulsory","KAPPA","0.0","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
  keys.add("compulsory","AT","the position of the restraint");
  keys.add("hidden","STRIDE","1","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.addOutputComponent("bias","default","scalar","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","scalar/vector","the instantaneous value of the squared force due to this bias potential");
  keys.addActionNameSuffix("_SCALAR");
  keys.needsAction("COMBINE");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
  keys.needsAction("BIASVALUE");
}

RestraintShortcut::RestraintShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the args
  std::vector<std::string> args;
  parseVector("ARG",args);
  if( args.size()==0 ) {
    error("found no input arguments");
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
    std::vector<std::string> slope(args.size());
    parseVector("SLOPE",slope);
    // all the slopes are passed: checking only the first one dropped e.g. SLOPE=0.0,1.0
    std::string slopestr="SLOPE=" + slope[0];
    for(unsigned i=1; i<slope.size(); ++i) {
      slopestr += "," + slope[i];
    }
    std::string allargs=args[0];
    for(unsigned i=1; i<args.size(); ++i) {
      allargs += "," + args[i];
    }
    readInputLine( getShortcutLabel() + ": RESTRAINT_SCALAR ARG=" + allargs + " " + slopestr + " " + convertInputLineToString() );
    return;
  }

  std::string stride;
  parse("STRIDE",stride);
  std::vector<std::string> at;
  parseVector("AT",at);
  std::vector<std::string> slope(at.size());
  parseVector("SLOPE",slope);
  std::vector<std::string> kappa(at.size());
  parseVector("KAPPA",kappa);

  std::string biasargs, forceargs;
  for(unsigned i=0; i<args.size(); ++i) {
    std::string argn;
    std::size_t dot=args[i].find_first_of(".");
    if(dot!=std::string::npos) {
      argn = args[i].substr(0,dot) + "_" + args[i].substr(dot+1);
    } else {
      argn = args[i];
    }
    double kap, slo;
    Tools::convert( kappa[i], kap );
    Tools::convert( slope[i], slo );
    if( kap==0 && slo==0 ) {
      continue;
    }
    // distance from AT: COMBINE applies the periodic boundary conditions; for a non periodic argument
    // the difference is computed directly in the CUSTOM functions, which saves one action
    std::vector<Value*> argvals;
    ActionWithArguments::interpretArgumentList( std::vector<std::string>(1,args[i]), plumed.getActionSet(), this, argvals );
    bool periodic=false;
    for(const auto & v : argvals) {
      periodic = periodic || v->isPeriodic();
    }
    std::string cvarg, cv;
    if( periodic ) {
      readInputLine( getShortcutLabel() + "_cv_" + argn + ": COMBINE PERIODIC=NO ARG=" + args[i] + " PARAMETERS=" + at[i] );
      cvarg = getShortcutLabel() + "_cv_" + argn;
      cv = "x";
    } else {
      cvarg = args[i];
      cv = "(x-(" + at[i] + "))";
    }
    if( kap!=0 ) {
      readInputLine( getShortcutLabel() + "_harm_" + argn + ": CUSTOM PERIODIC=NO FUNC=0.5*" + kappa[i] + "*" + cv + "^2 ARG=" + cvarg );
      readInputLine( getShortcutLabel() + "_kap_" + argn + ": SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_harm_" + argn );
      biasargs += (biasargs.length()==0 ? "ARG=" : ",") + getShortcutLabel() + "_kap_" + argn;
    }
    if( slo!=0 ) {
      readInputLine( getShortcutLabel() + "_linear_" + argn + ": CUSTOM PERIODIC=NO FUNC=" + slope[i] + "*" + cv + " ARG=" + cvarg );
      readInputLine( getShortcutLabel() + "_slope_" + argn + ": SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_linear_" + argn );
      biasargs += (biasargs.length()==0 ? "ARG=" : ",") + getShortcutLabel() + "_slope_" + argn;
    }
    // squared force (k*x+m)^2, with x the distance from AT
    readInputLine( getShortcutLabel() + "_f2_" + argn + ": CUSTOM PERIODIC=NO FUNC=((" + kappa[i] + ")*" + cv + "+(" + slope[i] + "))^2 ARG=" + cvarg );
    forceargs += (forceargs.length()==0 ? "ARG=" : ",") + getShortcutLabel() + "_f2_" + argn;
  }
  if( biasargs.length()==0 ) {
    error("all the values of KAPPA and SLOPE are zero");
  }
  // This is the bias
  readInputLine( getShortcutLabel() + "_bias: COMBINE PERIODIC=NO " + biasargs );
  readInputLine( getShortcutLabel() + ": BIASVALUE ARG=" + getShortcutLabel() + "_bias STRIDE=" + stride );
  readInputLine( getShortcutLabel() + "_force2: COMBINE PERIODIC=NO " + forceargs );
}

}


}
