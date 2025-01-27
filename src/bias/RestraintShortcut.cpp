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
    std::string slopestr="";
    if( slope[0]!="0.0" ) {
      slopestr="SLOPE=" + slope[0];
      for(unsigned i=1; i<slope.size(); ++i) {
        slopestr += "," + slope[i];
      }
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
  bool non_constant_force=false;
  for(unsigned i=0; i<args.size(); ++i) {
    std::string argn;
    std::size_t dot=args[i].find_first_of(".");
    if(dot!=std::string::npos) {
      argn = args[i].substr(0,dot) + "_" + args[i].substr(dot+1);
    } else {
      argn = args[i];
    }
    readInputLine( getShortcutLabel() + "_cv_" + argn + ": COMBINE PERIODIC=NO ARG=" + args[i] + " PARAMETERS=" + at[i] );
    double kap;
    Tools::convert(  kappa[i], kap );
    if( fabs(kap)>0 ) {
      non_constant_force = true;
      readInputLine( getShortcutLabel() + "_harm_" + argn + ": CUSTOM PERIODIC=NO FUNC=0.5*" + kappa[i] + "*x^2 ARG=" + getShortcutLabel() + "_cv_" + argn );
      readInputLine( getShortcutLabel() + "_kap_" + argn + ": SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_harm_" + argn );
      readInputLine( getShortcutLabel() + "_f2_" + argn + ": CUSTOM PERIODIC=NO FUNC=" + kappa[i] + "*2*x+" + slope[i] + "*" + slope[i] + " ARG=" + getShortcutLabel() + "_harm_" + argn );
      if( i==0 ) {
        biasargs = "ARG=" + getShortcutLabel() + "_kap_" + argn;
      } else {
        biasargs += "," + getShortcutLabel() + "_kap_" + argn;
      }
      if( i==0 ) {
        forceargs = "ARG=" + getShortcutLabel() + "_f2_" + argn;
      } else {
        forceargs += "," + getShortcutLabel() + "_f2_" + argn;
      }
    }
    double slo;
    Tools::convert( slope[i], slo );
    if( fabs(slo)>0 ) {
      readInputLine( getShortcutLabel() + "_linear_" + argn + ": CUSTOM PERIODIC=NO FUNC=" + slope[i] + "*x ARG=" + getShortcutLabel() + "_cv_" + argn );
      readInputLine( getShortcutLabel() + "_slope_" + argn + ": SUM PERIODIC=NO ARG=" + getShortcutLabel() + "_linear_" + argn );
      if( biasargs.length()==0 ) {
        biasargs = "ARG=" + getShortcutLabel() + "_slope_" + argn;
      } else {
        biasargs += "," + getShortcutLabel() + "_slope_" + argn;
      }
    }
  }
  // This is the bias
  readInputLine( getShortcutLabel() + "_bias: COMBINE PERIODIC=NO " + biasargs );
  readInputLine( getShortcutLabel() + ": BIASVALUE ARG=" + getShortcutLabel() + "_bias STRIDE=" + stride );
  if( non_constant_force ) {
    readInputLine( getShortcutLabel() + "_force2: COMBINE PERIODIC=NO " + forceargs );
  }
}

}


}
