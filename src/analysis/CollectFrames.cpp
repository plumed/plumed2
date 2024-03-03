/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC ANALYSIS COLLECT_FRAMES
/*
This allows you to convert a trajectory and a dissimilarity matrix into a dissimilarity object

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class CollectFrames : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& ao ); 
};

PLUMED_REGISTER_ACTION(CollectFrames,"COLLECT_FRAMES")

void CollectFrames::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","STRIDE","the frequency with which data should be stored for analysis.  By default data is collected on every step");
  keys.add("compulsory","CLEAR","0","the frequency with which data should all be deleted and restarted");
  keys.add("compulsory","ARG","the arguments you would like to collect");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
  keys.addOutputComponent("data","default","the data that is being collected by this action");
  keys.addOutputComponent("logweights","default","the logarithms of the weights of the data points");
}

CollectFrames::CollectFrames( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  std::string stride, clearstride; parse("STRIDE",stride); parse("CLEAR",clearstride);
  std::vector<std::string> argn; parseVector("ARG",argn); std::vector<Value*> theargs;
  ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, theargs );
  if( theargs.size()==0 ) error("no arguments were specified for collection");

  // Create all the collect actions
  for(unsigned i=0; i<theargs.size(); ++i) {
      if( theargs[i]->getNumberOfValues()!=theargs[0]->getNumberOfValues() ) error("mismatch between number of arguments calculated by each collected argument");
      readInputLine( getShortcutLabel() + "_" + theargs[i]->getName() + ": COLLECT ARG=" + theargs[i]->getName() + " STRIDE=" + stride + " CLEAR=" + clearstride );
  }
  // Make a list of collect actions
  std::string allcol = getShortcutLabel() + "_" + theargs[0]->getName(); for(unsigned i=1; i<theargs.size(); ++i) allcol + "," + getShortcutLabel() + "_" + theargs[i]->getName();
  // And transfer everything to a matrix
  readInputLine( getShortcutLabel() + "_data: VSTACK ARG=" + allcol );

  // Now get the logweights
  std::vector<std::string> logw; parseVector("LOGWEIGHTS",logw); std::vector<Value*> thew;
  if( logw.size()>0 ) ActionWithArguments::interpretArgumentList( logw, plumed.getActionSet(), this, thew );
  if( logw.size()>1 ) error("maximum of one argument should be specified for logweights");

  if( logw.size()==0 ) {
      std::string zeros="0"; for(unsigned i=1; i<theargs[0]->getNumberOfValues(); ++i) zeros += ",0";
      readInputLine( getShortcutLabel() + "_cweight: CONSTANT VALUE=" + zeros );
      readInputLine( getShortcutLabel() + "_logweights: COLLECT ARG=" + getShortcutLabel() + "_cweight STRIDE=" + stride + " CLEAR=" + clearstride ); 
  } else {
      if( theargs[0]->getNumberOfValues()!=thew[0]->getNumberOfValues() ) error("mismatch between number of weights and number of collected arguments");
      readInputLine( getShortcutLabel() + "_logweights: COLLECT ARG=" + thew[0]->getName() + " STRIDE=" + stride + " CLEAR=" + clearstride );
  }
}

}
}
