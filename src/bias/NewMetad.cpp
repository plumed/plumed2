/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/IFile.h"
#include "ReweightBase.h"

namespace PLMD {
namespace bias {

class NewMetad : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit NewMetad(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(NewMetad,"NEW_METAD")

void NewMetad::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","PACE","the frequency with which hills should be added");
  keys.add("compulsory","SIGMA","the widths of the Gaussians that we are using");
  keys.add("compulsory","HEIGHT","the heights of the hills that should be added");
  keys.add("optional","BIASFACTOR","the bias factor for the well tempered metadynamics");
  keys.add("compulsory","GRID_MIN","the minimum value to use for the the grid");
  keys.add("compulsory","GRID_MAX","the maximum value to use for all the grids");
  keys.add("compulsory","GRID_BIN","the number of bins to use for all the grids");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

NewMetad::NewMetad(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Calculate the well-tempered weight
  std::string biasfactor="", height; parse("HEIGHT",height); parse("BIASFACTOR",biasfactor);
  if( biasfactor.length()>0 ){ biasfactor = " BIASFACTOR=" + biasfactor; } else { biasfactor=" FIXED_HEIGHT"; }
  double temp=0.0; parse("TEMP",temp); std::string tempstr=""; 
  if( temp>0.0 ) { std::string tstr; Tools::convert( temp, tstr); tempstr = " TEMP=" + tstr; }
  readInputLine( getShortcutLabel() + "_wtfact: REWEIGHT_WELLTEMPERED HEIGHT=" + height + biasfactor + tempstr);

  // This stores the history dependent bias on a grid
  std::string pacestr; parse("PACE",pacestr); std::vector<std::string> args; parseVector("ARG",args);
  std::vector<std::string> gmin( args.size() ), gmax( args.size() ), grid_nbins(args.size()), sigma(args.size()); 
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax); 
  parseVector("GRID_BIN",grid_nbins); parseVector("SIGMA",sigma); 
  // This bit we are going to rewrite in KDE
  double bw; std::string center=" CENTER=0.0", band=" SIGMA=" + sigma[0], argstr=" READ_ARG=cv1";
  for(unsigned i=0;i<sigma.size();++i) {
      if( !Tools::convert( sigma[i], bw ) ) error("could not convert input bandwidth to real number");
      if( i>0 ) { center += ",0.0"; band += "," + sigma[i]; std::string nn; Tools::convert(i+1,nn); argstr += ",cv" + nn; }
  } 
  readInputLine( getShortcutLabel() + "_ref: READ_CLUSTER " + argstr + center + band );
  readInputLine( getShortcutLabel() + "_icov: CALCULATE_REFERENCE CONFIG=" + getShortcutLabel() + "_ref INPUT={MATHEVAL ARG1=" + getShortcutLabel() + "_ref.variance FUNC=1/x PERIODIC=NO}" ); 
  // This bit is good again
  readInputLine( getShortcutLabel() + "_height: CONSTANT VALUE=1.0");
  std::string input = getShortcutLabel() + "_kde: KDE_CALC METRIC=" + getShortcutLabel() + "_icov ARG1=" + args[0] + " HEIGHTS=" + getShortcutLabel() + "_height";
  std::string gminstr=" GRID_MIN=" + gmin[0]; std::string gmaxstr=" GRID_MAX=" + gmax[0]; std::string gbinstr=" GRID_BIN=" + grid_nbins[0];
  for(unsigned i=1;i<args.size();++i) { 
    std::string num; Tools::convert( i+1, num ); input += " ARG" + num + "=" + args[i]; 
    gminstr += "," + gmin[i]; gmaxstr += "," + gmax[i]; gbinstr += "," + grid_nbins[i]; 

  }
  readInputLine( input + " " + gminstr + " " + gmaxstr + " " + " " + gbinstr );
  readInputLine( getShortcutLabel() + "_grid: AVERAGE ARG=" + getShortcutLabel() + "_kde NORMALIZATION=false STRIDE=" + pacestr + " LOGWEIGHTS=" + getShortcutLabel() + "_wtfact");
  // Evaluate the instantaneous value of the bias potential 
  readInputLine( getShortcutLabel() + "_bias: EVALUATE_FUNCTION_FROM_GRID ARG=" + getShortcutLabel() + "_grid" );
  // Bias the simulation using this bias potentital
  readInputLine("BIASVALUE ARG=" + getShortcutLabel() + "_bias" );
  // Complete setup of the well tempered weights
  if( biasfactor.length()>0 ) {
      args.resize(1); args[0] = getShortcutLabel() + "_bias";
      ReweightBase* rwb = plumed.getActionSet().selectWithLabel<ReweightBase*>( getShortcutLabel() + "_wtfact" );
      plumed_assert( rwb ); rwb->setArguments( args );
  }
}

}
}
