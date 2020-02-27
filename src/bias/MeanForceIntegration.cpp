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
#include "core/AverageBase.h"

namespace PLMD {
namespace bias {

class MeanForceIntegration : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MeanForceIntegration(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MeanForceIntegration,"MEAN_FORCE_INTEGRATION")

void MeanForceIntegration::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","BANDWIDTH","the bandwidth to use when computing the average foces");
  keys.add("compulsory","STRIDE","the frequency with respect to accumulate the average force grid");
}

MeanForceIntegration::MeanForceIntegration(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // This retrieves the information about the grid that we will use to setup the new grids
  std::string arg; parse("ARG",arg); 
  AverageBase* mybias = plumed.getActionSet().selectWithLabel<AverageBase*>( arg + "_grid" );
  if( !mybias ) error("cannot do mean force integration as I couldn't find the grid that stores the metadynamics bias");
  unsigned dim = mybias->copyOutput(0)->getRank();
  std::vector<bool> grid_pbc(dim); std::vector<double> grid_spacing(dim); 
  std::string gtype; std::vector<std::string> args(dim), gmin(dim), gmax(dim); std::vector<unsigned> grid_nbins(dim); 
  mybias->getInfoForGridHeader( gtype, args, gmin, gmax, grid_nbins, grid_spacing, grid_pbc, false );  
  // Get metadynamics bias forces on grids
  readInputLine( getShortcutLabel() + "_fmetad: GET_GRID_DERIVATIVES ARG=" + arg + "_grid"); 
  // This stores the force under the effect of the peturbation
  int pace = mybias->getStride(); std::string pacestr; Tools::convert(pace, pacestr);
  std::string stride; parse("STRIDE",stride); std::vector<std::string> sigma(args.size()); parseVector("BANDWIDTH",sigma); 
  std::string input = getShortcutLabel() + "_pforce: HISTOGRAM NORMALIZATION=true CLEAR=" + pacestr + " STRIDE=" + stride + " ARG1=" + args[0];
  std::string gminstr=" GRID_MIN=" + gmin[0]; std::string gmaxstr=" GRID_MAX=" + gmax[0]; std::string bandstr=" BANDWIDTH=" + sigma[0];
  std::string str_nbins; Tools::convert( grid_nbins[0], str_nbins ); std::string gbinstr=" GRID_BIN=" + str_nbins;
  for(unsigned i=1;i<args.size();++i) { 
    std::string num; Tools::convert( i+1, num ); input += " ARG" + num + "=" + args[i]; 
    gminstr += "," + gmin[i]; gmaxstr += "," + gmax[i]; bandstr += "," + sigma[i]; 
    Tools::convert( grid_nbins[i], str_nbins ); gbinstr += "," + str_nbins; 

  }
  readInputLine( input + " " + gminstr + " " + gmaxstr + " " + bandstr + " " + gbinstr ); std::vector<std::string> nargs(args);
  // Get kernel forces on grid
  readInputLine( getShortcutLabel() + "_fpforce: GET_GRID_DERIVATIVES ARG=" + getShortcutLabel() + "_pforce");
  for(unsigned i=0;i<args.size();++i) { 
      // Replace dots in argument names so they can be used as action names
      std::size_t dot=args[i].find("."); if( dot!=std::string::npos ) nargs[i]=args[i].substr(0,dot) + "_" + args[i].substr(dot+1);
      // Now compute the unperturbed mean force
      readInputLine( getShortcutLabel() + "_uforce_" + nargs[i] + ": MATHEVAL ARG1=" + getShortcutLabel() + "_fmetad." + args[i] + "_der" + 
                                          " ARG2=" + getShortcutLabel() + "_fpforce." + args[i] + "_der ARG3=" + getShortcutLabel() + "_pforce FUNC=-x-y/z PERIODIC=NO");
      // And compute the quantity that must be incremeted over the simulation
      readInputLine( getShortcutLabel() + "_mforce_" + nargs[i] + ": MATHEVAL PERIODIC=NO FUNC=x*y ARG1=" + getShortcutLabel() + "_uforce_" + nargs[i] + 
                                          " ARG2=" + getShortcutLabel() + "_pforce");  
      // Now the quantity in the numerator for the average force
      readInputLine( getShortcutLabel() + "_numer_" + nargs[i] + ": AVERAGE NORMALIZATION=false ARG=" +  getShortcutLabel() + "_mforce_" + nargs[i] + " STRIDE=" + pacestr );
  }
  // And the denominator for the average force
  std::string arg_input="";
  readInputLine( getShortcutLabel() + "_denom: AVERAGE NORMALIZATION=false ARG=" + getShortcutLabel() + "_pforce STRIDE=" + pacestr );
  // Now compute all the average forces
  for(unsigned i=0;i<args.size();++i) {
      std::string num; Tools::convert( i+1, num ); arg_input += "ARG" + num + "=" + getShortcutLabel() + "_aforce_" + nargs[i];
      readInputLine( getShortcutLabel() + "_aforce_" + nargs[i] + ": MATHEVAL PERIODIC=NO FUNC=x/y ARG1=" + getShortcutLabel() + "_numer_" + nargs[i] +
                                                                    " ARG2=" + getShortcutLabel() + "_denom");
  }
  // And now we need something to compute the integral of all these forces
  readInputLine( getShortcutLabel() + ": CUMULATIVE_INTEGRAL " + arg_input );
}

}
}
