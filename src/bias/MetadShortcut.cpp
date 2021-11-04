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
#include "MetadShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "tools/IFile.h"
#include "gridtools/KDEShortcut.h"
#include "core/ReweightBase.h"

namespace PLMD {
namespace bias {

PLUMED_REGISTER_ACTION(MetadShortcut,"NEW_METAD")

void MetadShortcut::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that should be used as input to this method");
  keys.add("compulsory","PACE","the frequency with which hills should be added");
  keys.add("compulsory","SIGMA","the widths of the Gaussians that we are using");
  keys.add("compulsory","HEIGHT","the heights of the hills that should be added");
  keys.add("optional","BIASFACTOR","the bias factor for the well tempered metadynamics");
  keys.add("optional","GRID_MIN","the minimum value to use for the the grid");
  keys.add("optional","GRID_MAX","the maximum value to use for all the grids");
  keys.add("optional","GRID_BIN","the number of bins to use for all the grids");
  keys.add("optional","TEMP","the system temperature - this is only needed if you are doing well-tempered metadynamics");
}

MetadShortcut::MetadShortcut(const ActionOptions& ao):
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
  parseVector("GRID_MIN",gmin); parseVector("GRID_MAX",gmax); parseVector("GRID_BIN",grid_nbins); parseVector("SIGMA",sigma);
  // Convert the bandwidth to something constant actions
  gridtools::KDEShortcut::convertBandwiths( getShortcutLabel(), sigma, this );
  // Create a metadynamics bias
  createMetadBias( getShortcutLabel(), pacestr, args, gmin, gmax, grid_nbins, getShortcutLabel() + "_wtfact", "", "", this );
  // Bias the simulation using this bias potentital
  readInputLine("BIASVALUE ARG=" + getShortcutLabel() + "_bias" );
  // Complete setup of the well tempered weights
  if( biasfactor.length()>0 ) {
      args.resize(1); args[0] = getShortcutLabel() + "_bias";
      ReweightBase* rwb = plumed.getActionSet().selectWithLabel<ReweightBase*>( getShortcutLabel() + "_wtfact" );
      plumed_assert( rwb ); rwb->setArguments( args );
  }
}

void MetadShortcut::createMetadBias( const std::string& lab, const std::string& pacestr, const std::vector<std::string>& args, 
                                     const std::vector<std::string>& gmin, const std::vector<std::string>& gmax, const std::vector<std::string>& grid_nbins,
                                     const std::string& weight_str, const std::string& truncflag1, const std::string& truncflag2, ActionShortcut* act ) {
  if( gmin.size()>0 ) {
      // Now create the height
      act->readInputLine( lab + "_height: CONSTANT VALUE=1.0");
      std::string input = lab + "_kde: KDE_CALC METRIC=" + lab + "_icov ARG1=" + args[0] + " HEIGHTS=" + lab + "_height";
      std::string gminstr=" GRID_MIN=" + gmin[0]; std::string gmaxstr=" GRID_MAX=" + gmax[0]; std::string gbinstr=" GRID_BIN=" + grid_nbins[0];
      for(unsigned i=1;i<args.size();++i) { 
        std::string num; Tools::convert( i+1, num ); input += " ARG" + num + "=" + args[i]; 
        gminstr += "," + gmin[i]; gmaxstr += "," + gmax[i]; gbinstr += "," + grid_nbins[i]; 

      }
      act->readInputLine( input + " " + gminstr + " " + gmaxstr + " " + " " + gbinstr + " " + truncflag1 );
      act->readInputLine( lab + "_grid: AVERAGE ARG=" + lab + "_kde NORMALIZATION=false STRIDE=" + pacestr + " LOGWEIGHTS=" + weight_str );
      // Evaluate the instantaneous value of the bias potential 
      act->readInputLine( lab + "_bias: EVALUATE_FUNCTION_FROM_GRID ARG=" + lab + "_grid " + truncflag2 );
  } else {
      if( args.size()>1 ) {
          for(unsigned i=0; i<args.size();++i) {
              std::string num; Tools::convert( i+1, num ); std::size_t com=args[i].find_first_of(","); std::string thisarg=args[i].substr(0,com); 
              act->readInputLine( lab + "_sigma_" + thisarg + ": SELECT_COMPONENTS ARG=" + lab + "_sigma COMPONENTS=" + num );
          }
      }
      std::string store_args; for(unsigned i=0;i<args.size();++i) { std::string num; Tools::convert( i+1, num ); store_args += " ARG" + num + "=" + args[i]; }
      // Create a store to hold the list of Gaussians
      act->readInputLine( lab + "_store: COLLECT_FRAMES STRIDE=" + pacestr + store_args + " LOGWEIGHTS=" + weight_str );
      std::string names=" VAR=v1", func="v1*v1", func_args;
      for(unsigned i=1;i<args.size();++i) {  std::string num; Tools::convert( i+1, num ); names += ",v" + num; func += "+v" + num + "*v" + num; }
      for(unsigned i=0;i<args.size();++i) {
          std::string num; Tools::convert( i+1, num ); std::size_t com=args[i].find_first_of(","); std::string thisarg=args[i].substr(0,com);
          act->readInputLine( lab + "_sub_" + thisarg + ": DIFFERENCE MIX_HISTORY_DEPENDENCE ARG2=" + thisarg + " ARG1=" + lab + "_store." + thisarg ); 
          if( args.size()==1 ) act->readInputLine( lab  + "_scaled_" + thisarg + ": CUSTOM FUNC=x/y PERIODIC=NO ARG1=" + lab + "_sub_" + thisarg + " ARG2=" + lab + "_sigma");
          else act->readInputLine( lab  + "_scaled_" + thisarg + ": CUSTOM FUNC=x/y PERIODIC=NO ARG1=" + lab + "_sub_" + thisarg + " ARG2=" + lab + "_sigma_" + thisarg );
          func_args += " ARG" + num + "=" + lab + "_scaled_" + thisarg;
      }
      // Calculate the potential
      std::string num; Tools::convert( args.size()+1, num );
      act->readInputLine( lab + "_heights: CUSTOM PERIODIC=NO FUNC=exp(x) ARG1=" + lab + "_store.logweights" );
      act->readInputLine( lab + "_gvals: CUSTOM PERIODIC=NO FUNC=vh*exp(-(" + func + ")/2)" + names + ",vh " + func_args + " ARG" + num + "=" + lab  + "_heights");
      act->readInputLine( lab + "_bias: SUM ARG=" + lab + "_gvals PERIODIC=NO");
  }
}

}
}
