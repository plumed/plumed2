/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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

namespace PLMD {
namespace analysis {

//+PLUMEDOC REWEIGHTING WHAM_HISTOGRAM
/*
This can be used to output the a histogram using the weighted histogram technique

\par Examples

*/
//+ENDPLUMEDOC

class WhamHistogram : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  WhamHistogram( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(WhamHistogram,"WHAM_HISTOGRAM")

void WhamHistogram::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the arguments that you would like to make the histogram for");
  keys.add("compulsory","BIAS","*.bias","the value of the biases to use when performing WHAM");
  keys.add("compulsory","TEMP","the temperature at which the simulation was run");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be stored to perform WHAM");
  keys.add("compulsory","GRID_MIN","the minimum to use for the grid");
  keys.add("compulsory","GRID_MAX","the maximum to use for the grid");
  keys.add("compulsory","GRID_BIN","the number of bins to use for the grid");
  keys.add("optional","BANDWIDTH","the bandwidth for kernel density estimation");
}


WhamHistogram::WhamHistogram( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao)
{
  // Input for REWEIGHT_WHAM
  std::string rew_line = getShortcutLabel() + "_weights: REWEIGHT_WHAM";
  std::string bias; parse("BIAS",bias); rew_line += " ARG=" + bias;
  std::string temp; parse("TEMP",temp); rew_line += " TEMP=" + temp;
  readInputLine( rew_line );
  // Input for COLLECT_FRAMES
  std::string col_line = getShortcutLabel() + "_collect: COLLECT_FRAMES LOGWEIGHTS=" + getShortcutLabel() + "_weights";
  std::string stride; parse("STRIDE",stride); col_line += " STRIDE=" + stride;
  std::string arg; parse("ARG",arg); col_line += " ARG=" + arg;
  readInputLine( col_line );
  // Input for HISTOGRAM
  std::string histo_line = getShortcutLabel() + ": HISTOGRAM ARG=" + getShortcutLabel() + "_collect.*";
  std::string min; parse("GRID_MIN",min); histo_line += " GRID_MIN=" + min;
  std::string max; parse("GRID_MAX",max); histo_line += " GRID_MAX=" + max;
  std::string bin; parse("GRID_BIN",bin); histo_line += " GRID_BIN=" + bin;
  std::string bw=""; parse("BANDWIDTH",bw);
  if( bw!="" ) histo_line += " BANDWIDTH=" + bw;
  else histo_line += " KERNEL=DISCRETE";
  readInputLine( histo_line );
}

}
}
