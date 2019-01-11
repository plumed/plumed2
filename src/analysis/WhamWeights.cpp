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

//+PLUMEDOC REWEIGHTING WHAM_WEIGHTS
/*
This can be used to output the data that has been stored in an Analysis object.

\par Examples

*/
//+ENDPLUMEDOC

class WhamWeights : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  WhamWeights( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(WhamWeights,"WHAM_WEIGHTS")

void WhamWeights::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); keys.remove("LABEL");
  keys.add("compulsory","BIAS","*.bias","the value of the biases to use when performing WHAM");
  keys.add("compulsory","TEMP","the temperature at which the simulation was run");
  keys.add("compulsory","STRIDE","1","the frequency with which the bias should be stored to perform WHAM");
  keys.add("compulsory","FILE","the file on which to output the WHAM weights");
  keys.add("optional","FMT","the format to use for the real numbers in the output file");
}

WhamWeights::WhamWeights( const ActionOptions& ao ) :
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
  readInputLine( col_line );
  // Input for line to output data
  std::string out_line="OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=" + getShortcutLabel() + "_collect";
  std::string file; parse("FILE",file); out_line += " FILE=" + file;
  std::string fmt="%f"; parse("FMT",fmt); out_line += " FMT=" + fmt;
  readInputLine( out_line );
}

}
}
