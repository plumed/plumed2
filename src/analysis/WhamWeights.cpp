/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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
Calculate and output weights for configurations using the weighted histogram analysis method.

This shortcut action allows you to calculate and output weights computed using the weighted histogram
analysis technique.  For more detail on how this technique works see \ref REWEIGHT_WHAM

\par Examples

The following input can be used to analyze the output from a series of umbrella sampling calculations.
The trajectory from each of the simulations run with the different biases should be concatenated into a
single trajectory before running the following analysis script on the concatenated trajectory using PLUMED
driver.  The umbrella sampling simulations that will be analyzed using the script below applied a harmonic
restraint that restrained the torsional angle involving atoms 5, 7, 9 and 15 to particular values.  The script
below calculates the reweighting weights for each of the trajectories and then applies the binless WHAM algorithm
to determine a weight for each configuration in the concatenated trajectory.

\plumedfile
#SETTINGS NREPLICAS=4
phi: TORSION ATOMS=5,7,9,15
rp: RESTRAINT ARG=phi KAPPA=50.0 ...
  AT=@replicas:{
        -3.00000000000000000000
        -1.45161290322580645168
        .09677419354838709664
        1.64516129032258064496
     }
...

WHAM_WEIGHTS BIAS=rp.bias TEMP=300 FILE=wham-weights
\endplumedfile

The script above must be run with multiple replicas using the following command:

\verbatim
mpirun -np 4 plumed driver --mf_xtc alltraj.xtc --multi 4
\endverbatim

*/
//+ENDPLUMEDOC

class WhamWeights : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit WhamWeights( const ActionOptions& );
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
