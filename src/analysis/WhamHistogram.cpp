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

This shortcut action allows you to calculate a histogram using the weighted histogram
analysis technique.  For more detail on how this the weights for configurations are
computed see \ref REWEIGHT_WHAM

\par Examples

The following input can be used to analyze the output from a series of umbrella sampling calculations.
The trajectory from each of the simulations run with the different biases should be concatenated into a
single trajectory before running the following analysis script on the concatenated trajectory using PLUMED
driver.  The umbrella sampling simulations that will be analyzed using the script below applied a harmonic
restraint that restrained the torsional angle involving atoms 5, 7, 9 and 15 to particular values.  The script
below calculates the reweighting weights for each of the trajectories and then applies the binless WHAM algorithm
to determine a weight for each configuration in the concatenated trajectory.  A histogram is then constructed from
the configurations visited and their weights.  This histogram is then converted into a free energy surface and output
to a file called fes.dat

\plumedfile
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
rp: RESTRAINT ARG=phi KAPPA=50.0 ...
  AT=@replicas:{
        -3.00000000000000000000
        -2.22580645161290322584
        -1.45161290322580645168
        -.67741935483870967752
        .09677419354838709664
        .87096774193548387080
        1.64516129032258064496
        2.41935483870967741912
     }
...

hh: WHAM_HISTOGRAM ARG=phi BIAS=rp.bias TEMP=300 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=50
fes: CONVERT_TO_FES GRID=hh TEMP=300
DUMPGRID GRID=fes FILE=fes.dat
\endplumedfile

The script above must be run with multiple replicas using the following command:

\verbatim
mpirun -np 6 plumed driver --mf_xtc alltraj.xtc --multi 6
\endverbatim

*/
//+ENDPLUMEDOC

class WhamHistogram : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit WhamHistogram( const ActionOptions& );
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
