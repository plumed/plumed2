/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2023 The plumed team
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
namespace wham {

//+PLUMEDOC REWEIGHTING WHAM_WEIGHTS
/*
Calculate and output weights for configurations using the weighted histogram analysis method.

This shortcut action allows you to calculate and output weights computed using the weighted histogram
analysis technique.  The following input demonstrates how this works in practise by showing you how this
action can be used to analyze the output from a series of umbrella sampling calculations.
The trajectory from each of the simulations run with the different biases should be concatenated into a
single trajectory before running the following analysis script on the concatenated trajectory using PLUMED
driver.  The umbrella sampling simulations that will be analyzed using the script below applied a harmonic
restraint that restrained the torsional angle involving atoms 5, 7, 9 and 15 to particular values.  The script
below calculates the reweighting weights for each of the trajectories and then applies the binless WHAM algorithm
to determine a weight for each configuration in the concatenated trajectory.

```plumed
#SETTINGS NREPLICAS=4
phi: TORSION ATOMS=5,7,9,15
rp: RESTRAINT ARG=phi KAPPA=50.0 AT=@replicas:{-3.00,-1.45,0.10,1.65}
WHAM_WEIGHTS BIAS=rp.bias TEMP=300 FILE=wham-weights
```

The script above must be run with multiple replicas using the following command:

````
mpirun -np 4 plumed driver --mf_xtc alltraj.xtc --multi 4
````

For more details on how the weights for configurations are determined using the wham method see the documentation
for the [WHAM](WHAM.md) action.

*/
//+ENDPLUMEDOC

class WhamWeights : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit WhamWeights( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(WhamWeights,"WHAM_WEIGHTS")

void WhamWeights::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.remove("LABEL");
  keys.add("compulsory","BIAS","*.bias","the value of the biases to use when performing WHAM");
  keys.add("optional","TEMP","the temperature at which the simulation was run");
  keys.add("compulsory","STRIDE","1","the frequency with which the bias should be stored to perform WHAM");
  keys.add("compulsory","FILE","the file on which to output the WHAM weights");
  keys.add("optional","FMT","the format to use for the real numbers in the output file");
  keys.setValueDescription("vector","the weights that were calculated using WHAM");
  keys.needsAction("GATHER_REPLICAS");
  keys.needsAction("CONCATENATE");
  keys.needsAction("COLLECT");
  keys.needsAction("WHAM");
  keys.needsAction("DUMPVECTOR");
}

WhamWeights::WhamWeights( const ActionOptions& ao ) :
  Action(ao),
  ActionShortcut(ao) {
  // Input for collection of weights for WHAM
  std::string bias;
  parse("BIAS",bias);
  std::string stride;
  parse("STRIDE",stride);
  // Input for GATHER_REPLICAS
  readInputLine( getShortcutLabel() + "_gather: GATHER_REPLICAS ARG=" + bias );
  // Put all the replicas in a single vector
  readInputLine( getShortcutLabel() + "_gatherv: CONCATENATE ARG=" + getShortcutLabel() + "_gather.*");
  // Input for COLLECT_FRAMES
  readInputLine( getShortcutLabel() + "_collect: COLLECT TYPE=vector ARG=" + getShortcutLabel() + "_gatherv STRIDE=" + stride);
  // Input for WHAM
  std::string temp, tempstr="";
  parse("TEMP",temp);
  if( temp.length()>0 ) {
    tempstr="TEMP=" + temp;
  }
  readInputLine( getShortcutLabel() + ": WHAM ARG=" + getShortcutLabel() + "_collect " + tempstr );
  // Input for PRINT (will just output at end of calc
  std::string filename, fmt;
  parse("FILE",filename);
  parse("FMT",fmt);
  if(fmt.length()>0) {
    fmt=" FMT=" + fmt;
  }
  readInputLine( "DUMPVECTOR STRIDE=0 ARG=" + getShortcutLabel() + " FILE=" + filename + fmt );
}

}
}
