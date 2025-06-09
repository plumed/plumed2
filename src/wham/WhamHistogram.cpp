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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace wham {

//+PLUMEDOC REWEIGHTING WHAM_HISTOGRAM
/*
This can be used to output the a histogram using the weighted histogram technique

This shortcut action allows you to calculate a histogram using the weighted histogram analysis technique.
The following input illustrates how this is used in practise to analyze the output from a series of umbrella sampling calculations.
The trajectory from each of the simulations run with the different biases should be concatenated into a
single trajectory before running the following analysis script on the concatenated trajectory using PLUMED
driver.  The umbrella sampling simulations that will be analyzed using the script below applied a harmonic
restraint that restrained the torsional angle involving atoms 5, 7, 9 and 15 to particular values.  The script
below calculates the reweighting weights for each of the trajectories and then applies the binless WHAM algorithm
to determine a weight for each configuration in the concatenated trajectory.  A histogram is then constructed from
the configurations visited and their weights.  This histogram is then converted into a free energy surface and output
to a file called fes.dat

```plumed
#SETTINGS NREPLICAS=4
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
rp: RESTRAINT ARG=phi KAPPA=50.0 AT=@replicas:{-3.00,-1.45,0.10,1.65}
hh: WHAM_HISTOGRAM ARG=phi BIAS=rp.bias TEMP=300 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=50
fes: CONVERT_TO_FES ARG=hh TEMP=300
DUMPGRID ARG=fes FILE=fes.dat
```

The script above must be run with multiple replicas using the following command:

````
mpirun -np 4 plumed driver --mf_xtc alltraj.xtc --multi 4
````

Notice that if you use the BANDWIDTH keyword, as in the example below, PLUMED will estimate the histogram
using [kernel density estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation).

```plumed
#SETTINGS NREPLICAS=4
phi: TORSION ATOMS=5,7,9,15
psi: TORSION ATOMS=7,9,15,17
rp: RESTRAINT ARG=phi KAPPA=50.0 AT=@replicas:{-3.00,-1.45,0.10,1.65}
hh: WHAM_HISTOGRAM ARG=phi BIAS=rp.bias TEMP=300 GRID_MIN=-pi GRID_MAX=pi GRID_BIN=50 BANDWIDTH=0.1
fes: CONVERT_TO_FES ARG=hh TEMP=300
DUMPGRID ARG=fes FILE=fes.dat
```


For more details on how the weights for configurations are determined using the wham method see the documentation
for the [WHAM](WHAM.md) action.

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
  keys.setValueDescription("grid","the histogram that was generated using the WHAM weights");
  keys.needsAction("GATHER_REPLICAS");
  keys.needsAction("CONCATENATE");
  keys.needsAction("COLLECT");
  keys.needsAction("WHAM");
  keys.needsAction("KDE");
}


WhamHistogram::WhamHistogram( const ActionOptions& ao ) :
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
  readInputLine( getShortcutLabel() + "_wham: WHAM ARG=" + getShortcutLabel() + "_collect " + tempstr );
  // Input for COLLECT_FRAMES
  std::vector<std::string> args;
  parseVector("ARG",args);
  std::string argstr;
  for(unsigned i=0; i<args.size(); ++i) {
    readInputLine( getShortcutLabel() + "_data_" + args[i] + ": COLLECT ARG=" + args[i] );
    if( i==0 ) {
      argstr = " ARG=";
    } else {
      argstr += ",";
    }
    argstr += getShortcutLabel() + "_data_" + args[i];
  }
  // Input for HISTOGRAM
  std::string histo_line, bw="";
  parse("BANDWIDTH",bw);
  if( bw!="" ) {
    histo_line += " BANDWIDTH=" + bw;
  } else {
    histo_line += " KERNEL=DISCRETE";
  }
  std::string min;
  parse("GRID_MIN",min);
  histo_line += " GRID_MIN=" + min;
  std::string max;
  parse("GRID_MAX",max);
  histo_line += " GRID_MAX=" + max;
  std::string bin;
  parse("GRID_BIN",bin);
  histo_line += " GRID_BIN=" + bin;
  readInputLine( getShortcutLabel() + ": KDE " + argstr + " HEIGHTS=" + getShortcutLabel() + "_wham" + histo_line );
}

}
}
