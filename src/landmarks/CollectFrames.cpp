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
Collect atomic positions or argument values from the trajectory for later analysis

This shortcut uses [COLLECT](COLLECT.md) actions to collect data from the trajectory that is amenable for later analysis
using PLUMED's landmark selection actions or dimensionality reduction methods.  You can use this method to collect atomic position
data as shown in the following example:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL

# This should output the atomic positions for the frames that were collected to a pdb file called traj.pdb
DUMPPDB ATOMS=cc_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb
```

If you look at the expanded version of the shortcut above you can see how the position data is stored in a matrix. It is also worth
noting that all the stored structured are aligned to the first frame before being stored. When we take the difference between the
stored configurations we are thus excluding any rotation or the reference frame or translation of the center of mass that has taken place.

Instead of storing the positions of the atoms you you can store a time series of argument values as shown below:

```plumed
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=5,6,7,8
t3: TORSION ATOMS=9,10,11,12

# This collects the three torsion angles
cc: COLLECT_FRAMES ARG=t1,t2,t3 STRIDE=10 CLEAR=1000
DUMPVECTOR ARG=cc.* FILE=timeseries STRIDE=1000
```

Notice that the stored data is in a matrix once again here.  Further note that we are storing the three torsions on every tenth step.  In addition we also delete all the stored data on every 1000th step of trajectory.
The time series for each consecutive 1000 step segment of trajectory will be output to separate files in the above input.  These output files also contain weights for each of the stored frames.  For the above input
these weights are all equal to one.  We can, however, pass the weight as an input argument.  This functionality is useful if a bias is acting upon the system and we want to keep track of how much bias is acting upon
each frame of the traejctory as has been done in the input below.

```plumed
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=5,6,7,8
t3: TORSION ATOMS=9,10,11,12

r: RESTRAINT ARG=t1 AT=pi/2 KAPPA=10
bw: REWEIGHT_BIAS TEMP=300

# This collects the three torsion angles
cc: COLLECT_FRAMES ARG=t1,t2,t3 LOGWEIGHTS=bw
DUMPVECTOR ARG=cc.* FILE=timeseries
```

You can learn how to use COLLECT_FRAMES to perform dimensionality reduction calculations by working through [this tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html)

## Calculating dissimilary matrices

One of the simplest things that we can do if we have stored a set of trajectory frames using this action is we can calculate the dissimilarity between
every pair of frames we have stored.  This is the first step for many of the dimensionality reduction algorithms described in [the dimred module](module_dimred.md).
To calculate these dissimilarities you use the [DISSIMILARITIES](DISSIMILARITIES.md) command.

Notice that instead of reading in a trajectory and computing dissimilarities you can read in the dissimilarity matrix directly using [CONSTANT](CONSTANT.md) and then perform all analyses with PLUMED.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace landmarks {

class CollectFrames : public ActionShortcut {
private:
  std::string fixArgumentName( const std::string& argin );
public:
  static void registerKeywords( Keywords& keys );
  explicit CollectFrames( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(CollectFrames,"COLLECT_FRAMES")

void CollectFrames::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("optional","ARG","scalar/vector","the labels of the values whose time series you would like to collect for later analysis");
  keys.add("compulsory","STRIDE","1","the frequency with which data should be stored for analysis.  By default data is collected on every step");
  keys.add("compulsory","CLEAR","0","the frequency with which data should all be deleted and restarted");
  keys.add("compulsory","ALIGN","OPTIMAL","if storing atoms how would you like the alignment to be done can be SIMPLE/OPTIMAL");
  keys.add("optional","ATOMS","list of atomic positions that you would like to collect and store for later analysis");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
  keys.addOutputComponent("data","default","matrix","the data that is being collected by this action");
  keys.addOutputComponent("logweights","default","vector","the logarithms of the weights of the data points");
  keys.needsAction("POSITION");
  keys.needsAction("CONCATENATE");
  keys.needsAction("MEAN");
  keys.needsAction("CUSTOM");
  keys.needsAction("CONCATENATE");
  keys.needsAction("COLLECT");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("RMSD_VECTOR");
  keys.needsAction("COMBINE");
  keys.needsAction("VSTACK");
  keys.needsAction("CONSTANT");
}

std::string CollectFrames::fixArgumentName( const std::string& argin ) {
  std::string argout = argin;
  std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) {
    argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  }
  return argout;
}

CollectFrames::CollectFrames( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  std::string stride, clearstride;
  parse("STRIDE",stride);
  parse("CLEAR",clearstride);
  std::vector<std::string> argn;
  parseVector("ARG",argn);
  std::vector<Value*> theargs;
  ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, theargs );
  std::string indices;
  parse("ATOMS",indices);
  if( theargs.size()==0 && indices.length()==0 ) {
    error("no arguments or atoms were specified for collection");
  }

  // Create the values to collect the atomic positions
  if( indices.length()>0 ) {
    // Collect reference position
    readInputLine( getShortcutLabel() + "_getposx: POSITION ATOMS=" + indices );
    std::string align;
    parse("ALIGN",align);
    readInputLine( getShortcutLabel() + "_getpos: CONCATENATE ARG=" + getShortcutLabel() + "_getposx.x," + getShortcutLabel() + "_getposx.y," + getShortcutLabel() + "_getposx.z");
    // Find atomic center
    readInputLine( getShortcutLabel() + "_cposx: MEAN ARG=" + getShortcutLabel() + "_getposx.x PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cposy: MEAN ARG=" + getShortcutLabel() + "_getposx.y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cposz: MEAN ARG=" + getShortcutLabel() + "_getposx.z PERIODIC=NO");
    // Subtract atomimc center
    readInputLine( getShortcutLabel() + "_refx: CUSTOM ARG=" + getShortcutLabel() + "_getposx.x," + getShortcutLabel() + "_cposx FUNC=x-y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_refy: CUSTOM ARG=" + getShortcutLabel() + "_getposx.y," + getShortcutLabel() + "_cposy FUNC=x-y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_refz: CUSTOM ARG=" + getShortcutLabel() + "_getposx.z," + getShortcutLabel() + "_cposz FUNC=x-y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_ref: CONCATENATE ARG=" + getShortcutLabel() + "_refx," + getShortcutLabel() + "_refy," + getShortcutLabel() + "_refz");
    // Store the reference position in a collect action
    readInputLine( getShortcutLabel() + "_refpos: COLLECT TYPE=matrix ARG=" + getShortcutLabel() + "_ref STRIDE=" + clearstride + " CLEAR=" + clearstride );
    readInputLine( getShortcutLabel() + "_refposT: TRANSPOSE ARG=" + getShortcutLabel() + "_refpos");
    // Calculate the RMSD between the instaneous position and the reference position
    readInputLine( getShortcutLabel() + "_rmsd: RMSD_VECTOR ARG=" + getShortcutLabel() + "_getpos," + getShortcutLabel() + "_refpos DISPLACEMENT SQUARED TYPE=" + align );
    // Add the reference position to the RMSD displacement
    readInputLine( getShortcutLabel() + "_fpos: COMBINE ARG=" + getShortcutLabel() + "_refposT," + getShortcutLabel() + "_rmsd.disp PERIODIC=NO");
    // Store the reference data
    std::string suffix = "_atomdata";
    if( theargs.size()==0 ) {
      suffix = "_data";
    }
    readInputLine( getShortcutLabel() + suffix + ": COLLECT TYPE=matrix ARG=" + getShortcutLabel() + "_fpos STRIDE=" + stride + " CLEAR=" + clearstride );
  }

  // Create all the collect actions for arguments
  for(unsigned i=0; i<theargs.size(); ++i) {
    if( theargs[i]->getNumberOfValues()!=theargs[0]->getNumberOfValues() ) {
      error("mismatch between number of arguments calculated by each collected argument");
    }
    readInputLine( getShortcutLabel() + "_" + fixArgumentName( theargs[i]->getName() ) + ": COLLECT ARG=" + theargs[i]->getName() + " STRIDE=" + stride + " CLEAR=" + clearstride );
  }
  // Make a list of collect actions
  if( theargs.size()>0 ) {
    std::string allcol = getShortcutLabel() + "_" + fixArgumentName( theargs[0]->getName() );
    for(unsigned i=1; i<theargs.size(); ++i) {
      allcol += "," + getShortcutLabel() + "_" + fixArgumentName( theargs[i]->getName() );
    }
    // And transfer everything to a matrix
    std::string suffix = "_argdata";
    if( indices.length()==0 ) {
      suffix = "_data";
    }
    readInputLine( getShortcutLabel() + suffix + ": VSTACK ARG=" + allcol );
  }
  // Merge all the collected data together into a single matrix
  if( theargs.size()>0 && indices.length()>0 ) {
    readInputLine( getShortcutLabel() + "_data: CONCATENATE MATRIX11=" + getShortcutLabel() + "_atomdata MATRIX12=" + getShortcutLabel() + "_argdata");
  }

  // Now get the logweights
  std::vector<std::string> logw;
  parseVector("LOGWEIGHTS",logw);
  std::vector<Value*> thew;
  if( logw.size()>0 ) {
    ActionWithArguments::interpretArgumentList( logw, plumed.getActionSet(), this, thew );
  }
  if( logw.size()>1 ) {
    error("maximum of one argument should be specified for logweights");
  }

  if( logw.size()==0 ) {
    std::string zeros="0";
    if( theargs.size()>0 ) {
      for(unsigned i=1; i<theargs[0]->getNumberOfValues(); ++i) {
        zeros += ",0";
      }
    }
    readInputLine( getShortcutLabel() + "_cweight: CONSTANT VALUE=" + zeros );
    readInputLine( getShortcutLabel() + "_logweights: COLLECT ARG=" + getShortcutLabel() + "_cweight STRIDE=" + stride + " CLEAR=" + clearstride );
  } else {
    if( theargs[0]->getNumberOfValues()!=thew[0]->getNumberOfValues() ) {
      error("mismatch between number of weights and number of collected arguments");
    }
    readInputLine( getShortcutLabel() + "_logweights: COLLECT ARG=" + thew[0]->getName() + " STRIDE=" + stride + " CLEAR=" + clearstride );
  }
  // And finally create a value that contains as many ones as there are data points (this is used if we want to do Classical MDS
  readInputLine( getShortcutLabel() + "_one: CONSTANT VALUE=1");
  readInputLine( getShortcutLabel() + "_ones: COLLECT ARG=" + getShortcutLabel() + "_one STRIDE=" + stride + " CLEAR=" + clearstride );
}

}
}
