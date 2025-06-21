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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_STRIDE
/*
Select every ith frame from the stored set of configurations

If you have collected a set of trajectory frames using [COLLECT_FRAMES](COLLECT_FRAMES.md) you can use this action to
select a subset you have collected.  This particular method for landmark selection reduces the number of frames by selecting taking every
$i$th frame. So, for example, if you use the input below every 10th frame of the stored trajectory is transferred to the `ll_data` Value
that is output which is output in the PDB file.  This happens because we are collecting 1000 trajectory frames in total but only taking
100 landmarks from this data.

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_STRIDE ARG=cc NLANDMARKS=100

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If you expand the shortcuts in the input above you will notice that the LANDMARK_SELECT_STRIDE shortcut creates a [DISSIMILARITIES](DISSIMILARITIES.md) action
that calculates the distances between the input frames. We need to calculate these dissimilarities here because the LANDMARK_SELECT_STRIDE shortcut computes the
weights of the landmarks by doing a [VORONOI](VORONOI.md) analysis.  If you would like to turn this and the computing of dissimilarities off you can use the
NODISSIMILARITIES flag as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_STRIDE ARG=cc NLANDMARKS=100 NODISSIMILARITIES

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If the NODISSIMILARITIES flag is __not__ present the weights of the landmark points are calculated by using VORONOI procedure as is illustrated in the expanded
version of the first example input above.  If you wish to turn this off you can use the NOVORONOI flag as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_STRIDE ARG=cc NLANDMARKS=100 NOVORONOI

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

Be aware, however, that dissimilarities are still computed if you only the the NOVORONOI flag.  If you use NODISSIMILARITIES the voronoi weights are not computed
and the dissimilarities are not computed.

If you have already computed the dissimilarities between the collected frames you can pass them in input to the LANDMARK_SELECT_STRIDE funtion as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# This calculates the dissimilarities between the stored frames
cc_dataT: TRANSPOSE ARG=cc_data
dd: DISSIMILARITIES ARG=cc_data,cc_dataT

# Select landmarks
ll: LANDMARK_SELECT_STRIDE ARG=cc DISSIMILARITIES=dd NLANDMARKS=100

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

Notice that you can also read in dissimilarities from a file using a [CONSTANT](CONSTANT.md) action and pass these directly to the LANDMARK_SELECT_STRIDE action and avoid using COLLECT_FRAMES.

You can learn how to use landmark selection for dimensionality reduction calculations by working through [this tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html)

*/
//+ENDPLUMEDOC

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_RANDOM
/*
Select a random set of landmarks from a large set of configurations.

If you have collected a set of trajectory frames using [COLLECT_FRAMES](COLLECT_FRAMES.md) you can use this action to
select a subset of the configurations you have collected.  This particular method for landmark selection reduces the number of frames by
chooseing NLANDMARKS points from the data collected by COLLECT_FRAMES at random by using a quasi random number generator for which you may way to specify a seed using
the SEED keyword.  So, for example, if you use the input 100 randomly-selected
points from the 1000 trajectory frames that were by collected by the COLLECT_FRAMES action are transferred to the `ll_data` Value that is output which is output in the PDB file.

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_RANDOM ARG=cc SEED=23 NLANDMARKS=100

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If you expand the shortcuts in the input above you will notice that the LANDMARK_SELECT_RANDOM shortcut creates a [DISSIMILARITIES](DISSIMILARITIES.md) action
that calculates the distances between the input frames. We need to calculate these dissimilarities here because the LANDMARK_SELECT_RANDOM shortcut computes the
weights of the landmarks by doing a [VORONOI](VORONOI.md) analysis.  If you would like to turn this and the computing of dissimilarities off you can use the
NODISSIMILARITIES flag as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_RANDOM ARG=cc NLANDMARKS=100 NODISSIMILARITIES

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If the NODISSIMILARITIES flag is __not__ present the weights of the landmark points are calculated by using VORONOI procedure as is illustrated in the expanded
version of the first example input above.  If you wish to turn this off you can use the NOVORONOI flag as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_RANDOM ARG=cc NLANDMARKS=100 NOVORONOI

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

Be aware, however, that dissimilarities are still computed if you only the the NOVORONOI flag.  If you use NODISSIMILARITIES the voronoi weights are not computed
and the dissimilarities are not computed.

If you have already computed the dissimilarities between the collected frames you can pass them in input to the LANDMARK_SELECT_RANDOM funtion as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# This calculates the dissimilarities between the stored frames
cc_dataT: TRANSPOSE ARG=cc_data
dd: DISSIMILARITIES ARG=cc_data,cc_dataT

# Select landmarks
ll: LANDMARK_SELECT_RANDOM ARG=cc DISSIMILARITIES=dd NLANDMARKS=100

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

Notice that you can also read in dissimilarities from a file using a [CONSTANT](CONSTANT.md) action and pass these directly to the LANDMARK_SELECT_RANDOM shortcut and avoid using COLLECT_FRAMES.

You can learn how to use landmark selection for dimensionality reduction calculations by working through [this tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html)

*/
//+ENDPLUMEDOC

//+PLUMEDOC LANDMARKS LANDMARK_SELECT_FPS
/*
Select a of landmarks from a large set of configurations using farthest point sampling.

If you have collected a set of trajectory frames using [COLLECT_FRAMES](COLLECT_FRAMES.md) you can use this action to
select a subset of the configurations you have collected. This shortcut does this using [FARTHEST_POINT_SAMPLING](FARTHEST_POINT_SAMPLING.md)
the first point is thus selected at random, which is why you may want to set a random SEED using the `SEED` keyword.  The remaining points are then selected by taking the unselected point in the input data set that is the furthest
from all the points that have been selected thus far.  The following input demonstrates how you can use this method:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_FPS ARG=cc NLANDMARKS=100 SEED=23

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If you expand the shortcuts in the input above you will notice that the LANDMARK_SELECT_RANDOM shortcut creates a [DISSIMILARITIES](DISSIMILARITIES.md) action
that calculates the distances between the input frames. We have to compute these dissimilarities in order to perform the farthest point sampling here so you cannot use the
NODISSIMILARITIES flag with this action.  However, we also need the dissimilarities to compute the weights of the landmarks as this is done by performing a [VORONOI](VORONOI.md) analysis.
If you would like to turn off the computation of the VORONOI weights you can use the NOVORONOI flag as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# Select landmarks
ll: LANDMARK_SELECT_FPS ARG=cc NLANDMARKS=100 NOVORONOI

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

If you have already computed the dissimilarities between the collected frames you can pass them in input to the LANDMARK_SELECT_FPS funtion as shown below:

```plumed
# This stores the positions of all the first 10 atoms in the system for later analysis
cc: COLLECT_FRAMES ATOMS=1,2,3,4,5,6,7,8,9,10 ALIGN=OPTIMAL STRIDE=1 CLEAR=1000

# This calculates the dissimilarities between the stored frames
cc_dataT: TRANSPOSE ARG=cc_data
dd: DISSIMILARITIES ARG=cc_data,cc_dataT

# Select landmarks
ll: LANDMARK_SELECT_FPS ARG=cc DISSIMILARITIES=dd NLANDMARKS=100

# Output the data to a file
DUMPPDB ATOMS=ll_data ATOM_INDICES=1,2,3,4,5,6,7,8,9,10 FILE=traj.pdb STRIDE=1000
```

Notice that you can also read in dissimilarities from a file using a [CONSTANT](CONSTANT.md) action and pass these directly to the LANDMARK_SELECT_FPS shortcut and avoid using COLLECT_FRAMES.

You can learn how to use landmark selection for dimensionality reduction calculations by working through [this tutorial](https://www.plumed-tutorials.org/lessons/21/006/data/DIMENSIONALITY.html)

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace landmarks {

class LandmarkSelection : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LandmarkSelection( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(LandmarkSelection,"LANDMARK_SELECT_STRIDE")
PLUMED_REGISTER_ACTION(LandmarkSelection,"LANDMARK_SELECT_RANDOM")
PLUMED_REGISTER_ACTION(LandmarkSelection,"LANDMARK_SELECT_FPS")

void LandmarkSelection::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","ARG","the COLLECT_FRAMES action that you used to get the data");
  keys.add("optional","DISSIMILARITIES","the matrix of dissimilarities if this is not provided the squared dissimilarities are calculated");
  keys.add("compulsory","NLANDMARKS","the numbe rof landmarks you would like to create");
  if( keys.getDisplayName()!="LANDMARK_SELECT_STRIDE" ) {
    keys.add("optional","SEED","a random number seed");
  }
  keys.addFlag("NOVORONOI",false,"do not do a Voronoi analysis of the data to determine weights of final points");
  if( keys.getDisplayName()!="LANDMARK_SELECT_FPS" ) {
    keys.addFlag("NODISSIMILARITIES",false,"do not calculate the dissimilarities");
  }
  keys.addOutputComponent("data","ARG","matrix","the data that is being collected by this action");
  keys.addOutputComponent("logweights","ARG","vector","the logarithms of the weights of the data points");
  keys.addOutputComponent("rectdissims","DISSIMILARITIES","matrix","a rectangular matrix containing the distances between the landmark points and the rest of the points");
  keys.addOutputComponent("sqrdissims","DISSIMILARITIES","matrix","a square matrix containing the distances between each pair of landmark points");
  keys.needsAction("LOGSUMEXP");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("DISSIMILARITIES");
  keys.needsAction("ONES");
  keys.needsAction("CREATE_MASK");
  keys.needsAction("FARTHEST_POINT_SAMPLING");
  keys.needsAction("SELECT_WITH_MASK");
  keys.needsAction("COMBINE");
  keys.needsAction("VORONOI");
  keys.needsAction("MATRIX_PRODUCT");
  keys.needsAction("CUSTOM");
}

LandmarkSelection::LandmarkSelection( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  std::string nlandmarks;
  parse("NLANDMARKS",nlandmarks);
  bool novoronoi;
  parseFlag("NOVORONOI",novoronoi);

  bool nodissims=false;
  if( keywords.exists("NODISSIMILARITIES") ) {
    parseFlag("NODISSIMILARITIES",nodissims);
  }
  std::string argn, dissims;
  parse("ARG",argn);
  parse("DISSIMILARITIES",dissims);
  if( argn.length()>0 ) {
    ActionShortcut* as = plumed.getActionSet().getShortcutActionWithLabel( argn );
    if( !as || as->getName()!="COLLECT_FRAMES" ) {
      error("found no COLLECT_FRAMES action with label " + argn );
    }
    // Get the weights
    readInputLine( getShortcutLabel() + "_allweights: LOGSUMEXP ARG=" + argn + "_logweights");
  }
  if( dissims.length()>0 ) {
    ActionWithValue* ds = plumed.getActionSet().selectWithLabel<ActionWithValue*>( dissims );
    if( (ds->copyOutput(0))->getRank()!=2 ) {
      error("input for dissimilarities shoudl be a matrix");
    }
    // Calculate the dissimilarities if the user didn't specify them
  } else if( !nodissims ) {
    readInputLine( getShortcutLabel() + "_" + argn + "_dataT: TRANSPOSE ARG=" + argn + "_data");
    dissims = getShortcutLabel() + "_dissims";
    readInputLine( getShortcutLabel() + "_dissims: DISSIMILARITIES SQUARED ARG=" + argn + "_data," + getShortcutLabel() + "_" + argn + "_dataT");
  }
  // This deals with a corner case whereby users have a matrix of dissimilarities but no corresponding coordinates for these frames
  if( argn.length()==0 && dissims.size()>0 ) {
    ActionWithValue* ds = plumed.getActionSet().selectWithLabel<ActionWithValue*>( dissims );
    if( ds->getName()!="CONSTANT" || (ds->copyOutput(0))->getRank()!=2 ) {
      error("set ARG as well as DISSIMILARITIES");
    }
    std::string size;
    Tools::convert(  (ds->copyOutput(0))->getShape()[0], size );
    readInputLine( getShortcutLabel() + "_allweights: ONES SIZE=" + size );
  }

  if( getName()=="LANDMARK_SELECT_STRIDE" ) {
    readInputLine( getShortcutLabel() + "_mask: CREATE_MASK ARG=" + getShortcutLabel() + "_allweights TYPE=stride NZEROS=" + nlandmarks );
  } else if( getName()=="LANDMARK_SELECT_RANDOM" ) {
    if( argn.length()==0 ) {
      error("must set COLLECT_FRAMES object for landmark selection using ARG keyword");
    }
    std::string seed;
    parse("SEED",seed);
    if( seed.length()>0 ) {
      seed = " SEED=" + seed;
    }
    readInputLine( getShortcutLabel() + "_mask: CREATE_MASK ARG=" + getShortcutLabel() + "_allweights TYPE=random NZEROS=" + nlandmarks + seed );
  } else if( getName()=="LANDMARK_SELECT_FPS" ) {
    if( dissims.length()==0 ) {
      error("dissimiarities must be defined to use FPS sampling");
    }
    std::string seed;
    parse("SEED",seed);
    if( seed.length()>0 ) {
      seed = " SEED=" + seed;
    }
    readInputLine( getShortcutLabel() + "_mask: FARTHEST_POINT_SAMPLING ARG=" + dissims + " NZEROS=" + nlandmarks + seed );
  }

  if( argn.length()>0 ) {
    readInputLine( getShortcutLabel() + "_data: SELECT_WITH_MASK ARG=" + argn + "_data ROW_MASK=" + getShortcutLabel() + "_mask");
  }

  unsigned nland;
  Tools::convert( nlandmarks, nland );
  if( dissims.length()>0 ) {
    ActionWithValue* ds = plumed.getActionSet().selectWithLabel<ActionWithValue*>( dissims );
    if( (ds->copyOutput(0))->getShape()[0]==nland ) {
      if( !novoronoi ) {
        warning("cannot use voronoi procedure to give weights as not all distances between points are known");
        novoronoi=true;
      }
      readInputLine( getShortcutLabel() + "_sqrdissims: COMBINE ARG=" + dissims + " PERIODIC=NO");
    } else {
      readInputLine( getShortcutLabel() + "_rmask: CREATE_MASK ARG=" + getShortcutLabel() + "_allweights TYPE=nomask");
      readInputLine( getShortcutLabel() + "_rectdissims: SELECT_WITH_MASK ARG=" + dissims + " COLUMN_MASK=" + getShortcutLabel() + "_mask ROW_MASK=" + getShortcutLabel() + "_rmask");
      readInputLine( getShortcutLabel() + "_sqrdissims: SELECT_WITH_MASK ARG=" + dissims + " ROW_MASK=" + getShortcutLabel() + "_mask COLUMN_MASK=" + getShortcutLabel() + "_mask");
    }
  }

  if( !novoronoi && argn.length()>0 && dissims.length()>0 ) {
    readInputLine( getShortcutLabel() + "_voronoi: VORONOI ARG=" + getShortcutLabel() + "_rectdissims");
    readInputLine( getShortcutLabel() + "_allweightsT: TRANSPOSE ARG=" + getShortcutLabel() + "_allweights");
    readInputLine( getShortcutLabel() + "_weightsT: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_allweightsT," + getShortcutLabel() + "_voronoi");
    readInputLine( getShortcutLabel() + "_weights: TRANSPOSE ARG=" + getShortcutLabel() + "_weightsT");
    readInputLine( getShortcutLabel() + "_logweights: CUSTOM ARG=" + getShortcutLabel() + "_weights FUNC=log(x) PERIODIC=NO");
  } else if( argn.length()>0 ) {
    if( !novoronoi ) {
      warning("cannot use voronoi procedure to give weights to landmark points as DISSIMILARITIES was not set");
    }
    readInputLine( getShortcutLabel() + "_logweights: SELECT_WITH_MASK ARG=" + argn + "_logweights MASK=" + getShortcutLabel() + "_mask");
  }
  // Create the vector of ones that is needed by Classical MDS
  if( argn.length()>0 ) {
    readInputLine( getShortcutLabel() + "_ones: SELECT_WITH_MASK ARG=" + argn + "_ones MASK=" + getShortcutLabel() + "_mask");
  }
}

}
}
