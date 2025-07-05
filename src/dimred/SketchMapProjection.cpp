/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "mapping/Path.h"

//+PLUMEDOC DIMRED SKETCHMAP_PROJECTION
/*
Read in a sketch-map projection

This shortcut can be used to read in a projection of a trajectory that was generated using [SKETCHMAP](SKETCHMAP.md).
You can then use the coordinates that were read in to generate projections of other configurations.  Examples of how
this tool might be used are given in the last three papers cited below.  In these papers, a sketch-map projection is
computed from one particular set of data points.  This sketch-map projection that was output was then used to analyse
a second completely different data set.

An example input that illustrates how the sketch-map projection shortcut is used is shown below:

```plumed
#SETTING INPUTFILES=regtest/dimred/rt-smap-read/smap.pdb
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6

smap: SKETCHMAP_PROJECTION ...
    ARG=d1,d2,d3 REFERENCE=regtest/dimred/rt-smap-read/smap.pdb
    PROPERTY=smap_coord-1,smap_coord-2 CGTOL=1E-3
    WEIGHT=WEIGHT HIGH_DIM_FUNCTION={SMAP R_0=4 A=3 B=2} LOW_DIM_FUNCTION={SMAP R_0=4 A=1 B=2}
...

PRINT ARG=smap.* FILE=colvar
```

Each frame in your input trajectory generates the three distances so one set of sketch-map coordinates are generated from each frame.
Notice that if you want to generate projections of multiple input points at once you need to use [PROJECT_POINTS](PROJECT_POINTS.md)
directly rather than this wrapper.

The configurations of the landmark points in your sketch-map input file can also be defined in terms of a set of atomic positions.
In this case the distance between the reference configurations and the instantaneous coordinates are calculated using [RMSD](RMSD.md).
An input for this type of calculation is as follows:

```plumed
#SETTING INPUTFILES=regtest/dimred/rt-mds-rmsd/embed.pdb.reference

smap: SKETCHMAP_PROJECTION ...
   REFERENCE=regtest/dimred/rt-mds-rmsd/embed.pdb.reference
   PROPERTY=mds-1,mds-2 CGTOL=1E-3
   WEIGHT=weights HIGH_DIM_FUNCTION={SMAP R_0=4 A=3 B=2} LOW_DIM_FUNCTION={SMAP R_0=4 A=1 B=2}
...
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapProjection : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapProjection( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(SketchMapProjection,"SKETCHMAP_PROJECTION")

void SketchMapProjection::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  mapping::Path::registerInputFileKeywords( keys );
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","WEIGHT","the weight of each individual landmark in the stress fucntion that is to be optimised");
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","CGTOL","1E-6","The tolerance for the conjugate gradient minimization that finds the out of sample projections");
  keys.setValueDescription("scalar/vector","the out-of-sample projections of the input arguments using the input sketch-map projection");
  keys.addDOI("10.1073/pnas.1108486108");
  keys.addDOI("10.1073/pnas.1201152109");
  keys.addDOI("10.1021/ct3010563");
  keys.addDOI("10.1021/ct500950z");
  keys.addDOI("10.1021/acs.jctc.5b00714");
  keys.needsAction("RMSD");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("CONSTANT");
  keys.needsAction("CUSTOM");
  keys.needsAction("EUCLIDEAN_DISTANCE");
  keys.needsAction("NORMALIZED_EUCLIDEAN_DISTANCE");
  keys.needsAction("SUM");
  keys.needsAction("MORE_THAN");
  keys.needsAction("PROJECT_POINTS");
}

SketchMapProjection::SketchMapProjection( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Use path to read in the projections
  std::string refname, refactions, metric;
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  std::string type, reference_data, reference;
  parse("REFERENCE",reference);
  parse("TYPE",type);
  mapping::Path::readInputFrames( reference, type, argnames, false, this, reference_data );
  // And read in the data that we want on the projections
  std::vector<std::string> pnames;
  parseVector("PROPERTY",pnames);
  std::string weights;
  parse("WEIGHT",weights);
  pnames.push_back( weights );
  // Now create fixed vectors using some sort of reference action
  mapping::Path::readPropertyInformation( pnames, getShortcutLabel(), reference, this );
  // Normalise the vector of weights
  readInputLine( getShortcutLabel() + "_wsum: SUM PERIODIC=NO ARG=" + weights + "_ref");
  readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG=" + weights + "_ref," + getShortcutLabel() + "_wsum FUNC=x/y PERIODIC=NO");
  // Transform the high dimensional distances
  std::string hdfunc;
  parse("HIGH_DIM_FUNCTION",hdfunc);
  readInputLine( getShortcutLabel() + "_targ: MORE_THAN ARG=" + getShortcutLabel() + "_data SQUARED SWITCH={" + hdfunc + "}");
  // Create the projection object
  std::string ldfunc, cgtol;
  parse("LOW_DIM_FUNCTION",ldfunc);
  parse("CGTOL",cgtol);
  std::string argstr="ARG=" + pnames[0] + "_ref";
  for(unsigned i=1; i<pnames.size()-1; ++i) {
    argstr += "," + pnames[i] + "_ref";
  }
  readInputLine( getShortcutLabel() + ": PROJECT_POINTS " + argstr + " TARGET1=" + getShortcutLabel() + "_targ " +
                 "FUNC1={" + ldfunc + "} WEIGHTS1=" + getShortcutLabel() + "_weights CGTOL=" + cgtol );
}

}
}
