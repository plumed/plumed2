/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2018 The plumed team
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
#include "Path.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"

//+PLUMEDOC COLVAR GPATH
/*
Distance along and from a path calculated using geometric formulas

The Path Collective Variables developed by Branduardi and that are described in the first paper that is cited below alow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional path.  The method introduced in that
paper is implemented in the shortcut [PATH](PATH.md). This action provides an implementation of the alternative method for calculating
the position on and distance from the path that was proposed by Diaz Leines and Ensing in the second paper that is cited below.  In their
method, the progress along the path $s$ is calculated using:

$$
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
$$

where $\mathbf{v}_1$ and $\mathbf{v}_3$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and $i_1$ and $i_2$ are the projections of the closest and second closest frames of the path. $\mathbf{v}_2$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, $z$ is calculated using:

$$
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
$$

The symbols here are as they were for $s$.

## Examples

The example input below shows how to use this shortcut.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-adapt/mypath.pdb
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: GPATH ...
   ARG=d1.x,d1.y PROPERTY=path
   REFERENCE=regtest/mapping/rt-adapt/mypath.pdb
...
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
```

The curved path here is defined using a series of points in a two dimensional space.  Furthermore, the coordinates that these points should be projected at on the low dimensional path
(the $i$ values in the equation for $s$ above) are specified in the `path` property of the input pdb file.

The input below illustrates how to the positions of atoms to define the path:

```plumed
##SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1b: GPATH REFERENCE=regtest/trajectories/path_msd/all.pdb
PRINT ARG=p1b.* FILE=colvar_b STRIDE=1
```

When an input like the one above is used the vectors $\mathbf{v}_1$, $\mathbf{v}_2$ and $\mathbf{v}_3$ from the expression above are computed using an [RMSD](RMSD.md)
action with the DISPLACEMENT flag enabled.  The instaneous structure is thus aligned with the reference structures so as to remove motions due to translation of the center
of mass and rotation of the reference frame.  Furthermore, because we have not used the PROPERTY keyword in this input, the $i$ values in the equation for $s$
above are set equal to 1, 2, 3...

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class GeometricPathShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPathShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(GeometricPathShortcut,"GPATH")

void GeometricPathShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
  keys.addOutputComponent("s","default","scalar","the position on the path");
  keys.addOutputComponent("z","default","scalar","the distance from the path");
  keys.needsAction("DISPLACEMENT");
  keys.needsAction("GEOMETRIC_PATH");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("CONSTANT");
  keys.addDOI("10.1063/1.2432340");
  keys.addDOI("10.1103/PhysRevLett.109.020601");
}

GeometricPathShortcut::GeometricPathShortcut( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  std::string mtype, reference_data;
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  parse("TYPE", mtype);
  // Create list of reference configurations that PLUMED will use
  std::string reference;
  parse("REFERENCE",reference);
  FILE* fp=std::fopen(reference.c_str(),"r");
  PDB mypdb;
  if(!fp) {
    error("could not open reference file " + reference );
  }
  bool do_read=mypdb.readFromFilepointer(fp,false,0.1);
  if( !do_read ) {
    error("missing file " + reference );
  }
  Path::readInputFrames( reference, mtype, argnames, true, this, reference_data );
  // Now get coordinates on spath
  std::vector<std::string> pnames;
  parseVector("PROPERTY",pnames);
  Path::readPropertyInformation( pnames, getShortcutLabel(), reference, this );
  // Create action that computes the geometric path variablesa
  std::string propstr = getShortcutLabel() + "_ind";
  if( pnames.size()>0 ) {
    propstr = pnames[0] + "_ref";
  }
  if( argnames.size()>0 ) {
    readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data " + " PROPERTY=" + propstr + " REFERENCE=" + reference_data + " METRIC={DIFFERENCE}");
  } else {
    std::string num, align_str, displace_str;
    Tools::convert( mypdb.getOccupancy()[0], align_str );
    Tools::convert( mypdb.getBeta()[0], displace_str );
    for(unsigned j=1; j<mypdb.getAtomNumbers().size(); ++j ) {
      Tools::convert( mypdb.getOccupancy()[j], num );
      align_str += "," + num;
      Tools::convert( mypdb.getBeta()[0], num );
      displace_str += "," + num;
    }
    std::string metric = "RMSD DISPLACEMENT TYPE=" + mtype + " ALIGN=" + align_str + " DISPLACE=" + displace_str;
    readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data.disp " + " PROPERTY=" +  propstr + " REFERENCE=" + reference_data + " METRIC={" + metric + "} METRIC_COMPONENT=disp");
  }
}


}
}


