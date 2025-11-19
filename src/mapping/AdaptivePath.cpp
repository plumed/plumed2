/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "tools/PDB.h"
#include "tools/FileTools.h"
#include "Path.h"

//+PLUMEDOC COLVAR ADAPTIVE_PATH
/*
Compute path collective variables that adapt to the lowest free energy path connecting states A and B.

This shortcut can be used to implement the method discussed in the second paper cited below. This method allows you to
run a simulation with [PATH](PATH.md) collective variables that forces the system to transition between state A and state B.
If you run with ADAPTIVE_PATH rather than [PATH](PATH.md) the CV adapts as the simulation progresses and should find the
lowest free energy path that connects the two states.

The Path Collective Variables developed by Branduardi that are described in the first paper cited below allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  Instad of computing progress along the path ($s$) using Branduardi's method we use the expressions that are used in
[GEOMETRIC_PATH](GEOMETRIC_PATH.md) and that are introduced in the second paper that is cited below.  The progress along the path is is thus computed using:

$$
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
$$

In this expression $\mathbf{v}_1$ and $\mathbf{v}_3$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and $i_1$ and $i_2$ are the projections of the closest and second closest frames of the path. $\mathbf{v}_2$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, $z$ is calculated using:

$$
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
$$

Notice that these are the definitions of $s$ and $z$ from the [GEOMETRIC_PATH](GEOMETRIC_PATH.md).  The reason for this is that
the adaptive path method implemented in this action was inspired by the work of Diaz and Ensing that was introduced in the second paper cited below.
To learn more about how the path is adapted we strongly recommend reading this paper.

The example input below shows how to use the adaptive path.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-adapt/mypath.pdb
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: ADAPTIVE_PATH ...
   ARG=d1.x,d1.y UPDATE=50 FIXED=5,15
   WFILE=out-path.pdb WSTRIDE=50
   REFERENCE=regtest/mapping/rt-adapt/mypath.pdb
   PROPERTY=path
...
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
```

The curved path here is defined using a series of points in a two dimensional space.  Furthermore, the coordinates that these points should be projected at
on the low dimensional path (the $i$ values in the equation for $s$ above) are specified in the `path` property of the input pdb file.

If you expand the shortcuts in the input above you can see how
we use a [GEOMETRIC_PATH](GEOMETRIC_PATH.md) shortcut to calculate our position along the path and distance from the path.  In addition, we also collect
information on the average displacement from the path using an [AVERAGE_PATH_DISPLACEMENT](AVERAGE_PATH_DISPLACEMENT.md) so that we can update the path every step.
During these update steps the reference points on the path are shifted by the values that are stored in the [AVERAGE_PATH_DISPLACEMENT](AVERAGE_PATH_DISPLACEMENT.md)
action.  We then use the [REPARAMETERIZE_PATH](REPARAMETERIZE_PATH.md) to ensure that the new set of reference points on the path are all equally spaced.  After the adaption
steps the new reference points are used when calculating are position on and distance from the reference path.

The input below illustrates how to the positions of atoms to define the adaptive path:

```plumed
##SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1b: ADAPTIVE_PATH REFERENCE=regtest/trajectories/path_msd/all.pdb FIXED=1,42 UPDATE=1000 WFILE=out-path.pdb WSTRIDE=20000
PRINT ARG=p1b.* FILE=colvar_b STRIDE=1
```

As we have not used the PROPERTY keyword in the input the $i$ values in the equation for $s$ above are set equal to 1, 2, 3...
Furthermore, when an input like the one above is used the vectors $\mathbf{v}_1$, $\mathbf{v}_2$ and $\mathbf{v}_3$ from the expression above are computed using an [RMSD](RMSD.md)
action with the DISPLACEMENT flag enabled.  The instaneous structure is thus aligned with the reference structures so as to remove motions due to translation of the center
of mass and rotation of the reference frame. These vector of atomic displacements that do not include the translation and rotation are also used in the
[AVERAGE_PATH_DISPLACEMENT](AVERAGE_PATH_DISPLACEMENT.md) action that is used to update the path.

As we have not used the PROPERTY keyword in the input the $i$ values in the equation for $s$ above are set equal to 1, 2, 3...


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class AdaptivePath : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit AdaptivePath(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AdaptivePath,"ADAPTIVE_PATH")

void AdaptivePath::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
  keys.add("compulsory","FIXED","the positions in the list of input frames of the two path nodes whose positions remain fixed during the path optimization");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50 percent in the average. This option may increase convergence by allowing to forget the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","UPDATE","the frequency with which the path should be updated");
  keys.add("optional","WFILE","file on which to write out the path");
  keys.add("compulsory","FMT","%f","the format to use for output files");
  keys.add("compulsory","WSTRIDE","0,","frequency with which to write out the path");
  keys.setValueDescription("scalar","the position along and from the adaptive path");
  keys.addDOI("10.1063/1.2432340");
  keys.addDOI("10.1103/PhysRevLett.109.020601");
  keys.needsAction("GEOMETRIC_PATH");
  keys.needsAction("AVERAGE_PATH_DISPLACEMENT");
  keys.needsAction("REPARAMETERIZE_PATH");
  keys.needsAction("DUMPPDB");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("DISPLACEMENT");
  keys.needsAction("CONSTANT");
}

AdaptivePath::AdaptivePath(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the arguments
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  std::string reference_data, metric, mtype;
  parse("TYPE", mtype);
  std::string reference;
  parse("REFERENCE",reference);
  {
    //extra scope for fp and mypdb
    unique_FILE fp{std::fopen(reference.c_str(),"r")};
    PDB mypdb;
    if(!fp.get()) {
      error("could not open reference file " + reference );
    }
    bool do_read=mypdb.readFromFilepointer(fp.get(),false,0.1);
    if( !do_read ) {
      error("missing file " + reference );
    }
    // Create list of reference configurations that PLUMED will use
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
      metric = "RMSD DISPLACEMENT TYPE=" + mtype + " ALIGN=" + align_str + " DISPLACE=" + displace_str;
      readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data.disp " + " PROPERTY=" +  propstr + " REFERENCE=" + reference_data + " METRIC={" + metric + "} METRIC_COMPONENT=disp");
    }
  }
  // Create the object to accumulate the average path displacements
  std::string update, halflife;
  parse("HALFLIFE",halflife);
  parse("UPDATE",update);
  std::string refframes = " REFERENCE=" + getShortcutLabel() + "_pos";
  if( argnames.size()>0 ) {
    readInputLine( getShortcutLabel() + "_disp: AVERAGE_PATH_DISPLACEMENT ARG=" + getShortcutLabel() + "_data HALFLIFE=" + halflife + " CLEAR=" + update + " METRIC={DIFFERENCE} REFERENCE=" + reference_data );
  } else {
    readInputLine( getShortcutLabel() + "_disp: AVERAGE_PATH_DISPLACEMENT ARG=" + getShortcutLabel() + "_data.disp HALFLIFE=" + halflife + " CLEAR=" + update + " METRIC={" + metric + "} METRIC_COMPONENT=disp REFERENCE=" + reference_data );
  }

  // Create the object to update the path
  std::string fixedn;
  parse("FIXED",fixedn);
  std::string component="METRIC_COMPONENT=disp";
  if( argnames.size()>0 ) {
    metric="DIFFERENCE";
    component="";
  }
  readInputLine("REPARAMETERIZE_PATH DISPLACE_FRAMES=" + getShortcutLabel() + "_disp FIXED=" + fixedn + " STRIDE=" + update + " METRIC={" + metric + "} " + component + " REFERENCE=" + reference_data );

  // Information for write out
  std::string wfilename;
  parse("WFILE",wfilename);
  if( wfilename.length()>0 ) {
    // This just gets the atom numbers for output
    std::string atomstr;
    if( argnames.size()==0 ) {
      unique_FILE fp{std::fopen(reference.c_str(),"r")};
      double fake_unit=0.1;
      PDB mypdb;
      bool do_read=mypdb.readFromFilepointer(fp.get(),false,fake_unit);
      if( !do_read ) {
        error("missing file " + reference );
      }
      std::string num;
      Tools::convert( mypdb.getAtomNumbers()[0].serial(), atomstr );
      for(unsigned j=1; j<mypdb.getAtomNumbers().size(); ++j ) {
        Tools::convert( mypdb.getAtomNumbers()[j].serial(), num );
        atomstr += "," + num;
      }
      Tools::convert( mypdb.getResidueNumber( mypdb.getAtomNumbers()[0] ), num );
      atomstr += " RESIDUE_INDICES=" + num;
      for(unsigned i=1; i<mypdb.getAtomNumbers().size(); ++i ) {
        Tools::convert( mypdb.getResidueNumber( mypdb.getAtomNumbers()[i] ), num );
        atomstr += "," + num;
      }
      std::string anum, dnum;
      Tools::convert( mypdb.getOccupancy()[0], anum );
      Tools::convert( mypdb.getBeta()[0], dnum );
      atomstr += " OCCUPANCY=" + anum;
      for(unsigned i=1; i<mypdb.getOccupancy().size(); ++i) {
        Tools::convert( mypdb.getOccupancy()[i], anum );
        atomstr += "," + anum;
      }
      atomstr += " BETA=" + dnum;
      for(unsigned i=1; i<mypdb.getBeta().size(); ++i) {
        Tools::convert( mypdb.getBeta()[i], dnum );
        atomstr += "," + dnum;
      }
    }

    if( wfilename.find(".pdb")==std::string::npos ) {
      error("output must be to a pdb file");
    }
    std::string ofmt, pframes, wstride;
    parse("WSTRIDE",wstride);
    parse("FMT",ofmt);
    if( argnames.size()>0 ) {
      std::string argstr = argnames[0];
      for(unsigned i=1; i<argnames.size(); ++i) {
        argstr += "," + argnames[i];
      }
      readInputLine("DUMPPDB DESCRIPTION=PATH STRIDE=" + wstride + " FMT=" + ofmt + " FILE=" + wfilename + " ARG=" + reference_data );
    } else {
      readInputLine("DUMPPDB DESCRIPTION=PATH STRIDE=" + wstride + " FMT=" + ofmt + " FILE=" + wfilename + " ATOMS=" + reference_data + " ATOM_INDICES=" + atomstr );
    }
  }
  log<<"  Bibliography "<<plumed.cite("Diaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
}

}
}
