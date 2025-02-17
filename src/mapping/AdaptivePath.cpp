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
#include "Path.h"

//+PLUMEDOC COLVAR ADAPTIVE_PATH
/*
Compute path collective variables that adapt to the lowest free energy path connecting states A and B.

The Path Collective Variables developed by Branduardi and co-workers \cite brand07 allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

\f[
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
\f]

In this expression \f$\mathbf{v}_1\f$ and \f$\mathbf{v}_3\f$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and \f$i_1\f$ and \f$i_2\f$ are the projections of the closest and second closest frames of the path. \f$\mathbf{v}_2\f$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, \f$z\f$ is calculated using:

\f[
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
\f]

Notice that these are the definitions of \f$s\f$ and \f$z\f$ that are used by \ref PATH when the GPATH option is employed.  The reason for this is that
the adaptive path method implemented in this action was inspired by the work of Diaz and Ensing in which these formula were used \cite BerndAdaptivePath.
To learn more about how the path is adapted we strongly recommend reading this paper.

\par Examples

The input below provides an example that shows how the adaptive path works. The path is updated every 50 steps of
MD based on the data accumulated during the preceding 50 time steps.

\plumedfile
d1: DISTANCE ATOMS=1,2 COMPONENTS
pp: ADAPTIVE_PATH TYPE=EUCLIDEAN FIXED=2,5 UPDATE=50 WFILE=out-path.pdb WSTRIDE=50 REFERENCE=mypath.pdb
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
\endplumedfile

In the case above the distance between frames is calculated based on the \f$x\f$ and \f$y\f$ components of the vector connecting
atoms 1 and 2.  As such an extract from the input reference path (mypath.pdb) would look as follows:

\auxfile{mypath.pdb}
REMARK ARG=d1.x,d1.y d1.x=1.12 d1.y=-.60
END
REMARK ARG=d1.x,d1.y d1.x=.99 d1.y=-.45
END
REMARK ARG=d1.x,d1.y d1.x=.86 d1.y=-.30
END
REMARK ARG=d1.x,d1.y d1.x=.73 d1.y=-.15
END
REMARK ARG=d1.x,d1.y d1.x=.60 d1.y=0
END
REMARK ARG=d1.x,d1.y d1.x=.47 d1.y=.15
END
\endauxfile

Notice that one can also use RMSD frames in place of arguments like those above.

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
  keys.add("compulsory","TOLERANCE","1E-6","the tolerance to use for the path updating algorithm that makes all frames equidistant");
  keys.add("optional","WFILE","file on which to write out the path");
  keys.add("compulsory","FMT","%f","the format to use for output files");
  keys.add("compulsory","WSTRIDE","0,","frequency with which to write out the path");
  keys.setValueDescription("scalar","the position along and from the adaptive path");
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
  FILE* fp=std::fopen(reference.c_str(),"r");
  PDB mypdb;
  if(!fp) {
    error("could not open reference file " + reference );
  }
  bool do_read=mypdb.readFromFilepointer(fp,false,0.1);
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
    metric = "RMSD_VECTOR DISPLACEMENT TYPE=" + mtype + " ALIGN=" + align_str + " DISPLACE=" + displace_str;
    readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data.disp " + " PROPERTY=" +  propstr + " REFERENCE=" + reference_data + " METRIC={" + metric + "} METRIC_COMPONENT=disp");
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
      FILE* fp=std::fopen(reference.c_str(),"r");
      double fake_unit=0.1;
      PDB mypdb;
      bool do_read=mypdb.readFromFilepointer(fp,false,fake_unit);
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
