/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
pp: ADAPTIVE_PATH TYPE=EUCLIDEAN FIXED=5,15 UPDATE=50 WFILE=out-path.pdb WSTRIDE=50 REFERENCE=mypath.pdb
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
\endplumedfile

In the case above the distance between frames is calculated based on the \f$x\f$ and \f$y\f$ components of the vector connecting
atoms 1 and 2.  As such an extract from the input reference path (mypath.pdb) would look as follows:

\verbatim
REMARK ARG=d1.x,d1.y d1.x=1.12 d1.y=-.60
END
REMARK ARG=d1.x,d1.y d1.x=.99 d1.y=-.45
END
\endverbatim

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
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","ARG","the list of arguments you would like to use in your definition of the path");
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
  keys.add("compulsory","FIXED","the positions in the list of input frames of the two path nodes whose positions remain fixed during the path optimization");
  keys.add("compulsory","HALFLIFE","-1","the number of MD steps after which a previously measured path distance weighs only 50% in the average. This option may increase convergence by allowing to \"forget\" the memory of a bad initial guess path. The default is to set this to infinity");
  keys.add("compulsory","UPDATE","the frequency with which the path should be updated");
  keys.add("compulsory","TOLERANCE","1E-6","the tolerance to use for the path updating algorithm that makes all frames equidistant");
  keys.add("optional","WFILE","file on which to write out the path");
  keys.add("compulsory","FMT","%f","the format to use for output files");
  keys.add("compulsory","WSTRIDE","frequency with which to write out the path");
}

AdaptivePath::AdaptivePath(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Create the geometric path object to compute distances
  std::string reffile, mtype; parse("REFERENCE",reffile); parse("TYPE",mtype);
  std::vector<std::string> argnames; parseVector("ARG",argnames); std::string additional_input="";
  if( argnames.size()>0 ) {
      additional_input=" ARG=" + argnames[0]; for(unsigned i=1;i<argnames.size();++i) additional_input += "," + argnames[i];
      if( mtype=="OPTIMAL_FAST" ) mtype="EUCLIDEAN";
  }
  std::string pname; parse("PROPERTY",pname);
  if( pname.length()>0 ) additional_input += " PROPERTY=" + pname;
  readInputLine( getShortcutLabel() + ": GPATH REFERENCE=" + reffile + " TYPE=" + mtype + additional_input );

  // Get the number of frames in the path
  std::string metric; unsigned nframes = Path::getNumberOfFramesAndMetric( mtype, reffile, metric );

  // Create the object to accumulate the average path displacements
  std::string update, halflife; parse("HALFLIFE",halflife); parse("UPDATE",update);
  std::string refframes = " REFFRAMES=" + getShortcutLabel() + "_ref1"; 
  for(unsigned i=1;i<nframes;++i){ std::string num; Tools::convert(i+1,num); refframes += "," + getShortcutLabel() + "_ref" + num; }
  readInputLine( getShortcutLabel() + "_disp: AVERAGE_PATH_DISPLACEMENT ARG=" + getShortcutLabel() + "_data HALFLIFE=" + halflife + " CLEAR=" + update + metric + refframes );

  // Create the object to update the path
  std::string fixedn; parse("FIXED",fixedn);
  if( fixedn.length()>0 ) readInputLine("REPARAMETERIZE_PATH DISPLACE_FRAMES=" + getShortcutLabel() + "_disp FIXED=" + fixedn + " STRIDE=" + update + metric + refframes );
  else readInputLine("REPARAMETERIZE_PATH DISPLACE_FRAMES=" + getShortcutLabel() + "_disp STRIDE=" + update + metric + refframes );

  // Information for write out
  std::string wfilename; parse("WFILE",wfilename);
  if( wfilename.length()>0 ) {
      if( wfilename.find(".pdb")==std::string::npos ) error("output must be to a pdb file");
      std::string ofmt, pframes, wstride; parse("WSTRIDE",wstride); parse("FMT",ofmt); 
      for(unsigned i=0;i<nframes;++i) { std::string num; Tools::convert( i+1, num ); pframes += " CONFIG" + num + "=" + getShortcutLabel() + "_ref" + num; }
      readInputLine("PRINT STRIDE=" + wstride + " FMT=" + ofmt + " FILE=" + wfilename + pframes );
  }
  log<<"  Bibliography "<<plumed.cite("Diaz Leines and Ensing, Phys. Rev. Lett. 109, 020601 (2012)")<<"\n";
}

}
}
