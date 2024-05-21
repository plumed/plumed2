/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/PDB.h"

//+PLUMEDOC COLVAR PATH
/*
Path collective variables with a more flexible framework for the distance metric being used.

The Path Collective Variables developed by Branduardi and co-workers \cite brand07 allow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

\f[
s = \frac{ \sum_{i=1}^N i \exp( -\lambda R[X - X_i] ) }{ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) }
\f]

while the distance from the path (z) is measured using:

\f[
z = -\frac{1}{\lambda} \ln\left[ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) \right]
\f]

In these expressions \f$N\f$ high-dimensional frames (\f$X_i\f$) are used to describe the path in the high-dimensional
space. The two expressions above are then functions of the distances from each of the high-dimensional frames \f$R[X - X_i]\f$.
Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration.  You could calculate
the RMSD distance or you could calculate the amount by which a set of collective variables change.  As such this implementation
of the path CV allows one to use all the difference distance metrics that are discussed in \ref dists. This is as opposed to
the alternative implementation of path (\ref PATHMSD) which is a bit faster but which only allows one to use the RMSD distance.

The \f$s\f$ and \f$z\f$ variables are calculated using the above formulas by default.  However, there is an alternative method
of calculating these collective variables, which is detailed in \cite bernd-path.  This alternative method uses the tools of
geometry (as opposed to algebra, which is used in the equations above).  In this alternative formula the progress along the path
\f$s\f$ is calculated using:

\f[
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
\f]

where \f$\mathbf{v}_1\f$ and \f$\mathbf{v}_3\f$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and \f$i_1\f$ and \f$i_2\f$ are the projections of the closest and second closest frames of the path. \f$\mathbf{v}_2\f$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, \f$z\f$ is calculated using:

\f[
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
\f]

The symbols here are as they were for \f$s\f$.  If you would like to use these equations to calculate \f$s\f$ and \f$z\f$ then you should use the GPATH flag.
The values of \f$s\f$ and \f$z\f$ can then be referenced using the gspath and gzpath labels.

\par Examples

In the example below the path is defined using RMSD distance from frames.

\plumedfile
p1: PATH REFERENCE=file.pdb TYPE=OPTIMAL LAMBDA=500.0
PRINT ARG=p1.spath,p1.zpath STRIDE=1 FILE=colvar FMT=%8.4f
\endplumedfile

The reference frames in the path are defined in the pdb file shown below.  In this frame
each configuration in the path is separated by a line containing just the word END.

\auxfile{file.pdb}
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
END
ATOM      1  CL  ALA     1      -3.175   0.365   2.024  1.00  1.00
ATOM      5  CLP ALA     1      -1.814  -0.106   1.685  1.00  1.00
ATOM      6  OL  ALA     1      -1.201  -0.849   2.425  1.00  1.00
ATOM      7  NL  ALA     1      -1.296   0.337   0.534  1.00  1.00
END
ATOM      1  CL  ALA     1      -2.990   0.383   2.277  1.00  1.00
ATOM      5  CLP ALA     1      -1.664  -0.085   1.831  1.00  1.00
ATOM      6  OL  ALA     1      -0.987  -0.835   2.533  1.00  1.00
ATOM      7  NL  ALA     1      -1.227   0.364   0.646  1.00  1.00
END
\endauxfile

In the example below the path is defined using the values of two torsional angles (t1 and t2).
In addition, the \f$s\f$ and \f$z\f$ are calculated using the geometric expressions described
above rather than the algebraic expressions that are used by default.

\plumedfile
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: PATH TYPE=EUCLIDEAN REFERENCE=epath.pdb GPATH NOSPATH NOZPATH
PRINT ARG=pp.* FILE=colvar
\endplumedfile

Notice that the LAMBDA parameter is not required here as we are not calculating \f$s\f$ and \f$s\f$
using the algebraic formulas defined earlier.  The positions of the frames in the path are defined
in the file epath.pdb.  An extract from this file looks as shown below.

\auxfile{epath.pdb}
REMARK ARG=t1,t2 t1=-4.25053  t2=3.88053
END
REMARK ARG=t1,t2 t1=-4.11     t2=3.75
END
REMARK ARG=t1,t2 t1=-3.96947  t2=3.61947
END
\endauxfile

The remarks in this pdb file tell PLUMED the labels that are being used to define the position in the
high dimensional space and the values that these arguments have at each point on the path.

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR GPROPERTYMAP
/*
Property maps but with a more flexible framework for the distance metric being used.

This colvar calculates a property map using the formalism developed by Spiwok \cite Spiwok:2011ce.
In essence if you have the value of some property, \f$X_i\f$, that it takes at a set of high-dimensional
positions then you calculate the value of the property at some arbitrary point in the high-dimensional space
using:

\f[
X=\frac{\sum_i X_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))}
\f]

Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration, \f$D_i\f$.  You could calculate
the RMSD distance or you could calculate the amount by which a set of collective variables change.  As such this implementation
of the property map allows one to use all the different distance metric that are discussed in \ref dists. This is as opposed to
the alternative implementation \ref PROPERTYMAP which is a bit faster but which only allows one to use the RMSD distance.

\par Examples

The input shown below can be used to calculate the interpolated values of two properties called X and Y based on the values
that these properties take at a set of reference configurations and using the formula above.  For this input the distances
between the reference configurations and the instantaneous configurations are calculated using the OPTIMAL metric that is
discussed at length in the manual pages on \ref RMSD.

\plumedfile
p2: GPROPERTYMAP REFERENCE=allv.pdb PROPERTY=X,Y LAMBDA=69087
PRINT ARG=p2.X,p2.Y,p2.zpath STRIDE=1 FILE=colvar
\endplumedfile

The additional input file for this calculation, which contains the reference frames and the values of X and Y at these reference
points has the following format.

\auxfile{allv.pdb}
REMARK X=1 Y=2
ATOM      1  CL  ALA     1      -3.171   0.295   2.045  1.00  1.00
ATOM      5  CLP ALA     1      -1.819  -0.143   1.679  1.00  1.00
ATOM      6  OL  ALA     1      -1.177  -0.889   2.401  1.00  1.00
ATOM      7  NL  ALA     1      -1.313   0.341   0.529  1.00  1.00
ATOM      8  HL  ALA     1      -1.845   0.961  -0.011  1.00  1.00
ATOM      9  CA  ALA     1      -0.003  -0.019   0.021  1.00  1.00
ATOM     10  HA  ALA     1       0.205  -1.051   0.259  1.00  1.00
ATOM     11  CB  ALA     1       0.009   0.135  -1.509  1.00  1.00
ATOM     15  CRP ALA     1       1.121   0.799   0.663  1.00  1.00
ATOM     16  OR  ALA     1       1.723   1.669   0.043  1.00  1.00
ATOM     17  NR  ALA     1       1.423   0.519   1.941  1.00  1.00
ATOM     18  HR  ALA     1       0.873  -0.161   2.413  1.00  1.00
ATOM     19  CR  ALA     1       2.477   1.187   2.675  1.00  1.00
END
FIXED
REMARK X=2 Y=3
ATOM      1  CL  ALA     1      -3.175   0.365   2.024  1.00  1.00
ATOM      5  CLP ALA     1      -1.814  -0.106   1.685  1.00  1.00
ATOM      6  OL  ALA     1      -1.201  -0.849   2.425  1.00  1.00
ATOM      7  NL  ALA     1      -1.296   0.337   0.534  1.00  1.00
ATOM      8  HL  ALA     1      -1.807   0.951  -0.044  1.00  1.00
ATOM      9  CA  ALA     1       0.009  -0.067   0.033  1.00  1.00
ATOM     10  HA  ALA     1       0.175  -1.105   0.283  1.00  1.00
ATOM     11  CB  ALA     1       0.027   0.046  -1.501  1.00  1.00
ATOM     15  CRP ALA     1       1.149   0.725   0.654  1.00  1.00
ATOM     16  OR  ALA     1       1.835   1.491  -0.011  1.00  1.00
ATOM     17  NR  ALA     1       1.380   0.537   1.968  1.00  1.00
ATOM     18  HR  ALA     1       0.764  -0.060   2.461  1.00  1.00
ATOM     19  CR  ALA     1       2.431   1.195   2.683  1.00  1.00
END
\endauxfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_ACTION(Path,"GPROPERTYMAP")

void Path::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.addOutputComponent("gspath","GPATH","the position along the path calculated using the geometric formula");
  keys.addOutputComponent("gzpath","GPATH","the distance from the path calculated using the geometric formula");
  keys.addOutputComponent("spath","default","the position along the path calculated");
  keys.addOutputComponent("zpath","default","the distance from the path calculated");
}

void Path::registerInputFileKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","ARG","the list of arguments you would like to use in your definition of the path");
  keys.add("optional","COEFFICIENTS","the coefficients of the displacements along each argument that should be used when calculating the euclidean distance");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("NOSPATH",false,"do not calculate the spath CV");
  keys.addFlag("NOZPATH",false,"do not calculate the zpath CV");
  keys.addFlag("GPATH",false,"calculate the trigonometric path");
  keys.needsAction("DRMSD"); keys.needsAction("RMSD"); keys.needsAction("LOWEST"); keys.needsAction("GPATH");
  keys.needsAction("EUCLIDEAN_DISTANCE"); keys.needsAction("CUSTOM"); keys.needsAction("SUM"); keys.needsAction("COMBINE");
  keys.needsAction("NORMALIZED_EUCLIDEAN_DISTANCE"); keys.needsAction("PDB2CONSTANT"); keys.needsAction("CONSTANT");
}

Path::Path( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  bool nospath, nozpath, gpath; parseFlag("NOSPATH",nospath); parseFlag("NOZPATH",nozpath); parseFlag("GPATH",gpath);
  if( gpath ) {
    readInputLine( getShortcutLabel() + "_gpath: GPATH " + convertInputLineToString() );
    readInputLine( getShortcutLabel() + "_gspath: COMBINE ARG=" + getShortcutLabel() + "_gpath.s PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_gzpath: COMBINE ARG=" + getShortcutLabel() + "_gpath.z PERIODIC=NO");
  }
  if( nospath && nozpath ) return;
  // Setup the properties
  std::vector<std::string> properties, pnames;
  if( getName()=="PATH") { properties.resize(1); }
  else { parseVector("PROPERTY",pnames); properties.resize( pnames.size() ); }
  std::string type, reference_data, reference; parse("REFERENCE",reference);
  std::vector<std::string> argnames; parseVector("ARG",argnames); parse("TYPE",type);
  if( type.find("DRMSD")!=std::string::npos ) readInputLine( getShortcutLabel() + "_data: DRMSD SQUARED TYPE=" + type + " REFERENCE=" + reference );
  else readInputFrames( reference, type, argnames, false, this, reference_data );
  // Find the shortest distance to the frames
  readInputLine( getShortcutLabel() + "_mindist: LOWEST ARG=" + getShortcutLabel() + "_data");
  // Now create all other parts of the calculation
  std::string lambda; parse("LAMBDA",lambda);
  // Now create MATHEVAL object to compute exponential functions
  readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG=" + getShortcutLabel() + "_data," + getShortcutLabel() + "_mindist FUNC=exp(-(x-y)*" + lambda + ") PERIODIC=NO" );
  // Create denominator
  readInputLine( getShortcutLabel() + "_denom: SUM ARG=" + getShortcutLabel() + "_weights PERIODIC=NO");
  // Now compte zpath variable
  if( !nozpath ) {
    readInputLine( getShortcutLabel() + "_z: CUSTOM ARG=" + getShortcutLabel() + "_denom," + getShortcutLabel() + "_mindist FUNC=y-log(x)/" + lambda + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zpath: COMBINE ARG=" + getShortcutLabel() + "_z PERIODIC=NO");
  }
  // Now get coefficients for properies for spath
  readPropertyInformation( pnames, getShortcutLabel(), reference, this );
  // Now create COMBINE objects to compute numerator of path
  for(unsigned i=0; i<properties.size(); ++i) {
    if( pnames.size()>0 ) {
      readInputLine( pnames[i] + "_numer_prod: CUSTOM ARG=" + getShortcutLabel() + "_weights," + pnames[i] + "_ref FUNC=x*y PERIODIC=NO");
      readInputLine( pnames[i] + "_numer: SUM ARG=" + pnames[i]  + "_numer_prod PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_" + pnames[i] + ": CUSTOM ARG=" + pnames[i]  + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    } else if( !nospath ) {
      readInputLine( getShortcutLabel() + "_s_prod: CUSTOM ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_ind FUNC=x*y PERIODIC=NO");
      readInputLine( getShortcutLabel()  + "_numer: SUM ARG=" + getShortcutLabel() + "_s_prod PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_s: CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_spath: COMBINE ARG=" + getShortcutLabel() + "_s PERIODIC=NO");
    }
  }
}

std::string Path::fixArgumentName( const std::string& argin ) {
  std::string argout = argin; std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  return argout;
}

void Path::readInputFrames( const std::string& reference, const std::string& type, std::vector<std::string>& argnames, const bool& displacements, ActionShortcut* action, std::string& reference_data ) {
  FILE* fp=std::fopen(reference.c_str(),"r"); PDB pdb; if(!fp) action->error("could not open reference file " + reference );
  bool do_read=pdb.readFromFilepointer(fp,false,0.1); if( !do_read ) action->error("missing file " + reference );
  if( pdb.getPositions().size()!=0 && argnames.size()==0 ) {
    reference_data = action->getShortcutLabel() + "_data_ref";
    if( displacements ) action->readInputLine( action->getShortcutLabel() + "_data: RMSD DISPLACEMENT SQUARED REFERENCE=" + reference + " TYPE=" + type );
    else action->readInputLine( action->getShortcutLabel() + "_data: RMSD SQUARED REFERENCE=" + reference + " TYPE=" + type );
  } else if( pdb.getPositions().size()!=0 ) {
    reference_data = action->getShortcutLabel() + "_atomdata_ref";
    if( displacements ) action->readInputLine( action->getShortcutLabel() + "_atomdata: RMSD DISPLACEMENT SQUARED REFERENCE=" + reference + " TYPE=" + type );
    else action->readInputLine( action->getShortcutLabel() + "_atomdata: RMSD SQUARED REFERENCE=" + reference + " TYPE=" + type );
  } else if( argnames.size()==0 ) {
    argnames.resize( pdb.getArgumentNames().size() );
    for(unsigned i=0; i<argnames.size(); ++i) argnames[i] = pdb.getArgumentNames()[i];
  }
  std::vector<Value*> theargs; if( argnames.size()>0 ) ActionWithArguments::interpretArgumentList( argnames, action->plumed.getActionSet(), action, theargs );

  if( theargs.size()>0 ) {
    std::string instargs, refargs;
    for(unsigned i=0; i<theargs.size(); ++i) {
      std::string iargn = fixArgumentName( theargs[i]->getName() );
      action->readInputLine( action->getShortcutLabel() + "_ref_" + iargn + ": PDB2CONSTANT REFERENCE=" + reference + " ARG=" + theargs[i]->getName() );
      if( i==0 ) { instargs=" ARG1=" + theargs[i]->getName(); refargs=" ARG2=" + action->getShortcutLabel() + "_ref_" + iargn; }
      else { instargs +="," + theargs[i]->getName(); refargs +="," + action->getShortcutLabel() + "_ref_" + iargn; }
      if( pdb.getPositions().size()==0 && i==0 ) reference_data = action->getShortcutLabel() + "_ref_" + iargn;
      else reference_data += "," + action->getShortcutLabel() + "_ref_" + iargn;
    }
    std::string comname="EUCLIDEAN_DISTANCE SQUARED"; std::string coeffstr; action->parse("COEFFICIENTS",coeffstr);
    if( coeffstr.length()>0 ) {
      if( displacements ) action->error("cannot use COEFFICIENTS arguments with GEOMETRIC PATH");
      action->readInputLine( action->getShortcutLabel() + "_coeff: CONSTANT VALUES=" + coeffstr );
      action->readInputLine( action->getShortcutLabel() + "_coeff2: CUSTOM ARG=" + action->getShortcutLabel() + "_coeff FUNC=x*x PERIODIC=NO");
      comname = "NORMALIZED_EUCLIDEAN_DISTANCE SQUARED METRIC=" + action->getShortcutLabel() + "_coeff2";
    } else if( displacements ) comname = "DISPLACEMENT";

    if( pdb.getPositions().size()==0 ) action->readInputLine( action->getShortcutLabel() + "_data: " + comname + instargs + refargs );
    else action->readInputLine( action->getShortcutLabel() + "_argdata: " + comname + instargs + refargs );
  }
}

void Path::readPropertyInformation( const std::vector<std::string>& pnames, const std::string& lab, const std::string& refname, ActionShortcut* action ) {
  if( pnames.size()>0 ) {
    for(unsigned i=0; i<pnames.size(); ++i) {
      action->readInputLine( pnames[i] + ": CONSTANT VALUE=1");
      action->readInputLine( pnames[i] + "_ref: PDB2CONSTANT REFERENCE=" + refname + " ARG=" + pnames[i] );
    }
  } else {
    ActionWithValue* av=action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( lab + "_data" );
    unsigned nfram = av->copyOutput(0)->getShape()[0]; std::string indices = "VALUES=1";
    for(unsigned i=1; i<nfram; ++i) { std::string num; Tools::convert( i+1, num ); indices += "," + num; }
    action->readInputLine( lab + "_ind: CONSTANT " + indices );
  }
}

}
}
