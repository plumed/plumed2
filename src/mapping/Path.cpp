/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "colvar/RMSD.h"
#include "tools/PDB.h"

namespace PLMD {
namespace mapping {

//+PLUMEDOC FUNCTION PATH
/*
Calculate path collective variable given a set of distances from a collection of waymarkers.

This function calculates the Path Collective Variabels that were introduced in \cite brand07.
These variables calculate the system's progress along a curvilinear path ("s" component) and the
perpendicular distance from the curvilinear path ("z" component).

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

PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_ACTION(Path,"GPROPERTYMAP")

void Path::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys ); Path::registerInputFileKeywords( keys );
  keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
}

void Path::registerInputFileKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","ARG","the list of arguments you would like to use in your definition of the path");
  keys.add("optional","COEFFICIENTS","the coefficients of the displacements along each argument that should be used when calculating the euclidean distance");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

Path::Path( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  // Setup the properties
  std::vector<std::string> properties, pnames;
  if( getName()=="PATH") { properties.resize(1); }
  else { parseVector("PROPERTY",pnames); properties.resize( pnames.size() ); }
  // Create list of reference configurations that PLUMED will use
  std::string refname, refargs, metric; 
  std::vector<std::string> argnames; parseVector("ARG",argnames);
  readInputFrames( argnames, refname, false, this, refargs, metric );
  // Now create all other parts of the calculation
  std::string lambda; parse("LAMBDA",lambda);
  // Now create MATHEVAL object to compute exponential functions
  readInputLine( getShortcutLabel() + "_weights: MATHEVAL ARG1=" + getShortcutLabel() + "_data  FUNC=exp(-x*" + lambda + ") PERIODIC=NO" );
  // Create denominator
  readInputLine( getShortcutLabel() + "_denom: SUM ARG=" + getShortcutLabel() + "_weights PERIODIC=NO");
  // Now compte zpath variable
  readInputLine( getShortcutLabel() + "_z: MATHEVAL ARG=" + getShortcutLabel() + "_denom FUNC=-log(x)/" + lambda + " PERIODIC=NO");
  // Now get coefficients for properies for spath
  readPropertyInformation( pnames, getShortcutLabel(), refname, this );
  // Now create COMBINE objects to compute numerator of path
  for(unsigned i=0;i<properties.size();++i) {
      if( pnames.size()>0 ) {
          readInputLine( pnames[i] + "_numer_prod: CUSTOM ARG1=" + getShortcutLabel() + "_weights ARG2=" + pnames[i] + "_ref FUNC=x*y PERIODIC=NO");
          readInputLine( pnames[i] + "_numer: SUM ARG=" + pnames[i]  + "_numer_prod PERIODIC=NO");
          readInputLine( pnames[i] + ": CUSTOM ARG1=" + pnames[i]  + "_numer ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      } else {
          readInputLine( getShortcutLabel() + "_s_prod: CUSTOM ARG1=" + getShortcutLabel() + "_weights ARG2=" + getShortcutLabel() + "_ind FUNC=x*y PERIODIC=NO");
          readInputLine( getShortcutLabel()  + "_numer: SUM ARG=" + getShortcutLabel() + "_s_prod PERIODIC=NO");
          readInputLine( getShortcutLabel() + "_s: CUSTOM ARG1=" + getShortcutLabel() + "_numer ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      }
  }
}

std::string Path::fixArgumentName( const std::string& argin ) {
  std::string argout = argin; std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) argout = argin.substr(0,dot) + "_" + argin.substr(dot+1); 
  return argout;
}

void Path::readArgumentFromPDB( const std::string& argname, const std::string& lab, const std::string& fname, PlumedMain& plmd, const unsigned number ) {
  FILE* fp=std::fopen(fname.c_str(),"r"); bool do_read=true; double fake_unit=0.1; // N.B. units don't matter as we are not reading positions
  if(!fp) plumed_merror("could not open reference file " + fname); std::string strvals; unsigned nframes=0;

  while ( do_read ) {
     PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);
     if( !do_read ) break ; double val;
     if( !mypdb.getArgumentValue(argname,val) ) plumed_merror("did not find argument " + argname + " in file named " + fname );
     if( nframes+1==number || number==0 ) {
         if( strvals.length()==0 ) Tools::convert( val, strvals ); 
         else { std::string rnum; Tools::convert( val, rnum ); strvals += "," + rnum; }
     }
     nframes++;
  }
  plumed_assert( strvals.length()>0 ); plmd.readInputLine( lab + ": CONSTANT_VALUE VALUES=" + strvals );
}

void Path::readPropertyInformation( const std::vector<std::string>& pnames, const std::string& lab, const std::string& refname, ActionShortcut* action ) {
  if( pnames.size()>0 ) {
      for(unsigned i=0;i<pnames.size();++i) readArgumentFromPDB( pnames[i], pnames[i] + "_ref", refname, action->plumed ); 
  } else {
      ActionWithValue* av=action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( lab + "_data" );
      unsigned nfram = av->copyOutput(0)->getShape()[0]; std::string indices = "VALUES=1";
      for(unsigned i=1;i<nfram;++i) { std::string num; Tools::convert( i+1, num ); indices += "," + num; }
      action->readInputLine( lab + "_ind: CONSTANT_VALUE " + indices );
  }    
}

Value* Path::getValueWithLabel( ActionShortcut* action, const std::string& name ) {
  std::size_t dot=name.find("."); ActionWithValue* vv=action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( name.substr(0,dot) );
  if( !vv ) action->error("cannot find value with name " + name );
  if( dot==std::string::npos ) return vv->copyOutput(0);
  if( !vv->exists(name) ) action->error("cannot find value with name " + name );
  return vv->copyOutput( name );
}

void Path::readInputFrames( const std::vector<std::string>& argnames, std::string& refname, const bool& geometric, 
                            ActionShortcut* action, std::string& refargs, std::string& metric ) {
  action->parse("REFERENCE",refname); 

  if( argnames.size()>0 ) {
     // Check that all args are scalars
     for(unsigned i=0; i<argnames.size(); ++i) {
         if( getValueWithLabel( action, argnames[i] )->getRank()>0 ) action->error("arguments in path must be scalars and not vectors or matrices");
     }
     // Create the list of reference values for each argument
     for(unsigned i=0; i<argnames.size(); ++i) {
         if( i==0 ) refargs = fixArgumentName(argnames[i]) + "_ref"; else refargs += "," + fixArgumentName(argnames[i]) + "_ref";
         readArgumentFromPDB( argnames[i], fixArgumentName(argnames[i]) + "_ref", refname, action->plumed );
     }
     // Turn the vectors containing the arganems into comma separated lists
     std::string full_args=argnames[0], full_ref=fixArgumentName(argnames[0]) + "_ref";
     for(unsigned i=1;i<argnames.size(); ++i) { full_args += "," + argnames[i]; full_ref += "," + fixArgumentName(argnames[i]) + "_ref"; }
     std::string comname="EUCLIDEAN_DISTANCE SQUARED"; std::string coeffstr; action->parse("COEFFICIENTS",coeffstr); 
     if( coeffstr.length()>0 ) {
         action->readInputLine( action->getShortcutLabel() + "_coeff: CONSTANT_VALUE VALUES=" + coeffstr );
         action->readInputLine( action->getShortcutLabel() + "_coeff2: CUSTOM ARG=" + action->getShortcutLabel() + "_coeff FUNC=x*x PERIODIC=NO");
         comname = "NORMALIZED_EUCLIDEAN_DISTANCE SQUARED METRIC=" + action->getShortcutLabel() + "_coeff2";
     }
     // We have to calculate the displacements if we are doing a geometric path
     if( geometric ) { metric = "DIFFERENCE"; comname = "DISPLACEMENT"; }
     // And evaluate the Euclidean distance from this set of reference points
     action->readInputLine( action->getShortcutLabel() + "_data: " + comname + " ARG1=" + full_args + " ARG2=" + full_ref ); 
     return;
  }
  std::string mtype; action->parse("TYPE",mtype);
  if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) { 
     // Read the reference positions in from the input file
     refargs = action->getShortcutLabel() + "_ref"; colvar::RMSD::createReferenceConfiguration( refargs, refname, action->plumed ); 
     // Now create the vector that contains the atoms
     FILE* fp=std::fopen(refname.c_str(),"r"); if(!fp) action->error("could not open reference file " + refname );
     double fake_unit=0.1; PDB mypdb; bool do_read=mypdb.readFromFilepointer(fp,false,fake_unit); 
     colvar::RMSD::createPosVector( action->getShortcutLabel() + "_pos", mypdb, action );
     // Create align and displace
     std::string num, align_str, displace_str; Tools::convert( mypdb.getOccupancy()[0], align_str ); Tools::convert( mypdb.getBeta()[0], displace_str );
     for(unsigned j=1; j<mypdb.getAtomNumbers().size(); ++j ) { Tools::convert( mypdb.getOccupancy()[j], num ); align_str += "," + num; Tools::convert( mypdb.getBeta()[0], num ); displace_str += "," + num; }
     // And create the RMSD object
     std::string unfix_str; if( action->getName()=="ADAPTIVE_PATH" ) unfix_str = " UNFIX"; 
     std::string comname = "RMSD_CALC SQUARED " + unfix_str; if( geometric ) { comname = "RMSD_CALC DISPLACEMENT " + unfix_str; metric = comname + " TYPE=" + mtype + " ALIGN=" + align_str + " DISPLACE=" + displace_str; }
     action->readInputLine( action->getShortcutLabel() + "_data: " + comname + " TYPE=" + mtype + " ARG1=" + action->getShortcutLabel() + "_pos ARG2=" + action->getShortcutLabel() + "_ref ALIGN=" + align_str + " DISPLACE=" + displace_str ); 
     return;
  }
  if( mtype.find("DRMSD")!=std::string::npos ) {
     if( geometric ) action->error("DRMSD not available with GPATH shortcut");
     action->readInputLine( action->getShortcutLabel() + "_data: DRMSD SQUARED TYPE=" + mtype + " REFERENCE=" + refname );
     return;
  }
  action->error("metric type " + mtype + " has not been implemented");
}

}
}


