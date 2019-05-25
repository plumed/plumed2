/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "setup/DRMSD.h"
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
The reference frames in the path are defined in the pdb file.  In this frame
each configuration in the path is separated by a line containing just the word END.

\plumedfile
p1: PATH REFERENCE=file.pdb TYPE=OPTIMAL LAMBDA=500.0
PRINT ARG=p1.sss,p1.zzz STRIDE=1 FILE=colvar FMT=%8.4f
\endplumedfile

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
using the algebraic formulae defined earlier.  The positions of the frames in the path are defined
in the file epath.pdb.  An extract from this file looks as shown below.

\verbatim
REMARK ARG=t1,t2 t1=-4.25053  t2=3.88053
END
REMARK ARG=t1,t2 t1=-4.11     t2=3.75
END
REMARK ARG=t1,t2 t1=-3.96947  t2=3.61947
END
\endverbatim

The remarks in this pdb file tell PLUMED the labels that are being used to define the position in the
high dimensional space and the values that these arguments have at each point on the path.

The following input instructs PLUMED to calculate the values of the path collective variables.  The frames that make up this
path are defined in the file all.pdb and all distances are measured using the OPTIMAL metric that is discussed in the manual
page on \ref RMSD.

\plumedfile
p2: PATH REFERENCE=all.pdb LAMBDA=69087
PRINT ARG=p2.spath,p2.zpath STRIDE=1 FILE=colvar
\endplumedfile

If you wish to use collective variable values in the definition of your path you would use an input file with something like this:

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4a
p2: PATH REFERENCE=mypath.pdb LAMBDA=2 TYPE=EUCLIDEAN
PRINT ARG=p2.spath,p2.zpath STRIDE=1 FILE=colvar
\endplumedfile

The corresponding pdb file containing the  definitions of the frames in the path would then look like this:

\verbatim
DESCRIPTION: a defintiion of a PATH
REMARK TYPE=EUCLIDEAN
REMARK ARG=d1,d2
REMARK d1=1.0 d2=1.0
END
REMARK TYPE=EUCLIDEAN
REMARK ARG=d1,d2
REMARK d1=2.0 d2=2.0
END
\endverbatim

For each frame in the path you must specify the arguments that should be used to calculate the distance between the instantaneous configuration
of the system and the reference configurations together with the values that these arguments take in each of the reference configurations.

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_ACTION(Path,"GPROPERTYMAP")

void Path::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","ARG","the list of arguments you would like to use in your definition of the path");
  keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
}

Path::Path( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  // Check if we need to read in properties from the reference file
  std::vector<std::string> properties, pnames;
  if( getName()=="GPROPERTYMAP") {
      parseVector("PROPERTY",pnames); properties.resize( pnames.size() );
  } else {
      plumed_assert(getName()=="PATH"); properties.resize( 1 );
  }
  std::vector<std::string> argnames; parseVector("ARG",argnames);
  std::string mtype; parse("TYPE",mtype);
  if( argnames.size()>0 && mtype=="OPTIMAL_FAST" ) mtype="EUCLIDEAN";

  std::string refname; parse("REFERENCE",refname); 
  // Create list of reference configurations that PLUMED will use
  std::vector<std::string> refactions;
  createActionsToComputeDistances( mtype, refname, false, this, argnames, refactions );
  // Now create all other parts of the calculation
  std::string lambda; parse("LAMBDA",lambda);
  // Now create MATHEVAL object to compute exponential functions
  readInputLine( getShortcutLabel() + "_weights: MATHEVAL ARG1=" + getShortcutLabel() + "_data  FUNC=exp(-x*" + lambda + ") PERIODIC=NO" );
  // Create denominator
  readInputLine( getShortcutLabel() + "_denom: COMBINE ARG=" + getShortcutLabel() + "_weights PERIODIC=NO");
  // Now compte zpath variable
  readInputLine( getShortcutLabel() + "_z: MATHEVAL ARG=" + getShortcutLabel() + "_denom FUNC=-log(x)/" + lambda + " PERIODIC=NO");
  // Now get coefficients for properies for spath
  FILE* fp=std::fopen(refname.c_str(),"r"); bool do_read=true; double fake_unit=0.1; unsigned nfram = 0;
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
      // Break if we are done
      if( !do_read ) break ;
      // This creates the coefficients
      if( pnames.size()>0 ) {
          double pval; std::string propstr;
          for(unsigned i=0; i<pnames.size(); ++i) {
              if( !mypdb.getArgumentValue(pnames[i], pval) ) plumed_merror("could not find property named " + pnames[i] + " in input file " + refname );
              Tools::convert( pval, propstr ); 
              if( nfram==0 ) { properties[i] = "COEFFICIENTS=" + propstr; } else { properties[i] += "," + propstr; }
          }
      } else {
          std::string propstr; Tools::convert( nfram+1, propstr );
          if( nfram==0 ) { properties[0] = "COEFFICIENTS=" + propstr; } else { properties[0] += "," + propstr; }
      }
      nfram++;
  }
  // Now create COMBINE objects to compute numerator of path
  for(unsigned i=0;i<properties.size();++i) {
      std::string numer_input, path_input; 
      if( pnames.size()>0 ) {
          numer_input = pnames[i] + "_numer:";
          path_input = pnames[i] + ": MATHEVAL ARG1=" + pnames[i] + "_numer";
      } else {
          numer_input = getShortcutLabel()  + "_numer:";
          path_input = getShortcutLabel() + "_s: MATHEVAL ARG1=" + getShortcutLabel() + "_numer";
      }
      // Create numerators for SPATH variables
      readInputLine( numer_input + " COMBINE ARG=" + getShortcutLabel() + "_weights PERIODIC=NO " + properties[i] );
      // Create final values of SPATH variables
      readInputLine( path_input + " ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  }
}

void Path::createActionsToComputeDistances( const std::string mtype, const std::string& refname, const bool& geometric, 
                                            ActionShortcut* action, const std::vector<std::string>& argnames, std::vector<std::string>& refactions ) {
  std::vector<AtomNumber> indices; std::vector<double> alig, disp; std::string distances_str;
  FILE* fp=std::fopen(refname.c_str(),"r"); std::string scut_lab = action->getShortcutLabel();
  if(!fp) action->error("could not open reference file " + refname );
  bool do_read=true; double fake_unit=0.1; unsigned nfram = 0; std::string argstr; 
  if( argnames.size()>0 ) {
      argstr=" ARG=" + argnames[0]; for(unsigned i=1;i<argnames.size();++i) argstr += "," + argnames[i];
  }
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
      // Break if we are done
      if( !do_read ) break ;
      std::string num, stri; Tools::convert( nfram+1, num );
      action->readInputLine( scut_lab + "_ref" + num + ": READ_CONFIG REFERENCE=" + refname  + " NUMBER=" + num  + argstr );

      if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) { 
          if( nfram==0 ) {
              indices.resize( mypdb.getAtomNumbers().size() );
              for(unsigned i=0;i<indices.size();++i) indices[i]=mypdb.getAtomNumbers()[i];
              alig.resize( mypdb.getOccupancy().size() );
              for(unsigned i=0;i<alig.size();++i) alig[i]=mypdb.getOccupancy()[i];
              disp.resize( mypdb.getBeta().size() );
              for(unsigned i=0;i<disp.size();++i) disp[i]=mypdb.getBeta()[i]; 
          } else {
              if( indices.size()!=mypdb.getAtomNumbers().size() ) plumed_merror("mismatch between numbers of atoms in frames of path");
              for(unsigned i=0;i<indices.size();++i) {
                  if( indices[i]!=mypdb.getAtomNumbers()[i] ) plumed_merror("mismatch between atom numbers in frames of path");
                  if( alig[i]!=mypdb.getOccupancy()[i] ) plumed_merror("mismatch between occupancies in frames of path");
                  if( disp[i]!=mypdb.getBeta()[i] ) plumed_merror("mismatch between beta values in frames of path");
              }
          }
          refactions.push_back( scut_lab + "_ref" + num );
      } else if( mtype.find("DRMSD")!=std::string::npos ) {
          distances_str = setup::DRMSD::getDistancesString( action->plumed, scut_lab + "_ref" + num, mtype );
          action->readInputLine( scut_lab + "_refv" + num + ": CALCULATE_REFERENCE CONFIG=" + scut_lab + "_ref" + num + " INPUT={DISTANCE " + distances_str + "}" );
          refactions.push_back( scut_lab + "_refv" + num );
      } else if( argnames.size()==0 ) {
          action->readInputLine( scut_lab + "_refv" + num + ": CALCULATE_REFERENCE CONFIG=" + scut_lab + "_ref" + num + " INPUT=" + mtype );
          refactions.push_back( scut_lab + "_refv" + num );
      } else {
          refactions.push_back( scut_lab + "_ref" + num ); 
      }
      nfram++; if( refactions.size()!=nfram ) action->error("mismatch between number of reference action labels and number of frames");
  }
  unsigned nquantities=0;
  if( mtype!="OPTIMAL-FAST" && mtype!="OPTIMAL" && mtype!="SIMPLE" && mtype!="EUCLIDEAN" ) { 
      if( mtype.find("DRMSD")!=std::string::npos ) action->readInputLine( scut_lab + "_instantaneous: DISTANCE " + distances_str );
      else action->readInputLine( scut_lab + "_instantaneous: " + mtype );
      ActionWithValue* aval = action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( scut_lab + "_instantaneous" );
      nquantities = aval->copyOutput(0)->getNumberOfValues( scut_lab + "_instantaneous" );
  }
  // Now create PLUMED object that computes all distances
  std::string ref_line =  scut_lab + "_data: PLUMED_VECTOR ";
  for(unsigned i=0;i<nfram;++i) {
      std::string num; Tools::convert(i+1, num );
      if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
          ref_line += " INPUT" + num + "={RMSD REFERENCE_ATOMS=" + scut_lab + "_ref" + num;
          if( geometric ) ref_line += " DISPLACEMENT";
          std::string atnum; Tools::convert( indices[0].serial(), atnum ); ref_line += " ATOMS=" + atnum;
          for(unsigned i=1;i<indices.size();++i){ Tools::convert( indices[i].serial(), atnum ); ref_line += "," + atnum; } 
          // Get the align values 
          std::string anum; Tools::convert( alig[0], anum ); ref_line += " ALIGN=" + anum;
          for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); ref_line += "," + anum; }
          // Get the displace values
          std::string dnum; Tools::convert( disp[0], dnum ); ref_line += " DISPLACE=" + dnum;
          for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); ref_line += "," + dnum; }
          // Set the type
          ref_line += " TYPE=" + mtype + " SQUARED}";
      } else {
          ref_line += "INPUT" + num + "={" + scut_lab + "_diff" + num + ": DIFFERENCE ARG2=" + refactions[i]; 
          if( mtype!="EUCLIDEAN" ) {
            ref_line += " ARG1=" + scut_lab + "_instantaneous";
          } else {
            ActionWithValue* av = action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( refactions[i] );
            plumed_assert( av ); nquantities = av->copyOutput(0)->getNumberOfValues( av->getLabel() ); 
            ref_line +=" ARG1=" + argnames[0]; for(unsigned i=1;i<argnames.size();++i) ref_line += "," + argnames[i]; 
          }
          std::string powstr = "POWERS=2"; for(unsigned i=1;i<nquantities;++i) powstr += ",2";
          if( mtype=="DRMSD" ) powstr += " NORMALIZE"; 
          if( !geometric) ref_line += "; COMBINE ARG=" + scut_lab + "_diff" + num + " PERIODIC=NO " + powstr + "} ";
          else ref_line += "} ";
      } 
  }
  action->readInputLine( ref_line ); 
}

unsigned Path::getNumberOfFramesAndMetric( const std::string& mtype, const std::string& reffile, std::string& metric ) {
  std::vector<AtomNumber> indices; std::vector<double> alig, disp; 
  FILE* fp=std::fopen(reffile.c_str(),"r"); bool do_read=true; double fake_unit=0.1; unsigned nframes = 0;
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
      if( !do_read ) break ;
      if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
          indices.resize( mypdb.getAtomNumbers().size() );
          for(unsigned i=0;i<indices.size();++i) indices[i]=mypdb.getAtomNumbers()[i];
          alig.resize( mypdb.getOccupancy().size() );
          for(unsigned i=0;i<alig.size();++i) alig[i]=mypdb.getOccupancy()[i]; 
          disp.resize( mypdb.getBeta().size() );
          for(unsigned i=0;i<disp.size();++i) disp[i]=mypdb.getBeta()[i];
      }
      nframes++;
  }   
  // Now setup action to compute distances between configurations
  if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
      std::string atnum; Tools::convert( indices[0].serial(), atnum ); metric  = " METRIC={RMSD REFERENCE_ATOMS=" + atnum;
      for(unsigned i=1;i<alig.size();++i){ Tools::convert(indices[i].serial(), atnum); metric += "," + atnum; }
      unsigned natoms=indices[0].serial();
      for(unsigned i=1;i<indices.size();++i) {
          if( indices[i].serial()>natoms ) natoms = indices[i].serial();
      }
      Tools::convert( natoms+indices[0].serial(), atnum ); metric += " ATOMS=" + atnum;
      for(unsigned i=1;i<alig.size();++i){ Tools::convert(natoms+indices[i].serial(), atnum); metric += "," + atnum; }
      std::string anum; Tools::convert( alig[0], anum ); metric += " ALIGN=" + anum;
      for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); metric += "," + anum; }
      // Get the displace values
      std::string dnum; Tools::convert( disp[0], dnum ); metric += " DISPLACE=" + dnum;
      for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); metric += "," + dnum; }
      metric += " TYPE=" + mtype + " DISPLACEMENT}";
  } else {
      metric = " METRIC={DIFFERENCE ARG1=arg2 ARG2=arg1}";
  }
  return nframes;
}

}
}


