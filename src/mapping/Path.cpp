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
#include "core/ActionSetup.h"
#include "core/ActionRegister.h"
#include "tools/PDB.h"

using namespace std;

namespace PLMD {
namespace mapping {

//+PLUMEDOC FUNCTION PATH
/*
Calculate path collective variable given a set of distances from a collection of waymarkers.

This function calculates the Path Collective Variabels that were introduced in \cite brand07.
These variables calculate the system's progress along a curvilinear path ("s" component) and the
perpendicular distance from the curvilinear path ("z" component).

\par Examples

*/
//+ENDPLUMEDOC


class Path : public ActionSetup {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords(Keywords& keys);
  explicit Path(const ActionOptions&);
};


PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_SHORTCUT(Path,"PATH")
PLUMED_REGISTER_SHORTCUT(Path,"GPROPERTYMAP")

void Path::shortcutKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
}

void Path::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                           const std::map<std::string,std::string>& keys,
                           std::vector<std::vector<std::string> >& actions ) {
  // Check if we need to read in properties from the reference file
  std::vector<std::string> properties, pnames;
  if( words[0]=="GPROPERTYMAP") {
      pnames=Tools::getWords( keys.find("PROPERTY")->second, "\t\n ,");
      properties.resize( pnames.size() );
  } else {
      plumed_assert(words[0]=="PATH"); properties.resize( 1 );
  }
  // Create list of reference configurations that PLUMED will use
  std::vector<AtomNumber> indices; std::vector<double> alig, disp;
  FILE* fp=std::fopen(const_cast<char*>(keys.find("REFERENCE")->second.c_str()),"r");
  if(!fp) plumed_merror("could not open reference file " + keys.find("REFERENCE")->second );
  bool do_read=true; double fake_unit=0.1; unsigned nfram = 0; 
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
      // Break if we are done
      if( !do_read ) break ;
      std::string num; Tools::convert( nfram+1, num );
      std::vector<std::string> ref_input; ref_input.push_back( lab + "_ref" + num + ":" );
      ref_input.push_back("READ_ATOMS"); ref_input.push_back("REFERENCE=" + keys.find("REFERENCE")->second );
      ref_input.push_back("NUMBER=" + num ); actions.push_back( ref_input ); 
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
      if( pnames.size()>0 ) {
         std::vector<std::string> remarks( mypdb.getRemark() );
         for(unsigned i=0; i<pnames.size(); ++i) {
           std::string propstr; bool found=Tools::parse( remarks, pnames[i], propstr );
           if( !found ) plumed_merror("could not find property named " + pnames[i] + " in input file " + keys.find("REFERENCE")->second );
           if( nfram==0 ) { properties[i] = "COEFFICIENTS=" + propstr; } else { properties[i] += "," + propstr; }
         }
      } else {
         std::string propstr; Tools::convert( nfram+1, propstr );
         if( nfram==0 ) { properties[0] = "COEFFICIENTS=" + propstr; } else { properties[0] += "," + propstr; }
      }
      nfram++;   
  }
  // Now create PLUMED object that computes all distances
  std::vector<std::string> ref_line; ref_line.push_back( lab + "_data:" );
  ref_line.push_back("PLUMED_VECTOR");
  for(unsigned i=0;i<nfram;++i) {
      std::string num; Tools::convert(i+1, num );
      if( keys.find("TYPE")->second=="OPTIMAL-FAST" || keys.find("TYPE")->second=="OPTIMAL" || keys.find("TYPE")->second=="SIMPLE" ) {
          std::string inp_line = "INPUT" + num + "= RMSD REFERENCE_ATOMS=" + lab + "_ref" + num;
          std::string atnum; Tools::convert( indices[0].serial(), atnum ); inp_line += " ATOMS=" + atnum;
          for(unsigned i=1;i<indices.size();++i){ Tools::convert( indices[i].serial(), atnum ); inp_line += "," + atnum; } 
          // Get the align values 
          std::string anum; Tools::convert( alig[0], anum ); inp_line += " ALIGN=" + anum;
          for(unsigned i=1;i<alig.size();++i){ Tools::convert( alig[i], anum ); inp_line += "," + anum; }
          // Get the displace values
          std::string dnum; Tools::convert( disp[0], dnum ); inp_line += " DISPLACE=" + dnum;
          for(unsigned i=1;i<disp.size();++i){ Tools::convert( disp[i], dnum ); inp_line += "," + dnum; }
          // Set the type
          inp_line += " TYPE=" + keys.find("TYPE")->second + " SQUARED ";
          ref_line.push_back( inp_line );
      } else {
           plumed_merror("TYPE " + keys.find("TYPE")->second + " is not defined in shortcut for PATH");
      } 
  }
  actions.push_back( ref_line );
  // Now create MATHEVAL object to compute exponential functions
  std::vector<std::string> exp_line; exp_line.push_back(lab + "_weights:");
  exp_line.push_back("MATHEVAL"); exp_line.push_back("ARG1=" + lab + "_data");
  exp_line.push_back("FUNC=exp(-" + keys.find("LAMBDA")->second + "*x)" ); exp_line.push_back("PERIODIC=NO");
  actions.push_back( exp_line );
  // Create denominator
  std::vector<std::string> denom_line; denom_line.push_back( lab + "_denom:");
  denom_line.push_back("COMBINE"); denom_line.push_back("ARG=" + lab + "_weights");
  denom_line.push_back("PERIODIC=NO"); actions.push_back( denom_line );
  // Now compte zpath variable
  std::vector<std::string> zpath_line; zpath_line.push_back( lab + "_z:");
  zpath_line.push_back("MATHEVAL"); zpath_line.push_back("ARG=" + lab + "_denom");
  zpath_line.push_back("FUNC=-log(x) /" + keys.find("LAMBDA")->second );
  zpath_line.push_back("PERIODIC=NO"); actions.push_back( zpath_line );
  // Now create COMBINE objects to compute numerator of path
  for(unsigned i=0;i<properties.size();++i) {
      std::vector<std::string> numer_input, path_input;  
      if( pnames.size()>0 ) { numer_input.push_back( pnames[i] + "_numer:"); path_input.push_back( pnames[i] + ":"); }
      else { numer_input.push_back( lab + "_numer:"); path_input.push_back( lab + "_s:"); }
      // Create numerators for SPATH variables
      numer_input.push_back("COMBINE"); numer_input.push_back("ARG=" + lab + "_weights");
      numer_input.push_back( properties[i] ); numer_input.push_back("PERIODIC=NO");
      actions.push_back( numer_input );
      // Now create spath variables
      path_input.push_back("MATHEVAL"); 
      if( pnames.size()>0 ) path_input.push_back("ARG1=" + pnames[i] + "_numer");
      else path_input.push_back("ARG1=" + lab + "_numer");
      path_input.push_back("ARG2=" + lab + "_denom");
      path_input.push_back("FUNC=x/y"); path_input.push_back("PERIODIC=NO");
      actions.push_back( path_input );
  }
}

void Path::registerKeywords(Keywords& keys) { }

Path::Path(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao)
{
  plumed_error();
}

}
}


