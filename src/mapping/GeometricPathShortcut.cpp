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

using namespace std;

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
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.add("optional","PROPERTY","read in path coordinates by finding option with this label in remark of pdb frames");
}

GeometricPathShortcut::GeometricPathShortcut( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  std::vector<std::string> refactions;
  std::string mtype; parse("TYPE",mtype);
  std::string refname; parse("REFERENCE",refname); 
  // Create list of reference configurations that PLUMED will use
  Path::createActionsToComputeDistances( mtype, refname, true, this, refactions );
  // Now get coordinates on spath
  std::string pname, ref_str, coord_str; parse("PROPERTY",pname);
  std::vector<AtomNumber> indices; std::vector<double> alig, disp; 
  FILE* fp=std::fopen(refname.c_str(),"r"); bool do_read=true; double fake_unit=0.1; unsigned nfram = 0;
  while (do_read ) {
      PDB mypdb; do_read=mypdb.readFromFilepointer(fp,false,fake_unit);  // Units don't matter here
      // Break if we are done
      if( !do_read ) break ;
      if( mtype=="OPTIMAL-FAST" || mtype=="OPTIMAL" || mtype=="SIMPLE" ) {
          indices.resize( mypdb.getAtomNumbers().size() );
          for(unsigned i=0;i<indices.size();++i) indices[i]=mypdb.getAtomNumbers()[i];
          alig.resize( mypdb.getOccupancy().size() );
          for(unsigned i=0;i<alig.size();++i) alig[i]=mypdb.getOccupancy()[i];
          disp.resize( mypdb.getBeta().size() );
          for(unsigned i=0;i<disp.size();++i) disp[i]=mypdb.getBeta()[i]; 
      }
      // This creates the coefficients
      if( nfram==0 ) { ref_str = " REFFRAMES=" + refactions[nfram]; } else { ref_str += "," + refactions[nfram]; }
      if( pname.length()>0 ) {
          std::vector<std::string> remarks( mypdb.getRemark() );
          std::string propstr; bool found=Tools::parse( remarks, pname, propstr );
          if( !found ) plumed_merror("could not find property named " + pname + " in input file " + refname );
          if( nfram==0 ) { coord_str = "COORDINATES=" + propstr; } else { coord_str += "," + propstr; }
      } else {
          std::string propstr; Tools::convert( nfram+1, propstr );
          if( nfram==0 ) { coord_str = "COORDINATES=" + propstr; } else { coord_str += "," + propstr; }
      }
      nfram++;
  }
  // Now setup action to compute distances between configurations
  std::string metric;
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
      metric = " METRIC={DIFFERENCE ARG1=arg1 ARG2=arg2}";
  }
  // Create action that computes the geometric path variables
  readInputLine( getShortcutLabel() + ": GEOMETRIC_PATH ARG=" + getShortcutLabel() + "_data " + coord_str + ref_str + metric ); 
}


}
}


