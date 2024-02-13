/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "RMSDShortcut.h"
#include "core/ActionRegister.h"
#include "tools/Pdb.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

PLUMED_REGISTER_ACTION(RMSDShortcut,"RMSD")

void RMSDShortcut::registerRMSD(Keywords& keys) {
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

void RMSDShortcut::registerKeywords(Keywords& keys){
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  registerRMSD( keys );
}

unsigned RMSDShortcut::createReferenceConfiguration( const std::string& lab, const std::string& input, PlumedMain& plumed, const unsigned& number, const bool& disp ) {
  FILE* fp=std::fopen(input.c_str(),"r"); bool do_read=true; std::vector<double> vals;
  if(!fp) plumed_merror("could not open reference file " + input); unsigned natoms=0, nframes=0;

  while ( do_read ) {
     PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
     if( !do_read && nframes>0 ) break ;

     if( natoms==0 ) natoms = mypdb.getPositions().size();
     else if( mypdb.getPositions().size()!=natoms ) plumed_merror("mismatch between sizes of reference configurations");

     if( nframes+1==number || number==0 ) {
         std::vector<double> align( mypdb.getOccupancy() );
         double asum=0; for(unsigned i=0;i<align.size();++i) asum += align[i];
         if( asum>epsilon ) {
             double iasum = 1 / asum; for(unsigned i=0;i<align.size();++i) align[i] *= iasum;
         } else {
             double iasum = 1 / mypdb.size(); for(unsigned i=0;i<align.size();++i) align[i] = iasum;
         }
         Vector center; center.zero(); for(unsigned i=0;i<mypdb.getPositions().size();++i) center += align[i]*mypdb.getPositions()[i];

         for(unsigned j=0; j<3; ++j) {
             for(unsigned i=0; i<mypdb.getPositions().size(); ++i) vals.push_back( mypdb.getPositions()[i][j] - center[j] );
         }
     }
     nframes++;
  }
  if( disp ) {
      std::fclose(fp); std::string rnum; plumed_assert( vals.size()>0 );
      Tools::convert( vals[0], rnum ); std::string valstr = " VALUES=" + rnum;
      for(unsigned i=1; i<vals.size();++i) { Tools::convert( vals[i], rnum ); valstr += "," + rnum; }
      if( number==0 && nframes>1 ) {
          std::string nc, nr; Tools::convert( nframes, nr ); Tools::convert( 3*natoms, nc );
          plumed.readInputLine( lab + ": CONSTANT NROWS=" + nr + " NCOLS=" + nc + valstr );
      } else plumed.readInputLine( lab + ": CONSTANT" + valstr );
  }
  return nframes;
}

RMSDShortcut::RMSDShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  bool disp; parseFlag("DISPLACEMENT",disp);
  std::string reference; parse("REFERENCE",reference);
  // Read the reference pdb file
  PDB pdb; if( !pdb.read(reference,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength()) ) plumed_merror("missing file " + reference );
  // Find out the position of the center of mass
  unsigned nf = createReferenceConfiguration( getShortcutLabel() + "_ref", reference, plumed, 1, disp );
  // Now create the RMSD object
  std::string num, rmsd_line = getShortcutLabel() + ": ";
  if( nf==1 && !disp ) {
      rmsd_line += "RMSD_SCALAR REFERENCE=" + reference;
  } else rmsd_line = "RMSD_VECTOR ARG1=" + getShortcutLabel() + "_pos";
  // Now align
  // std::vector<double> align( pdb.getOccupancy() ); Tools::convert( align[0], num ); rmsd_line += " ALIGN=" + num;
  // for(unsigned i=1; i<align.size(); ++i) { Tools::convert( align[i], num ); rmsd_line += "," + num; }
  // // And displace
  // std::vector<double> displace( pdb.getBeta() ); Tools::convert( displace[0], num ); rmsd_line += " DISPLACE=" + num;
  // for(unsigned i=1; i<displace.size(); ++i) { Tools::convert( displace[i], num ); rmsd_line += "," + num; }
  // And create the RMSD object
  bool numder; parseFlag("NUMERICAL_DERIVATIVES",numder); 
  if(numder && nf==1 && !disp ) rmsd_line += " NUMERICAL_DERIVATIVES"; else if( numder ) error("can only use NUMERICAL_DERIVATIVES flag when RMSD is calculating a single scalar value");
  bool squared; parseFlag("SQUARED",squared); if(squared) rmsd_line += " SQUARED";
  bool nopbc; parseFlag("NOPBC",nopbc); if(nopbc) rmsd_line += " NOPBC";
  std::string tt; parse("TYPE",tt); readInputLine( rmsd_line + " TYPE=" + tt );
}

}
}
