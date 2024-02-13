/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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

namespace PLMD {
namespace setup {

//+PLUMEDOC COLVAR PDB2CONSTANT
/*
Create a constant value from a PDB input file

\par Examples

*/
//+ENDPLUMEDOC

class PDB2Constant : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit PDB2Constant(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PDB2Constant,"PDB2CONSTANT")

void PDB2Constant::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure");
  keys.add("compulsory","NUMBER","0","if there are multiple structures in the pdb file you can specify that you want the RMSD from a specific structure by specifying its place in the file here. If NUMBER=0 then the RMSD from all structures are computed");
}

PDB2Constant::PDB2Constant(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string input; parse("REFERENCE",input);
  unsigned frame; parse("NUMBER",frame);

  FILE* fp=std::fopen(input.c_str(),"r"); bool do_read=true; std::vector<double> vals;
  if(!fp) plumed_merror("could not open reference file " + input); unsigned natoms=0, nframes=0;

  while ( do_read ) {
     PDB mypdb; do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
     if( !do_read && nframes>0 ) break ;

     if( natoms==0 ) natoms = mypdb.getPositions().size();
     else if( mypdb.getPositions().size()!=natoms ) plumed_merror("mismatch between sizes of reference configurations");

     if( nframes+1==frame || frame==0 ) {
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
  std::fclose(fp); std::string rnum; plumed_assert( vals.size()>0 );
  Tools::convert( vals[0], rnum ); std::string valstr = " VALUES=" + rnum;
  for(unsigned i=1; i<vals.size();++i) { Tools::convert( vals[i], rnum ); valstr += "," + rnum; }
  if( frame==0 && nframes>1 ) {
      std::string nc, nr; Tools::convert( nframes, nr ); Tools::convert( 3*natoms, nc );
      plumed.readInputLine( getShortcutLabel() + ": CONSTANT NROWS=" + nr + " NCOLS=" + nc + valstr );
  } else plumed.readInputLine( getShortcutLabel() + ": CONSTANT" + valstr );
}

}
}
