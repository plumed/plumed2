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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "tools/PDB.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

class RMSDShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit RMSDShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(RMSDShortcut,"RMSD")

void RMSDShortcut::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  keys.add("compulsory","ALIGN","1.0","the weights to use when aligning to the reference structure");
  keys.add("compulsory","DISPLACE","1.0","the weights to use when calculating the displacement from the reference structure");
  keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("DISPLACEMENT",false,"Calculate the vector of displacements instead of the length of this vector");
  keys.addInputKeyword("optional","ARG","vector/matrix","instead of using the REFERENCE option you can use this action to specify the labels of two actions that you are calculating the RMSD between");
  keys.add("compulsory","NUMBER","0","if there are multiple structures in the pdb file you can specify that you want the RMSD from a specific structure by specifying its place in the file here. If NUMBER=0 then the RMSD from all structures are computed");
  keys.addOutputComponent("disp","DISPLACEMENT","vector/matrix","the vector of displacements for the atoms");
  keys.addOutputComponent("dist","DISPLACEMENT","scalar/vector","the RMSD distance the atoms have moved");
  keys.setValueDescription("scalar/vector","the RMSD distance between the instaneous structure and the reference structure/s that were input");
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_VECTOR");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("WHOLEMOLECULES");
  keys.needsAction("POSITION");
  keys.needsAction("CONCATENATE");
  keys.needsAction("RMSD");
  keys.addDOI("10.1107/S0108767388010128");
}

RMSDShortcut::RMSDShortcut(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string argn;
  parse("ARG",argn);
  if( argn.length()>0 ) {
    readInputLine( getShortcutLabel() + ": RMSD_VECTOR ARG=" + argn + " " + convertInputLineToString() );
    return;
  }
  bool disp;
  parseFlag("DISPLACEMENT",disp);
  std::string reference;
  parse("REFERENCE",reference);
  // Read the reference pdb file
  PDB pdb;
  if( !pdb.read(reference,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength()) ) {
    plumed_merror("missing file " + reference );
  }
  unsigned frame;
  parse("NUMBER",frame);
  unsigned nf=1;
  if( frame==0 ) {
    FILE* fp=std::fopen(reference.c_str(),"r");
    bool do_read=true;
    nf=0;
    while ( do_read ) {
      PDB mypdb;
      do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
      if( !do_read && nf>0 ) {
        break ;
      }
      nf++;
    }
  }
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  // Now create the RMSD object
  std::string rmsd_line = getShortcutLabel() + ": ";
  if( nf==1 && !disp ) {
    rmsd_line += "RMSD_SCALAR REFERENCE=" + reference;
    if(nopbc) {
      rmsd_line += " NOPBC";
    }
  } else {
    std::string ffnum;
    Tools::convert( frame, ffnum );
    readInputLine( getShortcutLabel() + "_ref: PDB2CONSTANT REFERENCE=" + reference + " NUMBER=" + ffnum );
    std::vector<AtomNumber> anum( pdb.getAtomNumbers() );
    if( !nopbc ) {
      std::string num;
      Tools::convert( anum[0].serial(), num );
      std::string wm_line = "WHOLEMOLECULES ENTITY0=" + num;
      for(unsigned i=1; i<anum.size(); ++i) {
        Tools::convert( anum[i].serial(), num );
        wm_line += "," + num;
      }
      readInputLine( wm_line );
    }
    std::string num;
    Tools::convert( anum[0].serial(), num );
    std::string pos_line = getShortcutLabel() + "_cpos: POSITION NOPBC ATOMS=" + num;
    for(unsigned i=1; i<anum.size(); ++i) {
      Tools::convert( anum[i].serial(), num );
      pos_line += "," + num;
    }
    readInputLine( pos_line );
    // Concatenate the three positions together
    readInputLine( getShortcutLabel() + "_pos: CONCATENATE ARG=" + getShortcutLabel() + "_cpos.x," + getShortcutLabel() + "_cpos.y," + getShortcutLabel() + "_cpos.z");
    rmsd_line += "RMSD ARG=" + getShortcutLabel() + "_pos," + getShortcutLabel() + "_ref";
    if( disp ) {
      rmsd_line += " DISPLACEMENT";
    }
    // Now align
    std::vector<double> align( pdb.getOccupancy() );
    Tools::convert( align[0], num );
    rmsd_line += " ALIGN=" + num;
    for(unsigned i=1; i<align.size(); ++i) {
      Tools::convert( align[i], num );
      rmsd_line += "," + num;
    }
    // And displace
    std::vector<double> displace( pdb.getBeta() );
    Tools::convert( displace[0], num );
    rmsd_line += " DISPLACE=" + num;
    for(unsigned i=1; i<displace.size(); ++i) {
      Tools::convert( displace[i], num );
      rmsd_line += "," + num;
    }
  }
  // And create the RMSD object
  bool numder;
  parseFlag("NUMERICAL_DERIVATIVES",numder);
  if(numder && nf==1 && !disp ) {
    rmsd_line += " NUMERICAL_DERIVATIVES";
  } else if( numder ) {
    error("can only use NUMERICAL_DERIVATIVES flag when RMSD is calculating a single scalar value");
  }
  bool squared;
  parseFlag("SQUARED",squared);
  if(squared) {
    rmsd_line += " SQUARED";
  }
  std::string tt;
  parse("TYPE",tt);
  readInputLine( rmsd_line + " TYPE=" + tt );
}

}
}
