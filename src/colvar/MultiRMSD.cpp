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

//+PLUMEDOC DCOLVAR MULTI_RMSD
/*
Calculate RMSD distances for different domains and combine them.

This action is largely depracated.  In previous versions of PLUMED a more complex version of this method was implemented.
You can see an example in [example from the plumed nest](https://www.plumed-nest.org/eggs/20/026/).

We felt that the input syntax for the method was not very transparant.  We have thus provided this minimal action
that creates the input for calculating the MultiDomain RMSD for simple cases.  This action is a shortcut.  If you look at the inputs in the
egg that is linked above you can see how we
use the various actions that are in PLUMED to calculate the final quantity.  If you would like to implement some
complicated CVs that are linear combinations of multiple RMSD calculations looking at how this shortcut works will help you start.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

class MultiRMSD : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MultiRMSD(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MultiRMSD,"MULTI_RMSD")

void MultiRMSD::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be MULTI-OPTIMAL, MULTI-OPTIMAL-FAST,  MULTI-SIMPLE or MULTI-DRMSD.");
  keys.addFlag("SQUARED",false," This should be set if you want the mean squared displacement instead of the root mean squared displacement");
  keys.addFlag("NOPBC",false,"don't use periodic boundary conditions");
  keys.setValueDescription("scalar","the sum of the multiple RMSD distances");
  keys.setDeprecated("RMSD");
  keys.needsAction("CONSTANT");
  keys.needsAction("WHOLEMOLECULES");
  keys.needsAction("POSITION");
  keys.needsAction("CONCATENATE");
  keys.needsAction("RMSD_VECTOR");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
}

MultiRMSD::MultiRMSD(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  warning("this action is depracated.  look at the log to see how it is implemented using the new syntax");
  std::string type;
  parse("TYPE",type);
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  std::size_t dash=type.find_first_of("-");
  if( dash!=std::string::npos ) {
    if( type.substr(0,dash)=="MULTI" ) {
      warning("MULTI is deprecated.  You can just use OPTIMAL/SIMPLE");
    } else {
      error("cannot understand type " + type );
    }
    type = type.substr(dash+1);
  }
  std::string reference;
  parse("REFERENCE",reference);
  PDB pdb;
  if( !pdb.read(reference,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + reference );
  }

  unsigned nblocks =  pdb.getNumberOfAtomBlocks();
  if( nblocks<2 ) {
    error("multidomain RMSD only has one block of atoms");
  }
  std::string num;
  std::vector<unsigned> blocks( nblocks+1 );
  blocks[0]=0;
  for(unsigned i=0; i<nblocks; ++i) {
    blocks[i+1]=pdb.getAtomBlockEnds()[i];
  }

  for(unsigned i=1; i<=nblocks; ++i) {
    // Setup a constant
    double asum=0;
    std::string bnum;
    Tools::convert( i, bnum );
    for(unsigned j=blocks[i-1]; j<blocks[i]; ++j) {
      asum += pdb.getOccupancy()[j];
    }
    Vector center;
    center.zero();
    for(unsigned j=blocks[i-1]; j<blocks[i]; ++j) {
      center += ( pdb.getOccupancy()[j] / asum )*pdb.getPositions()[j];
    }
    std::vector<double> vals;
    for(unsigned k=0; k<3; ++k) {
      for(unsigned j=blocks[i-1]; j<blocks[i]; ++j) {
        vals.push_back( pdb.getPositions()[j][k] - center[k] );
      }
    }
    std::string valstr;
    Tools::convert( vals[0], valstr );
    for(unsigned ii=1; ii<vals.size(); ++ii) {
      std::string rnum;
      Tools::convert( vals[ii], rnum );
      valstr += "," + rnum;
    }
    // Create the reference value
    readInputLine( getShortcutLabel() + "_ref" + bnum + ": CONSTANT VALUES=" + valstr );
    // Do whole molecules
    if( !nopbc ) {
      Tools::convert( pdb.getAtomNumbers()[blocks[i-1]].serial(), num );
      std::string wm_line = "WHOLEMOLECULES ENTITY0=" + num;
      for(unsigned j=blocks[i-1]+1; j<blocks[i]; ++j) {
        Tools::convert( pdb.getAtomNumbers()[j].serial(), num );
        wm_line += "," + num;
      }
      readInputLine( wm_line );
    }
    // Get the positions of the atoms in this block
    Tools::convert( pdb.getAtomNumbers()[blocks[i-1]].serial(), num );
    std::string pos_line = getShortcutLabel() + "_cpos" + bnum + ": POSITION NOPBC ATOMS=" + num;
    for(unsigned j=blocks[i-1]+1; j<blocks[i]; ++j) {
      Tools::convert( pdb.getAtomNumbers()[j].serial(), num );
      pos_line += "," + num;
    }
    readInputLine( pos_line );
    // Concatenate the positiosn together
    readInputLine( getShortcutLabel() + "_pos" + bnum + ": CONCATENATE ARG=" + getShortcutLabel() + "_cpos" + bnum + ".x," + getShortcutLabel() + "_cpos" + bnum + ".y," + getShortcutLabel() + "_cpos" + bnum + ".z");
    // Computer the RMSD for this block
    std::string rmsd_line = getShortcutLabel() + "_rmsd" + bnum + ": RMSD_VECTOR SQUARED ARG=" + getShortcutLabel() + "_pos" + bnum + "," + getShortcutLabel() + "_ref" + bnum;
    // Now align
    Tools::convert( pdb.getOccupancy()[blocks[i-1]], num );
    rmsd_line += " ALIGN=" + num;
    for(unsigned j=blocks[i-1]+1; j<blocks[i]; ++j) {
      Tools::convert( pdb.getOccupancy()[j], num );
      rmsd_line += "," + num;
    }
    // And displace
    Tools::convert( pdb.getBeta()[blocks[i-1]], num );
    rmsd_line += " DISPLACE=" + num;
    for(unsigned j=blocks[i-1]+1; j<blocks[i]; ++j) {
      Tools::convert( pdb.getBeta()[j], num );
      rmsd_line += "," + num;
    }
    readInputLine( rmsd_line + " TYPE=" + type );
  }
  std::string argstr = getShortcutLabel() + "_rmsd1";
  for(unsigned i=1; i<nblocks; ++i) {
    std::string bnum;
    Tools::convert( i+1, bnum);
    argstr += "," + getShortcutLabel() + "_rmsd" + bnum;
  }
  bool squared;
  parseFlag("SQUARED",squared);
  if( !squared ) {
    readInputLine( getShortcutLabel() + "_2: COMBINE ARG=" + argstr + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + ": COMBINE ARG=" + argstr + " PERIODIC=NO");
  }
}

}
}
