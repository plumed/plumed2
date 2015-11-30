/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "MultiColvar.h"
#include "HBPammObject.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR HBPAMM_PERH
/*
Numberof HBPAMM hydrogen bonds formed by each hydrogen atom in the system

\par Examples

*/
//+ENDPLUMEDOC


class HBPammHydrogens : public MultiColvar {
private:
  double rcut2;
  HBPammObject hbpamm_obj;
public:
  static void registerKeywords( Keywords& keys );
  explicit HBPammHydrogens(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(HBPammHydrogens,"HBPAMM_PERH")

void HBPammHydrogens::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.add("atoms-1","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list.");
  keys.add("atoms-1","SITES","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
                             "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified"
                             "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
                             "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
                             "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
                             "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms-2","DONORS","The list of atoms which can donate a hydrogen bond.  The atoms involved must be specified "
                              "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
                              "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
                              "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
                              "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("atoms-2","ACCEPTORS","The list of atoms which can accept a hydrogen bond.  The atoms involved must be specified "
                                 "as a list of labels of \\ref mcolv or labels of a \\ref multicolvarfunction actions.  If you would just like to use "
                                 "the atomic positions you can use a \\ref DENSITY command to specify a group of atoms.  Specifying your atomic positions using labels of "
                                 "other \\ref mcolv or \\ref multicolvarfunction commands is useful, however, as you can then exploit a much wider "
                                 "variety of functions of the contact matrix as described in \\ref contactmatrix");
  keys.add("numbered","CLUSTERS","the name of the file that contains the definitions of all the kernels for PAMM");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.reset_style("CLUSTERS","compulsory");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("SUM");
}

HBPammHydrogens::HBPammHydrogens(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  usespecies=true; weightHasDerivatives=false;
  int natoms=2; std::vector<AtomNumber> all_atoms;
  readSpeciesKeyword("HYDROGENS","SITES",natoms,all_atoms);

  double reg; parse("REGULARISE",reg);
  if( getNumberOfInputAtomTypes()==1 ){  
      std::string errormsg, desc; parse("CLUSTERS",desc);
      hbpamm_obj.setup(desc, reg, this, errormsg );
      if( errormsg.length()>0 ) error( errormsg );
  } else {
      plumed_error();
  }

  // Set the link cell cutoff
  double sfmax=hbpamm_obj.get_cutoff();
  setLinkCellCutoff( sfmax );
  rcut2 = sfmax*sfmax;

  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms, true );
  // And setup the ActionWithVessel
  checkRead();
}

double HBPammHydrogens::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   double value=0, md_da ; 

   // Calculate the coordination number
   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      for(unsigned j=1;j<myatoms.getNumberOfAtoms();++j){
          if( i==j ) continue ;
          Vector d_da=getSeparation( myatoms.getPosition(i), myatoms.getPosition(j) );
          if ( (md_da=d_da[0]*d_da[0])<rcut2 && 
               (md_da+=d_da[1]*d_da[1])<rcut2 &&
               (md_da+=d_da[2]*d_da[2])<rcut2) value += hbpamm_obj.evaluate( i, j, 0, d_da, sqrt(md_da), myatoms ); 
      }
   }

   return value;
}

}
}

