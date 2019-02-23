/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "multicolvar/MultiColvarBase.h"
#include "HBPammObject.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace pamm {

//+PLUMEDOC MCOLVAR HBPAMM_SH
/*
Number of HBPAMM hydrogen bonds formed by each hydrogen atom in the system

\par Examples

*/
//+ENDPLUMEDOC


class HBPammHydrogens : public multicolvar::MultiColvarBase {
private:
  double rcut2;
  unsigned block1upper,block2lower;
  Matrix<HBPammObject> hbpamm_obj;
public:
  static void registerKeywords( Keywords& keys );
  explicit HBPammHydrogens(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(HBPammHydrogens,"HBPAMM_SH")

void HBPammHydrogens::registerKeywords( Keywords& keys ) {
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.add("atoms-1","HYDROGENS","The list of hydrogen atoms that can form part of a hydrogen bond.  The atoms must be specified using a comma separated list.");
  keys.add("atoms-1","SITES","The list of atoms which can be part of a hydrogen bond.  When this command is used the set of atoms that can donate a "
           "hydrogen bond is assumed to be the same as the set of atoms that can form hydrogen bonds.  The atoms involved must be specified "
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
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the atoms
  usespecies=true; weightHasDerivatives=false;
  // Read in hydrogen atom indicees
  std::vector<AtomNumber> all_atoms; parseMultiColvarAtomList("HYDROGENS",-1,all_atoms);
  if( atom_lab.size()==0 ) error("no hydrogens specified in input file");
  // Now create a task list - one task per hydrogen
  unsigned nH = atom_lab.size(); for(unsigned i=0; i<nH; ++i) addTaskToList(i);
  // Read in other atoms in hydrogen bond
  ablocks.resize(1); parseMultiColvarAtomList("SITES",-1,all_atoms);
  if( atom_lab.size()>nH ) {
    block1upper=atom_lab.size() - nH + 1; block2lower=0; ablocks[0].resize( atom_lab.size() - nH );
    for(unsigned i=nH; i<atom_lab.size(); ++i) ablocks[0][i-nH]=i;
  } else {
    parseMultiColvarAtomList("DONORS",-1,all_atoms);
    block1upper=block2lower=atom_lab.size() - nH + 1;
    for(unsigned i=nH; i<atom_lab.size(); ++i) ablocks[0].push_back(i);
    parseMultiColvarAtomList("ACCEPTORS",-1,all_atoms);
    if( atom_lab.size()>(block1upper+nH-1) || (block1upper-1)>0 ) error("no acceptor donor pairs in input specified therefore no hydrogen bonds");
    for(unsigned i=nH+block2lower-1; i<atom_lab.size(); ++i) ablocks[0].push_back(i);
  }
  setupMultiColvarBase( all_atoms );

  double reg; parse("REGULARISE",reg);
  unsigned nnode_t=mybasemulticolvars.size(); if( nnode_t==0 ) nnode_t=1;
  if( nnode_t==1 ) {
    std::string errormsg, desc; parse("CLUSTERS",desc);
    hbpamm_obj.resize(1,1);
    hbpamm_obj(0,0).setup(desc, reg, this, errormsg );
    if( errormsg.length()>0 ) error( errormsg );
  } else {
    unsigned nr=nnode_t, nc=nnode_t;
    hbpamm_obj.resize( nr, nc );
    for(unsigned i=0; i<nr; ++i) {
      // Retrieve the base number
      unsigned ibase;
      if( nc<10 ) {
        ibase=(i+1)*10;
      } else if ( nc<100 ) {
        ibase=(i+1)*100;
      } else {
        error("wow this is an error I never would have expected");
      }

      for(unsigned j=i; j<nc; ++j) {
        std::string errormsg, desc; parseNumbered("CLUSTERS",ibase+j+1,desc);
        if( i==j ) {
          hbpamm_obj(i,j).setup( desc, reg, this, errormsg );
          if( errormsg.length()>0 ) error( errormsg );
        } else {
          hbpamm_obj(i,j).setup( desc, reg, this, errormsg );
          hbpamm_obj(j,i).setup( desc, reg, this, errormsg );
          if( errormsg.length()>0 ) error( errormsg );
        }
      }
    }
  }

  // Set the link cell cutoff
  double sfmax=0;
  for(unsigned i=0; i<hbpamm_obj.ncols(); ++i) {
    for(unsigned j=i; j<hbpamm_obj.nrows(); ++j) {
      double rcut=hbpamm_obj(i,j).get_cutoff();
      if( rcut>sfmax ) { sfmax=rcut; }
    }
  }
  setLinkCellCutoff( sfmax );
  rcut2 = sfmax*sfmax;

  // Setup the multicolvar base
  setupMultiColvarBase( all_atoms );
  // And setup the ActionWithVessel
  checkRead();
}

double HBPammHydrogens::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  double value=0, md_da ;

  // Calculate the coordination number
  for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
    if( i>block1upper ) continue;
    for(unsigned j=1; j<myatoms.getNumberOfAtoms(); ++j) {
      if( i==j || j<block2lower ) continue ;
      // Get the base colvar numbers
      unsigned dno = atom_lab[myatoms.getIndex(i)].first;
      unsigned ano = atom_lab[myatoms.getIndex(j)].first;
      Vector d_da=getSeparation( myatoms.getPosition(i), myatoms.getPosition(j) );
      if ( (md_da=d_da[0]*d_da[0])<rcut2 &&
           (md_da+=d_da[1]*d_da[1])<rcut2 &&
           (md_da+=d_da[2]*d_da[2])<rcut2) value += hbpamm_obj(dno,ano).evaluate( i, j, 0, d_da, sqrt(md_da), myatoms );
    }
  }

  return value;
}

}
}

