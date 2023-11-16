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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "tools/OpenMP.h"
#include "tools/Tree.h"

#include <vector>
#include <string>

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC WHOLEMOLECULES
/*
This action is used to rebuild molecules that can become split by the periodic boundary conditions.

It is similar to the ALIGN_ATOMS keyword of plumed1, and is needed since some
MD dynamics code (e.g. GROMACS) can break molecules during the calculation.

Running some CVs without this command can cause there to be discontinuities changes
in the CV value and artifacts in the calculations.  This command can be applied
more than once.  To see what effect is has use a variable without pbc or use
the \ref DUMPATOMS directive to output the atomic positions.

\attention
This directive modifies the stored position at the precise moment
it is executed. This means that only collective variables
which are below it in the input script will see the corrected positions.
As a general rule, put it at the top of the input file. Also, unless you
know exactly what you are doing, leave the default stride (1), so that
this action is performed at every MD step.

The way WHOLEMOLECULES modifies each of the listed entities is this:
- First atom of the list is left in place
- Each atom of the list is shifted by a lattice vectors so that it becomes as close as possible
  to the previous one, iteratively.

In this way, if an entity consists of a list of atoms such that consecutive atoms in the
list are always closer than half a box side the entity will become whole.
This can be usually achieved selecting consecutive atoms (1-100), but it is also possible
to skip some atoms, provided consecutive chosen atoms are close enough.

\par Examples

This command instructs plumed to reconstruct the molecule containing atoms 1-20
at every step of the calculation and dump them on a file.

\plumedfile
# to see the effect, one could dump the atoms as they were before molecule reconstruction:
# DUMPATOMS FILE=dump-broken.xyz ATOMS=1-20
WHOLEMOLECULES ENTITY0=1-20
DUMPATOMS FILE=dump.xyz ATOMS=1-20
\endplumedfile

This command instructs plumed to reconstruct two molecules containing atoms 1-20 and 30-40

\plumedfile
WHOLEMOLECULES ENTITY0=1-20 ENTITY1=30-40
DUMPATOMS FILE=dump.xyz ATOMS=1-20,30-40
\endplumedfile

This command instructs plumed to reconstruct the chain of backbone atoms in a
protein

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO STRUCTURE=helix.pdb
WHOLEMOLECULES RESIDUES=all MOLTYPE=protein
\endplumedfile

*/
//+ENDPLUMEDOC


class WholeMolecules:
  public ActionPilot,
  public ActionAtomistic
{
  std::vector<std::vector<AtomNumber> > groups;
  std::vector<std::vector<AtomNumber> > roots;
  std::vector<Vector> refs;
  bool doemst, addref;
public:
  explicit WholeMolecules(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate() override;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(WholeMolecules,"WHOLEMOLECULES")

void WholeMolecules::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("numbered","ENTITY","the atoms that make up a molecule that you wish to align. To specify multiple molecules use a list of ENTITY keywords: ENTITY0, ENTITY1,...");
  keys.reset_style("ENTITY","atoms");
  keys.add("residues","RESIDUES","this command specifies that the backbone atoms in a set of residues all must be aligned. It must be used in tandem with the \\ref MOLINFO "
           "action and the MOLTYPE keyword. If you wish to use all the residues from all the chains in your system you can do so by "
           "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
           "you are interested in as a list of numbers");
  keys.add("optional","MOLTYPE","the type of molecule that is under study.  This is used to define the backbone atoms");
  keys.addFlag("EMST", false, "Define atoms sequence in entities using an Euclidean minimum spanning tree");
  keys.addFlag("ADDREFERENCE", false, "Define the reference position of the first atom of each entity using a PDB file");
}

WholeMolecules::WholeMolecules(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  doemst(false), addref(false)
{
  // parse optional flags
  parseFlag("EMST", doemst);
  parseFlag("ADDREFERENCE", addref);

  // create groups from ENTITY
  for(int i=0;; i++) {
    std::vector<AtomNumber> group;
    parseAtomList("ENTITY",i,group);
    if( group.empty() ) break;
    groups.push_back(group);
  }

  // Read residues to align from MOLINFO
  std::vector<std::string> resstrings; parseVector("RESIDUES",resstrings);
  if( resstrings.size()>0 ) {
    if( resstrings.size()==1 ) {
      if( resstrings[0]=="all" ) resstrings[0]="all-ter";   // Include terminal groups in alignment
    }
    std::string moltype; parse("MOLTYPE",moltype);
    if(moltype.length()==0) error("Found RESIDUES keyword without specification of the molecule - use MOLTYPE");
    auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
    if( !moldat ) error("MOLINFO is required to use RESIDUES");
    std::vector< std::vector<AtomNumber> > backatoms;
    moldat->getBackbone( resstrings, moltype, backatoms );
    for(unsigned i=0; i<backatoms.size(); ++i) {
      groups.push_back( backatoms[i] );
    }
  }

  // check number of groups
  if(groups.size()==0) error("no atoms found for WHOLEMOLECULES!");

  // if using PDBs reorder atoms in groups based on proximity in PDB file
  if(doemst) {
    auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
    if( !moldat ) error("MOLINFO is required to use EMST");
    // initialize tree
    Tree tree = Tree(moldat);
    // cycle on groups and reorder atoms
    for(unsigned i=0; i<groups.size(); ++i) {
      groups[i] = tree.getTree(groups[i]);
      // store root atoms
      roots.push_back(tree.getRoot());
    }
  } else {
    // fill root vector with previous atom in groups
    for(unsigned i=0; i<groups.size(); ++i) {
      std::vector<AtomNumber> root;
      for(unsigned j=0; j<groups[i].size()-1; ++j) root.push_back(groups[i][j]);
      // store root atoms
      roots.push_back(root);
    }
  }

  // adding reference if needed
  if(addref) {
    auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
    if( !moldat ) error("MOLINFO is required to use ADDREFERENCE");
    for(unsigned i=0; i<groups.size(); ++i) {
      // add reference position of first atom in entity
      refs.push_back(moldat->getPosition(groups[i][0]));
    }
  }

  // print out info
  for(unsigned i=0; i<groups.size(); ++i) {
    log.printf("  atoms in entity %d : ",i);
    for(unsigned j=0; j<groups[i].size(); ++j) log.printf("%d ",groups[i][j].serial() );
    log.printf("\n");
    if(addref) log.printf("     with reference position : %lf %lf %lf\n",refs[i][0],refs[i][1],refs[i][2]);
  }

  // collect all atoms
  std::vector<AtomNumber> merge;
  for(unsigned i=0; i<groups.size(); ++i) {
    merge.insert(merge.end(),groups[i].begin(),groups[i].end());
  }

  checkRead();
  Tools::removeDuplicates(merge);
  requestAtoms(merge);
  doNotRetrieve();
  doNotForce();
}

void WholeMolecules::calculate() {
  for(unsigned i=0; i<groups.size(); ++i) {
    if(addref) {
      Vector & first (modifyGlobalPosition(groups[i][0]));
      first = refs[i]+pbcDistance(refs[i],first);
    }
    for(unsigned j=0; j<groups[i].size()-1; ++j) {
      const Vector & first (getGlobalPosition(roots[i][j]));
      Vector & second (modifyGlobalPosition(groups[i][j+1]));
      second=first+pbcDistance(first,second);
    }
  }
}


}
}
