/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "core/SetupMolInfo.h"
#include "tools/OpenMP.h"

#include <vector>
#include <string>

using namespace std;

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
  vector<vector<AtomNumber> > groups;
  bool doref;
  vector<Vector> refs;
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
  keys.addFlag("ADDREFERENCE", false, "Set this flag if you want to define a reference position for the first atom of each entity");
  keys.add("numbered", "REF", "Add reference position for first atom of each entity");
}

WholeMolecules::WholeMolecules(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  doref(false)
{
  vector<AtomNumber> merge;
  for(int i=0;; i++) {
    vector<AtomNumber> group;
    parseAtomList("ENTITY",i,group);
    if( group.empty() ) break;
    log.printf("  atoms in entity %d : ",i);
    for(unsigned j=0; j<group.size(); ++j) log.printf("%d ",group[j].serial() );
    log.printf("\n");
    groups.push_back(group);
    merge.insert(merge.end(),group.begin(),group.end());
  }
  // read reference position of first atom of each entity
  parseFlag("ADDREFERENCE", doref);
  if(doref) {
    for(int i=0; i<groups.size(); ++i) {
      vector<double> ref;
      parseNumberedVector("REF",i,ref);
      refs.push_back(Vector(ref[0],ref[1],ref[2]));
      log.printf("  reference position in entity %d : %lf %lf %lf\n",i,ref[0],ref[1],ref[2]);
    }
  }

  // Read residues to align from MOLINFO
  vector<string> resstrings; parseVector("RESIDUES",resstrings);
  if( resstrings.size()>0 ) {
    if( resstrings.size()==1 ) {
      if( resstrings[0]=="all" ) resstrings[0]="all-ter";   // Include terminal groups in alignment
    }
    string moltype; parse("MOLTYPE",moltype);
    if(moltype.length()==0) error("Found RESIDUES keyword without specification of the moleclue - use MOLTYPE");
    std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
    if( moldat.size()==0 ) error("Unable to find MOLINFO in input");
    std::vector< std::vector<AtomNumber> > backatoms;
    moldat[0]->getBackbone( resstrings, moltype, backatoms );
    for(unsigned i=0; i<backatoms.size(); ++i) {
      log.printf("  atoms in entity %u : ", static_cast<unsigned>(groups.size()+1));
      for(unsigned j=0; j<backatoms[i].size(); ++j) log.printf("%d ",backatoms[i][j].serial() );
      log.printf("\n");
      groups.push_back( backatoms[i] );
      merge.insert(merge.end(),backatoms[i].begin(),backatoms[i].end());
    }
  }

  if(groups.size()==0) error("no atom found for WHOLEMOLECULES!");

  checkRead();
  Tools::removeDuplicates(merge);
  requestAtoms(merge);
  doNotRetrieve();
  doNotForce();
}

void WholeMolecules::calculate() {
  if(doref) {
    #pragma omp parallel num_threads(OpenMP::getNumThreads())
    {
      #pragma omp for nowait
      for(unsigned i=0; i<groups.size(); ++i) {
        Vector & first (modifyGlobalPosition(groups[i][0]));
        first = refs[i]+pbcDistance(refs[i],first);
        for(unsigned j=0; j<groups[i].size()-1; ++j) {
          const Vector & first (getGlobalPosition(groups[i][j]));
          Vector & second (modifyGlobalPosition(groups[i][j+1]));
          second=first+pbcDistance(first,second);
        }
      }
    }
  } else {
    for(unsigned i=0; i<groups.size(); ++i) {
      for(unsigned j=0; j<groups[i].size()-1; ++j) {
        const Vector & first (getGlobalPosition(groups[i][j]));
        Vector & second (modifyGlobalPosition(groups[i][j+1]));
        second=first+pbcDistance(first,second);
      }
    }
  }
}



}

}
