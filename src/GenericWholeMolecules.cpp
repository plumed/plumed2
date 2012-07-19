/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionAtomistic.h"
#include "ActionPilot.h"
#include "ActionRegister.h"
#include "Vector.h"
#include "AtomNumber.h"
#include "Tools.h"
#include "Atoms.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "MolInfo.h"

#include <vector>
#include <string>

using namespace std;
using namespace PLMD;

namespace PLMD {

//+PLUMEDOC GENERIC WHOLEMOLECULES
/*
This action is used to rebuild molecules that can become split by the periodic
boundary conditions in a manner similar to the ALIGN_ATOMS keyword of plumed1.

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

\par Examples
This command instructs plumed to reconstruct the molecule containing atoms 1-20
at every step of the calculation.

\verbatim
WHOLEMOLECULES STRIDE=1 MOLECULE=1-20
\endverbatim

This command instructs plumed to reconstruct the chain of backbone atoms in a 
protein

\verbatim
MOLINFO REFERENCE=helix.pdb
WHOLEMOLECULES STRIDE=1 RESIDUES=ALL RES_ATOMS=N,CA,CB,C,O
\endverbatim
(See also \ref MOLINFO)

*/
//+ENDPLUMEDOC


class GenericWholeMolecules:
  public ActionPilot,
  public ActionAtomistic
{
  vector<vector<AtomNumber> > groups;
  Vector & modifyPosition(AtomNumber);
public:
  GenericWholeMolecules(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate();
  void apply(){};
};

PLUMED_REGISTER_ACTION(GenericWholeMolecules,"WHOLEMOLECULES")

void GenericWholeMolecules::registerKeywords( Keywords& keys ){
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("numbered","ENTITY","the atoms that make up a molecule that you wish to align. To specify multiple molecules use a list of ENTITY keywords: ENTITY1, ENTITY2,...");
  keys.reset_style("ENTITY","atoms");
  keys.add("atoms","RESIDUES","this command specifies a set of residues which all must be aligned. It must be used in tandem with the \\ref MOLINFO "
                              "action and the RES_ATOMS keyword. If you wish to use all the residues from all the chains in your system you can do so by "
                              "specifying all. Alternatively, if you wish to use a subset of the residues you can specify the particular residues "
                              "you are interested in as a list of numbers"); 
  keys.add("optional","RES_ATOMS","this command tells plumed what atoms should be aligned in each of the residues that are being aligned");
}

inline
Vector & GenericWholeMolecules::modifyPosition(AtomNumber i){
  return atoms.positions[i.index()];
}

GenericWholeMolecules::GenericWholeMolecules(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao)
{
  vector<AtomNumber> merge;
  for(int i=0;;i++){
    vector<AtomNumber> group;
    parseAtomList("ENTITY",i,group); 
    if( group.empty() ) break;
    log.printf("  atoms in entity %d : ",i);
    for(unsigned j=0;j<group.size();++j) log.printf("%d ",group[j].serial() );
    log.printf("\n");
    groups.push_back(group);
    merge.insert(merge.end(),group.begin(),group.end());
  }

  // Read residues to align from MOLINFO
  vector<string> resstrings; parseVector("RESIDUES",resstrings);
  if( resstrings.size()>0 ){
      vector<string> backnames; parseVector("RES_ATOMS",backnames);
      if(backnames.size()==0) error("Found RESIDUES keyword without any specification of the atoms that should be in a residue - use RES_ATOMS");
      std::vector<MolInfo*> moldat=plumed.getActionSet().select<MolInfo*>();
      if( moldat.size()==0 ) error("Unable to find MOLINFO in input");
      std::vector< std::vector<AtomNumber> > backatoms;
      moldat[0]->getBackbone( resstrings, backnames, backatoms );
      for(unsigned i=0;i<backatoms.size();++i){
          log.printf("  atoms in entity %d : ", groups.size()+1 );
          for(unsigned j=0;j<backatoms[i].size();++j) log.printf("%d ",backatoms[i][j].serial() );
          log.printf("\n");
          groups.push_back( backatoms[i] );
          merge.insert(merge.end(),backatoms[i].begin(),backatoms[i].end()); 
      }
  }

  checkRead();
  Tools::removeDuplicates(merge);
  requestAtoms(merge);
}

void GenericWholeMolecules::calculate(){
  for(unsigned i=0;i<groups.size();++i){
    for(unsigned j=0;j<groups[i].size()-1;++j){
      Vector & first (modifyPosition(groups[i][j]));
      Vector & second (modifyPosition(groups[i][j+1]));
      second=first+pbcDistance(first,second);
    }
  }
}



}

