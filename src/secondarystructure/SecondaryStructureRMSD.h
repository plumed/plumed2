/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#ifndef __PLUMED_secondarystructure_SecondaryStructureRMSD_h
#define __PLUMED_secondarystructure_SecondaryStructureRMSD_h

#include "core/ActionWithVector.h"
#include "tools/RMSD.h"
#include <vector>

namespace PLMD {

class ActionShortcut;

namespace secondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class SecondaryStructureRMSD : public ActionWithVector {
private:
/// Are we operating without periodic boundary conditions
  bool nopbc;
/// The type of rmsd we are calculating
  std::string alignType;
/// List of all the atoms we require
  std::vector<AtomNumber> all_atoms;
/// The list of tasks that we need to do on the round
  std::vector<unsigned> ss_tasks;
/// The atoms involved in each of the secondary structure segments
  std::vector< std::vector<unsigned> > colvar_atoms;
/// The list of reference configurations
  std::vector<RMSD> myrmsd;
  std::vector<std::map<std::pair<unsigned,unsigned>, double> > drmsd_targets;
/// Variables for strands cutoff
  bool align_strands;
  double s_cutoff2;
  unsigned align_atom_1, align_atom_2;
  bool verbose_output;
/// Get the index of an atom
  unsigned getAtomIndex( const unsigned& current, const unsigned& iatom ) const ;
public:
  static void registerKeywords( Keywords& keys );
  static void readBackboneAtoms( ActionShortcut* action, PlumedMain& plumed, const std::string& backnames, std::vector<unsigned>& chain_lengths, std::string& all_atoms );
  static bool readShortcutWords( std::string& ltmap, ActionShortcut* action );
  static void expandShortcut( const bool& uselessthan, const std::string& labout, const std::string& labin, const std::string& ltmap, ActionShortcut* action );
  explicit SecondaryStructureRMSD(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) override;
  int checkTaskStatus( const unsigned& taskno, int& flag ) const override;
  void calculate() override;
  void performTask( const unsigned&, MultiValue& ) const override;
};

inline
unsigned SecondaryStructureRMSD::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned SecondaryStructureRMSD::getAtomIndex( const unsigned& current, const unsigned& iatom ) const {
  return colvar_atoms[current][iatom];
}

}
}

#endif
