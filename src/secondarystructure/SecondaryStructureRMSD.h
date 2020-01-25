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

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "tools/RMSD.h"
#include <vector>

namespace PLMD {

class ActionShortcut;

namespace secondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class SecondaryStructureRMSD :
  public ActionAtomistic,
  public ActionWithValue
{
private:
/// Are we operating without periodic boundary conditions
  bool nopbc;
/// The type of rmsd we are calculating
  std::string alignType;
/// List of all the atoms we require
  std::vector<AtomNumber> all_atoms;
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
/// Tempory variables for getting positions of atoms and applying forces
  std::vector<double> forcesToApply;
/// Get the index of an atom
  unsigned getAtomIndex( const unsigned& current, const unsigned& iatom ) const ;
protected:
/// Get the atoms in the backbone
  void readBackboneAtoms( const std::string& backnames, std::vector<unsigned>& chain_lengths );
/// Add a set of atoms to calculat ethe rmsd from
  void addColvar( const std::vector<unsigned>& newatoms );
/// Set a reference configuration
  void setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units );
/// Setup a pair of atoms to use for strands cutoff
  void setAtomsFromStrands( const unsigned& atom1, const unsigned& atom2 );
/// Finish by setting up the values that will hold what is calculated by this action
  void setupValues();
public:
  static void registerKeywords( Keywords& keys );
  static void readShortcutWords( std::string& ltmap, ActionShortcut* action );
  static void expandShortcut( const std::string& labout, const std::string& labin, const std::string& ltmap, ActionShortcut* action );
  explicit SecondaryStructureRMSD(const ActionOptions&);
  virtual ~SecondaryStructureRMSD();
  unsigned getNumberOfDerivatives() const override ;
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) override;
  void calculate() override;
  void performTask( const unsigned&, MultiValue& ) const override;
  void apply() override;
};

inline
unsigned SecondaryStructureRMSD::getNumberOfDerivatives() const {
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned SecondaryStructureRMSD::getAtomIndex( const unsigned& current, const unsigned& iatom ) const {
  return colvar_atoms[current][iatom];
}

}
}

#endif
