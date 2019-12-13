/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "vesselbase/ActionWithVessel.h"
#include <vector>

namespace PLMD {

class SingleDomainRMSD;

namespace secondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class SecondaryStructureRMSD :
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
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
  std::vector<std::unique_ptr<SingleDomainRMSD>> references;
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
public:
  static void registerKeywords( Keywords& keys );
  explicit SecondaryStructureRMSD(const ActionOptions&);
  virtual ~SecondaryStructureRMSD();
  unsigned getNumberOfFunctionsInAction();
  unsigned getNumberOfDerivatives() override;
  unsigned getNumberOfQuantities() const override;
  void turnOnDerivatives() override;
  void calculate() override;
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override;
  void apply() override;
  bool isPeriodic() override { return false; }
};

inline
unsigned SecondaryStructureRMSD::getNumberOfQuantities() const {
  return 1 + references.size();
}


inline
unsigned SecondaryStructureRMSD::getNumberOfFunctionsInAction() {
  return colvar_atoms.size();
}

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
