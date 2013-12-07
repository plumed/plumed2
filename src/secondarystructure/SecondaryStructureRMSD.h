/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_secondarystructure_SecondaryStructureRMSD_h
#define __PLUMED_secondarystructure_SecondaryStructureRMSD_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "vesselbase/ActionWithVessel.h"
#include <vector>

namespace PLMD {

class SingleDomainRMSD;
class DRMSD;
class RMSD;

namespace secondarystructure {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class SecondaryStructureRMSD : 
  public ActionAtomistic,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
{
private:
/// Are we using pbc
  bool pbcon;
/// Tempory integer to say which refernce configuration is the closest
  unsigned closest;
/// The type of rmsd we are calculating
  std::string alignType;
/// List of all the atoms we require
  DynamicList<AtomNumber> all_atoms;
/// The atoms involved in each of the secondary structure segments
  std::vector< std::vector<unsigned> > colvar_atoms;
/// The reference configurations
  std::vector<RMSD*> secondary_rmsd;
  std::vector<DRMSD*> secondary_drmsd;
/// Stuff for derivatives
  std::vector< std::vector<Vector> > der;
  std::vector<Tensor> vir;
/// Everything for controlling the updating of neighbor lists
  int updateFreq;
  bool firsttime;
/// Variables for strands cutoff
  bool align_strands;
  double s_cutoff;
  unsigned align_atom_1, align_atom_2;
  bool verbose_output;
/// Tempory variables for getting positions of atoms and applying forces
  std::vector<Vector> pos;
  std::vector<double> forcesToApply;
/// Get the index of an atom
  unsigned getAtomIndex( const unsigned& iatom );
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
  SecondaryStructureRMSD(const ActionOptions&);
  virtual ~SecondaryStructureRMSD();
  unsigned getNumberOfFunctionsInAction();
  unsigned getNumberOfDerivatives();
  void prepare();
  void finishTaskListUpdate();
  void calculate();
  void performTask();
  void clearDerivativesAfterTask( const unsigned& );
  void apply();
  void mergeDerivatives( const unsigned& , const double& );
  bool isPeriodic(){ return false; }
};

inline
unsigned SecondaryStructureRMSD::getNumberOfFunctionsInAction(){
  return colvar_atoms.size();
}

inline
unsigned SecondaryStructureRMSD::getNumberOfDerivatives(){
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned SecondaryStructureRMSD::getAtomIndex( const unsigned& iatom ){
  return all_atoms.linkIndex( colvar_atoms[getCurrentTask()][iatom] );
}

}
}

#endif
