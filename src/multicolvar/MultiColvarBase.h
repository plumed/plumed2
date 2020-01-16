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
#ifndef __PLUMED_multicolvar_MultiColvarBase_h
#define __PLUMED_multicolvar_MultiColvarBase_h

#include <vector>
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"

namespace PLMD {

class ActionSet;
class ActionShortcut;

namespace multicolvar {

class MultiColvarBase :
  public ActionAtomistic,
  public ActionWithValue
{
private:
/// Use periodic boundary conditions
  bool usepbc;
/// The atoms indices of the centers of the group
  std::vector<AtomNumber> mygroup;
/// Vector of forces
  std::vector<double> forcesToApply;
/// Blocks of atom numbers
  std::vector< std::vector<unsigned> > ablocks;
/// Atom numbers for centers
  std::vector<unsigned> catom_indices;
/// Forces on virtual atoms if required
  std::vector<Vector> vatom_forces;
protected:
/// Add a value to the action
  void addValue();
  void addValueWithDerivatives();
/// Add a component to the action
  void addComponent( const std::string& name );
  void addComponentWithDerivatives( const std::string& name );
/// Get the number of atoms involved in each CV
  unsigned getNumberOfAtomsInEachCV() const ;
/// This is used in angles and planes to make three atoms in input into four during execution
  void useFourAtomsForEachCV();
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Set the value of the colvar
  void setValue( const unsigned&, const double&, MultiValue& ) const ;
/// Add some atom derivatives
  void addAtomsDerivatives( const unsigned& ival, const unsigned& jder, const Vector& der, MultiValue& myvals ) const ;
/// Add some box derivatives
  void addBoxDerivatives( const unsigned& ival, const Tensor& vir, MultiValue& myvals ) const ;
/// Calculate the box derivatives using the virial formula
  void setBoxDerivativesNoPbc( const unsigned& ival, const std::vector<Vector>& fpositions, MultiValue& myvals ) const ;
public:
  static void shortcutKeywords( Keywords& keys );
  static void readShortcutKeywords( std::map<std::string,std::string>& keymap, ActionShortcut* action );
  static void expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights,
                               const std::map<std::string,std::string>& keymap, ActionShortcut* action );
  static void expandFunctions( const std::string& labout, const std::string& argin, const std::string& weights, ActionShortcut* action );
  static void registerKeywords( Keywords& keys );
  explicit MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase();
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const ;
/// Do the calculation
  void calculate();
/// Perform one of the tasks
  void performTask( const unsigned&, MultiValue& ) const ;
/// Compute the value of the CV
  virtual void compute( const std::vector<Vector>& pos, MultiValue& myvals ) const = 0 ;
/// Apply the forces from this action
  virtual void apply();
/// The number of virtual atoms that are calculated by this action
  unsigned getNumberOfVirtualAtoms() const ;
};

inline
unsigned MultiColvarBase::getNumberOfVirtualAtoms() const {
  return getFullNumberOfTasks();
}

inline
unsigned MultiColvarBase::getNumberOfDerivatives() const {
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned MultiColvarBase::getNumberOfAtomsInEachCV() const {
  return ablocks.size();
}

inline
void MultiColvarBase::setValue( const unsigned& ival, const double& vv, MultiValue& myvals ) const {
  myvals.setValue( getPntrToOutput(ival)->getPositionInStream(), vv );
}

inline
void MultiColvarBase::addAtomsDerivatives( const unsigned& ival, const unsigned& jder, const Vector& der, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned itask=myvals.getTaskIndex(), jval=getPntrToOutput(ival)->getPositionInStream();
  myvals.addDerivative( jval, 3*ablocks[jder][itask] + 0, der[0] );
  myvals.addDerivative( jval, 3*ablocks[jder][itask] + 1, der[1] );
  myvals.addDerivative( jval, 3*ablocks[jder][itask] + 2, der[2] );
}

inline
void MultiColvarBase::addBoxDerivatives( const unsigned& ival, const Tensor& vir, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned nvir=3*getNumberOfAtoms(), jval=getPntrToOutput(ival)->getPositionInStream();
  for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) myvals.addDerivative( jval, nvir + 3*i+j, vir(i,j) );
}

}
}

#endif
