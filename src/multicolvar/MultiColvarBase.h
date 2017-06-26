/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "tools/DynamicList.h"

namespace PLMD {
namespace multicolvar {

class AtomValuePack;

class MultiColvarBase :
  public ActionAtomistic,
  public ActionWithValue
{
  friend class AtomValuePack;
private:
/// Use periodic boundary conditions
  bool usepbc;
/// Vector of forces
  std::vector<double> forcesToApply;
/// Blocks of atom numbers
  std::vector< std::vector<unsigned> > ablocks;
protected:
/// Add a value to the action
  void addValue();
/// Add a component to the action
  void addComponent( const std::string& name );
/// Get the number of atoms involved in each CV
  unsigned getNumberOfAtomsInEachCV() const ;
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit MultiColvarBase(const ActionOptions&);
  ~MultiColvarBase() {}
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();
/// Buld the current lists of tags
  void buildCurrentTaskList( std::vector<unsigned>& tflags ) const ;
/// Do the calculation
  void calculate();
/// Perform one of the tasks
  void performTask( const unsigned&, MultiValue& ) const ;
/// Compute the value of the CV
  virtual void compute( const unsigned& index, AtomValuePack& myatoms ) const = 0 ;
/// Apply the forces from this action
  virtual void apply();
};

inline
unsigned MultiColvarBase::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms()+9;
}

inline
unsigned MultiColvarBase::getNumberOfAtomsInEachCV() const {
  return ablocks.size();
}

}
}

#endif
