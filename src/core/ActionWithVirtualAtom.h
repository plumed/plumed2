/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_ActionWithVirtualAtom_h
#define __PLUMED_core_ActionWithVirtualAtom_h

#include "ActionAtomistic.h"
#include "tools/AtomNumber.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "Atoms.h"

namespace PLMD {

/**
\ingroup INHERIT
Inherit from here if you are calculating the position of a virtual atom (eg a center of mass)
*/

/// Class to add a single virtual atom to the system.
/// (it might be extended to add multiple virtual atoms).
class ActionWithVirtualAtom:
  public ActionAtomistic
{
  AtomNumber index;
  std::vector<Tensor> derivatives;
  std::vector<Tensor> boxDerivatives;
  std::map<AtomNumber,Tensor> gradients;
  void apply() override;
protected:
/// Set position of the virtual atom
  void setPosition(const Vector &);
/// Set its mass
  void setMass(double);
/// Set its charge
  void setCharge(double);
/// Request atoms on which the calculation depends
  void requestAtoms(const std::vector<AtomNumber> & a);
/// Set the derivatives of virtual atom coordinate wrt atoms on which it dependes
  void setAtomsDerivatives(const std::vector<Tensor> &d);
/// Set the box derivatives.
/// This should be a vector of size 3. First index corresponds
/// to the components of the virtual atom.
/// Notice that this routine subtract the trivial term coming from cell deformation
/// since this term is already implicitly included. Indeed, if the vatom
/// position is a linear function of atomic coordinates it is not necessary
/// to call this function (implicit term is fine) (e.g. vatom::COM and vatom::Center).
/// On the other hand if the vatom position is a non-linear function of atomic coordinates this
/// should be called (see vatom::Ghost).
  void setBoxDerivatives(const std::vector<Tensor> &d);
/// Set box derivatives automatically.
/// It should be called after the settomsDerivatives has been used for all
/// single atoms.
/// \warning It only works for virtual atoms NOT using PBCs!
///          This implies that all atoms used + the new virtual atom should be
///          in the same periodic image.
  void setBoxDerivativesNoPbc();
public:
  void setGradients();
  const std::map<AtomNumber,Tensor> & getGradients()const;
/// Return the atom id of the corresponding virtual atom
  AtomNumber getIndex()const;
  explicit ActionWithVirtualAtom(const ActionOptions&ao);
  ~ActionWithVirtualAtom();
  static void registerKeywords(Keywords& keys);
  void setGradientsIfNeeded();
};

inline
AtomNumber ActionWithVirtualAtom::getIndex()const {
  return index;
}

inline
void ActionWithVirtualAtom::setPosition(const Vector & pos) {
  atoms.positions[index.index()]=pos;
}

inline
void ActionWithVirtualAtom::setMass(double m) {
  atoms.masses[index.index()]=m;
}

inline
void ActionWithVirtualAtom::setCharge(double c) {
  atoms.charges[index.index()]=c;
}

inline
void ActionWithVirtualAtom::setAtomsDerivatives(const std::vector<Tensor> &d) {
  derivatives=d;
}

inline
const std::map<AtomNumber,Tensor> & ActionWithVirtualAtom::getGradients()const {
  return gradients;
}

}

#endif
