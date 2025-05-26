/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 of Alexander Humeniuk.

   This file is part of the liquid_crystal plumed module.

   The liquid_crystal plumed module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The liquid_crystal plumed module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR FERRONEMATIC_ORDER
/*
Calculate the ferronematic order parameter.

The ferronamtic order parameter P depends on the relative orientations of the molecular
axes. If the axes all point into the same direction, giving rise to a net polarization if
the molecules have permanent dipole moments, P is close to 1. If the molecular axes are
oriented isotropically or are aligned antiparallel so that there is no net polarization,
P is close 0.

The nematic and ferronematic order parameters can be used to distinguish the isotropic,
nematic and ferronematic phases of liquid crystals.

$P$ is length of the average of the molecular axes ($\hat{u}_i$ for $i=1,\ldots,N$)
$$
P = \vert \frac{1}{N} \sum_{i=1}^N \hat{u}_i \vert
$$
Since the molecular axes are unit vectors, P ranges from a minimum of 0 to a maximum of 1.

By adding a bias to the ferronematic order parameter, one can drive a liquid crystal from the
isotropic to the ferronematic phase.

The axis of a rod-like molecule is defined as the distance vector between two atoms,
it points from the tail atom to the head atom.

```plumed
# Assume there are three molecules with 20 atoms each.
# In the first molecule the molecular axis vector points from atom 1 to atom 20,
# in the second molecule it points from atom 21 to atom 40
# and in the third from atom 41 to atom 60.
GROUP LABEL=tails ATOMS=1,21,41
GROUP LABEL=heads ATOMS=20,40,60

# Compute ferronematic order parameter for the three molecules.
P: FERRONEMATIC_ORDER ATOMS=tails,heads
PRINT FILE=colvar ARG=P

# Add a bias to the ferronematic order parameter P.
BIASVALUE ARG=P
```

*/
//+ENDPLUMEDOC

class FerroNematicOrder : public Colvar {
  bool pbc;
  size_t num_molecules;

public:
  explicit FerroNematicOrder(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FerroNematicOrder,"FERRONEMATIC_ORDER")

void FerroNematicOrder::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS",
    "The molecular axes are specified by pairs of atoms. "
    "For N molecules, therefore, 2*N atom indices have to be provided. "
    "The first half of the atom list contains the head and the second half the tail atoms.");
  keys.setValueDescription("scalar",
    "The ferronematic order parameter P, P=0 for the isotropic or nematic(antiparallel) phases and P=1 for the ferronematic(parallel) phase)");
}

FerroNematicOrder::FerroNematicOrder(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size() % 2 != 0) {
    error("Number of atoms must be multiple of 2");
  }
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  num_molecules = atoms.size() / 2;
  log.printf("  %d molecules\n", num_molecules);
  for (size_t i = 0; i < num_molecules; i++) {
    log.printf("  molecular axis for molecule %d points from atom %d to atom %d\n",
      i, atoms[i].serial(), atoms[num_molecules+i].serial());
  }
  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atoms);
}


// calculator
void FerroNematicOrder::calculate() {

  // Polarization vector, average over molecular axes u(i)
  //  polarization_a = 1/N ∑_i u_a(i)
  // where N is the number of molecules and i runs over 1,...,N
  Vector polarization;

  for (size_t i = 0; i < num_molecules; i++) {
    // The axis of a molecule is defined by two atoms at opposite ends of the molecule.
    size_t head = i;
    size_t tail = num_molecules + i;
    // The vector `distance` defines the molecular axis.
    Vector distance;
    if(pbc) {
      distance=pbcDistance(getPosition(tail),getPosition(head));
    } else {
      distance=delta(getPosition(tail),getPosition(head));
    }
    // normalize vector defining the molecular axis
    Vector u = distance / distance.modulo();

    // Add contribution from molecule i to the polarization vector.
    polarization += u;
  }
  // compute the average, polarization_a = 1/N ∑_i u_a(i)
  polarization *= 1.0/num_molecules;
  // The ferronematic order parameter is the length of the polarization vector,
  // P = |polarization|
  double order_parameter = polarization.modulo();
  setValue(order_parameter);

  // Now compute the gradients of the order parameter with respect to the atomic positions.

  // virial = - ∑_i (r_head(i) - r_tail(i)) dP/d(r_head(i))
  Tensor virial;

  // If there is no net polarization (P=0), either because the molecules are perfectly antiparallel
  // or because the orientation is isotropic, the polarization vector is not defined, in this
  // case we set it to the 0-vector.
  Vector director; // = (0,0,0)
  if (order_parameter > 0.0) {
    // unit vector in direction of net polarization
    director = polarization / polarization.modulo();
  }

  for (size_t i = 0; i < num_molecules; i++) {
    // The axis of a molecule is defined by two atoms at opposite ends of the molecule.
    size_t head = i;
    size_t tail = num_molecules + i;
    // The vector `distance` defines the molecular axis.
    Vector distance;
    if(pbc) {
      distance=pbcDistance(getPosition(tail),getPosition(head));
    } else {
      distance=delta(getPosition(tail),getPosition(head));
    }
    // normalize vector defining the molecular axis
    double length = distance.modulo();
    Vector u = distance / length;

    // Compute the scalar product between the ferronematic director (net polarization)
    // and the molecular axis
    // cos(angle) = <director,u>
    double cos = dotProduct(director, u);

    // gradient on the head atom of the molecular axis, dP/d(r_head(i))
    Vector deriv = 1.0/num_molecules * (1.0/length) * (director - cos * u);
    setAtomsDerivatives(head, deriv);
    // gradient on the tail atom of the molecular axis,
    // dP/d(r_tail(i)) = - dP/d(r_head(i))
    setAtomsDerivatives(tail, -deriv);

    // contribution to virial
    virial -= Tensor(deriv, distance);
  }
  setBoxDerivatives(virial);
}

}
}
