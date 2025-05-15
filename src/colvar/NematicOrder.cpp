/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR NEMATIC_ORDER
/*
This file provides a template for if you want to introduce a new CV.

The descrition of the CV that you include here will appear in the further details and examples
section of the manual

You can (and should) include examples showing how to use your CV as follows:

```plumed
# This should be a sample input.
t: NEMATIC_ORDER ATOMS=1,2
PRINT ARG=t STRIDE=100 FILE=COLVAR
```

*/
//+ENDPLUMEDOC

class NematicOrder : public Colvar {
  bool pbc;
  size_t num_molecules;

public:
  explicit NematicOrder(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(NematicOrder,"NEMATIC_ORDER")

void NematicOrder::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS",
    "The molecular axes (arrows) are specified by pairs of atoms. "
    "For n molecules, therefore, 2*n atom indices have to be provided. "
    "The first half of the atom list contains the head and the second half the tail atoms.");
  keys.setValueDescription("scalar","The nematic order parameter S, S=0 for the isotropic phase and S=1 for the nematic phase)");
}

NematicOrder::NematicOrder(const ActionOptions&ao):
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
void NematicOrder::calculate() {

  // build the 3x3 tensor order parameter
  //  Q_ab = 1/N ∑_i [ 3/2 u_a(i) u_b(i) - 1/2 delta_ab ]
  //       = 1/N ∑_i Q_ab(i)
  // where N is the number of molecules and i runs over 1,...,N
  Matrix<double> Q(3,3);
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

    // Add contribution from molecule i to the nematic order tensor.
    for (int a=0; a<3; a++) {
      for (int b=0; b<3; b++) {
        Q(a,b) += 3.0/2.0 * u(a) * u(b);
      }
      Q(a,a) -= 0.5;
    }
  }
  // compute the average, Q = 1/N ∑_i Q(i)
  Q *= 1.0/num_molecules;

  // diagonalize Q to get the nematic director and the order parameter
  Matrix<double> eigenvecs(3,3);
  std::vector<double> eigenvals(3);

  //// DEBUG
  log.printf("nematic order tensor Q\n");
  log.printf(" %f %f %f\n", Q(0,0), Q(0,1), Q(0,2));
  log.printf(" %f %f %f\n", Q(1,0), Q(1,1), Q(1,2));
  log.printf(" %f %f %f\n", Q(2,0), Q(2,1), Q(2,2));
  log.printf("\n");
  ////

  if(diagMat( Q, eigenvals, eigenvecs )!=0) {
    plumed_merror("diagonalization in NematicOrder failed! This matrix is weird\n");
  };

  //// DEBUG
  log.printf("eigenvectors\n");
  log.printf(" %f %f %f\n", eigenvecs(0,0), eigenvecs(0,1), eigenvecs(0,2));
  log.printf(" %f %f %f\n", eigenvecs(1,0), eigenvecs(1,1), eigenvecs(1,2));
  log.printf(" %f %f %f\n", eigenvecs(2,0), eigenvecs(2,1), eigenvecs(2,2));
  log.printf("\n");
  ////

  // The tensor order parameter Q has three eigenvalues lambda_- <= lambda_0 <= lambda_+
  // and the order parameter is defined as S=lambda_+ or equivalently as S=-2*lambda_0
  // see section 2.5 "Nematic Order Parameter" in
  // R. Eppenga & D. Frenkel,
  // "Monte Carlo study of the isotropic and nematic phases of infinitely thin hard plates",
  // Molecular Physics, 52(6), 1303–1334.
  // https://doi.org/10.1080/00268978400101951
  double order_parameter = eigenvals[2];

  setValue(order_parameter);

  // Now compute the gradients of the order parameter S with respect to the atomic positions.

  // The nematic director n is the eigenvalue belonging to
  // the largest eigenvalue lambda_+ (with index 2).
  Vector director;
  for (size_t a = 0; a < 3; a++) {
    // The eigenvectors are stored as rows in `eigenvecs`.
    director(a) = eigenvecs(2,a);
  }

  // virial = - ∑_i (r_head(i) - r_tail(i)) dS/d(r_head(i))
  Tensor virial;

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

    // Compute the scalar product between the nematic director and the molecular axis
    // cos(angle) = <director,u>
    double cos = dotProduct(director, u);

    //// DEBUG
    log.printf("eigenvalues= %f %f %f \n", eigenvals[0], eigenvals[1], eigenvals[2]);
    log.printf("cos= %f\n", cos);
    log.printf("director       = %f %f %f\n", director(0), director(1), director(2));
    log.printf("molecular axis = %f %f %f\n", u(0), u(1), (2));
    ////

    // gradient on the head atom of the molecular axis, dS/d(r_head(i))
    Vector deriv = 1.0/num_molecules * (3.0/length) * (cos * director - pow(cos,2) * u);
    setAtomsDerivatives(head, deriv);
    // gradient on the tail atom of the molecular axis,
    // dS/d(r_tail(i)) = - dS/d(r_head(i))
    setAtomsDerivatives(tail, -deriv);

    // contribution to virial
    virial -= Tensor(deriv, distance);
  }
  setBoxDerivatives(virial);
}

}
}
