/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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

/*
 This vast majority of the source code in this file was written by
 Sandro Bottaro with some help from Giovanni Bussi
*/

#include "Colvar.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/PDB.h"
#include "tools/ERMSD.h"
#include <cmath>

namespace PLMD {
namespace colvar {


//+PLUMEDOC COLVAR ERMSD
/*
Calculate eRMSD with respect to a reference structure.

eRMSD is a metric developed for measuring distances between three-dimensional RNA structures.
The standard RMSD measure is highly inaccurate when measuring distances among three-dimensional
structures of nucleic acids.
It is not unusual, for example, that two RNA structures with low RMSD (i.e. less than 0.4nm) display a completely different network of base-base interactions.

eRMSD measures the distance between structures by considering only the relative positions and orientations of nucleobases. The eRMSD can be considered as a vectorial version of contact maps and it is calculated as follows:

1. Set up a local reference system in the center of the six-membered ring of each nucleobase in a molecule.
   The xy plane lies on the plane of the nucleobase, and it is oriented such that the Watson-Crick interaction is always at $\theta\approx 60^{\circ}$.

2. Calculate all pairwise distance vectors $\vec{r}_{i,j}$ among base centers.

3. Rescale distance vectors as $\tilde{\vec{r}}_{i,j}=(r_x/a,r_y/a,r_z/b)$, where a=5 Angstrom and b=3 Angstrom. This rescaling has the effect of weighting more deviations on the z-axis with respect to the x/y directions.

4. Calculate the G vectors

$$
\vec{G}(\tilde{\vec{r}}) = (\sin(\gamma \tilde{r}) \tilde{r}_x/\tilde{r},\sin(\gamma \tilde{r}) \tilde{r}_y/\tilde{r},\sin(\gamma \tilde{r}) \tilde{r}_z/\tilde{r}, 1+\cos(\gamma \tilde{r})) \times
\frac{\Theta(\tilde{r}_{cutoff}-\tilde{r})}{\gamma}
$$

Here, $ \gamma = \pi/\tilde{r}_{cutoff}$ and $ \Theta $ is the Heaviside step function. The default cutoff is set to 2.4.

5. The eRMSD between two structures $ \alpha $ and $ \beta $ reads

$$
eRMSD = \sqrt{\frac{1}{N} \sum_{j,k} \vert \vec{G}(\tilde{\vec{r}}_{jk}^{\alpha}) - \vec{G}(\tilde{\vec{r}}_{jk}^{\beta}) \vert^2 }
$$

Using the default cutoff, two structures with eRMSD of 0.7 or lower can be considered as significantly similar. A full description of the eRMSD can be found in the first paper in the bibliography below.

ERMSD is computed using the position of three atoms on the 6-membered ring of each involved nucleobase. The atoms should be:
- C2,C4,C6 for pyrimdines
- C2,C6,C4 for purines

The different order for purines and pyrimidines is fundamental and allows you to compute ERMSD between structures with different
sequences as well! Notice that the simplest way to avoid mistakes in choosing these atoms is to use the `@lcs-#` strings
as shown in the examples (see also [MOLINFO](MOLINFO.md)).

!!! caution "ERMSD paths not possible"

    Notice that the ERMSD implemented here is not integrated with the other metrics in plumed. As a consequence, it is not (yet) possible
    to e.g. build path collective variables using ERMSD

!!! caution "ERMSD with multiple molecules"

    Notice that ERMSD expect a single molecule and makes coordinates whole before anything else. As such, results might be unexpected
    for a multi molecular system.

## Pair search

By default, a buffered neighbor list is built using PLUMED linked cells.
The search sphere encloses the anisotropic cutoff ellipsoid: its radius is
`CUTOFF/2`, corresponding to the longest semiaxis. The default buffer is one sixth
of this radius; `NL_BUFFER`, expressed in the selected length unit (nm by default),
can override it. The list is rebuilt whenever any
base center moves more than half the buffer from its position at the last build.
Centers are checked on every evaluation, so no update stride is required.
Base rotations do not invalidate the enclosing sphere. Reference contacts are
retained even when the corresponding instantaneous contact disappears.

Use `PAIR_SEARCH=CELLS` to rebuild an unbuffered linked-cell list each evaluation,
or `PAIR_SEARCH=DIRECT` to visit every pair. These give the same metric and
analytical derivatives. Pair evaluation uses OpenMP and the action MPI communicator.

## Examples

Calculate the eRMSD from reference structure reference.pdb using the default cutoff (2.4). The list of residues involved in the calculation has to be specified. In this example, the eRMSD is calculated
considering residues 1,2,3,4,5,6.

```plumed
#SETTINGS MOLFILE=regtest/basic/rt-ermsd/ref.pdb
MOLINFO STRUCTURE=regtest/basic/rt-ermsd/ref.pdb
eRMSD1: ERMSD REFERENCE=regtest/basic/rt-ermsd/ref.pdb ATOMS=@lcs-1,@lcs-2,@lcs-3,@lcs-4,@lcs-5,@lcs-6
```

*/
//+ENDPLUMEDOC


class ERMSD : public Colvar {


  std::vector<Vector> derivs;
  PLMD::ERMSD ermsd;
  bool pbc;

public:
  explicit ERMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ERMSD,"ERMSD")

void ERMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);

  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","CUTOFF","2.4","dimensionless ellipsoidal cutoff for pairs of base centers.");
  keys.add("compulsory","PAIR_SEARCH","VERLET","pair search: DIRECT (all pairs), CELLS (linked cells rebuilt each step), or VERLET (buffered linked-cell list with automatic displacement updates).");
  keys.add("optional","NL_BUFFER","buffer in the selected length unit added to the enclosing sphere radius for VERLET; defaults to one sixth of that radius. No update stride is needed.");
  keys.add("atoms","ATOMS","the list of atoms (use lcs)");
  keys.add("optional","PAIRS","List of pairs considered. All pairs are considered if this value is not specified.");
  keys.setValueDescription("scalar","the eRMSD between the instantaneous structure and the reference structure that was input");
  keys.addDOI("10.1093/nar/gku972");
}

ERMSD::ERMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao), pbc(true) {
  std::string reference;
  parse("REFERENCE",reference);
  double cutoff=2.4;
  parse("CUTOFF",cutoff);


  bool nopbc(false);
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  std::vector<AtomNumber> atoms_;
  parseAtomList("ATOMS",atoms_);

  std::vector<unsigned> pairs_;
  parseVector("PAIRS",pairs_);
  std::string pairSearch="VERLET";
  parse("PAIR_SEARCH",pairSearch);
  double buffer=-1.;
  std::string bufferInput;
  parse("NL_BUFFER",bufferInput);
  const bool haveBuffer=!bufferInput.empty();
  if(haveBuffer && !Tools::convertNoexcept(bufferInput,buffer)) {
    error("NL_BUFFER must be a number");
  }
  if(cutoff<=0. || !std::isfinite(cutoff)) {
    error("CUTOFF must be finite and positive");
  }
  if(haveBuffer && (buffer<=0. || !std::isfinite(buffer))) {
    error("NL_BUFFER must be finite and positive");
  }
  PLMD::ERMSD::PairSearch method;
  if(pairSearch=="DIRECT") {
    method=PLMD::ERMSD::PairSearch::Direct;
  } else if(pairSearch=="CELLS") {
    method=PLMD::ERMSD::PairSearch::Cells;
  } else if(pairSearch=="VERLET") {
    method=PLMD::ERMSD::PairSearch::Verlet;
  } else {
    error("PAIR_SEARCH must be DIRECT, CELLS, or VERLET");
  }
  if(haveBuffer && method!=PLMD::ERMSD::PairSearch::Verlet) {
    error("NL_BUFFER requires PAIR_SEARCH=VERLET");
  }
  ermsd.setPairSearch(method,haveBuffer ? buffer : -1.);
  checkRead();

  addValueWithDerivatives();
  setNotPeriodic();

  if(atoms_.size()<6) {
    error("at least six atoms should be specified");
  }
  if(atoms_.size()%3!=0) {
    error("Atoms are not multiple of 3");
  }
  if(pairs_.size()%2!=0) {
    error("pairs are not multiple of 2");
  }



  // read everything in ang and transform to nm if we are not in natural units
  PDB pdb;
  if( !pdb.read(reference,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + reference );
  }
  // store target_ distance
  std::vector <Vector> reference_positions;
  unsigned natoms = atoms_.size();
  log.printf("Read %u atoms\n",natoms);

  reference_positions.resize(natoms);
  for(unsigned i=0; i<natoms; i++) {
    reference_positions[i] = pdb.getPosition(atoms_[i]);
  }

// shift to count from zero (pairs are residue indices, 1-based in input)
  for(unsigned i=0; i<pairs_.size(); ++i) {
    if(pairs_[i]<1 || pairs_[i]>natoms/3) {
      error("PAIRS should contain residue indices between 1 and the number of residues (" + std::to_string(natoms/3) + ")");
    }
    pairs_[i]--;
  }

  ermsd.setCommunicator(comm);
  ermsd.setReference(reference_positions,pairs_,cutoff/getUnits().getLength());

  requestAtoms(atoms_);
  derivs.resize(natoms);

  log.printf("  pair search %s; enclosing sphere radius %g\n",pairSearch.c_str(),cutoff/(2.*getUnits().getLength()));
  if(method==PLMD::ERMSD::PairSearch::Verlet) {
    log.printf("  neighbor buffer %g; rebuild when a center moves more than half the buffer\n",
               haveBuffer ? buffer : cutoff/(12.*getUnits().getLength()));
  }
  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %u atoms\n",natoms);

  log<<"  Bibliography "
     <<plumed.cite("Bottaro, Di Palma, and Bussi, Nucleic Acids Res. 42, 13306 (2014)")
     <<plumed.cite("Bottaro, Banas, Sponer, and Bussi, J. Phys. Chem. Lett. 7, 4032 (2016)")<<"\n";

}

// calculator
void ERMSD::calculate() {
// set derivatives to zero
  Tools::set_to_zero(derivs);
  double ermsdist;
  Tensor virial;
// This is a trick to avoid explicit virial calculation
// 1. we make the molecule whole (unless NOPBC was requested)
  if(pbc) {
    makeWhole();
  }
// 2. we ignore pbcs
  Pbc fake_pbc;
// Notice that this might have problems when having 2 RNA molecules (hybridization).

  ermsdist=ermsd.calculate(getPositions(),fake_pbc,derivs,virial);
  const double scale=getUnits().getLength();
  setValue(ermsdist*scale);

  for(unsigned i=0; i<derivs.size(); ++i) {
    setAtomsDerivatives(i,derivs[i]*scale);
  }

  setBoxDerivativesNoPbc();
}

}
}
