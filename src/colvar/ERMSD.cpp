/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
 This vast majority of the source code in this file was writting by
 Sandro Bottaro with some help from Giovanni Bussi
*/

#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "tools/ERMSD.h"
#include "core/Atoms.h"
#include <iostream>

using namespace std;

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
   The xy plane lies on the plane of the nucleobase, and it is oriented such that the Watson-Crick interaction is always at \f$\theta\approx 60^{\circ}\f$.

2. Calculate all pairwise distance vectors \f$\vec{r}_{i,j}\f$ among base centers.

3. Rescale distance vectors as \f$\tilde{\vec{r}}_{i,j}=(r_x/a,r_y/a,r_z/b)\f$, where  a=b=5 \f$\r{A}\f$, c=3 \f$\r{A}\f$. This rescaling has the effect of weighting more deviations on the z-axis with respect to the x/y directions.

4. Calculate the G vectors

\f[
\vec{G}(\tilde{\vec{r}}) = (\sin(\gamma \tilde{r}) \tilde{r}_x/\tilde{r},\sin(\gamma \tilde{r}) \tilde{r}_y/\tilde{r},\sin(\gamma \tilde{r}) \tilde{r}_z/\tilde{r}, 1+\cos(\gamma \tilde{r})) \times
\frac{\Theta(\tilde{r}_{cutoff}-\tilde{r})}{\gamma}
\f]

Here, \f$ \gamma = \pi/\tilde{r}_{cutoff}\f$ and \f$ \Theta \f$ is the Heaviside step function. The default cutoff is set to 2.4.

5. The eRMSD between two structures \f$ \alpha \f$ and \f$ \beta \f$ reads

\f[
eRMSD = \sqrt{\frac{1}{N} \sum_{j,k} \vert \vec{G}(\tilde{\vec{r}}_{jk}^{\alpha}) - \vec{G}(\tilde{\vec{r}}_{jk}^{\beta}) \vert^2 }
\f]

Using the default cutoff, two structures with eRMSD of 0.7 or lower can be considered as significantly similar. A full description of the eRMSD can be found in \cite bott14

ERMSD is computed using the position of three atoms on the 6-membered ring of each involved nucleobase. The atoms should be:
- C2,C4,C6 for pyrimdines
- C2,C6,C4 for purines

The different order for purines and pyrimidines is fundamental and allows you to compute ERMSD between structures with different
sequences as well! Notice that the simplest way to avoid mistakes in choosing these atoms is to use the `@lcs-#` strings
as shown in the examples (see also \ref MOLINFO).

\warning Notice that the ERMSD implemented here is not integrated with the other metrics in plumed. As a consequence, it is not (yet) possible
to e.g. build path collective variables using ERMSD

\warning Notice that ERMSD expect a single molecule and makes coordinate whole before anything else. As such, results might be unexpected
for a multi molecular system.

\par Examples

Calculate the eRMSD from reference structure reference.pdb using the default cutoff (2.4). The list of residues involved in the calculation has to be specified. In this example, the eRMSD is calculated
considering residues 1,2,3,4,5,6.

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt-ermsd/ref.pdb
MOLINFO STRUCTURE=reference.pdb
eRMSD1: ERMSD REFERENCE=reference.pdb ATOMS=@lcs-1,@lcs-2,@lcs-3,@lcs-4,@lcs-5,@lcs-6
\endplumedfile

*/
//+ENDPLUMEDOC


class ERMSD : public Colvar {


  vector<Vector> derivs;
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
  keys.add("compulsory","CUTOFF","2.4","only pairs of atoms closer than CUTOFF are considered in the calculation.");
  keys.add("atoms","ATOMS","the list of atoms (use lcs)");
  keys.add("optional","PAIRS","List of pairs considered. All pairs are considered if this value is not specified.");

}

ERMSD::ERMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao), pbc(true)
{
  string reference;
  parse("REFERENCE",reference);
  double cutoff=2.4;
  parse("CUTOFF",cutoff);


  bool nopbc(false);
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  vector<AtomNumber> atoms_;
  parseAtomList("ATOMS",atoms_);

  vector<unsigned> pairs_;
  parseVector("PAIRS",pairs_);
  checkRead();

  addValueWithDerivatives(); setNotPeriodic();

  if(atoms_.size()<6) error("at least six atoms should be specified");
  if(atoms_.size()%3!=0) error("Atoms are not multiple of 3");
  if(pairs_.size()%2!=0) error("pairs are not multiple of 2");


  //checkRead();
  //log.printf("  of atoms");
  //for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  //requestAtoms(atoms);

  // read everything in ang and transform to nm if we are not in natural units
  PDB pdb;
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + reference );
  // store target_ distance
  vector <Vector> reference_positions;
  unsigned natoms = atoms_.size();
  log.printf("Read %u atoms\n",natoms);

  reference_positions.resize(natoms);
  for(unsigned i=0; i<natoms; i++) {
    reference_positions[i] = pdb.getPosition(atoms_[i]);
    //log.printf("%f %f %f \n",reference_positions[i][0],reference_positions[i][1],reference_positions[i][2]);
  }

// shift to count from zero
  for(unsigned i=0; i<pairs_.size(); ++i) pairs_[i]--;

  ermsd.setReference(reference_positions,pairs_,cutoff/atoms.getUnits().getLength());

  requestAtoms(atoms_);
  derivs.resize(natoms);

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %u atoms\n",natoms);

  log<<"  Bibliography "
     <<plumed.cite("Bottaro, Di Palma, and Bussi, Nucleic Acids Res. 42, 13306 (2014)")
     <<plumed.cite("Bottaro, Banas, Sponer, and Bussi, J. Phys. Chem. Lett. 7, 4032 (2016)")<<"\n";

}

// calculator
void ERMSD::calculate() {
// set derivatives to zero
  for(unsigned i=0; i<derivs.size(); ++i) {derivs[i].zero();}
  double ermsdist;
  Tensor virial;
// This is a trick to avoid explicit virial calculation
// 1. we make the molecule whole
  makeWhole();
// 2. we ignore pbcs
  Pbc fake_pbc;
// Notice that this might have problems when having 2 RNA molecules (hybridization).

  ermsdist=ermsd.calculate(getPositions(),fake_pbc,derivs,virial);
  const double scale=atoms.getUnits().getLength();
  setValue(ermsdist*scale);

  for(unsigned i=0; i<derivs.size(); ++i) {setAtomsDerivatives(i,derivs[i]*scale);}

  setBoxDerivativesNoPbc();

//setBoxDerivatives(virial);

}

}
}
