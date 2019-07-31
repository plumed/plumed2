/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "reference/DRMSD.h"
#include "reference/MetricRegister.h"
#include "core/Atoms.h"
#include <memory>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC DCOLVAR DRMSD
/*
Calculate the distance RMSD with respect to a reference structure.

To calculate the root-mean-square deviation between the atoms in two configurations
you must first superimpose the two structures in some ways.  Obviously, it is the internal vibrational
motions of the structure - i.e. not the translations and rotations - that are interesting. However,
aligning two structures by removing the translational and rotational motions is not easy.  Furthermore,
in some cases there can be alignment issues caused by so-called frame-fitting problems. It is thus
often cheaper and easier to calculate the distances between all the pairs of atoms.  The distance
between the two structures, \f$\mathbf{X}^a\f$ and \f$\mathbf{X}^b\f$ can then be measured as:

\f[
d(\mathbf{X}^A, \mathbf{X}^B) = \sqrt{\frac{1}{N(N-1)} \sum_{i \ne j} [ d(\mathbf{x}_i^a,\mathbf{x}_j^a) - d(\mathbf{x}_i^b,\mathbf{x}_j^b) ]^2}
\f]

where \f$N\f$ is the number of atoms and \f$d(\mathbf{x}_i,\mathbf{x}_j)\f$ represents the distance between
atoms \f$i\f$ and \f$j\f$.  Clearly, this representation of the configuration is invariant to translation and rotation.
However, it can become expensive to calculate when the number of atoms is large.  This can be resolved
within the DRMSD colvar by setting LOWER_CUTOFF and UPPER_CUTOFF.  These keywords ensure that only
pairs of atoms that are within a certain range are incorporated into the above sum.

In PDB files the atomic coordinates and box lengths should be in Angstroms unless
you are working with natural units.  If you are working with natural units then the coordinates
should be in your natural length unit.  For more details on the PDB file format visit http://www.wwpdb.org/docs.html

\par Examples

The following tells plumed to calculate the distance RMSD between
the positions of the atoms in the reference file and their instantaneous
position. Only pairs of atoms whose distance in the reference structure is within
0.1 and 0.8 nm are considered.

\plumedfile
DRMSD REFERENCE=file1.pdb LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8
\endplumedfile

The reference file is a PDB file that looks like this

\auxfile{file1.pdb}
ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
END
\endauxfile

The following tells plumed to calculate a DRMSD value for a pair of molecules.

\plumedfile
DRMSD REFERENCE=file2.pdb LOWER_CUTOFF=0.1 UPPER_CUTOFF=0.8 TYPE=INTER-DRMSD
\endplumedfile

In the input reference file (file.pdb) the atoms in each of the two molecules are separated by a TER
command as shown below.

\auxfile{file2.pdb}
ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
TER
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
END
\endauxfile

In this example the INTER-DRMSD type ensures that the set of distances from which the final
quantity is computed involve one atom from each of the two molecules.  If this is replaced
by INTRA-DRMSD then only those distances involving pairs of atoms that are both in the same
molecule are computed.

*/
//+ENDPLUMEDOC


class DRMSD : public Colvar {

  bool pbc_;
  MultiValue myvals;
  ReferenceValuePack mypack;
  std::unique_ptr<PLMD::DRMSD> drmsd_;

public:
  explicit DRMSD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(DRMSD,"DRMSD")

void DRMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","LOWER_CUTOFF","only pairs of atoms further than LOWER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","UPPER_CUTOFF","only pairs of atoms closer than UPPER_CUTOFF are considered in the calculation.");
  keys.add("compulsory","TYPE","DRMSD","what kind of DRMSD would you like to calculate.  You can use either the normal DRMSD involving all the distances between "
           "the atoms in your molecule.  Alternatively, if you have multiple molecules you can use the type INTER-DRMSD "
           "to compute DRMSD values involving only those distances between the atoms at least two molecules or the type INTRA-DRMSD "
           "to compute DRMSD values involving only those distances between atoms in the same molecule");
}

DRMSD::DRMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao), pbc_(true), myvals(1,0), mypack(0,0,myvals)
{
  string reference;
  parse("REFERENCE",reference);
  double lcutoff;
  parse("LOWER_CUTOFF",lcutoff);
  double ucutoff;
  parse("UPPER_CUTOFF",ucutoff);
  bool nopbc(false);
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  addValueWithDerivatives(); setNotPeriodic();

  // read everything in ang and transform to nm if we are not in natural units
  PDB pdb;
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + reference );

  // store target_ distance
  std::string type; parse("TYPE",type);
  drmsd_=metricRegister().create<PLMD::DRMSD>( type );
  drmsd_->setBoundsOnDistances( !nopbc, lcutoff, ucutoff );
  drmsd_->read( pdb );
  checkRead();

  std::vector<AtomNumber> atoms;
  drmsd_->getAtomRequests( atoms );
//   drmsd_->setNumberOfAtoms( atoms.size() );
  requestAtoms( atoms );

  // Setup the derivative pack
  myvals.resize( 1, 3*atoms.size()+9 ); mypack.resize( 0, atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) mypack.setAtomIndex( i, i );

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  with indices : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");
}

// calculator
void DRMSD::calculate() {

  double drmsd; Tensor virial; mypack.clear();
  drmsd=drmsd_->calculate(getPositions(), getPbc(), mypack, false);

  setValue(drmsd);
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) { if( myvals.isActive(3*i) ) setAtomsDerivatives( i, mypack.getAtomDerivative(i) ); }
  setBoxDerivatives( mypack.getBoxDerivatives() );
}

}
}
