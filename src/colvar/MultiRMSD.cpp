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
#include "Colvar.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "reference/MultiDomainRMSD.h"
#include "reference/MetricRegister.h"
#include "core/Atoms.h"
#include <memory>


using namespace std;

namespace PLMD {
namespace colvar {

class MultiRMSD : public Colvar {

  std::unique_ptr<PLMD::MultiDomainRMSD> rmsd;
  bool squared;
  MultiValue myvals;
  ReferenceValuePack mypack;
  bool nopbc;

public:
  explicit MultiRMSD(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};


using namespace std;

//+PLUMEDOC DCOLVAR MULTI_RMSD
/*
Calculate the RMSD distance moved by a number of separated domains from their positions in a reference structure.


When you have large proteins the calculation of the root mean squared deviation between all the atoms in a reference
structure and the instantaneous configuration becomes prohibitively expensive.  You may thus instead want to calculate
the RMSD between the atoms in a set of domains of your protein and your reference structure.  That is to say:

\f[
d(X,X_r) = \sqrt{ \sum_{i} w_i\vert X_i - X_i' \vert^2 }
\f]

where here the sum is over the domains of the protein, \f$X_i\f$ represents the positions of the atoms in domain \f$i\f$
in the instantaneous configuration and \f$X_i'\f$ is the positions of the atoms in domain \f$i\f$ in the reference
configuration.  \f$w_i\f$ is an optional weight.

The distances for each of the domains in the above sum can be calculated using the \ref DRMSD or \ref RMSD measures or
using a combination of these distance.  The reference configuration is specified in a pdb file like the one below:

\verbatim
ATOM      2  O   ALA     2      -0.926  -2.447  -0.497  1.00  1.00      DIA  O
ATOM      4  HNT ALA     2       0.533  -0.396   1.184  1.00  1.00      DIA  H
ATOM      6  HT1 ALA     2      -0.216  -2.590   1.371  1.00  1.00      DIA  H
ATOM      7  HT2 ALA     2      -0.309  -1.255   2.315  1.00  1.00      DIA  H
ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
TER
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
ATOM     16  HN  ALA     2       1.713   1.021  -0.873  1.00  1.00      DIA  H
ATOM     18  HA  ALA     2       0.099  -0.774  -2.218  1.00  1.00      DIA  H
ATOM     19  CB  ALA     2       2.063  -1.223  -1.276  1.00  1.00      DIA  C
ATOM     20  HB1 ALA     2       2.670  -0.716  -2.057  1.00  1.00      DIA  H
ATOM     21  HB2 ALA     2       2.556  -1.051  -0.295  1.00  1.00      DIA  H
ATOM     22  HB3 ALA     2       2.070  -2.314  -1.490  1.00  1.00      DIA  H
END
\endverbatim

with the TER keyword being used to separate the various domains in you protein.


\par Examples

The following tells plumed to calculate the RMSD distance between
the positions of the atoms in the reference file and their instantaneous
position.  The Kearsley algorithm for each of the domains.

\plumedfile
MULTI_RMSD REFERENCE=file.pdb TYPE=MULTI-OPTIMAL
\endplumedfile

The following tells plumed to calculate the RMSD distance between the positions of
the atoms in the domains of reference the reference structure and their instantaneous
positions.  Here distances are calculated using the \ref DRMSD measure.

\plumedfile
MULTI_RMSD REFERENCE=file.pdb TYPE=MULTI-DRMSD
\endplumedfile

in this case it is possible to use the following DRMSD options in the pdb file using the REMARK syntax:
\verbatim
NOPBC to calculate distances without PBC
LOWER_CUTOFF=# only pairs of atoms further than LOWER_CUTOFF are considered in the calculation
UPPER_CUTOFF=# only pairs of atoms further than UPPER_CUTOFF are considered in the calculation
\endverbatim
as shown in the following example

\verbatim
REMARK NOPBC
REMARK LOWER_CUTOFF=0.1
REMARK UPPER_CUTOFF=0.8
ATOM      2  O   ALA     2      -0.926  -2.447  -0.497  1.00  1.00      DIA  O
ATOM      4  HNT ALA     2       0.533  -0.396   1.184  1.00  1.00      DIA  H
ATOM      6  HT1 ALA     2      -0.216  -2.590   1.371  1.00  1.00      DIA  H
ATOM      7  HT2 ALA     2      -0.309  -1.255   2.315  1.00  1.00      DIA  H
ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
TER
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
ATOM     16  HN  ALA     2       1.713   1.021  -0.873  1.00  1.00      DIA  H
ATOM     18  HA  ALA     2       0.099  -0.774  -2.218  1.00  1.00      DIA  H
ATOM     19  CB  ALA     2       2.063  -1.223  -1.276  1.00  1.00      DIA  C
ATOM     20  HB1 ALA     2       2.670  -0.716  -2.057  1.00  1.00      DIA  H
ATOM     21  HB2 ALA     2       2.556  -1.051  -0.295  1.00  1.00      DIA  H
ATOM     22  HB3 ALA     2       2.070  -2.314  -1.490  1.00  1.00      DIA  H
END
\endverbatim


*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(MultiRMSD,"MULTI_RMSD")

void MultiRMSD::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","TYPE","MULTI-SIMPLE","the manner in which RMSD alignment is performed.  Should be MULTI-OPTIMAL, MULTI-OPTIMAL-FAST,  MULTI-SIMPLE or MULTI-DRMSD.");
  keys.addFlag("SQUARED",false," This should be set if you want the mean squared displacement instead of the root mean squared displacement");
}

MultiRMSD::MultiRMSD(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),squared(false),myvals(1,0), mypack(0,0,myvals),nopbc(false)
{
  string reference;
  parse("REFERENCE",reference);
  string type;
  type.assign("SIMPLE");
  parse("TYPE",type);
  parseFlag("SQUARED",squared);
  parseFlag("NOPBC",nopbc);
  checkRead();

  addValueWithDerivatives(); setNotPeriodic();
  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
    error("missing input file " + reference );

  rmsd=metricRegister().create<MultiDomainRMSD>(type,pdb);
  // Do not align molecule if we are doing DRMSD for domains and NOPBC has been specified in input
  if( pdb.hasFlag("NOPBC") ) nopbc=true;

  std::vector<AtomNumber> atoms;
  rmsd->getAtomRequests( atoms );
  requestAtoms( atoms );

  myvals.resize( 1, 3*atoms.size()+9 ); mypack.resize( 0, atoms.size() );
  for(unsigned i=0; i<atoms.size(); ++i) mypack.setAtomIndex( i, i );

  log.printf("  reference from file %s\n",reference.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",type.c_str() );
  if(squared)log.printf("  chosen to use SQUARED option for MSD instead of RMSD\n");
}

// calculator
void MultiRMSD::calculate() {
  if(!nopbc) makeWhole();
  double r=rmsd->calculate( getPositions(), getPbc(), mypack, squared );

  setValue(r);
  for(unsigned i=0; i<getNumberOfAtoms(); i++) setAtomsDerivatives( i, mypack.getAtomDerivative(i) );

  if( !mypack.virialWasSet() ) setBoxDerivativesNoPbc();
  else setBoxDerivatives( mypack.getBoxDerivatives() );
}

}
}



