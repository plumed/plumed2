/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2022 The plumed team
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
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"
#include "tools/RMSD.h"
#include <string>

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SHADOW
/*
Communicate atoms positions among replicas and calculate the RMSD with respect to a mother (reference) simulation.

The option UPDATE allows to specify the stride for communication between mother and replica systems.
The flag REFERENCE needs to be specified in the input file of the mother replica.
This action must be run in a multi-replica framework (such as the -multi option in GROMACS).

## Examples

In this example, we perform a simulation of a RNA molecule using two replicas: a mother and a shadow replica.
The mother simulation communicates the coordinates of the RNA backbone to the replica every 100 steps.
The RMSD of the replica with respect to the mother is calculated on the RNA backbone atoms and an [UPPER_WALLS](UPPER_WALLS.md) is applied at 0.2 nm.
The mother replica contains also the [UPPER_WALLS](UPPER_WALLS.md) action. However, the forces on the RNA atoms of the mother replica are automatically set to zero
inside the [SHADOW](SHADOW.md) action.

The input file for the mother simulation looks as follows:

```plumed
# Reference PDB
MOLINFO STRUCTURE=conf_emin_PLUMED.pdb WHOLE
# Define RNA nucleic backbone
rna: GROUP ATOMS=1,2,5,6,33,36,37,40,41,67,70,71,74,75,98,101,102,105,106,131,134,135,138,139,165,168,169,172,173,198,201,202,205,206,228,231,232,235,236,259,262,263,266,267,289,292,293,296,297,323,326,327,330,331,356,359,360,363,364,390,393,394,397,398,421,424,425,428,429,452,455,456,459,460,482,485,486,489,490,516,519,520,523,524,550,553,554,557,558,584,587,588,591,592,617,620,621,624,625,651,654,655,658,659,682,685,686,689,690,712,715,716,719,720,743,746,747,750,751,773,776,777,780,781,804,807,808,811,812,834,837,838,841,842,868,871,872,875,876,899,902,903,906,907
# Reconstruct RNA PBC
WHOLEMOLECULES ENTITY0=rna EMST STRIDE=1

# Define shadow RMSD on RNA backbone
rmsd: SHADOW ATOMS=rna NOPBC UPDATE=100 REFERENCE
# Add upper wall - derivatives are set to zero inside SHADOW action
uws: UPPER_WALLS ARG=rmsd AT=0.2 KAPPA=10000.0 STRIDE=1

# Print useful info
PRINT FILE=COLVAR STRIDE=500 ARG=rmsd,uws.bias
```

while the input file for a shadow replica looks like:

```plumed
# Reference PDB
MOLINFO STRUCTURE=conf_emin_PLUMED.pdb WHOLE
# Define RNA nucleic backbone
rna: GROUP ATOMS=1,2,5,6,33,36,37,40,41,67,70,71,74,75,98,101,102,105,106,131,134,135,138,139,165,168,169,172,173,198,201,202,205,206,228,231,232,235,236,259,262,263,266,267,289,292,293,296,297,323,326,327,330,331,356,359,360,363,364,390,393,394,397,398,421,424,425,428,429,452,455,456,459,460,482,485,486,489,490,516,519,520,523,524,550,553,554,557,558,584,587,588,591,592,617,620,621,624,625,651,654,655,658,659,682,685,686,689,690,712,715,716,719,720,743,746,747,750,751,773,776,777,780,781,804,807,808,811,812,834,837,838,841,842,868,871,872,875,876,899,902,903,906,907
# Reconstruct RNA PBC
WHOLEMOLECULES ENTITY0=rna EMST STRIDE=1

# Define shadow RMSD on RNA backbone
rmsd: SHADOW ATOMS=rna NOPBC UPDATE=100
# Add upper wall
uws: UPPER_WALLS ARG=rmsd AT=0.2 KAPPA=10000.0 STRIDE=1

# Print useful info
PRINT FILE=COLVAR STRIDE=500 ARG=rmsd,uws.bias
```

*/
//+ENDPLUMEDOC

class Shadow : public Colvar {
// private stuff
  bool isreference_;
  unsigned nupdate_;
  bool pbc_;
  bool first_time_;
// RMSD object
  PLMD::RMSD rmsd_;
// parallel stuff
  unsigned size_;
  unsigned rank_;
// update reference
  void update_reference();

public:
  static void registerKeywords( Keywords& keys );
  explicit Shadow(const ActionOptions&);
// active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(Shadow,"SHADOW")

void Shadow::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","atoms for which we calculate the shadow RMSD");
  keys.add("compulsory","UPDATE","stride for updating reference coordinates");
  keys.addFlag("REFERENCE",false,"this is the reference replica");
  keys.setValueDescription("scalar","the value of the shadow RMSD");
}

Shadow::Shadow(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  isreference_(false), nupdate_(1), pbc_(true), first_time_(true) {
  // list of atoms
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  // update stride
  parse("UPDATE",nupdate_);

  // is this the reference replica
  parseFlag("REFERENCE",isreference_);

  // periodic boundary conditions
  bool nopbc=!pbc_;
  parseFlag("NOPBC",nopbc);
  pbc_=!nopbc;

  // set intra-replica (openmp) parallel stuff
  size_ = comm.Get_size();
  rank_ = comm.Get_rank();

  // get number of (MPI) replicas
  int nrep = 0;
  int replica = 0;
  // only if openmp master
  if(rank_==0) {
    nrep    = multi_sim_comm.Get_size();
    replica = multi_sim_comm.Get_rank();
  }
  comm.Sum(&nrep,1);
  comm.Sum(&replica,1);
  // check number of replicas
  //if(nrep<2) error("SHADOW must be used with at least two replicas");

  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf("%d ",atoms[i].serial());
  }
  log.printf("\n");
  log.printf("  stride for updating reference coordinates : %d\n", nupdate_);
  log.printf("  number of replicas : %d\n", nrep);
  log.printf("  replica id : %d\n", replica);
  if(isreference_) {
    log.printf("  this is the reference replica\n");
  }

  // add value and set periodicity
  addValueWithDerivatives();
  setNotPeriodic();

  // request atoms
  requestAtoms(atoms);
}

void Shadow::update_reference() {
// number of atoms
  unsigned natoms = getNumberOfAtoms();
// initialize rmsd variables
  std::vector<double> align(natoms,1.0);
  std::vector<double> displace(natoms,1.0);
  std::vector<Vector> reference(natoms);

// first get the reference coordinates
// if master openmp task
  if(rank_==0) {
    // if reference replica
    if(isreference_) {
      reference = getPositions();
    }
    // share coordinates
    multi_sim_comm.Sum(&reference[0][0], 3*natoms);
  }
// now intra replica (openmp) communication
  if(size_>1) {
    comm.Sum(&reference[0][0], 3*natoms);
  }

// clear the rmsd object
  rmsd_.clear();
// and initialize it
  rmsd_.set(align,displace,reference,"OPTIMAL");
}

void Shadow::calculate() {
  // make whole
  if(pbc_) {
    makeWhole();
  }

  // if it is time, update reference coordinates
  if(first_time_ || getStep()%nupdate_==0) {
    update_reference();
    first_time_ = false;
  }

  // calculate RMSD and derivatives
  std::vector<Vector> derivatives(getNumberOfAtoms());
  double rmsd = rmsd_.calculate(getPositions(), derivatives);

  // set RMSD value
  setValue(rmsd);
  // if this is not the reference replica, add derivatives
  if(!isreference_) {
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      setAtomsDerivatives(i, derivatives[i]);
    }
  }
  // set virial
  setBoxDerivativesNoPbc();
}

}
}
