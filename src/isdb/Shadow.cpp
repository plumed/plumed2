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
#include "tools/Pbc.h"
#include "tools/RMSD.h"
#include <string>

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SHADOW
/*


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
}

Shadow::Shadow(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  isreference_(false), nupdate_(1), pbc_(true), first_time_(true)
{
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
  if(nrep<2) error("SHADOW must be used with at least two replicas");

  checkRead();

  log.printf("  atoms involved : ");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf("%d ",atoms[i].serial());
  log.printf("\n");
  log.printf("  stride for updating reference coordinates : %d\n", nupdate_);
  log.printf("  number of replicas : %d\n", nrep);
  log.printf("  replica id : %d\n", replica);
  if(isreference_) log.printf("  this is the reference replica\n");

  // add value and set periodicity
  addValueWithDerivatives(); setNotPeriodic();

  // request atoms
  requestAtoms(atoms);
}

void Shadow::update_reference()
{
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
    if(isreference_) reference = getPositions();
    // share coordinates
    multi_sim_comm.Sum(&reference[0][0], 3*natoms);
  }
// now intra replica (openmp) communication
  if(size_>1) comm.Sum(&reference[0][0], 3*natoms);

// clear the rmsd object
  rmsd_.clear();
// and initialize it
  rmsd_.set(align,displace,reference,"OPTIMAL");
}

void Shadow::calculate()
{
  // make whole
  if(pbc_) makeWhole();

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
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) setAtomsDerivatives(i, derivatives[i]);
  }
  // set virial
  setBoxDerivativesNoPbc();
}

}
}
