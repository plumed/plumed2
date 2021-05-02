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
#include "Atoms.h"
#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include "ActionToPutData.h"
#include "ActionSet.h"
#include "tools/Pbc.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

namespace PLMD {

/// We assume that charges and masses are constant along the simulation
/// Set this to false if you want to revert to the original (expensive) behavior
// static const bool shareMassAndChargeOnlyAtFirstStep=true;

class PlumedMain;

Atoms::Atoms(PlumedMain&plumed):
  natoms(0),
  // dataCanBeSet(false),
  // positionsHaveBeenSet(0),
  // forcesHaveBeenSet(0),
  shuffledAtoms(0),
  // mdatoms(MDAtomsBase::create(sizeof(double))),
  plumed(plumed),
  naturalUnits(false),
  MDnaturalUnits(false),
  timestep(0.0),
  forceOnEnergy(0.0),
  zeroallforces(false),
  kbT(0.0),
  asyncSent(false),
  atomsNeeded(false),
  ddStep(0)
{
}

Atoms::~Atoms() {
  if(actions.size()>0) {
    std::cerr<<"WARNING: there is some inconsistency in action added to atoms, as some of them were not properly destroyed. This might indicate an internal bug!!\n";
  }
}

void Atoms::startStep() {
//   positionsHaveBeenSet=0; 
//   forcesHaveBeenSet=0; 
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions()); 
  for(const auto & ip : inputs) { ip.second->dataCanBeSet=true; }
//  dataCanBeSet=true;
}

// void Atoms::setPositions(void*p) {
//   plumed_massert( dataCanBeSet,"setPositions must be called after setStep in MD code interface");
//   plumed_massert( p || gatindex.size()==0, "NULL position pointer with non-zero local atoms");
//   mdatoms->setp(p); positionsHaveBeenSet=3;
// }
// 
// void Atoms::setForces(void*p) {
//   plumed_massert( dataCanBeSet,"setForces must be called after setStep in MD code interface");
//   plumed_massert( p || gatindex.size()==0, "NULL force pointer with non-zero local atoms");
//   forcesHaveBeenSet=3;
//   mdatoms->setf(p);
// }
// 
// void Atoms::setPositions(void*p,int i) {
//   plumed_massert( dataCanBeSet,"setPositions must be called after setStep in MD code interface");
//   plumed_massert( p || gatindex.size()==0, "NULL positions pointer with non-zero local atoms");
//   mdatoms->setp(p,i); positionsHaveBeenSet++;
// }
// 
// void Atoms::setForces(void*p,int i) {
//   plumed_massert( dataCanBeSet,"setForces must be called after setStep in MD code interface");
//   plumed_massert( p || gatindex.size()==0, "NULL force pointer with non-zero local atoms");
//   mdatoms->setf(p,i); forcesHaveBeenSet++;
// }

void Atoms::share() {
// At first step I scatter all the atoms so as to store their mass and charge
// Notice that this works with the assumption that charges and masses are
// not changing during the simulation!
  if( needsAllAtoms() ) {
    shareAll();
    return;
  }

  if(!(int(gatindex.size())==natoms && shuffledAtoms==0)) {
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        if(!actions[i]->getUnique().empty()) {
          atomsNeeded=true;
          // unique are the local atoms
          unique.insert(actions[i]->getUniqueLocal().begin(),actions[i]->getUniqueLocal().end());
        }
      }
    }
  } else {
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        if(!actions[i]->getUnique().empty()) {
          atomsNeeded=true;
        }
      }
    }

  }

  share(unique);
}

void Atoms::shareAll() {
  unique.clear();
  // keep in unique only those atoms that are local
  if(dd && shuffledAtoms>0) {
    for(int i=0; i<natoms; i++) if(g2l[i]>=0) unique.insert(AtomNumber::index(i));
  } else {
    for(int i=0; i<natoms; i++) unique.insert(AtomNumber::index(i));
  }
  atomsNeeded=true;
  share(unique);
}

void Atoms::share(const std::set<AtomNumber>& unique) {
  // plumed_assert( positionsHaveBeenSet==3 );   // && massesHaveBeenSet );

  // virial.zero();
  if(zeroallforces || int(gatindex.size())==natoms) {
    for(int i=0; i<natoms; i++) forces[i].zero();
  } else {
    for(const auto & p : unique) forces[p.index()].zero();
  }
  for(unsigned i=getNatoms(); i<positions.size(); i++) forces[i].zero(); // virtual atoms
  forceOnEnergy=0.0;
  //mdatoms->getBox(box);

  if(!atomsNeeded) return;
  atomsNeeded=false; 
  // This is how many double per atom should be scattered and the values we are fetching
  int ndata=0; std::vector<Value*> values_to_get; std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions());

  if(int(gatindex.size())==natoms && shuffledAtoms==0) {
// faster version, which retrieves all atoms
    for(const auto & ip : inputs) {
      if( ip.second->collectFromDomains() ) { ip.second->share( 0, natoms ); values_to_get.push_back(ip.second->copyOutput(0)); ndata++; }
    }
    // mdatoms->getPositions(0,natoms,positions);
  } else {
    uniq_index.clear();
    uniq_index.reserve(unique.size());
    if(shuffledAtoms>0) {
      for(const auto & p : unique) uniq_index.push_back(g2l[p.index()]);
    }
    for(const auto & ip : inputs) {
      if( ip.second->collectFromDomains() ) { ip.second->share( unique, uniq_index ); values_to_get.push_back(ip.second->copyOutput(0)); ndata++; }
    }
    // mdatoms->getPositions(unique,uniq_index,positions);
  }

// how many double per atom should be scattered:
  // int ndata=3; std::vector<Value*> values_to_get;
  // std::vector<ActionToPutData*> & inputs(plumed.getInputActions()); 
  // for(unsigned i=0;i<inputs.size();++i) { 
  //     if( inputs[i]->collectFromDomains() ) { inputs[i]->share( gatindex ); values_to_get.push_back(inputs[i]->copyOutput(0)); ndata++; }
  // }
  // if(!massAndChargeOK) {
  //   ndata=5;
  //   masses.assign(natoms,std::numeric_limits<double>::quiet_NaN());
  //   charges.assign(natoms,std::numeric_limits<double>::quiet_NaN());
  //   mdatoms->getCharges(gatindex,charges);
  //   mdatoms->getMasses(gatindex,masses);
  // }

  if(dd && shuffledAtoms>0) {
    if(dd.async) {
      for(unsigned i=0; i<dd.mpi_request_positions.size(); i++) dd.mpi_request_positions[i].wait();
      for(unsigned i=0; i<dd.mpi_request_index.size(); i++)     dd.mpi_request_index[i].wait();
    }
    int count=0;
    for(const auto & p : unique) {
      dd.indexToBeSent[count]=p.index(); int dpoint=0;
      for(unsigned i=0;i<values_to_get.size();++i) { 
          dd.positionsToBeSent[ndata*count+dpoint]=values_to_get[i]->get(p.index()); 
          dpoint++; 
      }
      count++;
    }
    if(dd.async) {
      asyncSent=true;
      dd.mpi_request_positions.resize(dd.Get_size());
      dd.mpi_request_index.resize(dd.Get_size());
      for(int i=0; i<dd.Get_size(); i++) {
        dd.mpi_request_index[i]=dd.Isend(&dd.indexToBeSent[0],count,i,666);
        dd.mpi_request_positions[i]=dd.Isend(&dd.positionsToBeSent[0],ndata*count,i,667);
      }
    } else {
      const int n=(dd.Get_size());
      std::vector<int> counts(n);
      std::vector<int> displ(n);
      std::vector<int> counts5(n);
      std::vector<int> displ5(n);
      dd.Allgather(count,counts);
      displ[0]=0;
      for(int i=1; i<n; ++i) displ[i]=displ[i-1]+counts[i-1];
      for(int i=0; i<n; ++i) counts5[i]=counts[i]*ndata;
      for(int i=0; i<n; ++i) displ5[i]=displ[i]*ndata;
      dd.Allgatherv(&dd.indexToBeSent[0],count,&dd.indexToBeReceived[0],&counts[0],&displ[0]);
      dd.Allgatherv(&dd.positionsToBeSent[0],ndata*count,&dd.positionsToBeReceived[0],&counts5[0],&displ5[0]);
      int tot=displ[n-1]+counts[n-1];
      for(int i=0; i<tot; i++) {
        int dpoint=0;
        for(unsigned j=0;j<values_to_get.size();++j) {
            values_to_get[j]->set( dd.indexToBeReceived[i], dd.positionsToBeReceived[ndata*i+dpoint] ); 
            dpoint++; 
        }
      }
    }
  }
}

void Atoms::setPbcFromBox() {
  ActionToPutData* apb=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Box");
  plumed_assert( apb ); Tensor box; Value* myval( apb->copyOutput(0) );
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) box(i,j) = myval->get(3*i+j);
  if(dd) dd.Bcast(box,0);
  pbc.setBox(box);
}

void Atoms::wait() {
  // dataCanBeSet=false; // Everything should be set by this stage
  if(getNatoms()==0) return;
// How many double per atom should be scattered
  int ndata=0; std::vector<Value*> values_to_set;
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions());
  for(const auto & ip : inputs) { 
      if( ip.second->collectFromDomains() ) { values_to_set.push_back(ip.second->copyOutput(0)); ndata++; }  
  }

  // if(dd) {
  //   dd.Bcast(box,0);
  // }
  setPbcFromBox();

  if(dd && shuffledAtoms>0) {
// receive toBeReceived
    if(asyncSent) {
      Communicator::Status status;
      std::size_t count=0;
      for(int i=0; i<dd.Get_size(); i++) {
        dd.Recv(&dd.indexToBeReceived[count],dd.indexToBeReceived.size()-count,i,666,status);
        int c=status.Get_count<int>();
        dd.Recv(&dd.positionsToBeReceived[ndata*count],dd.positionsToBeReceived.size()-ndata*count,i,667);
        count+=c;
      }
      for(int i=0; i<count; i++) {
        int dpoint=0;
        for(unsigned j=0;j<values_to_set.size();++j) {
            values_to_set[j]->set(dd.indexToBeReceived[i], dd.positionsToBeReceived[ndata*i+dpoint] ); 
            dpoint++;
        }
      }
      asyncSent=false;
    }
    for(const auto & ip : inputs) {
        if( ip.second->sumOverDomains() ) dd.Sum( (ip.second->copyOutput(0))->data  );
    }
  }
// I take note that masses and charges have been set once for all
// at the beginning of the simulation.
  // if(shareMassAndChargeOnlyAtFirstStep) massAndChargeOK=true;
}

// void Atoms::updateForces() {
//   if( getNatoms()==0 ) return;
//   // plumed_assert( forcesHaveBeenSet==3 );
// 
//   if(forceOnEnergy*forceOnEnergy>epsilon) {
//     double alpha=1.0-forceOnEnergy;
//     mdatoms->rescaleForces(gatindex,alpha);
//     mdatoms->updateForces(gatindex,forces);
//   } else {
//     if(int(gatindex.size())==natoms && shuffledAtoms==0) mdatoms->updateForces(gatindex,forces);
//     else mdatoms->updateForces(unique,uniq_index,forces);
//   }
//   // if( !plumed.novirial && dd.Get_rank()==0 ) {
//   //   // plumed_assert( virialHasBeenSet );
//   //   //mdatoms->updateVirial(virial);
//   // }
// }

void Atoms::setNatoms(int n) {
  natoms=n;
  positions.resize(n);
  forces.resize(n);
  masses.resize(n);
  charges.resize(n);
  gatindex.resize(n);
  for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=i;
}


void Atoms::add(ActionAtomistic*a) {
  actions.push_back(a);
}

void Atoms::remove(ActionAtomistic*a) {
  auto f=find(actions.begin(),actions.end(),a);
  plumed_massert(f!=actions.end(),"cannot remove an action registered to atoms");
  actions.erase(f);
}


void Atoms::DomainDecomposition::enable(Communicator& c) {
  on=true;
  Set_comm(c.Get_comm());
  async=Get_size()<10;
  if(std::getenv("PLUMED_ASYNC_SHARE")) {
    std::string s(std::getenv("PLUMED_ASYNC_SHARE"));
    if(s=="yes") async=true;
    else if(s=="no") async=false;
    else plumed_merror("PLUMED_ASYNC_SHARE variable is set to " + s + "; should be yes or no");
  }
}

void Atoms::setAtomsNlocal(int n) {
  gatindex.resize(n);
  g2l.resize(natoms,-1);
  if(dd) {
// Since these vectors are sent with MPI by using e.g.
// &dd.positionsToBeSent[0]
// we make sure they are non-zero-sized so as to
// avoid errors when doing boundary check
    if(n==0) n++;
    dd.positionsToBeSent.resize(n*5,0.0);
    dd.positionsToBeReceived.resize(natoms*5,0.0);
    dd.indexToBeSent.resize(n,0);
    dd.indexToBeReceived.resize(natoms,0);
  }
}

void Atoms::setAtomsGatindex(int*g,bool fortran) {
  plumed_massert( g || gatindex.size()==0, "NULL gatindex pointer with non-zero local atoms");
  ddStep=plumed.getStep();
  if(fortran) {
    for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=g[i]-1;
  } else {
    for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=g[i];
  }
  for(unsigned i=0; i<g2l.size(); i++) g2l[i]=-1;
  if( gatindex.size()==natoms ) {
    shuffledAtoms=0;
    for(unsigned i=0; i<gatindex.size(); i++) {
      if( gatindex[i]!=i ) { shuffledAtoms=1; break; }
    }
  } else {
    shuffledAtoms=1;
  }
  if(dd) {
    dd.Sum(shuffledAtoms);
  }
  for(unsigned i=0; i<gatindex.size(); i++) g2l[gatindex[i]]=i;

  for(unsigned i=0; i<actions.size(); i++) {
    // keep in unique only those atoms that are local
    actions[i]->updateUniqueLocal();
  }
  unique.clear();
}

void Atoms::setAtomsContiguous(int start) {
  ddStep=plumed.getStep();
  for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=start+i;
  for(unsigned i=0; i<g2l.size(); i++) g2l[i]=-1;
  for(unsigned i=0; i<gatindex.size(); i++) g2l[gatindex[i]]=i;
  if(gatindex.size()<natoms) shuffledAtoms=1;
  for(unsigned i=0; i<actions.size(); i++) {
    // keep in unique only those atoms that are local
    actions[i]->updateUniqueLocal();
  }
  unique.clear();
}

//void Atoms::setRealPrecision(int p) {
//  //mdatoms=MDAtomsBase::create(p);
//}

// int Atoms::getRealPrecision()const {
//   return mdatoms->getRealPrecision();
// }

// void Atoms::MD2double(const void*m,double&d)const {
//   plumed_assert(mdatoms); mdatoms->MD2double(m,d);
// }
// void Atoms::double2MD(const double&d,void*m)const {
//   plumed_assert(mdatoms); mdatoms->double2MD(d,m);
// }

void Atoms::updateUnits() {
  // mdatoms->setUnits(units,MDUnits);
  if( natoms==0 ) return;
  // Set the units of the energy
  ActionToPutData* ap=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Energy");
  plumed_assert( ap ); ap->setUnit( MDUnits.getEnergy()/units.getEnergy() );
  // Update the units of the box
  ActionToPutData* apb=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Box");
  plumed_assert( apb ); apb->setUnit( MDUnits.getLength()/units.getLength() ); 
  apb->setForceUnit( units.getEnergy()/MDUnits.getEnergy() );
  // Update the units of the positions
  ActionToPutData* apx=plumed.getActionSet().selectWithLabel<ActionToPutData*>("posx");
  plumed_assert( apx ); apx->setUnit( MDUnits.getLength()/units.getLength() );  
  apx->setForceUnit( (units.getEnergy()/MDUnits.getEnergy())/(units.getLength()/MDUnits.getLength()) );
  ActionToPutData* apy=plumed.getActionSet().selectWithLabel<ActionToPutData*>("posy");
  plumed_assert( apy ); apy->setUnit( MDUnits.getLength()/units.getLength() );  
  apy->setForceUnit( (units.getEnergy()/MDUnits.getEnergy())/(units.getLength()/MDUnits.getLength()) );
  ActionToPutData* apz=plumed.getActionSet().selectWithLabel<ActionToPutData*>("posz");
  plumed_assert( apz ); apz->setUnit( MDUnits.getLength()/units.getLength() );  
  apz->setForceUnit( (units.getEnergy()/MDUnits.getEnergy())/(units.getLength()/MDUnits.getLength()) );
  // Update the units of the mass
  ActionToPutData* apm=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Masses");
  plumed_assert( apm ); apm->setUnit( MDUnits.getMass()/units.getMass() );
  // Update the units of the charge
  ActionToPutData* apq=plumed.getActionSet().selectWithLabel<ActionToPutData*>("Charges");
  plumed_assert( apq ); apq->setUnit( MDUnits.getCharge()/units.getCharge() );
}

void Atoms::setTimeStep(const double tstep) {
  timestep=tstep;
// The following is to avoid extra digits in case the MD code uses floats
// e.g.: float f=0.002 when converted to double becomes 0.002000000094995
// To avoid this, we keep only up to 6 significant digits after first one
  double magnitude=std::pow(10,std::floor(std::log10(timestep)));
  timestep=std::floor(timestep/magnitude*1e6)/1e6*magnitude;
}

double Atoms::getTimeStep()const {
  return timestep/units.getTime()*MDUnits.getTime();
}

void Atoms::setKbT(const double t) {
  kbT=t;
}

double Atoms::getKbT()const {
  return kbT/units.getEnergy()*MDUnits.getEnergy();
}

bool Atoms::needsAllAtoms() const {
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions()); bool getall=false;
  for(const auto & ip : inputs) {
      if( ip.second->collectAllFromDomains() ) { getall=true; }
  } 
  return getall;
}

void Atoms::createFullList(int*n) {
  if( needsAllAtoms() ) {
    *n=natoms;
    fullList.resize(natoms);
    for(unsigned i=0; i<natoms; i++) fullList[i]=i;
  } else {
// We update here the unique list defined at Atoms::unique.
// This is not very clear, and probably should be coded differently.
// Hopefully this fix the longstanding issue with NAMD.
    unique.clear();
    for(unsigned i=0; i<actions.size(); i++) {
      if(actions[i]->isActive()) {
        if(!actions[i]->getUnique().empty()) {
          atomsNeeded=true;
          // unique are the local atoms
          unique.insert(actions[i]->getUnique().begin(),actions[i]->getUnique().end());
        }
      }
    }
    fullList.resize(0);
    fullList.reserve(unique.size());
    for(const auto & p : unique) fullList.push_back(p.index());
    *n=fullList.size();
  }
}

void Atoms::getFullList(int**x) {
  if(!fullList.empty()) *x=&fullList[0];
  else *x=NULL;
}

void Atoms::clearFullList() {
  fullList.resize(0);
}

void Atoms::init() {
// Default: set domain decomposition to NO-decomposition, waiting for
// further instruction
  if(dd) {
    setAtomsNlocal(natoms);
    setAtomsContiguous(0);
  }
}

void Atoms::setup() {
  if( natoms==0 ) return ;
  // Create holders for the positions
  std::string str_natoms; Tools::convert( natoms, str_natoms ); int dim=3;
  plumed.cmd("createVector pos: PUT SHAPE=" + str_natoms + " SCATTERED PERIODIC=NO",&dim);
  // Create holder for the cell
  std::string noforce=""; 
  if( plumed.novirial || dd.Get_rank()!=0 ) noforce = " NOFORCE";
  plumed.cmd("createValue Box: PUT PERIODIC=NO SHAPE=3,3 " + noforce );
  // Create holder for the energy
  plumed.cmd("createValue Energy: PUT SUM_OVER_DOMAINS FORCES_FOR_POTENTIAL=posx,posy,posz,Box PERIODIC=NO");
  // Create holder for the masses
  plumed.cmd("createValue Masses: PUT SHAPE=" + str_natoms + " SCATTERED CONSTANT PERIODIC=NO");
  // Create holder for the charges 
  plumed.cmd("createValue Charges: PUT SHAPE=" + str_natoms + " SCATTERED CONSTANT PERIODIC=NO");
}

void Atoms::setDomainDecomposition(Communicator& comm) {
  dd.enable(comm);
}

void Atoms::resizeVectors(unsigned n) {
  positions.resize(n);
  forces.resize(n);
  masses.resize(n);
  charges.resize(n);
}

AtomNumber Atoms::addVirtualAtom(ActionAtomistic*a) {
  bool found=false;
  for( const auto & p : virtualAtomsActions ) {
    if( a==p ) { found=true; }
  }
  if( !found ) virtualAtomsActions.push_back( a );

  unsigned n=positions.size(); resizeVectors(n+1);
  return AtomNumber::index(n);
}

void Atoms::removeVirtualAtom(ActionAtomistic*a) {
  unsigned n=positions.size();
  plumed_massert(a==virtualAtomsActions[virtualAtomsActions.size()-1],"virtual atoms should be destroyed in reverse creation order");
  resizeVectors(n-1); virtualAtomsActions.pop_back();
}

void Atoms::writeBinary(std::ostream&o)const {
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions()); 
  for(const auto & ip : inputs) ip.second->writeBinary(o);
}

void Atoms::readBinary(std::istream&i) {
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions());
  for(const auto & ip : inputs) ip.second->readBinary(i);
  setPbcFromBox();
}

double Atoms::getKBoltzmann()const {
  if(naturalUnits || MDnaturalUnits) return 1.0;
  else return kBoltzmann/units.getEnergy();
}

double Atoms::getMDKBoltzmann()const {
  if(naturalUnits || MDnaturalUnits) return 1.0;
  else return kBoltzmann/MDUnits.getEnergy();
}

// void Atoms::getLocalMasses(std::vector<double>& localMasses) {
//   localMasses.resize(gatindex.size());
//   for(unsigned i=0; i<gatindex.size(); i++) localMasses[i] = masses[gatindex[i]];
// }

// void Atoms::getLocalPositions(std::vector<Vector>& localPositions) {
//   localPositions.resize(gatindex.size());
//   mdatoms->getLocalPositions(localPositions);
// }

// void Atoms::getLocalForces(std::vector<Vector>& localForces) {
//   localForces.resize(gatindex.size());
//   for(unsigned i=0; i<gatindex.size(); i++) localForces[i] = forces[gatindex[i]];
// }
// 
// void Atoms::getLocalMDForces(std::vector<Vector>& localForces) {
//   localForces.resize(gatindex.size());
//   for(unsigned i=0; i<gatindex.size(); i++) {
//     localForces[i] = mdatoms->getMDforces(i);
//   }
// }

ActionAtomistic* Atoms::getVirtualAtomsAction(AtomNumber i)const {
  bool found=false; unsigned va_action_ind=0, nn=0, va_index = i.index() - getNatoms();
  for(unsigned i=0; i<virtualAtomsActions.size(); ++i) {
    nn += virtualAtomsActions[i]->getNumberOfVirtualAtoms();
    if( va_index<nn ) { found=true; va_action_ind = i; break; }
  }
  plumed_dbg_assert( found );
  return virtualAtomsActions[va_action_ind];
}

//void Atoms::setExtraCV(const std::string &name,void*p) {
//  mdatoms->setExtraCV(name,p);
//}
//
//void Atoms::setExtraCVForce(const std::string &name,void*p) {
//  mdatoms->setExtraCVForce(name,p);
//}

//double Atoms::getExtraCV(const std::string &name) {
//  return mdatoms->getExtraCV(name);
//}
//
//void Atoms::updateExtraCVForce(const std::string &name,double f) {
//  mdatoms->updateExtraCVForce(name,f);
//}

}
