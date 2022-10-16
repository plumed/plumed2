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
  shuffledAtoms(0),
  plumed(plumed),
  naturalUnits(false),
  MDnaturalUnits(false),
  timestep(0.0),
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

void Atoms::clearAtomValues() {
  names.resize(0); posx.resize(0); posy.resize(0); posz.resize(0); masses.resize(0); charges.resize(0);
}

void Atoms::addAtomValues( const std::string& n, Value* x, Value* y, Value* z, Value* m, Value* q ) {
  plumed_massert( x->getNumberOfValues()==y->getNumberOfValues(), "inconsistent number of values in atom y values" );
  plumed_massert( x->getNumberOfValues()==z->getNumberOfValues(), "inconsistent number of values in atom z values" );
  plumed_massert( x->getNumberOfValues()==m->getNumberOfValues(), "inconsistent number of values in atom masses" );
  plumed_massert( x->getNumberOfValues()==q->getNumberOfValues(), "inconsistent number of values in atom charges" );
  names.push_back(n); posx.push_back(x); posy.push_back(y); posz.push_back(z); masses.push_back(m); charges.push_back(q);
}

std::string Atoms::getAtomString( const AtomNumber& i ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); std::string ind;
  if( nn==0 ) { Tools::convert( i.serial(), ind ); return ind; }
  return names[nn];
}

void Atoms::getValueIndices( const AtomNumber& i, unsigned& valno, unsigned& k ) const {
  valno=0; k = i.index();
  for(unsigned j=0; j<posx.size(); ++j) {
      if( k<posx[j]->getNumberOfValues() ) { valno=j; break; }
      k = k - posx[j]->getNumberOfValues();
  }
}

Vector Atoms::getPosition( const AtomNumber& i ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); Vector pos; 
  pos[0]=posx[nn]->get(kk); pos[1]=posy[nn]->get(kk);  pos[2]=posz[nn]->get(kk);
  return pos; 
}

void Atoms::setPosition( const AtomNumber& i, const Vector& pos ) {
  unsigned nn, kk; getValueIndices( i, nn, kk );  
  posx[nn]->set(kk,pos[0]); posy[nn]->set(kk,pos[1]); posz[nn]->set(kk,pos[2]);
}

double Atoms::getMass( const AtomNumber& i ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); return masses[nn]->get(kk);
}

double Atoms::getCharge( const AtomNumber& i ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); return charges[nn]->get(kk);
}

bool Atoms::checkConstant( const AtomNumber& i, const std::string& name ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); 
  if( name=="MASSES" ) return masses[nn]->isConstant();
  else if( name=="CHARGES" ) return charges[nn]->isConstant();
  plumed_error();
}

void Atoms::addForce( const AtomNumber& i, Vector f ) {
  unsigned nn, kk; getValueIndices( i, nn, kk );
  posx[nn]->addForce( kk, f[0] ); posy[nn]->addForce( kk, f[1] ); posz[nn]->addForce( kk, f[2] );
}

void Atoms::getGradient( const AtomNumber& i, Vector& deriv, std::map<AtomNumber,Vector>& gradients ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk );
  if( nn==0 ) { gradients[i] += deriv; return; }
  posx[nn]->passGradients( deriv[0], gradients );
  posy[nn]->passGradients( deriv[1], gradients );
  posz[nn]->passGradients( deriv[2], gradients );
}

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
          for(auto pp=actions[i]->getUnique().begin(); pp!=actions[i]->getUnique().end(); ++pp) {
              if( g2l[pp->index()]>=0 ) unique.insert(*pp);
          }
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

  if(!atomsNeeded) return;
  atomsNeeded=false; 
  // This is how many double per atom should be scattered and the values we are fetching
  int ndata=0; std::vector<Value*> values_to_get; std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions());

  if(int(gatindex.size())==natoms && shuffledAtoms==0) {
// faster version, which retrieves all atoms
    for(const auto & ip : inputs) {
      if( ip.second->collectFromDomains() ) { ip.second->share( 0, natoms ); values_to_get.push_back(ip.second->copyOutput(0)); ndata++; }
    }
  } else {
    uniq_index.clear();
    uniq_index.reserve(unique.size());
    if(shuffledAtoms>0) {
      for(const auto & p : unique) uniq_index.push_back(g2l[p.index()]);
    }
    for(const auto & ip : inputs) {
      if( ip.second->collectFromDomains() ) { ip.second->share( unique, uniq_index ); values_to_get.push_back(ip.second->copyOutput(0)); ndata++; }
    }
  }

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
  pbc.setBox(box);
}

void Atoms::wait() {
  if(getNatoms()==0) return;
// How many double per atom should be scattered
  int ndata=0; std::vector<Value*> values_to_set;
  std::map<std::string,ActionToPutData*> & inputs(plumed.getInputActions());
  for(const auto & ip : inputs) { 
      if( ip.second->collectFromDomains() ) { values_to_set.push_back(ip.second->copyOutput(0)); ndata++; }  
      if( dd && ip.second->isActive() && ip.second->bcastToDomains() ) { dd.Bcast((ip.second->copyOutput(0))->data,0); }
  }

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
        if( ip.second->isActive() && ip.second->getName()=="ENERGY" ) dd.Sum( (ip.second->copyOutput(0))->data  );
    }
  }
}

void Atoms::setNatoms(int n) {
  natoms=n;
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

  unique.clear();
}

void Atoms::setAtomsContiguous(int start) {
  ddStep=plumed.getStep();
  for(unsigned i=0; i<gatindex.size(); i++) gatindex[i]=start+i;
  for(unsigned i=0; i<g2l.size(); i++) g2l[i]=-1;
  for(unsigned i=0; i<gatindex.size(); i++) g2l[gatindex[i]]=i;
  if(gatindex.size()<natoms) shuffledAtoms=1;
  unique.clear();
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

void Atoms::setDomainDecomposition(Communicator& comm) {
  dd.enable(comm);
}

double Atoms::getKBoltzmann()const {
  if(naturalUnits || MDnaturalUnits) return 1.0;
  else return kBoltzmann/units.getEnergy();
}

double Atoms::getMDKBoltzmann()const {
  if(naturalUnits || MDnaturalUnits) return 1.0;
  else return kBoltzmann/MDUnits.getEnergy();
}

}
