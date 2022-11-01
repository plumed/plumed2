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
#include "ActionForInterface.h"
#include "ActionToPutData.h"
#include "ActionSet.h"
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
  plumed(plumed),
  asyncSent(false)
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

void Atoms::setNatoms(int n) {
  natoms=n;
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

bool Atoms::needsAllAtoms() const {
  std::vector<ActionForInterface*> inputs(plumed.getActionSet().select<ActionForInterface*>()); bool getall=false;
  for(const auto & ip : inputs) {
      ActionToPutData* ap = dynamic_cast<ActionToPutData*>( ip );
      if( !ap && ip->firststep ) getall=true;
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

void Atoms::setDomainDecomposition(Communicator& comm) {
  dd.enable(comm);
}

}
