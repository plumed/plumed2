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
#include "ActionAtomistic.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "GenericMolInfo.h"
#include <vector>
#include <string>
#include "ActionWithValue.h"
#include "Colvar.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "AverageBase.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "ActionToPutData.h"
#include "Group.h"

namespace PLMD {

ActionAtomistic::~ActionAtomistic() {
// forget the pending request
  atoms.remove(this);
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
  Action(ao),
  lockRequestAtoms(false),
  donotretrieve(false),
  donotforce(false),
  atoms(plumed.getAtoms())
{
  atoms.add(this);
  if( atoms.getNatoms()>0 ) {
      ActionWithValue* bv = plumed.getActionSet().selectWithLabel<ActionWithValue*>("Box");
      boxValue=bv->copyOutput(0); addDependency(bv); bv->copyOutput(0)->userdata.insert( getLabel() );
  }
}

void ActionAtomistic::registerKeywords( Keywords& keys ) {
  (void) keys; // avoid warning
}


void ActionAtomistic::requestAtoms(const std::vector<AtomNumber> & a, const bool clearDep) {
  plumed_massert(!lockRequestAtoms,"requested atom list can only be changed in the prepare() method");
  int nat=a.size();
  indexes=a;
  positions.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
  int n=0; for(unsigned i=0;i<atoms.posx.size();++i) n += atoms.posx[i]->getNumberOfValues(); 
  if(clearDep) clearDependencies();
  unique.clear(); std::vector<bool> requirements( atoms.posx.size(), false ); 
  if(a.size()>0 ) addDependency( boxValue->getPntrToAction() );
  for(unsigned i=0; i<indexes.size(); i++) {
    if(indexes[i].index()>=n) { std::string num; Tools::convert( indexes[i].serial(),num ); error("atom " + num + " out of range"); }
    // Check from atoms which of the atom values are required
    unsigned valno, k; atoms.getValueIndices( indexes[i], valno, k ); requirements[valno] = true; 
    if( valno==0 ) unique.insert(indexes[i]); 
    else if( k>0 ) error("action atomistic is not set up to deal with multiple vectors in position input");
  }
  // Add the dependencies to the actions that we require
  for(unsigned i=0; i<requirements.size(); ++i ) {
      if( !requirements[i] ) continue;
      addDependency( atoms.posx[i]->getPntrToAction() ); atoms.posx[i]->userdata.insert( getLabel() );
      addDependency( atoms.posy[i]->getPntrToAction() ); atoms.posy[i]->userdata.insert( getLabel() );
      addDependency( atoms.posz[i]->getPntrToAction() ); atoms.posz[i]->userdata.insert( getLabel() );
      addDependency( atoms.masses[i]->getPntrToAction() ); atoms.masses[i]->userdata.insert( getLabel() );
      addDependency( atoms.charges[i]->getPntrToAction() ); atoms.charges[i]->userdata.insert( getLabel() );
  }
  updateUniqueLocal();
  atoms.unique.clear();
}

Vector ActionAtomistic::pbcDistance(const Vector &v1,const Vector &v2)const {
  return pbc.distance(v1,v2);
}

void ActionAtomistic::pbcApply(std::vector<Vector>& dlist, unsigned max_index)const {
  pbc.apply(dlist, max_index);
}

void ActionAtomistic::calculateNumericalDerivatives( ActionWithValue* a ) {
  calculateAtomicNumericalDerivatives( a, 0 );
}

void ActionAtomistic::changeBox( const Tensor& newbox ) {
  pbc.setBox( newbox );
}

void ActionAtomistic::calculateAtomicNumericalDerivatives( ActionWithValue* a, const unsigned& startnum ) {
  if(!a) {
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"only Actions with a value can be differentiated");
  }

  unsigned nargs=0; std::vector<Value*> myvals;
  a->retrieveAllScalarValuesInLoop( getLabel(), nargs, myvals );
  const size_t natoms=getNumberOfAtoms();
  std::vector<Vector> value(myvals.size()*natoms);
  std::vector<Tensor> valuebox(myvals.size());
  std::vector<Vector> savedPositions(natoms);
  const double delta=sqrt(epsilon);

  for(int i=0; i<natoms; i++) for(int k=0; k<3; k++) {
      savedPositions[i][k]=positions[i][k];
      positions[i][k]=positions[i][k]+delta;
      a->calculate();
      positions[i][k]=savedPositions[i][k];
      for(int j=0; j<myvals.size(); j++) {
        value[j*natoms+i][k]=myvals[j]->get();
      }
    }
  Tensor box(pbc.getBox());
  for(int i=0; i<3; i++) for(int k=0; k<3; k++) {
      double arg0=box(i,k);
      for(int j=0; j<natoms; j++) positions[j]=pbc.realToScaled(positions[j]);
      box(i,k)=box(i,k)+delta;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++) positions[j]=pbc.scaledToReal(positions[j]);
      a->calculate();
      box(i,k)=arg0;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++) positions[j]=savedPositions[j];
      for(int j=0; j<myvals.size(); j++) valuebox[j](i,k)=myvals[j]->get();
    }

  a->calculate();
  a->clearDerivatives();
  for(int j=0; j<myvals.size(); j++) {
    double ref=myvals[j]->get();
    if(myvals[j]->hasDerivatives()) {
      for(int i=0; i<natoms; i++) for(int k=0; k<3; k++) {
          double d=(value[j*natoms+i][k]-ref)/delta;
          myvals[j]->addDerivative(startnum+3*i+k,d);
        }
      Tensor virial;
      for(int i=0; i<3; i++) for(int k=0; k<3; k++)virial(i,k)= (valuebox[j](i,k)-ref)/delta;
// BE CAREFUL WITH NON ORTHOROMBIC CELL
      virial=-matmul(box.transpose(),virial);
      for(int i=0; i<3; i++) for(int k=0; k<3; k++) myvals[j]->addDerivative(startnum+3*natoms+3*k+i,virial(k,i));
    }
  }
}

void ActionAtomistic::parseAtomList(const std::string&key, std::vector<AtomNumber> &t) {
  parseAtomList(key,-1,t);
}

void ActionAtomistic::parseAtomList(const std::string&key,const int num, std::vector<AtomNumber> &t) {
  plumed_massert( keywords.style(key,"atoms") || keywords.style(key,"hidden"), "keyword " + key + " should be registered as atoms");
  std::vector<std::string> strings;
  if( num<0 ) {
    parseVector(key,strings);
    if(strings.empty()) return;
  } else {
    if ( !parseNumberedVector(key,num,strings) ) return;
  }
  t.resize(0); interpretAtomList( strings, t );
}

void ActionAtomistic::interpretAtomList(std::vector<std::string>& strings, std::vector<AtomNumber> &t) {
  Tools::interpretRanges(strings);
  for(unsigned i=0; i<strings.size(); ++i) {
    AtomNumber atom;
    bool ok=Tools::convert(strings[i],atom); // this is converting strings to AtomNumbers
    if(ok) t.push_back(atom);
// here we check if this is a special symbol for MOLINFO
    if( !ok && strings[i].compare(0,1,"@")==0 ) {
      std::string symbol=strings[i].substr(1);
      if(symbol=="allatoms") {
        auto n=0; for(unsigned i=0;i<atoms.posx.size();++i) n += atoms.posx[i]->getNumberOfValues(); 
        t.reserve(n); for(unsigned i=0; i<n; i++) t.push_back(AtomNumber::index(i));
        ok=true;
      } else if(symbol=="mdatoms") {
        const auto n=atoms.posx[0]->getNumberOfValues();  
        t.reserve(t.size()+n);
        for(unsigned i=0; i<n; i++) t.push_back(AtomNumber::index(i));
        ok=true;
      } else {
        auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
        if( moldat ) {
          std::vector<AtomNumber> atom_list; moldat->interpretSymbol( symbol, atom_list );
          if( atom_list.size()>0 ) { ok=true; t.insert(t.end(),atom_list.begin(),atom_list.end()); }
          else { error(strings[i] + " is not a label plumed knows"); }
        } else {
          error("atoms specified using @ symbol but no MOLINFO was available");
        }
      }
    }
// here we check if the atom name is the name of a group
    if(!ok) {
      Group* mygrp=plumed.getActionSet().selectWithLabel<Group*>(strings[i]);
      if(mygrp) {
        std::vector<std::string> grp_str( mygrp->getGroupAtoms() );
        interpretAtomList( grp_str, t ); ok=true;
      } else {
        Group* mygrp2=plumed.getActionSet().selectWithLabel<Group*>(strings[i]+"_grp");
        if(mygrp2) {
           std::vector<std::string> grp_str( mygrp2->getGroupAtoms() );
           interpretAtomList( grp_str, t ); ok=true;
        }
      }
    }
// here we check if the atom name is the name of an added virtual atom
    if(!ok) {
      for(unsigned j=0; j<atoms.names.size(); ++j) {
          if( atoms.names[j]==strings[i] ) { 
              t.push_back( AtomNumber::index(atoms.posx[0]->getNumberOfValues()+j-1) ); ok=true;
          }
      }
    }
  }
}

void ActionAtomistic::retrieveAtoms() {
  pbc=atoms.pbc;
  if(donotretrieve || indexes.size()==0 ) return;
  ActionToPutData* cv = dynamic_cast<ActionToPutData*>( atoms.charges[0]->getPntrToAction() );
  chargesWereSet=cv->hasBeenSet(); 
  for(unsigned j=0; j<indexes.size(); j++) {
     positions[j]=atoms.getPosition(indexes[j]);
     charges[j]=atoms.getCharge(indexes[j]);
     masses[j]=atoms.getMass(indexes[j]);
  } 
}

void ActionAtomistic::setForcesOnAtoms( const std::vector<double>& forcesToApply, unsigned& ind ) {
  if( donotforce || (indexes.size()==0 && getName()!="FIXEDATOM") ) return;
  for(unsigned i=0; i<indexes.size(); ++i) {
    Vector ff;
    for(unsigned k=0; k<3; ++k) { 
        plumed_dbg_massert( ind<forcesToApply.size(), "problem setting forces in " + getLabel() ); 
        ff[k]=forcesToApply[ind]; ind++; 
    }
    atoms.addForce( indexes[i], ff );
  }
  for(unsigned i=0;i<9;++i) { 
     plumed_dbg_massert( ind<forcesToApply.size(), "problem setting forces in " + getLabel() );
     boxValue->addForce( i, forcesToApply[ind] ); ind++; 
  }
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) virial(i,j) = forcesToApply[3*indexes.size()+3*i+j]; 
}

void ActionAtomistic::readAtomsFromPDB( const PDB& pdb ) {
  for(unsigned j=0; j<indexes.size(); j++) {
    if( indexes[j].index()>pdb.size() ) error("there are not enough atoms in the input pdb file");
    if( pdb.getAtomNumbers()[j].index()!=indexes[j].index() ) error("there are atoms missing in the pdb file");
    positions[j]=pdb.getPositions()[indexes[j].index()];
  }
  for(unsigned j=0; j<indexes.size(); j++) charges[j]=pdb.getBeta()[indexes[j].index()];
  for(unsigned j=0; j<indexes.size(); j++) masses[j]=pdb.getOccupancy()[indexes[j].index()];
}

unsigned ActionAtomistic::getTotAtoms()const {
  return atoms.masses[0]->getNumberOfValues();
}

Vector ActionAtomistic::getGlobalPosition(AtomNumber i) const {
  return atoms.getPosition(i);
}

void ActionAtomistic::setGlobalPosition(AtomNumber i, const Vector& pos ) {
  atoms.setPosition(i, pos );
}

void ActionAtomistic::makeWhole( const unsigned start, const unsigned end ) {
  plumed_dbg_assert( start>=0 && end<=positions.size() );
  unsigned ss=0, ff=positions.size();
  if( start>0 ) ss = start; if( end>0 ) ff = end; 
  for(unsigned j=ss; j<ff-1; ++j) {
    const Vector & first (positions[j]);
    Vector & second (positions[j+1]); 
    second=first+pbcDistance(first,second);
  }
}

void ActionAtomistic::updateUniqueLocal() {
  unique_local.clear();
  if(atoms.dd && atoms.shuffledAtoms>0) {
    for(auto pp=unique.begin(); pp!=unique.end(); ++pp) {
      if(atoms.g2l[pp->index()]>=0) unique_local.insert(*pp);
    }
  } else {
    unique_local.insert(unique.begin(),unique.end());
  }
}

void ActionAtomistic::addForce( const unsigned& i, const Vector& f ) {
  atoms.addForce( indexes[i], f );
}

void ActionAtomistic::addVirial( const Tensor& v ) {
  if( !boxValue ) return ;
  virial = v;
  for(unsigned j=0; j<3; ++j) for(unsigned k=0; k<3; ++k) boxValue->addForce( 3*j+k, v[j][k] );
}

}
