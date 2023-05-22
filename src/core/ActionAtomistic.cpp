/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "PbcAction.h"
#include <vector>
#include <string>
#include "ActionWithValue.h"
#include "Group.h"
#include "ActionWithVirtualAtom.h"
#include "tools/Exception.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"

namespace PLMD {

ActionAtomistic::~ActionAtomistic() {
// empty destructor to delete unique_ptr
}

ActionAtomistic::ActionAtomistic(const ActionOptions&ao):
  Action(ao),
  unique_local_needs_update(true),
  boxValue(NULL),
  lockRequestAtoms(false),
  donotretrieve(false),
  donotforce(false),
  chargesWereSet(false)
{
  ActionWithValue* bv = plumed.getActionSet().selectWithLabel<ActionWithValue*>("Box");
  if( bv ) boxValue=bv->copyOutput(0);
  // We now get all the information about atoms that are lying about
  std::vector<ActionWithValue*> vatoms = plumed.getActionSet().select<ActionWithValue*>();
  for(const auto & vv : vatoms ) {
    ActionToPutData* ap = dynamic_cast<ActionToPutData*>(vv);
    if( ap ) {
      if( ap->getRole()=="x" ) xpos.push_back( ap->copyOutput(0) );
      if( ap->getRole()=="y" ) ypos.push_back( ap->copyOutput(0) );
      if( ap->getRole()=="z" ) zpos.push_back( ap->copyOutput(0) );
      if( ap->getRole()=="m" ) masv.push_back( ap->copyOutput(0) );
      if( ap->getRole()=="q" ) chargev.push_back( ap->copyOutput(0) );
    }
    ActionWithVirtualAtom* av = dynamic_cast<ActionWithVirtualAtom*>(vv);
    if( av || vv->getName()=="ARGS2VATOM" ) {
      xpos.push_back( vv->copyOutput( vv->getLabel() + ".x") );
      ypos.push_back( vv->copyOutput( vv->getLabel() + ".y") );
      zpos.push_back( vv->copyOutput( vv->getLabel() + ".z") );
      masv.push_back( vv->copyOutput( vv->getLabel() + ".mass") );
      chargev.push_back( vv->copyOutput( vv->getLabel() + ".charge") );
    }
  }
  if( xpos.size()!=ypos.size() || xpos.size()!=zpos.size() || xpos.size()!=masv.size() || xpos.size()!=chargev.size() )
    error("mismatch between value arrays");
//  if(atoms.getNatoms()==0) error("Cannot perform calculations involving atoms without atoms");
}

void ActionAtomistic::registerKeywords( Keywords& keys ) {
  (void) keys; // avoid warning
}


void ActionAtomistic::requestAtoms(const std::vector<AtomNumber> & a, const bool clearDep) {
  plumed_massert(!lockRequestAtoms,"requested atom list can only be changed in the prepare() method");
  int nat=a.size();
  indexes=a;
  positions.resize(nat);
  forces.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
//  int n=0; for(unsigned i=0;i<xpos.size();++i) n += xpos[i]->getNumberOfValues();
  int n=getTotAtoms();
  if(clearDep) clearDependencies();
  unique.clear(); std::vector<bool> requirements( xpos.size(), false );
  if( boxValue ) addDependency( boxValue->getPntrToAction() );
  for(unsigned i=0; i<indexes.size(); i++) {
    if(indexes[i].index()>=n) { std::string num; Tools::convert( indexes[i].serial(),num ); error("atom " + num + " out of range"); }
    unsigned valno, k; getValueIndices( indexes[i], valno, k ); requirements[valno] = true;
    if( valno==0 ) unique.push_back(indexes[i]);
    else if( k>0 ) error("action atomistic is not set up to deal with multiple vectors in position input");
  }
  // Add the dependencies to the actions that we require
  Tools::removeDuplicates(unique);
  for(unsigned i=0; i<requirements.size(); ++i ) {
    if( !requirements[i] ) continue;
    addDependency( xpos[i]->getPntrToAction() );
    addDependency( ypos[i]->getPntrToAction() );
    addDependency( zpos[i]->getPntrToAction() );
    addDependency( masv[i]->getPntrToAction() );
    addDependency( chargev[i]->getPntrToAction() );
  }
  unique_local_needs_update=true;
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

  const size_t nval=a->getNumberOfComponents();
  const size_t natoms=getNumberOfAtoms();
  std::vector<Vector> value(nval*natoms);
  std::vector<Tensor> valuebox(nval);
  std::vector<Vector> savedPositions(natoms);
  const double delta=std::sqrt(epsilon);

  for(int i=0; i<natoms; i++) for(int k=0; k<3; k++) {
      savedPositions[i][k]=positions[i][k];
      positions[i][k]=positions[i][k]+delta;
      a->calculate();
      positions[i][k]=savedPositions[i][k];
      for(int j=0; j<nval; j++) {
        value[j*natoms+i][k]=a->getOutputQuantity(j);
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
      for(int j=0; j<nval; j++) valuebox[j](i,k)=a->getOutputQuantity(j);
    }

  a->calculate();
  a->clearDerivatives();
  for(int j=0; j<nval; j++) {
    Value* v=a->copyOutput(j);
    double ref=v->get();
    if(v->hasDerivatives()) {
      for(int i=0; i<natoms; i++) for(int k=0; k<3; k++) {
          double d=(value[j*natoms+i][k]-ref)/delta;
          v->addDerivative(startnum+3*i+k,d);
        }
      Tensor virial;
      for(int i=0; i<3; i++) for(int k=0; k<3; k++)virial(i,k)= (valuebox[j](i,k)-ref)/delta;
// BE CAREFUL WITH NON ORTHOROMBIC CELL
      virial=-matmul(box.transpose(),virial);
      for(int i=0; i<3; i++) for(int k=0; k<3; k++) v->addDerivative(startnum+3*natoms+3*k+i,virial(k,i));
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
    bool ok=Tools::convertNoexcept(strings[i],atom); // this is converting strings to AtomNumbers
    if(ok) t.push_back(atom);
// here we check if this is a special symbol for MOLINFO
    if( !ok && strings[i].compare(0,1,"@")==0 ) {
      std::string symbol=strings[i].substr(1);
      if(symbol=="allatoms") {
        auto n=0; for(unsigned i=0; i<xpos.size(); ++i) n += xpos[i]->getNumberOfValues();
        t.reserve(n); for(unsigned i=0; i<n; i++) t.push_back(AtomNumber::index(i));
        ok=true;
      } else if(symbol=="mdatoms") {
        const auto n=xpos[0]->getNumberOfValues();
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
      unsigned ind = 0;
      for(unsigned j=0; j<xpos.size(); ++j) {
        if( xpos[j]->getPntrToAction()->getLabel()==strings[i] ) {
          t.push_back( AtomNumber::index(ind) ); ok=true; break;
        }
        ind = ind + xpos[j]->getNumberOfValues();
      }
    }
    if(!ok) error("it was not possible to interpret atom name " + strings[i]);
    // plumed_massert(ok,"it was not possible to interpret atom name " + strings[i]);
  }
}

void ActionAtomistic::getValueIndices( const AtomNumber& i, unsigned& valno, unsigned& k ) const {
  valno=0; k = i.index();
  for(unsigned j=0; j<xpos.size(); ++j) {
    if( k<xpos[j]->getNumberOfValues() ) { valno=j; break; }
    k = k - xpos[j]->getNumberOfValues();
  }
}

void ActionAtomistic::retrieveAtoms() {
  if( boxValue ) {
    PbcAction* pbca = dynamic_cast<PbcAction*>( boxValue->getPntrToAction() );
    plumed_assert( pbca ); pbc=pbca->pbc;
  }
  if( donotretrieve || indexes.size()==0 ) return;
  ActionToPutData* cv = dynamic_cast<ActionToPutData*>( chargev[0]->getPntrToAction() );
  chargesWereSet=cv->hasBeenSet();
  for(unsigned j=0; j<indexes.size(); j++) {
    unsigned nn, kk; getValueIndices( indexes[j], nn, kk );
    positions[j][0] = xpos[nn]->get(kk);
    positions[j][1] = ypos[nn]->get(kk);
    positions[j][2] = zpos[nn]->get(kk);
    charges[j] = chargev[nn]->get(kk);
    masses[j] = masv[nn]->get(kk);
  }
}

void ActionAtomistic::setForcesOnAtoms(const std::vector<double>& forcesToApply, unsigned& ind) {
  if( donotforce || (indexes.size()==0 && getName()!="FIXEDATOM") ) return;
  for(unsigned i=0; i<indexes.size(); ++i) {
    Vector ff;
    for(unsigned k=0; k<3; ++k) {
      plumed_dbg_massert( ind<forcesToApply.size(), "problem setting forces in " + getLabel() );
      ff[k]=forcesToApply[ind]; ind++;
    }
    addForce( indexes[i], ff );
  }
  setForcesOnCell( forcesToApply, ind );
}

void ActionAtomistic::setForcesOnCell(const std::vector<double>& forcesToApply, unsigned& ind) {
  for(unsigned i=0; i<9; ++i) {
    plumed_dbg_massert( ind<forcesToApply.size(), "problem setting forces in " + getLabel() );
    boxValue->addForce( i, forcesToApply[ind] ); ind++;
  }
}

Tensor ActionAtomistic::getVirial() const {
  Tensor vir; for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) vir[i][j] = boxValue->getForce(3*i+j);
  return vir;
}

void ActionAtomistic::readAtomsFromPDB(const PDB& pdb) {

  for(unsigned j=0; j<indexes.size(); j++) {
    if( indexes[j].index()>pdb.size() ) error("there are not enough atoms in the input pdb file");
    if( pdb.getAtomNumbers()[j].index()!=indexes[j].index() ) error("there are atoms missing in the pdb file");
    positions[j]=pdb.getPositions()[indexes[j].index()];
  }
  for(unsigned j=0; j<indexes.size(); j++) charges[j]=pdb.getBeta()[indexes[j].index()];
  for(unsigned j=0; j<indexes.size(); j++) masses[j]=pdb.getOccupancy()[indexes[j].index()];
}

unsigned ActionAtomistic::getTotAtoms()const {
  unsigned natoms = 0;
  for(unsigned i=0; i<xpos.size(); ++i ) natoms += xpos[i]->getNumberOfValues();
  return natoms;
}

Vector ActionAtomistic::getGlobalPosition(const AtomNumber& i) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); Vector pos;
  pos[0]=xpos[nn]->get(kk);
  pos[1]=ypos[nn]->get(kk);
  pos[2]=zpos[nn]->get(kk);
  return pos;
}

void ActionAtomistic::setGlobalPosition(const AtomNumber& i, const Vector& pos ) {
  unsigned nn, kk; getValueIndices( i, nn, kk );
  xpos[nn]->set(kk,pos[0]);
  ypos[nn]->set(kk,pos[1]);
  zpos[nn]->set(kk,pos[2]);
}

void ActionAtomistic::makeWhole() {
  for(unsigned j=0; j<positions.size()-1; ++j) {
    const Vector & first (positions[j]);
    Vector & second (positions[j+1]);
    second=first+pbcDistance(first,second);
  }
}

Vector ActionAtomistic::getForce( const AtomNumber& i ) const {
  unsigned nn, kk; getValueIndices( i, nn, kk ); Vector f;
  f[0]=xpos[nn]->getForce(kk);
  f[1]=ypos[nn]->getForce(kk);
  f[2]=zpos[nn]->getForce(kk);
  return f;
}

void ActionAtomistic::addForce( const AtomNumber& i, const Vector& f ) {
  unsigned nn, kk; getValueIndices( i, nn, kk );
  xpos[nn]->addForce( kk, f[0] );
  ypos[nn]->addForce( kk, f[1] );
  zpos[nn]->addForce( kk, f[2] );
}

void ActionAtomistic::getGradient( const unsigned& ind, Vector& deriv, std::map<AtomNumber,Vector>& gradients ) const {
  unsigned nn, kk; getValueIndices( indexes[ind], nn, kk );
  if( nn==0 ) { gradients[indexes[ind]] += deriv; return; }
  xpos[nn]->passGradients( deriv[0], gradients );
  ypos[nn]->passGradients( deriv[1], gradients );
  zpos[nn]->passGradients( deriv[2], gradients );
}

const std::vector<AtomNumber> & ActionAtomistic::getUniqueLocal( const bool& useunique, const std::vector<int>& g2l ) {
  if( useunique ) return unique;
  if( !unique_local_needs_update ) return unique_local;
  // Update unique local if it needs an update
  unique_local_needs_update=false; unique_local.clear();
  for(auto pp=unique.begin(); pp!=unique.end(); ++pp) {
    if(g2l[pp->index()]>=0) unique_local.push_back(*pp); // already sorted
  }
  return unique_local;
}

}
