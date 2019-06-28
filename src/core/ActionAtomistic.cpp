/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "SetupMolInfo.h"
#include <vector>
#include <string>
#include "ActionWithValue.h"
#include "Colvar.h"
#include "ActionWithVirtualAtom.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"

using namespace std;

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
//  if(atoms.getNatoms()==0) error("Cannot perform calculations involving atoms without atoms");
}

void ActionAtomistic::registerKeywords( Keywords& keys ) {
  (void) keys; // avoid warning
}


void ActionAtomistic::requestAtoms(const vector<AtomNumber> & a, const bool clearDep) {
  plumed_massert(!lockRequestAtoms,"requested atom list can only be changed in the prepare() method");
  int nat=a.size();
  indexes=a;
  positions.resize(nat);
  forces.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
  int n=atoms.positions.size();
  if(clearDep) clearDependencies();
  unique.clear();
  for(unsigned i=0; i<indexes.size(); i++) {
    if(indexes[i].index()>=n) { std::string num; Tools::convert( indexes[i].serial(),num ); error("atom " + num + " out of range"); }
    if(atoms.isVirtualAtom(indexes[i])) addDependency(atoms.getVirtualAtomsAction(indexes[i]));
// only real atoms are requested to lower level Atoms class
    else unique.insert(indexes[i]);
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

  const int nval=a->getNumberOfComponents();
  const int natoms=getNumberOfAtoms();
  std::vector<Vector> value(nval*natoms);
  std::vector<Tensor> valuebox(nval);
  std::vector<Vector> savedPositions(natoms);
  const double delta=sqrt(epsilon);

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
  vector<string> strings;
  if( num<0 ) {
    parseVector(key,strings);
    if(strings.empty()) return;
  } else {
    if ( !parseNumberedVector(key,num,strings) ) return;
  }
  interpretAtomList( strings, t );
}

void ActionAtomistic::interpretAtomList( std::vector<std::string>& strings, std::vector<AtomNumber> &t) {
  Tools::interpretRanges(strings); t.resize(0);
  for(unsigned i=0; i<strings.size(); ++i) {
    AtomNumber atom;
    bool ok=Tools::convert(strings[i],atom); // this is converting strings to AtomNumbers
    if(ok) t.push_back(atom);
// here we check if this is a special symbol for MOLINFO
    if( !ok && strings[i].compare(0,1,"@")==0 ) {
      std::string symbol=strings[i].substr(1);
      if(symbol=="allatoms") {
        const auto n=plumed.getAtoms().getNatoms() + plumed.getAtoms().getNVirtualAtoms();
        t.reserve(t.size()+n);
        for(unsigned i=0; i<n; i++) t.push_back(AtomNumber::index(i));
        ok=true;
      } else if(symbol=="mdatoms") {
        const auto n=plumed.getAtoms().getNatoms();
        t.reserve(t.size()+n);
        for(unsigned i=0; i<n; i++) t.push_back(AtomNumber::index(i));
        ok=true;
      } else {
        vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
        if( moldat.size()>0 ) {
          vector<AtomNumber> atom_list; moldat[0]->interpretSymbol( symbol, atom_list );
          ok=true; t.insert(t.end(),atom_list.begin(),atom_list.end());
        } else {
          error("atoms specified using @ symbol but no MOLINFO was available");
        }
      }
    }
// here we check if the atom name is the name of a group
    if(!ok) {
      if(atoms.groups.count(strings[i])) {
        const auto m=atoms.groups.find(strings[i]);
        t.insert(t.end(),m->second.begin(),m->second.end());
        ok=true;
      }
    }
// here we check if the atom name is the name of an added virtual atom
    if(!ok) {
      const ActionSet&actionSet(plumed.getActionSet());
      for(const auto & a : actionSet) {
        ActionWithVirtualAtom* c=dynamic_cast<ActionWithVirtualAtom*>(a.get());
        if(c) if(c->getLabel()==strings[i]) {
            ok=true;
            t.push_back(c->getIndex());
            break;
          }
      }
    }
    if(!ok) error("it was not possible to interpret atom name " + strings[i]);
    // plumed_massert(ok,"it was not possible to interpret atom name " + strings[i]);
  }
}


void ActionAtomistic::retrieveAtoms() {
  pbc=atoms.pbc;
  Colvar*cc=dynamic_cast<Colvar*>(this);
  if(cc && cc->checkIsEnergy()) energy=atoms.getEnergy();
  if(donotretrieve) return;
  chargesWereSet=atoms.chargesWereSet();
  const vector<Vector> & p(atoms.positions);
  const vector<double> & c(atoms.charges);
  const vector<double> & m(atoms.masses);
  for(unsigned j=0; j<indexes.size(); j++) positions[j]=p[indexes[j].index()];
  for(unsigned j=0; j<indexes.size(); j++) charges[j]=c[indexes[j].index()];
  for(unsigned j=0; j<indexes.size(); j++) masses[j]=m[indexes[j].index()];
}

void ActionAtomistic::setForcesOnAtoms( const std::vector<double>& forcesToApply, unsigned ind ) {
  if(donotforce) return;
  for(unsigned i=0; i<indexes.size(); ++i) {
    forces[i][0]=forcesToApply[ind]; ind++;
    forces[i][1]=forcesToApply[ind]; ind++;
    forces[i][2]=forcesToApply[ind]; ind++;
  }
  virial(0,0)=forcesToApply[ind]; ind++;
  virial(0,1)=forcesToApply[ind]; ind++;
  virial(0,2)=forcesToApply[ind]; ind++;
  virial(1,0)=forcesToApply[ind]; ind++;
  virial(1,1)=forcesToApply[ind]; ind++;
  virial(1,2)=forcesToApply[ind]; ind++;
  virial(2,0)=forcesToApply[ind]; ind++;
  virial(2,1)=forcesToApply[ind]; ind++;
  virial(2,2)=forcesToApply[ind];
  plumed_dbg_assert( ind+1==forcesToApply.size());
}

void ActionAtomistic::applyForces() {
  if(donotforce) return;
  vector<Vector>   & f(atoms.forces);
  Tensor           & v(atoms.virial);
  for(unsigned j=0; j<indexes.size(); j++) f[indexes[j].index()]+=forces[j];
  v+=virial;
  atoms.forceOnEnergy+=forceOnEnergy;
  if(extraCV.length()>0) atoms.updateExtraCVForce(extraCV,forceOnExtraCV);
}

void ActionAtomistic::clearOutputForces() {
  virial.zero();
  if(donotforce) return;
  for(unsigned i=0; i<forces.size(); ++i)forces[i].zero();
  forceOnEnergy=0.0;
  forceOnExtraCV=0.0;
}


void ActionAtomistic::readAtomsFromPDB( const PDB& pdb ) {
  Colvar*cc=dynamic_cast<Colvar*>(this);
  if(cc && cc->checkIsEnergy()) error("can't read energies from pdb files");

  for(unsigned j=0; j<indexes.size(); j++) {
    if( indexes[j].index()>pdb.size() ) error("there are not enough atoms in the input pdb file");
    if( pdb.getAtomNumbers()[j].index()!=indexes[j].index() ) error("there are atoms missing in the pdb file");
    positions[j]=pdb.getPositions()[indexes[j].index()];
  }
  for(unsigned j=0; j<indexes.size(); j++) charges[j]=pdb.getBeta()[indexes[j].index()];
  for(unsigned j=0; j<indexes.size(); j++) masses[j]=pdb.getOccupancy()[indexes[j].index()];
}

void ActionAtomistic::makeWhole() {
  for(unsigned j=0; j<positions.size()-1; ++j) {
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

}
