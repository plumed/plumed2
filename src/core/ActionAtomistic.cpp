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
#include "ActionShortcut.h"
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
  massesWereSet(false),
  chargesWereSet(false) {
  ActionWithValue* bv = plumed.getActionSet().selectWithLabel<ActionWithValue*>("Box");
  if( bv ) {
    boxValue=bv->copyOutput(0);
  }
  // We now get all the information about atoms that are lying about
  getAtomValuesFromPlumedObject( plumed, xpos, ypos, zpos, masv, chargev );
  if( xpos.size()!=ypos.size() || xpos.size()!=zpos.size() || xpos.size()!=masv.size() || xpos.size()!=chargev.size() ) {
    error("mismatch between value arrays");
  }
}

void ActionAtomistic::getAtomValuesFromPlumedObject( const PlumedMain& plumed, std::vector<Value*>& xpos, std::vector<Value*>& ypos, std::vector<Value*>& zpos, std::vector<Value*>& masv, std::vector<Value*>& chargev ) {
  std::vector<ActionShortcut*> shortcuts = plumed.getActionSet().select<ActionShortcut*>();
  bool foundpdb=false;
  for(const auto & ss : shortcuts ) {
    if( ss->getName()=="READMASSCHARGE" ) {
      foundpdb=true;
      ActionWithValue* mv = plumed.getActionSet().selectWithLabel<ActionWithValue*>( ss->getShortcutLabel() + "_mass");
      plumed_assert( mv );
      masv.push_back( mv->copyOutput(0) );
      ActionWithValue* qv = plumed.getActionSet().selectWithLabel<ActionWithValue*>( ss->getShortcutLabel() + "_charges");
      plumed_assert( qv );
      chargev.push_back( qv->copyOutput(0) );
    }
  }
  std::vector<ActionWithValue*> vatoms = plumed.getActionSet().select<ActionWithValue*>();
  for(const auto & vv : vatoms ) {
    plumed_assert(vv); // needed for following calls, see #1046
    ActionToPutData* ap = vv->castToActionToPutData();
    if( ap ) {
      if( ap->getRole()=="x" ) {
        xpos.push_back( ap->copyOutput(0) );
      }
      if( ap->getRole()=="y" ) {
        ypos.push_back( ap->copyOutput(0) );
      }
      if( ap->getRole()=="z" ) {
        zpos.push_back( ap->copyOutput(0) );
      }
      if( !foundpdb && ap->getRole()=="m" ) {
        masv.push_back( ap->copyOutput(0) );
      }
      if( !foundpdb && ap->getRole()=="q" ) {
        chargev.push_back( ap->copyOutput(0) );
      }
    }
    ActionWithVirtualAtom* av = vv->castToActionWithVirtualAtom();
    if( av || vv->getName()=="ARGS2VATOM" ) {
      xpos.push_back( vv->copyOutput( vv->getLabel() + ".x") );
      ypos.push_back( vv->copyOutput( vv->getLabel() + ".y") );
      zpos.push_back( vv->copyOutput( vv->getLabel() + ".z") );
      masv.push_back( vv->copyOutput( vv->getLabel() + ".mass") );
      chargev.push_back( vv->copyOutput( vv->getLabel() + ".charge") );
    }
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
  forces.resize(nat);
  masses.resize(nat);
  charges.resize(nat);
  atom_value_ind.resize( a.size() );
  int n=getTotAtoms();
  if(clearDep) {
    clearDependencies();
  }
  unique.clear();
  std::vector<bool> requirements( xpos.size(), false );
  if( boxValue ) {
    addDependency( boxValue->getPntrToAction() );
  }
  for(unsigned i=0; i<indexes.size(); i++) {
    if(indexes[i].index()>=n) {
      std::string num;
      Tools::convert( indexes[i].serial(),num );
      error("atom " + num + " out of range");
    }
    atom_value_ind[i] = getValueIndices( indexes[i] );
    requirements[atom_value_ind[i].first] = true;
    if( atom_value_ind[i].first==0 ) {
      unique.push_back(indexes[i]);
    } else if( atom_value_ind[i].second>0 ) {
      error("action atomistic is not set up to deal with multiple vectors in position input");
    }
  }

  atom_value_ind_grouped.clear();

  if(atom_value_ind.size()>0) {
    auto nn = atom_value_ind[0].first;
    auto kk = atom_value_ind[0].second;
    atom_value_ind_grouped.push_back(std::pair<std::size_t,std::vector<std::size_t>>(nn, {}));
    atom_value_ind_grouped.back().second.push_back(kk);
    auto prev_nn=nn;
    for(unsigned i=1; i<atom_value_ind.size(); i++) {
      auto nn = atom_value_ind[i].first;
      auto kk = atom_value_ind[i].second;
      if(nn!=prev_nn)
        atom_value_ind_grouped.push_back(std::pair<std::size_t,std::vector<std::size_t>>(nn, {}));
      atom_value_ind_grouped.back().second.push_back(kk);
      prev_nn=nn;
    }
  }

  // Add the dependencies to the actions that we require
  Tools::removeDuplicates(unique);
  value_depends.resize(0);
  for(unsigned i=0; i<requirements.size(); ++i ) {
    if( !requirements[i] ) {
      continue;
    }
    value_depends.push_back( i );
    addDependency( xpos[i]->getPntrToAction() );
    addDependency( ypos[i]->getPntrToAction() );
    addDependency( zpos[i]->getPntrToAction() );
    addDependency( masv[i]->getPntrToAction() );
    addDependency( chargev[i]->getPntrToAction() );
  }
  unique_local_needs_update=true;
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
    a=castToActionWithValue();
    plumed_massert(a,"only Actions with a value can be differentiated");
  }

  const size_t nval=a->getNumberOfComponents();
  const size_t natoms=getNumberOfAtoms();
  std::vector<Vector> value(nval*natoms);
  std::vector<Tensor> valuebox(nval);
  std::vector<Vector> savedPositions(natoms);
  const double delta=std::sqrt(epsilon);

  for(int i=0; i<natoms; i++)
    for(int k=0; k<3; k++) {
      savedPositions[i][k]=positions[i][k];
      positions[i][k]=positions[i][k]+delta;
      a->calculate();
      positions[i][k]=savedPositions[i][k];
      for(int j=0; j<nval; j++) {
        value[j*natoms+i][k]=a->getOutputQuantity(j);
      }
    }
  Tensor box(pbc.getBox());
  for(int i=0; i<3; i++)
    for(int k=0; k<3; k++) {
      double arg0=box(i,k);
      for(int j=0; j<natoms; j++) {
        positions[j]=pbc.realToScaled(positions[j]);
      }
      box(i,k)=box(i,k)+delta;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++) {
        positions[j]=pbc.scaledToReal(positions[j]);
      }
      a->calculate();
      box(i,k)=arg0;
      pbc.setBox(box);
      for(int j=0; j<natoms; j++) {
        positions[j]=savedPositions[j];
      }
      for(int j=0; j<nval; j++) {
        valuebox[j](i,k)=a->getOutputQuantity(j);
      }
    }

  a->calculate();
  a->clearDerivatives();
  for(int j=0; j<nval; j++) {
    Value* v=a->copyOutput(j);
    double ref=v->get();
    if(v->hasDerivatives()) {
      for(int i=0; i<natoms; i++)
        for(int k=0; k<3; k++) {
          double d=(value[j*natoms+i][k]-ref)/delta;
          v->addDerivative(startnum+3*i+k,d);
        }
      Tensor virial;
      for(int i=0; i<3; i++)
        for(int k=0; k<3; k++) {
          virial(i,k)= (valuebox[j](i,k)-ref)/delta;
        }
// BE CAREFUL WITH NON ORTHOROMBIC CELL
      virial=-matmul(box.transpose(),virial);
      for(int i=0; i<3; i++)
        for(int k=0; k<3; k++) {
          v->addDerivative(startnum+3*natoms+3*k+i,virial(k,i));
        }
    }
  }
}

bool ActionAtomistic::actionHasForces() {
  ActionWithValue* av = castToActionWithValue();
  if( av ) {
    return !av->doNotCalculateDerivatives();
  }
  if( indexes.size()>0 ) {
    plumed_merror("you have to overwrite the function actionHasForce to tell plumed if you method applies forces");
  }
  return true;
}

void ActionAtomistic::parseAtomList(const std::string&key, std::vector<AtomNumber> &t) {
  parseAtomList(key,-1,t);
}

void ActionAtomistic::parseAtomList(const std::string&key,const int num, std::vector<AtomNumber> &t) {
  plumed_massert( keywords.style(key,"atoms") || keywords.style(key,"hidden"), "keyword " + key + " should be registered as atoms");
  std::vector<std::string> strings;
  if( num<0 ) {
    parseVector(key,strings);
    if(strings.empty()) {
      return;
    }
  } else {
    if ( !parseNumberedVector(key,num,strings) ) {
      return;
    }
  }
  t.resize(0);
  interpretAtomList( strings, xpos, this, t );
}

void ActionAtomistic::interpretAtomList(std::vector<std::string>& strings, std::vector<AtomNumber> &t) {
  interpretAtomList( strings, xpos, this, t );
}

void ActionAtomistic::interpretAtomList(std::vector<std::string>& strings, const std::vector<Value*>& xpos, Action* action, std::vector<AtomNumber> &t) {
  Tools::interpretRanges(strings);
  for(unsigned i=0; i<strings.size(); ++i) {
    AtomNumber atom;
    bool ok=Tools::convertNoexcept(strings[i],atom); // this is converting strings to AtomNumbers
    if(ok) {
      t.push_back(atom);
    }
// here we check if this is a special symbol for MOLINFO
    if( !ok && strings[i].compare(0,1,"@")==0 ) {
      std::string symbol=strings[i].substr(1);
      if(symbol=="allatoms") {
        auto n=0;
        for(unsigned i=0; i<xpos.size(); ++i) {
          n += xpos[i]->getNumberOfValues();
        }
        t.reserve(n);
        for(unsigned i=0; i<n; i++) {
          t.push_back(AtomNumber::index(i));
        }
        ok=true;
      } else if(symbol=="mdatoms") {
        const auto n=xpos[0]->getNumberOfValues();
        t.reserve(t.size()+n);
        for(unsigned i=0; i<n; i++) {
          t.push_back(AtomNumber::index(i));
        }
        ok=true;
      } else if(Tools::startWith(symbol,"ndx:")) {
        auto words=Tools::getWords(symbol.substr(4));
        std::string ndxfile,ndxgroup;
        if(words.size()==1) {
          ndxfile=words[0];
        } else if(words.size()==2) {
          ndxfile=words[0];
          ndxgroup=words[1];
        } else {
          plumed_error()<<"Cannot intepret selection "<<symbol;
        }

        if(ndxgroup.size()>0) {
          action->log<<"  importing group '"+ndxgroup+"'";
        } else {
          action->log<<"  importing first group";
        }
        action->log<<" from index file "<<ndxfile<<"\n";

        IFile ifile;
        ifile.open(ndxfile);
        std::string line;
        std::string groupname;
        bool firstgroup=true;
        bool groupfound=false;
        while(ifile.getline(line)) {
          std::vector<std::string> words=Tools::getWords(line);
          if(words.size()>=3 && words[0]=="[" && words[2]=="]") {
            if(groupname.length()>0) {
              firstgroup=false;
            }
            groupname=words[1];
            if(groupname==ndxgroup || ndxgroup.length()==0) {
              groupfound=true;
            }
          } else if(groupname==ndxgroup || (firstgroup && ndxgroup.length()==0)) {
            for(unsigned i=0; i<words.size(); i++) {
              AtomNumber at;
              Tools::convert(words[i],at);
              t.push_back(at);
            }
          }
        }
        if(!groupfound) {
          plumed_error()<<"group has not been found in index file";
        }
        ok=true;
      } else {
        auto* moldat=action->plumed.getActionSet().selectLatest<GenericMolInfo*>(action);
        if( moldat ) {
          std::vector<AtomNumber> atom_list;
          moldat->interpretSymbol( symbol, atom_list );
          if( atom_list.size()>0 ) {
            ok=true;
            t.insert(t.end(),atom_list.begin(),atom_list.end());
          } else {
            action->error(strings[i] + " is not a label plumed knows");
          }
        } else {
          action->error("atoms specified using @ symbol but no MOLINFO was available");
        }
      }
    }
// here we check if the atom name is the name of a group
    if(!ok) {
      Group* mygrp=action->plumed.getActionSet().selectWithLabel<Group*>(strings[i]);
      if(mygrp) {
        std::vector<std::string> grp_str( mygrp->getGroupAtoms() );
        interpretAtomList( grp_str, xpos, action, t );
        ok=true;
      } else {
        Group* mygrp2=action->plumed.getActionSet().selectWithLabel<Group*>(strings[i]+"_grp");
        if(mygrp2) {
          std::vector<std::string> grp_str( mygrp2->getGroupAtoms() );
          interpretAtomList( grp_str, xpos, action, t );
          ok=true;
        }
      }
    }
// here we check if the atom name is the name of an added virtual atom
    if(!ok) {
      unsigned ind = 0;
      for(unsigned j=0; j<xpos.size(); ++j) {
        if( xpos[j]->getPntrToAction()->getLabel()==strings[i] ) {
          t.push_back( AtomNumber::index(ind) );
          ok=true;
          break;
        }
        ind = ind + xpos[j]->getNumberOfValues();
      }
    }
    if(!ok) {
      action->error("it was not possible to interpret atom name " + strings[i]);
    }
    // plumed_massert(ok,"it was not possible to interpret atom name " + strings[i]);
  }
}

std::pair<std::size_t, std::size_t> ActionAtomistic::getValueIndices( const AtomNumber& i ) const {
  std::size_t valno=0, k = i.index();
  for(unsigned j=0; j<xpos.size(); ++j) {
    if( k<xpos[j]->getNumberOfValues() ) {
      valno=j;
      break;
    }
    k = k - xpos[j]->getNumberOfValues();
  }
  return std::pair<std::size_t, std::size_t>( valno, k );
}

void ActionAtomistic::retrieveAtoms( const bool& force ) {
  if( boxValue ) {
    auto* ptr=boxValue->getPntrToAction();
    plumed_assert(ptr); // needed for following calls, see #1046
    PbcAction* pbca = ptr->castToPbcAction();
    plumed_assert( pbca );
    pbc=pbca->pbc;
  }
  if( donotretrieve || indexes.size()==0 ) {
    return;
  }
  auto * mtr=masv[0]->getPntrToAction();
  plumed_assert(mtr); // needed for following calls, see #1046
  ActionToPutData* mv = mtr->castToActionToPutData();
  if(mv) {
    massesWereSet=mv->hasBeenSet();
  } else if( (masv[0]->getPntrToAction())->getName()=="CONSTANT" ) {
    massesWereSet=true;  // Read masses from PDB file
  }
  auto * ptr=chargev[0]->getPntrToAction();
  plumed_assert(ptr); // needed for following calls, see #1046
  ActionToPutData* cv = ptr->castToActionToPutData();
  if(cv) {
    chargesWereSet=cv->hasBeenSet();
  } else if( (chargev[0]->getPntrToAction())->getName()=="CONSTANT" ) {
    chargesWereSet=true;  // Read masses from PDB file
  }
  unsigned j = 0;

// for(const auto & a : atom_value_ind) {
//   std::size_t nn = a.first, kk = a.second;
//   positions[j][0] = xpos[nn]->data[kk];
//   positions[j][1] = ypos[nn]->data[kk];
//   positions[j][2] = zpos[nn]->data[kk];
//   charges[j] = chargev[nn]->data[kk];
//   masses[j] = masv[nn]->data[kk];
//   j++;
// }

  for(const auto & a : atom_value_ind_grouped) {
    const auto nn=a.first;
    auto & xp=xpos[nn]->data;
    auto & yp=ypos[nn]->data;
    auto & zp=zpos[nn]->data;
    auto & ch=chargev[nn]->data;
    auto & ma=masv[nn]->data;
    for(const auto & kk : a.second) {
      positions[j][0] = xp[kk];
      positions[j][1] = yp[kk];
      positions[j][2] = zp[kk];
      charges[j] = ch[kk];
      masses[j] = ma[kk];
      j++;
    }
  }

}

void ActionAtomistic::setForcesOnAtoms(const std::vector<double>& forcesToApply, unsigned& ind) {
  if( donotforce || (indexes.size()==0 && getName()!="FIXEDATOM") ) {
    return;
  }
  for(unsigned i=0; i<value_depends.size(); ++i) {
    xpos[value_depends[i]]->hasForce = true;
    ypos[value_depends[i]]->hasForce = true;
    zpos[value_depends[i]]->hasForce = true;
  }

// for(const auto & a : atom_value_ind) {
//   plumed_dbg_massert( ind<forcesToApply.size(), "problem setting forces in " + getLabel() );
//   std::size_t nn = a.first, kk = a.second;
//   xpos[nn]->inputForce[kk] += forcesToApply[ind]; ind++;
//   ypos[nn]->inputForce[kk] += forcesToApply[ind]; ind++;
//   zpos[nn]->inputForce[kk] += forcesToApply[ind]; ind++;
// }

  for(const auto & a : atom_value_ind_grouped) {
    const auto nn=a.first;
    plumed_dbg_assert(ind<forcesToApply.size()) << "problem setting forces in " << getLabel();
    auto & xp=xpos[nn]->inputForce;
    auto & yp=ypos[nn]->inputForce;
    auto & zp=zpos[nn]->inputForce;
    for(const auto & kk : a.second) {
      xp[kk] += forcesToApply[ind];
      ind++;
      yp[kk] += forcesToApply[ind];
      ind++;
      zp[kk] += forcesToApply[ind];
      ind++;
    }
  }

  setForcesOnCell( forcesToApply, ind );
}

void ActionAtomistic::setForcesOnCell(const std::vector<double>& forcesToApply, unsigned& ind) {
  setForcesOnCell(forcesToApply.data(),forcesToApply.size(),ind);
}

void ActionAtomistic::setForcesOnCell(const double* forcesToApply, std::size_t size, unsigned& ind) {
  for(unsigned i=0; i<9; ++i) {
    plumed_dbg_massert( ind<size, "problem setting forces in " + getLabel() );
    boxValue->addForce( i, forcesToApply[ind] );
    ind++;
  }
}

Tensor ActionAtomistic::getVirial() const {
  Tensor vir;
  for(unsigned i=0; i<3; ++i)
    for(unsigned j=0; j<3; ++j) {
      vir[i][j] = boxValue->getForce(3*i+j);
    }
  return vir;
}

void ActionAtomistic::readAtomsFromPDB(const PDB& pdb) {

  for(unsigned j=0; j<indexes.size(); j++) {
    if( indexes[j].index()>pdb.size() ) {
      error("there are not enough atoms in the input pdb file");
    }
    if( pdb.getAtomNumbers()[j].index()!=indexes[j].index() ) {
      error("there are atoms missing in the pdb file");
    }
    positions[j]=pdb.getPositions()[indexes[j].index()];
  }
  for(unsigned j=0; j<indexes.size(); j++) {
    charges[j]=pdb.getBeta()[indexes[j].index()];
  }
  for(unsigned j=0; j<indexes.size(); j++) {
    masses[j]=pdb.getOccupancy()[indexes[j].index()];
  }
}

unsigned ActionAtomistic::getTotAtoms()const {
  unsigned natoms = 0;
  for(unsigned i=0; i<xpos.size(); ++i ) {
    natoms += xpos[i]->getNumberOfValues();
  }
  return natoms;
}

void ActionAtomistic::makeWhole() {
  for(unsigned j=0; j<positions.size()-1; ++j) {
    const Vector & first (positions[j]);
    Vector & second (positions[j+1]);
    second=first+pbcDistance(first,second);
  }
}

void ActionAtomistic::getGradient( const unsigned& ind, Vector& deriv, std::map<AtomNumber,Vector>& gradients ) const {
  std::size_t nn = atom_value_ind[ind].first;
  if( nn==0 ) {
    gradients[indexes[ind]] += deriv;
    return;
  }
  xpos[nn]->passGradients( deriv[0], gradients );
  ypos[nn]->passGradients( deriv[1], gradients );
  zpos[nn]->passGradients( deriv[2], gradients );
}

void ActionAtomistic::updateUniqueLocal( const bool& useunique, const std::vector<int>& g2l ) {
  if( useunique ) {
    unique_local=unique;
    return;
  }
  // Update unique local if it needs an update
  unique_local_needs_update=false;
  unique_local.clear();
  for(auto pp=unique.begin(); pp!=unique.end(); ++pp) {
    if(g2l[pp->index()]>=0) {
      unique_local.push_back(*pp);  // already sorted
    }
  }
}

}
