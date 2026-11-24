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
#include "ActionSet.h"
#include "ActionForInterface.h"
#include "ActionToPutData.h"
#include "ActionWithValue.h"
#include "ActionWithVirtualAtom.h"
#include "GenericMolInfo.h"
#include <algorithm>

namespace PLMD {

ActionSet::ActionSet(PlumedMain&p):
  plumed(p) {
  (void) plumed; // to suppress warning about "unused plumed"
}

ActionSet::~ActionSet() {
// required in order to deallocate in reverse order:
  for(int i=size()-1; i>=0; i--) {
    (*this)[i].reset();
  }
}

void ActionSet::clearDelete() {
  for(int i=size()-1; i>=0; i--) {
    (*this)[i].reset();
  }
  clear();
}

ActionShortcut* ActionSet::getShortcutActionWithLabel( const std::string& s ) const {
  updateIndex();
  auto found=index.shortcuts.find(s);
  if(found==index.shortcuts.end()) {
    return NULL;
  }
  return found->second;
}

std::size_t ActionSet::getNumberOfActionsForInterface() const {
  updateIndex();
  return index.nforinterface;
}

GenericMolInfo* ActionSet::selectLatestMolInfo(const Action*action) const {
  updateIndex();
  // position of action in the set (actions that are not in the set yet, e.g. the one
  // being constructed, come after all the others)
  std::size_t pos=size();
  auto found=index.positions.find(action);
  if(found!=index.positions.end()) {
    pos=found->second;
  }
  GenericMolInfo* latest=nullptr;
  for(const auto & m : index.molinfos) {
    if(m.first>=pos) {
      break;
    }
    latest=m.second;
  }
  return latest;
}

void ActionSet::getAtomValues( std::vector<Value*>& xpos, std::vector<Value*>& ypos, std::vector<Value*>& zpos, std::vector<Value*>& masv, std::vector<Value*>& chargev ) const {
  updateIndex();
  xpos.insert(xpos.end(),index.xpos.begin(),index.xpos.end());
  ypos.insert(ypos.end(),index.ypos.begin(),index.ypos.end());
  zpos.insert(zpos.end(),index.zpos.begin(),index.zpos.end());
  if( index.foundpdb ) {
    // masses and charges read from a PDB file replace the ones passed by the MD code
    masv.insert(masv.end(),index.pdbmass.begin(),index.pdbmass.end());
    masv.insert(masv.end(),index.vatommass.begin(),index.vatommass.end());
    chargev.insert(chargev.end(),index.pdbcharge.begin(),index.pdbcharge.end());
    chargev.insert(chargev.end(),index.vatomcharge.begin(),index.vatomcharge.end());
  } else {
    masv.insert(masv.end(),index.mass.begin(),index.mass.end());
    chargev.insert(chargev.end(),index.charge.begin(),index.charge.end());
  }
}

void ActionSet::updateIndex() const {
  if( size()<index.nindexed || (index.nindexed>0 && (*this)[index.nindexed-1].get()!=index.last) ) {
    // the set has been cleared or modified: start again
    index=Index();
  }
  // look up a label among the actions indexed so far (selectWithLabel would call updateIndex again)
  auto indexedWithLabel=[this](const std::string& l) -> ActionWithValue* {
    auto found=index.labels.find(l);
    if(found!=index.labels.end()) {
      for(const auto & a : found->second) {
        if(auto* av=dynamic_cast<ActionWithValue*>(a)) {
          return av;
        }
      }
    }
    return nullptr;
  };
  for(std::size_t i=index.nindexed; i<size(); ++i) {
    Action* p=(*this)[i].get();
    index.labels[p->getLabel()].push_back(p);
    index.positions[p]=i;
    if( auto* sc=dynamic_cast<ActionShortcut*>(p) ) {
      // the last shortcut with a given label is returned
      index.shortcuts[sc->shortcutlabel]=sc;
      if( sc->getName()=="READMASSCHARGE" ) {
        index.foundpdb=true;
        ActionWithValue* mv=indexedWithLabel( sc->getShortcutLabel() + "_mass");
        plumed_assert( mv );
        index.pdbmass.push_back( mv->copyOutput(0) );
        ActionWithValue* qv=indexedWithLabel( sc->getShortcutLabel() + "_charges");
        plumed_assert( qv );
        index.pdbcharge.push_back( qv->copyOutput(0) );
      }
    }
    if( auto* mi=dynamic_cast<GenericMolInfo*>(p) ) {
      index.molinfos.emplace_back(i,mi);
    }
    if( dynamic_cast<ActionForInterface*>(p) ) {
      index.nforinterface++;
    }
    if( auto* vv=dynamic_cast<ActionWithValue*>(p) ) {
      if( ActionToPutData* ap=vv->castToActionToPutData() ) {
        const std::string role=ap->getRole();
        if( role=="x" ) {
          index.xpos.push_back( ap->copyOutput(0) );
        }
        if( role=="y" ) {
          index.ypos.push_back( ap->copyOutput(0) );
        }
        if( role=="z" ) {
          index.zpos.push_back( ap->copyOutput(0) );
        }
        if( role=="m" ) {
          index.mass.push_back( ap->copyOutput(0) );
        }
        if( role=="q" ) {
          index.charge.push_back( ap->copyOutput(0) );
        }
      }
      if( vv->castToActionWithVirtualAtom() || vv->getName()=="ARGS2VATOM" ) {
        const std::string & l=vv->getLabel();
        index.xpos.push_back( vv->copyOutput( l + ".x") );
        index.ypos.push_back( vv->copyOutput( l + ".y") );
        index.zpos.push_back( vv->copyOutput( l + ".z") );
        index.mass.push_back( vv->copyOutput( l + ".mass") );
        index.charge.push_back( vv->copyOutput( l + ".charge") );
        index.vatommass.push_back( vv->copyOutput( l + ".mass") );
        index.vatomcharge.push_back( vv->copyOutput( l + ".charge") );
      }
    }
    index.nindexed=i+1;
    index.last=p;
  }
}

}
