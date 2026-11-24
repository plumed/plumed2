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
#ifndef __PLUMED_core_ActionSet_h
#define __PLUMED_core_ActionSet_h

#include "Action.h"
#include "ActionShortcut.h"
#include <memory>
#include <type_traits>
#include <unordered_map>

namespace PLMD {

class PlumedMain;
class GenericMolInfo;
class Value;

/// std::vector containing the sequence of Action to be done.
/// It is a vector of Action*, and as such it has the entire
/// std::vector interface. Moreover, it implements methods to extract
/// Acion* of a given type (select<T>()), NOT of a given type (selectNot<T>())
/// or to find an Action with a given label (selectWithLabel())
/// Finally, since it holds pointers, there is a clearDelete() function
/// which deletes the pointers before deleting the vector
class ActionSet:
  public std::vector<std::unique_ptr<Action>> {
  PlumedMain& plumed;
/// Index of the actions, used to avoid scanning the whole set (with a dynamic_cast
/// per action) every time a new action looks for labels, MOLINFO or atom values.
/// Actions are only appended to the set, so the index is extended lazily with the
/// actions added since the last lookup. It is rebuilt if the set has been cleared.
  struct Index {
    std::size_t nindexed=0;
    const Action* last=nullptr;
    std::unordered_map<std::string,std::vector<Action*>> labels;
    std::unordered_map<std::string,ActionShortcut*> shortcuts;
    std::unordered_map<const Action*,std::size_t> positions;
    std::vector<std::pair<std::size_t,GenericMolInfo*>> molinfos;
    std::size_t nforinterface=0;
/// atom values, see ActionAtomistic::getAtomValuesFromPlumedObject
    bool foundpdb=false;
    std::vector<Value*> xpos, ypos, zpos;
    std::vector<Value*> pdbmass, pdbcharge;
    std::vector<Value*> mass, charge;
    std::vector<Value*> vatommass, vatomcharge;
  };
  mutable Index index;
  void updateIndex() const;
  GenericMolInfo* selectLatestMolInfo(const Action*action) const;
public:
  explicit ActionSet(PlumedMain&p);
  ~ActionSet();
/// Clear and deletes all the included pointers.
  void clearDelete();

/// Extract pointers to all Action's of type T
/// To extract all Colvar , use select<Colvar*>();
  template <class T>
  std::vector<T> select()const;
/// Extract pointers to all Action's which are not of type T
/// E.g., to extract all noncolvars, use
///    selectNot<Colvar*>();
  template <class T>
  std::vector<Action*> selectNot()const;
/// Extract pointer to an action labeled s, only if it is of
/// type T. E.g., to extract an action labeled "pippo", use selectWithLabel<Action*>("pippo")
/// If you want it to be a Colvar, use selectWithLabel<Colvar*>(pippo). If it is
/// not found, it returns NULL
  template <class T>
  T selectWithLabel(const std::string&s)const;
/// get the labels in the list of actions in form of a string (useful to debug)
/// Only classes that can be dynamic casted to T are reported
  template <class T>
  std::string getLabelList() const;
/// get the labels in the form of a vector of strings
/// Only classes that can be dynamic casted to T are reported
  template <class T>
  std::vector<std::string> getLabelVector() const;
/// Extract pointer to last action of type T located before
/// action. If action is not in ActionSet or if there is no action of type T
/// returns NULL.
/// Typically to be used as selectLatest<Type>(this);
  template <class T>
  T selectLatest(const Action*action)const;
/// Get any shortcuts with this shortcut label
  ActionShortcut* getShortcutActionWithLabel( const std::string& s ) const ;
/// Number of actions that are ActionForInterface
  std::size_t getNumberOfActionsForInterface() const ;
/// Get the values that hold positions, masses and charges of all the atoms
/// (the MD atoms and the virtual atoms), appending them to the vectors
  void getAtomValues( std::vector<Value*>& xpos, std::vector<Value*>& ypos, std::vector<Value*>& zpos, std::vector<Value*>& masv, std::vector<Value*>& chargev ) const ;
};

/////
// INLINE IMPLEMENTATIONS:

template <class T>
std::vector<T> ActionSet::select()const {
  std::vector<T> ret;
  for(const auto & p : (*this)) {
    T t=dynamic_cast<T>(p.get());
    if(t) {
      ret.push_back(t);
    }
  };
  return ret;
}

template <class T>
T ActionSet::selectWithLabel(const std::string&s)const {
  updateIndex();
  auto found=index.labels.find(s);
  if(found==index.labels.end()) {
    return NULL;
  }
  // several actions can have the same label (a shortcut and the action it creates):
  // the first one of the requested type is returned, as when scanning the set
  for(const auto & p : found->second) {
    T t=dynamic_cast<T>(p);
    if(t) {
      return t;
    }
  }
  return NULL;
}

template <class T>
std::vector<Action*> ActionSet::selectNot()const {
  std::vector<Action*> ret;
  for(const auto & p : (*this)) {
    T t=dynamic_cast<T>(p);
    if(!t) {
      ret.push_back(p.get());
    }
  };
  return ret;
}

template <class T>
std::string ActionSet::getLabelList() const {
  std::string outlist;
  for(const auto & p : (*this)) {
    if(dynamic_cast<T>(p.get())) {
      outlist+=p->getLabel()+" ";
    }
  };
  return  outlist;
}


template <class T>
std::vector<std::string> ActionSet::getLabelVector() const {
  std::vector<std::string> outlist;
  for(const auto & p : (*this)) {
    if(dynamic_cast<T>(p.get())) {
      outlist.push_back(p->getLabel());
    }
  };
  return  outlist;
}

template <class T>
T ActionSet::selectLatest(const Action*action) const {
  if constexpr (std::is_same_v<T,GenericMolInfo*>) {
    // indexed, as this is looked up by every ActionAtomistic
    return selectLatestMolInfo(action);
  } else {
    T t=nullptr;
    for(const auto & p : (*this)) {
      if(p.get()==action) {
        return t;
      }
      T r=dynamic_cast<T>(p.get());
      if(r) {
        t=r;
      }
    }
    return t;
  }
}



}

#endif

