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
#ifndef __PLUMED_core_ActionSet_h
#define __PLUMED_core_ActionSet_h

#include "Action.h"
#include <memory>

namespace PLMD {

class PlumedMain;

/// std::vector containing the sequence of Action to be done.
/// It is a vector of Action*, and as such it has the entire
/// std::vector interface. Moreover, it implements methods to extract
/// Acion* of a given type (select<T>()), NOT of a given type (selectNot<T>())
/// or to find an Action with a given label (selectWithLabel())
/// Finally, since it holds pointers, there is a clearDelete() function
/// which deletes the pointers before deleting the vector
class ActionSet:
  public std::vector<std::unique_ptr<Action>>
{
  PlumedMain& plumed;
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
};

/////
// INLINE IMPLEMENTATIONS:

template <class T>
std::vector<T> ActionSet::select()const {
  std::vector<T> ret;
  for(const auto & p : (*this)) {
    T t=dynamic_cast<T>(p.get());
    if(t) ret.push_back(t);
  };
  return ret;
}

template <class T>
T ActionSet::selectWithLabel(const std::string&s)const {
  for(const auto & p : (*this)) {
    T t=dynamic_cast<T>(p.get());
    if(t && dynamic_cast<Action*>(t)->getLabel()==s) return t;
  };
  return NULL;
}

template <class T>
std::vector<Action*> ActionSet::selectNot()const {
  std::vector<Action*> ret;
  for(const auto & p : (*this)) {
    T t=dynamic_cast<T>(p);
    if(!t) ret.push_back(p.get());
  };
  return ret;
}

template <class T>
std::string ActionSet::getLabelList() const {
  std::string outlist;
  for(const auto & p : (*this)) {
    if(dynamic_cast<T>(p.get())) outlist+=p->getLabel()+" ";
  };
  return  outlist;
}


template <class T>
std::vector<std::string> ActionSet::getLabelVector() const {
  std::vector<std::string> outlist;
  for(const auto & p : (*this)) {
    if(dynamic_cast<T>(p.get())) outlist.push_back(p->getLabel());
  };
  return  outlist;
}



}

#endif

