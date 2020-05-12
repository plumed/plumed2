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
#ifndef __PLUMED_core_ActionRegister_h
#define __PLUMED_core_ActionRegister_h

#include <string>
#include <map>
#include <set>
#include <iosfwd>
#include "tools/Keywords.h"
#include <memory>

namespace PLMD {

class Action;
class ActionOptions;

/// Register holding all the allowed keywords.
/// This is a register which holds a map between strings (directives) and function pointers.
/// The function pointers are pointing to functions which create an object of
/// the corresponding class given the corresponding options (ActionOptions).
/// There should be only one of there objects allocated.
/// Actions should be registered here at the beginning of execution
/// If the same directive is used for different classes, it is automatically disabled
/// to avoid random results.
///
class ActionRegister {
/// Write on a stream the list of registered directives
  friend std::ostream &operator<<(std::ostream &,const ActionRegister&);
/// Pointer to a function which, given the options, create an Action
  typedef std::unique_ptr<Action>(*creator_pointer)(const ActionOptions&);
/// Pointer to a function which, returns the keywords allowed
  typedef void(*keywords_pointer)(Keywords&);
/// Map action to a function which creates the related object
  std::map<std::string,creator_pointer> m;
/// Map action to a function which documents the related object
  std::map<std::string,keywords_pointer> mk;
/// Set of disabled actions (which were registered more than once)
  std::set<std::string> disabled;
public:
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param cp A pointer to a function which creates an object of that class
/// \param kp A pointer to a function which returns the allowed keywords
  void add(std::string key,creator_pointer cp,keywords_pointer kp);
/// Verify if a directive is present in the register
  bool check(std::string action);
/// Create an Action of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  std::unique_ptr<Action> create(const ActionOptions&ao);
/// Print out the keywords for an action in html/vim ready for input into the manual
  bool printManual(const std::string& action, const bool& vimout, const bool& spellout);
/// Print out a template command for an action
  bool printTemplate(const std::string& action, bool include_optional);
  void remove(creator_pointer);
  ~ActionRegister();
};

/// Function returning a reference to the ActionRegister.
/// \relates ActionRegister
/// To avoid problems with order of initialization, this function contains
/// a static ActionRegister which is built the first time the function is called.
/// In this manner, it is always initialized before it's used
ActionRegister& actionRegister();

std::ostream & operator<<(std::ostream &log,const ActionRegister&ar);

}

#define PLUMED_CONCATENATE_DIRECT(s1, s2) s1##s2
#define PLUMED_CONCATENATE(s1, s2) PLUMED_CONCATENATE_DIRECT(s1, s2)
#define PLUMED_UNIQUENAME(str) PLUMED_CONCATENATE(str, __LINE__)

/// Shortcut for Action registration
/// \relates PLMD::ActionRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_ACTION.
/// \param classname the name of the class to be registered
/// \param directive a string containing the corresponding directive
/// This macro should be used in the .cpp file of the corresponding class
#define PLUMED_REGISTER_ACTION(classname,directive) \
  static class  PLUMED_UNIQUENAME(classname##RegisterMe){ \
    static std::unique_ptr<PLMD::Action> create(const PLMD::ActionOptions&ao){return std::unique_ptr<classname>(new classname(ao));} \
  public: \
    PLUMED_UNIQUENAME(classname##RegisterMe)(){PLMD::actionRegister().add(directive,create,classname::registerKeywords);} \
    ~PLUMED_UNIQUENAME(classname##RegisterMe)(){PLMD::actionRegister().remove(create);} \
  } PLUMED_UNIQUENAME(classname##RegisterMe);


#endif

