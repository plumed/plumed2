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
#ifndef __PLUMED_core_ActionRegister_h
#define __PLUMED_core_ActionRegister_h

#include "RegisterBase.h"

#include "tools/Keywords.h"

namespace PLMD {

class Action;
class ActionOptions;

struct ActionRegisterPointers {
/// Pointer to a function which, given the options, create an Action
  typedef std::unique_ptr<Action>(*creator_pointer)(const ActionOptions&);
/// Pointer to a function which, returns the keywords allowed
  typedef void(*keywords_pointer)(Keywords&);
  creator_pointer create;
  keywords_pointer keys;
};

/// Register holding all the allowed keywords.
/// This is a register which holds a map between strings (directives) and function pointers.
/// The function pointers are pointing to functions which create an object of
/// the corresponding class given the corresponding options (ActionOptions).
/// There should be only one of there objects allocated.
/// Actions should be registered here at the beginning of execution
///
class ActionRegister:
  public RegisterBase<ActionRegisterPointers> {

  typedef ActionRegisterPointers::creator_pointer creator_pointer;
  typedef ActionRegisterPointers::keywords_pointer keywords_pointer;
  typedef ActionRegisterPointers Pointers;

public:
  ID add(std::string key,creator_pointer cp,keywords_pointer kp);
/// Create an Action of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  std::unique_ptr<Action> create(const ActionOptions&ao);
  std::unique_ptr<Action> create(const std::vector<void*> & images,const ActionOptions&ao);
/// Print out the keywords for an action in html/vim ready for input into the manual
  bool printManual(const std::string& action, const bool& vimout, const bool& spellout);
/// Retrieve a keywords object for a particular action
  bool getKeywords( const std::string& action, Keywords& keys );
  void getKeywords(const std::vector<void*> & images, const std::string& action, Keywords& keys);
/// Print out a template command for an action
  bool printTemplate(const std::string& action, bool include_optional);
  std::vector<std::string> getActionNames() const;
};

/// Function returning a reference to the ActionRegister.
/// \relates ActionRegister
/// To avoid problems with order of initialization, this function contains
/// a static ActionRegister which is built the first time the function is called.
/// In this manner, it is always initialized before it's used
ActionRegister& actionRegister();

template<typename T>
inline constexpr bool isActionType = std::is_base_of<Action, T>::value;
//in C++20 you we'll make this a concept
//template<typename T>
//concept ActionType = std::is_base_of<::PLMD::Action, T>::value;
//so the template will be template<ActionType ActionType>class ActionRegistration{...}
//without the explicit need of the static assert

///Each instance of this specialized class represents an action that can be called
///with  the specified directive.
///As soon it goes out of scope it will deregister the directive from the singleton ActionRegister
template<typename ActionClass>
class ActionRegistration {
  ActionRegister::ID id;
  static std::unique_ptr<Action> create(const ActionOptions&ao) {
    return std::make_unique<ActionClass>(ao);
  }
public:
  ///On construction register the ActionClass with the wanted directive
  ActionRegistration(std::string_view directive):
    id(actionRegister().add(directive.data(),create,ActionClass::registerKeywords)) {
    static_assert(isActionType<ActionClass>,
                  "ActionRegistration accepts only class that inherit from Action");
  }
  ///On destruction deregister the ActionClass (useful when you unload a shared object)
  ~ActionRegistration() {
    actionRegister().remove(id);
  }
};
} //PLMD

#define PLUMED_CONCATENATE_DIRECT(s1, s2) s1##s2
#define PLUMED_CONCATENATE(s1, s2) PLUMED_CONCATENATE_DIRECT(s1, s2)

/// Shortcut for Action registration
/// \relates PLMD::ActionRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_ACTION.
/// \param classname the name of the class to be registered
/// \param directive a string containing the corresponding directive
/// This macro should be used in the .cpp file of the corresponding class
#define PLUMED_REGISTER_ACTION(classname,directive) \
  namespace {::PLMD::ActionRegistration<classname> \
             PLUMED_CONCATENATE(classname##Registerer,__LINE__)(directive);}
#endif

