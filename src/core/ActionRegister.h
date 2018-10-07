/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
/// Action is a friend so that we can access the keywords for all actions so we know what to do with UPDATE_FROM/UPDATE_UNTIL and RESTART
  friend class Action;
/// Write on a stream the list of registered directives
  friend std::ostream &operator<<(std::ostream &,const ActionRegister&);
/// Pointer to a function which, given the options, create an Action
  typedef std::unique_ptr<Action>(*creator_pointer)(const ActionOptions&);
/// Pointer to a function which, returns the keywords allowed
  typedef void(*keywords_pointer)(Keywords&);
/// Pointer to a function which expands a shortcut
  typedef void(*shortcut_pointer)(const std::string&, const std::vector<std::string>&, const std::map<std::string,std::string>&, std::vector<std::vector<std::string> >&);
/// Map action to a function which creates the related object
  std::map<std::string,creator_pointer> m;
/// Map action to a function which documents the related object
  std::map<std::string,keywords_pointer> mk;
/// Map action name to a function that documents the shortcut
  std::map<std::string,keywords_pointer> sk;
/// Map action name to a function that expands the shortcut
  std::map<std::string,shortcut_pointer> s;
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
/// Add a shortcut that expands a simpler action into something more complex
  void addShortcut(std::string key, keywords_pointer kp, shortcut_pointer sp );
/// Get the list of shortcuts for a particular action
  void getShortcutKeywords( const std::string name, Keywords& keys );
/// Remove a shortcut
  void removeShortcut(std::string key);
/// Check if an action has a shortcut associated with it
  bool checkForShortcut(std::string action);
/// And expand any shortcuts in action input
  std::vector<std::vector<std::string> > expandShortcuts( const unsigned& replica_index, std::vector<std::string>& words );
/// Create an Action of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  std::unique_ptr<Action> create(const ActionOptions&ao);
/// Print out the keywords for an action in html/vim ready for input into the manual
  bool printManual(const std::string& action, const bool& vimout);
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

/// Shortcut for registering shortcuts
/// \relates PLMD::ActionRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_SHORTCUT.
/// \param classname the name of the class for which a shortcut is to registered
/// \param directive a string containing the corresponding directive for this action
/// This macro should be used in the .cpp file of the class for which the shortcut is created
#define PLUMED_REGISTER_SHORTCUT(classname,directive) \
  static class  PLUMED_UNIQUENAME(classname##ShortcutMe){ \
  public: \
    PLUMED_UNIQUENAME(classname##ShortcutMe)(){PLMD::actionRegister().addShortcut(directive,classname::shortcutKeywords,classname::expandShortcut);} \
    ~PLUMED_UNIQUENAME(classname##ShortcutMe)(){PLMD::actionRegister().removeShortcut(directive);} \
  } PLUMED_UNIQUENAME(classname##ShortcutMe);

#endif

