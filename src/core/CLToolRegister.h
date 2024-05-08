/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_core_CLToolRegister_h
#define __PLUMED_core_CLToolRegister_h

#include "RegisterBase.h"
#include "tools/Keywords.h"

namespace PLMD {

class CLTool;
class CLToolOptions;

struct CLToolRegisterPointers {
/// Pointer to a function which, given the options, create an CLTool
  typedef std::unique_ptr<CLTool>(*creator_pointer)(const CLToolOptions&);
/// Pointer to a function which, returns the keywords allowed
  typedef void(*keywords_pointer)(Keywords&);
  creator_pointer create;
  Keywords keys;
};


/// Same as ActionRegister, but for CLTools
class CLToolRegister :
  public RegisterBase<CLToolRegisterPointers> {
  typedef CLToolRegisterPointers::creator_pointer creator_pointer;
  typedef CLToolRegisterPointers::keywords_pointer keywords_pointer;
  typedef CLToolRegisterPointers Pointers;

  // this is necessary to avoid warnings on getKeys overload.
  // notice that RegisterBase<CLToolRegisterPointers>::getKeys returns
  // the list of keys (here: CLTools). conversely
  // CLToolRegister::getKeys(std::string name) returns the options
  // associated to CLTool name.
  // We hide the former here, which is actually implemented by list()
  // not to break existing code.
  using RegisterBase<CLToolRegisterPointers>::getKeys;

public:
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param cp  A pointer to a function which creates an object of that class
/// \param kp  A pointer to a function which returns the allowed keywords
  ID add(std::string key,creator_pointer cp,keywords_pointer kp);
/// Create an CLTool of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  std::unique_ptr<CLTool> create(const CLToolOptions&ao);
  std::unique_ptr<CLTool> create(const std::vector<void*> & images,const CLToolOptions&ao);
/// Returns a list of the allowed CLTools
  std::vector<std::string> list()const;
/// Print out the instructions for using the tool in html ready for input into the manual
  bool printManual(const std::string& cltool,const bool& spelling);
/// Return all the keys of this cltool
  std::vector<std::string> getKeys(const std::string& cltool)const;
};

/// Function returning a reference to the CLToolRegister.
/// \relates CLToolRegister
/// To avoid problems with order of initialization, this function contains
/// a static CLToolRegister which is built the first time the function is called.
/// In this manner, it is always initialized before it's used
CLToolRegister& cltoolRegister();

}


/// Shortcut for CLTool registration
/// \relates PLMD::CLToolRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_CLTOOL.
/// \param classname the name of the class to be registered
/// \param directive a string containing the corresponding directive
/// This macro should be used in the .cpp file of the corresponding class
#define PLUMED_REGISTER_CLTOOL(classname,directive) \
  namespace { class classname##RegisterMe{ \
    CLToolRegister::ID id; \
    static std::unique_ptr<CLTool> create(const CLToolOptions&ao) { \
      return std::make_unique<classname>(ao); \
    } \
  public: \
    classname##RegisterMe() : \
      id(PLMD::cltoolRegister().add(directive,create,classname::registerKeywords)) \
    {} \
    ~classname##RegisterMe(){PLMD::cltoolRegister().remove(id);} \
  } classname##RegisterMeObject; }


#endif

