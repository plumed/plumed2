/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

#include <string>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>
#include "tools/Keywords.h"
#include <memory>

namespace PLMD {

class CLTool;
class CLToolOptions;

/// Same as ActionRegister, but for CLTools
class CLToolRegister {
/// Write on a stream the list of registered directives
  friend std::ostream &operator<<(std::ostream &,const CLToolRegister&);
/// Pointer to a function which, given the options, create an CLTool
  typedef std::unique_ptr<CLTool>(*creator_pointer)(const CLToolOptions&);
/// Pointer to a function which, returns the keywords allowed
  typedef void(*keywords_pointer)(Keywords&);
/// Map cltool to a function which creates the related object
  std::map<std::string,creator_pointer> m;
/// Map cltool name to the keywords for this function
  std::map<std::string,Keywords> mk;
/// Set of disabled cltools (which were registered more than once)
  std::set<std::string> disabled;
public:
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param cp  A pointer to a function which creates an object of that class
/// \param kp  A pointer to a function which returns the allowed keywords
  void add(std::string key,creator_pointer cp,keywords_pointer kp);
/// Verify if a directive is present in the register
  bool check(std::string cltool)const;
/// Create an CLTool of the type indicated in the options
/// \param ao object containing information for initialization, such as the full input line, a pointer to PlumedMain, etc
  std::unique_ptr<CLTool> create(const CLToolOptions&ao);
  void remove(creator_pointer);
  ~CLToolRegister();
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

std::ostream & operator<<(std::ostream &log,const CLToolRegister&ar);

}


/// Shortcut for CLTool registration
/// \relates PLMD::CLToolRegister
/// For easier registration, this file also provides a macro PLUMED_REGISTER_CLTOOL.
/// \param classname the name of the class to be registered
/// \param directive a string containing the corresponding directive
/// This macro should be used in the .cpp file of the corresponding class
#define PLUMED_REGISTER_CLTOOL(classname,directive) \
  static class classname##RegisterMe{ \
    static std::unique_ptr<PLMD::CLTool> create(const PLMD::CLToolOptions&ao){return std::unique_ptr<classname>(new classname(ao));} \
  public: \
    classname##RegisterMe(){PLMD::cltoolRegister().add(directive,create,classname::registerKeywords);} \
    ~classname##RegisterMe(){PLMD::cltoolRegister().remove(create);} \
  } classname##RegisterMeObject;


#endif

