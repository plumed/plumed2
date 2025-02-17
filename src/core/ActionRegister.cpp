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
#include "ActionRegister.h"
#include "ModuleMap.h"
#include "Action.h"
#include <algorithm>

namespace PLMD {

ActionRegister& actionRegister() {
  static ActionRegister ans;
  return ans;
}

std::unique_ptr<Action> ActionRegister::create(const ActionOptions&ao) {
  std::vector<void*> images; // empty vector
  return create(images,ao);
}

std::unique_ptr<Action> ActionRegister::create(const std::vector<void*> & images,const ActionOptions&ao) try {
  if(ao.line.size()<1) {
    return nullptr;
  }

  auto content=get(images,ao.line[0]);
  Keywords keys;
  //the name of the function is not clear but does `keys.thisactname = ao.line[0];`
  keys.setDisplayName(ao.line[0]);
  content.keys(keys);
  ActionOptions nao( ao,keys );
  auto fullPath=getFullPath(images,ao.line[0]);
  nao.setFullPath(fullPath);
  return content.create(nao);
} catch (PLMD::ExceptionRegisterError &e ) {
  auto& actionName = e.getMissingKey();
  e <<"Action \"" << actionName << "\" is not known.";
  if (getModuleMap().count(actionName)>0) {
    e << "\nAn Action named \""
      <<actionName
      <<"\" is available in module \""
      << getModuleMap().at(actionName)
      << "\".\nPlease consider installing PLUMED with that module enabled.";
  }
  throw e;
}

bool ActionRegister::printManual(const std::string& action, const bool& vimout, const bool& spellout) {
  if ( check(action) ) {
    Keywords keys;
    getKeywords( action, keys );
    if( vimout ) {
      printf("%s",action.c_str());
      keys.print_vim();
      printf("\n");
    } else if( spellout ) {
      keys.print_spelling();
    } else {
      keys.print_html();
    }
    return true;
  } else {
    return false;
  }
}

bool ActionRegister::printTemplate(const std::string& action, bool include_optional) {
  //no need to insert the try/catch block: check will ensure that action is known
  if( check(action) ) {
    Keywords keys;
    //the name of the function is not clear but does `keys.thisactname = action;`
    keys.setDisplayName(action);
    get(action).keys(keys);
    keys.print_template(action, include_optional);
    return true;
  } else {
    return false;
  }
}

std::vector<std::string> ActionRegister::getActionNames() const {
  return getKeys();
}

ActionRegister::ID ActionRegister::add(std::string key,creator_pointer cp,keywords_pointer kp) {
  // this force each action to be registered as an uppercase string
  if ( std::any_of( std::begin( key ), std::end( key ), []( char c ) {
  return ( std::islower( c ) )
           ;
  } ) ) plumed_error() << "Action: " + key + " cannot be registered, use only UPPERCASE characters";
  return RegisterBase::add(key,Pointers{cp,kp});
}

bool ActionRegister::getKeywords(const std::string& action, Keywords& keys) {
  //no need to insert the try/catch block: check will ensure that action is known
  if(check(action)) {
    //the name of the function is not clear but does `keys.thisactname = action;`
    keys.setDisplayName(action);
    get(action).keys(keys);
    return true;
  }
  return false;
}

void ActionRegister::getKeywords(const std::vector<void*> & images, const std::string& action, Keywords& keys) {
  auto content=get(images,action);
  //the name of the function is not clear but does `keys.thisactname = action;`
  keys.setDisplayName(action);
  content.keys(keys);
}

}
