/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "tools/Tools.h"
#include "Action.h"
#include <algorithm>
#include <iostream>


using namespace std;
namespace PLMD {

ActionRegister::~ActionRegister() {
  if(m.size()>0) {
    string names="";
    for(const auto & p : m) names+=p.first+" ";
    std::cerr<<"WARNING: Directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
  if(mk.size()>0) {
    string names="";
    for(const auto & p : mk) names+=p.first+" ";
    std::cerr<<"WARNING: Directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
  if(sk.size()>0) {
    string names="";
    for(const auto & p : sk) names+=p.first+" ";
    std::cerr<<"WARNING: Shortcut for directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
  if(s.size()>0) {
    string names="";
    for(const auto & p : s) names+=p.first+" ";
    std::cerr<<"WARNING: Shorcut for directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

ActionRegister& actionRegister() {
  static ActionRegister ans;
  return ans;
}

void ActionRegister::addShortcut(std::string key, keywords_pointer kp, shortcut_pointer sp ) {
  plumed_massert( !sk.count(key) && !s.count(key), "should only be one shortcut registered per action");
  sk.insert(pair<string,keywords_pointer>(key,kp)); s.insert(pair<string,shortcut_pointer>(key,sp));
}

void ActionRegister::removeShortcut(std::string key){
  for(auto p=sk.begin(); p!=sk.end(); ++p) {
    if((*p).first==key) { sk.erase(p); break; }
  }
  for(auto p=s.begin(); p!=s.end(); ++p) {
    if((*p).first==key) { s.erase(p); break; }
  }
}

bool ActionRegister::checkForShortcut(std::string action){
  if(s.count(action)>0 && sk.count(action)>0) return true;
  return false;
}

std::vector<std::vector<std::string> > ActionRegister::expandShortcuts( const unsigned& replica_index, std::vector<std::string>& words ){
  std::vector<std::vector<std::string> > actions;
  if( s.count(words[0])>0 ){
      plumed_assert( sk.count(words[0])>0 );
      Keywords keys; sk[words[0]](keys); std::map<std::string,std::string> keymap;
      for(unsigned i=0;i<keys.size();++i){
          std::string t, def, keyname = keys.get(i);
          if( keys.style( keyname, "compulsory") && !keys.getDefaultValue( keyname, def ) ){
              bool found=Tools::parse(words,keyname,t,replica_index);
              if( found ) keymap.insert(pair<std::string,std::string>(keyname,t));
          }
      } 
      if( keymap.size()>0 ){
         for(unsigned i=0;i<keys.size();++i){
             std::string t, def, keyname = keys.get(i);
             if( keys.getDefaultValue( keyname, def ) ){ 
                 bool found=Tools::parse(words,keyname,t,replica_index);
                 if( !found ) keymap.insert(pair<std::string,std::string>(keyname,def));
                 else keymap.insert(pair<std::string,std::string>(keyname,t));
             }
         }
         std::string lab; bool found=Tools::parse( words, "LABEL", lab, replica_index);
         plumed_assert( found ); s[words[0]](lab, words, keymap, actions );
         return actions;
      }
  }
  actions.push_back( words );
  return actions;
}

void ActionRegister::remove(creator_pointer f) {
  std::string directive;
  for(auto p=m.begin(); p!=m.end(); ++p) {
    if((*p).second==f) {
      directive=(*p).first; m.erase(p); break;
    }
  }
  for(auto p=mk.begin(); p!=mk.end(); ++p) {
    if((*p).first==directive) { mk.erase(p); break; }
  }
}

void ActionRegister::add(string key,creator_pointer f,keywords_pointer k) {
  if(m.count(key)) {
    m.erase(key);
    disabled.insert(key);
  } else {
    m.insert(pair<string,creator_pointer>(key,f));
    // Store a pointer to the function that creates keywords
    // A pointer is stored and not the keywords because all
    // Vessels must be dynamically loaded before the actions.
    mk.insert(pair<string,keywords_pointer>(key,k));
  };
}

bool ActionRegister::check(string key) {
  if(m.count(key)>0 && mk.count(key)>0) return true;
  return false;
}

std::unique_ptr<Action> ActionRegister::create(const ActionOptions&ao) {
  if(ao.line.size()<1)return NULL;
  // Create a copy of the manual locally. The manual is
  // then added to the ActionOptions. This allows us to
  // ensure during construction that all the keywords for
  // the action have been documented. In addition, we can
  // generate the documentation when the user makes an error
  // in the input.
  std::unique_ptr<Action> action;
  if( check(ao.line[0]) ) {
    Keywords keys; mk[ao.line[0]](keys);
    ActionOptions nao( ao,keys );
    action=m[ao.line[0]](nao);
    keys.destroyData();
  }
  return action;
}

bool ActionRegister::printManual( const std::string& action, const bool& vimout ) {
  if ( check(action) ) {
    Keywords keys; mk[action](keys);
    if( vimout ) {
      printf("%s",action.c_str()); keys.print_vim(); printf("\n");
    } else {
      keys.print_html();
    }
    keys.destroyData();
    return true;
  } else {
    return false;
  }
}

bool ActionRegister::printTemplate( const std::string& action, bool include_optional ) {
  if( check(action) ) {
    Keywords keys; mk[action](keys);
    keys.print_template(action, include_optional); keys.destroyData();
    return true;
  } else {
    return false;
  }
}

std::ostream & operator<<(std::ostream &log,const ActionRegister&ar) {
  vector<string> s;
  for(const auto & it : ar.m) s.push_back(it.first);
  sort(s.begin(),s.end());
  for(unsigned i=0; i<s.size(); i++) log<<"  "<<s[i]<<"\n";
  if(!ar.disabled.empty()) {
    s.assign(ar.disabled.size(),"");
    copy(ar.disabled.begin(),ar.disabled.end(),s.begin());
    sort(s.begin(),s.end());
    log<<"+++++++ WARNING +++++++\n";
    log<<"The following keywords have been registered more than once and will be disabled:\n";
    for(unsigned i=0; i<s.size(); i++) log<<"  - "<<s[i]<<"\n";
    log<<"+++++++ END WARNING +++++++\n";
  };
  return log;
}


}
