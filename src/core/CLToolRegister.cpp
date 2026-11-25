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
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "CLTool.h"
#include <algorithm>
#include <iostream>


namespace PLMD {

CLToolRegister& cltoolRegister() {
  static CLToolRegister ans;
  return ans;
}

std::unique_ptr<CLTool> CLToolRegister::create(const CLToolOptions&ao) {
  std::vector<void*> images; // empty vector
  return create(images,ao);
}

std::unique_ptr<CLTool> CLToolRegister::create(const std::vector<void*> & images,const CLToolOptions&ao) try {
  if(ao.line.size()<1) {
    return nullptr;
  }
  auto & content=get(images,ao.line[0]);
  CLToolOptions nao( ao,content.keys );
  return content.create(nao);
} catch (PLMD::ExceptionRegisterError &e ) {
  auto& toolName = e.getMissingKey();
  throw e <<"CL tool \"" << toolName << "\" is not known.";
}

CLToolRegister::ID CLToolRegister::add(std::string key,creator_pointer cp,keywords_pointer kp) {
  // this force each action to be registered as an uppercase string
  if ( std::any_of( std::begin( key ), std::end( key ), []( char c ) {
  return ( std::isupper( c ) )
           ;
  } ) ) plumed_error() << "CLTool: " + key + " cannot be registered, use only LOWERCASE characters";

  Keywords keys;
  kp(keys);
  return RegisterBase::add(key,Pointers{cp,keys});
}

bool CLToolRegister::printManual( const std::string& cltool, const bool& spelling ) {
  if( spelling && check(cltool) ) {
    auto cl=get(cltool);
    cl.keys.print_spelling();
    return true;
  } else if ( check(cltool) ) {
    auto cl=get(cltool);
    cl.keys.print_html();
    return true;
  } else {
    return false;
  }
}

std::vector<std::string> CLToolRegister::getKeys(const std::string& cltool)const {
  if ( check(cltool) ) {
    auto cl=get(cltool);
    auto k=cl.keys.getKeys();
    std::cerr<<k.size()<<"\n";
    for(unsigned i=0; i<k.size(); i++) {
      std::cerr<<k[i]<<"\n";
    }
    return k;
  } else {
    std::vector<std::string> empty;
    return empty;
  }
}


std::vector<std::string> CLToolRegister::list()const {
  return RegisterBase::getKeys();
}

}
