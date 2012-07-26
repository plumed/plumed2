/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ActionWithDistribution.h"
#include "Vessel.h"
#include "PlumedException.h"
#include "Value.h"
#include "Log.h"

namespace PLMD {

VesselRegister::~VesselRegister(){
  if(m.size()>0){
    std::string names="";
    for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p) names+=p->first+" ";
    plumed_merror("Directive "+ names +" has not been properly unregistered");
  }
}

VesselRegister& vesselRegister(){
  static VesselRegister ans;
  return ans;
}

void VesselRegister::remove(creator_pointer f){
  for(std::map<std::string,creator_pointer>::iterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void VesselRegister::add(std::string keyword,creator_pointer f,keyword_pointer k){
  plumed_massert(m.count(keyword)==0,"keyword has already been registered");
  m.insert(std::pair<std::string,creator_pointer>(keyword,f));
  k( keywords );   // Store the keywords for all the things
}

bool VesselRegister::check(std::string key){
  if( m.count(key)>0 ) return true;
  return false;
}

Vessel* VesselRegister::create(std::string keyword, const VesselOptions&da){
  Vessel* df;
  if(check(keyword)) df=m[keyword](da);
  else df=NULL;
  return df;
}

Keywords VesselRegister::getKeywords(){
  return keywords;
}

VesselOptions::VesselOptions(const std::string& thisname, const std::string& params, ActionWithDistribution* aa ):
myname(thisname),
action(aa),
parameters(params)
{
} 

Vessel::Vessel( const VesselOptions& da ):
myname(da.myname),
label("unset"),
action(da.action),
log((da.action)->log),
comm((da.action)->comm),
serial((da.action)->serial)
{
}

void Vessel::setLabel( const std::string& mylab ){
  plumed_massert( label=="unset", "label has already been set");
  label=mylab;
}

bool Vessel::getLabel( std::string& mylab ) const {
  if( label=="unset" ) return false;
  mylab=label;
  return true;
}

void Vessel::error( const std::string& msg ){
  action->log.printf("ERROR for keyword %s in action %s with label %s : %s \n \n",myname.c_str(), (action->getName()).c_str(), (action->getLabel()).c_str(), msg.c_str() );
  printKeywords();
  plumed_merror("ERROR for keyword " + myname + " in action "  + action->getName() + " with label " + action->getLabel() + " : " + msg );
  action->exit(1);
}

}
