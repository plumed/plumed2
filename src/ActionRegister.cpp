#include "ActionRegister.h"
#include "Log.h"
#include "Tools.h"
#include <algorithm>
#include <iostream>


using namespace std;
using namespace PLMD;

ActionRegister::~ActionRegister(){
  if(m.size()>0)
    for(mIterator p=m.begin();p!=m.end();++p)
       std::cerr<<"WARNING: directive "<<p->first<<" has not been properly unregistered\n";
}

ActionRegister& PLMD::actionRegister(){
  static ActionRegister ans;
  return ans;
}

void ActionRegister::remove(creator_pointer f){
  for(mIterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void ActionRegister::add(string key,creator_pointer f){
  if(m.count(key)){
    m.erase(key);
    disabled.insert(key);
  }else{
    m.insert(pair<string,creator_pointer>(key,f));
  };
}

bool ActionRegister::check(string key){
  if(m.count(key)>0) return true;
  return false;
}

Action* ActionRegister::create(const ActionOptions&ao){
  if(ao.line.size()<1)return NULL;
  Action* action;
  if(check(ao.line[0])) action=m[ao.line[0]](ao);
  else action=NULL;
  return action;
}


std::ostream & PLMD::operator<<(std::ostream &log,const ActionRegister&ar){
  vector<string> s;
  for(ActionRegister::const_mIterator it=ar.m.begin();it!=ar.m.end();++it)
    s.push_back((*it).first);
  sort(s.begin(),s.end());
  for(unsigned i=0;i<s.size();i++) log<<"  "<<s[i]<<"\n";
  if(ar.disabled.size()>0){
    s.assign(ar.disabled.size(),"");
    copy(ar.disabled.begin(),ar.disabled.end(),s.begin());
    sort(s.begin(),s.end());
    log<<"+++++++ WARNING +++++++\n";
    log<<"The following keywords have been registered more than once and will be disabled:\n";
    for(unsigned i=0;i<s.size();i++) log<<"  - "<<s[i]<<"\n";
    log<<"+++++++ END WARNING +++++++\n";
  };
  return log;
}


