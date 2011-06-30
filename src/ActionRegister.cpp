#include "ActionRegister.h"
#include "Log.h"
#include "Tools.h"
#include <algorithm>
#include <cassert>


using namespace std;
using namespace PLMD;

ActionRegister& PLMD::actionRegister(){
  static ActionRegister ans;
  return ans;
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
  if(m[ao.line[0]]) action=m[ao.line[0]](ao);
  else action=NULL;
  return action;
}

void ActionRegister::log(Log & log){
  vector<string> s;
  for(mIterator it=m.begin();it!=m.end();++it)
    s.push_back((*it).first);
  sort(s.begin(),s.end());
  for(unsigned i=0;i<s.size();i++){
    log.printf("  %s\n",s[i].c_str());
  };
  if(disabled.size()>0){
    s.assign(disabled.size(),"");
    copy(disabled.begin(),disabled.end(),s.begin());
    sort(s.begin(),s.end());
    log.printf("+++++++ WARNING +++++++\n");
    log.printf("The following keywords have been registered more than once and will be disabled:\n");
    for(unsigned i=0;i<s.size();i++){
      log.printf("  - %s\n",s[i].c_str());
    };
    log.printf("+++++++ END WARNING +++++++\n");
  };
}


