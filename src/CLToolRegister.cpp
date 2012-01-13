#include "CLToolRegister.h"
#include "Tools.h"
#include "CLTool.h"
#include <algorithm>
#include <iostream>


using namespace std;
using namespace PLMD;

CLToolRegister::~CLToolRegister(){
  if(m.size()>0){
    string names="";
    for(mIterator p=m.begin();p!=m.end();++p)names+=p->first+" ";
    plumed_merror("Directive "+ names +" has not been properly unregistered");
  }
}

CLToolRegister& PLMD::cltoolRegister(){
  static CLToolRegister ans;
  return ans;
}

void CLToolRegister::remove(creator_pointer f){
  for(mIterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void CLToolRegister::add(string key,creator_pointer f){
  if(m.count(key)){
    m.erase(key);
    disabled.insert(key);
  }else{
    m.insert(pair<string,creator_pointer>(key,f));
  };
}

bool CLToolRegister::check(string key){
  if(m.count(key)>0) return true;
  return false;
}

CLTool* CLToolRegister::create(const CLToolOptions&ao){
  if(ao.line.size()<1)return NULL;
  CLTool* cltool;
  if(check(ao.line[0])) cltool=m[ao.line[0]](ao);
  else cltool=NULL;
  return cltool;
}


std::ostream & PLMD::operator<<(std::ostream &log,const CLToolRegister&ar){
  vector<string> s(ar.list());
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

vector<string> CLToolRegister::list()const{
  vector<string> s;
  for(const_mIterator it=m.begin();it!=m.end();++it)
    s.push_back((*it).first);
  sort(s.begin(),s.end());
  return s;
}



