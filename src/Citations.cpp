#include "Citations.h"
#include "PlumedException.h"
#include "Tools.h"
#include <iostream>

using namespace std;
using namespace PLMD;

std::string Citations::cite(const std::string & item){
  unsigned i;
  for(i=0;i<items.size();++i) if(items[i]==item) break;
  if(i==items.size()) items.push_back(item);
  plumed_assert(i<items.size());
  string ret;
  Tools::convert(i+1,ret);
  ret="["+ret+"]";
  return ret;
}

std::ostream & PLMD::operator<<(std::ostream &log,const Citations&cit){
  for(unsigned i=0;i<cit.items.size();++i)
    log<<"  ["<<i+1<<"] "<<cit.items[i]<<"\n";
  return log;
}


