#include "Tools.h"
#include <cstring>
#include <assert.h>
#include <sstream>

using namespace PLMD;
using namespace std;

#define IMPLEMENT(T) \
bool Tools::convert(const string & str,T & t){ \
        istringstream istr(str.c_str()); \
        bool ok=istr>>t; \
        if(!ok) return false; \
        string remaining; \
        istr>>remaining; \
        return remaining.length()==0; \
}

IMPLEMENT(double)
IMPLEMENT(int)

bool Tools::convert(const string & str,string & t){
        t=str;
        return true;
}

vector<string> Tools::getWords(const string & line,const char* separators){
  char* ww;
  char* s3;
  char* copy;
  copy= new char[strlen(line.c_str())+1];
  strcpy(copy,line.c_str());
  vector<string> words;
  ww=strtok_r(copy,separators,&s3);
  if(ww){
    words.push_back(string(ww));
    while((ww=strtok_r(NULL,separators,&s3))) words.push_back(string(ww));
  }
  delete [] copy;
  return words;
}

bool Tools::getParsedLine(FILE* fp,vector<string> & words){
  string line("");
  words.clear();
  bool stat;
  bool inside=false;
  while(stat=getline(fp,line)){
    trimComments(line);
    trim(line);
    if(line.length()==0) continue;
    vector<string> w=getWords(line);
    if(w.size()==0) continue;
    if(inside && *(w.begin())=="..."){
      inside=false;
      if(w.size()==2) assert(w[1]==words[0]);
      assert(w.size()<=2);
      w.clear();
    }else if(*(w.end()-1)=="..."){
      inside=true;
      w.erase(w.end()-1);
    };
    for(unsigned i=0;i<w.size();++i) words.push_back(w[i]);
    if(!inside)break;
  }
  if(words.size()>0) return true;
  return stat;
}


bool Tools::getline(FILE* fp,string & line){
  line="";
  const int bufferlength=5;
  char buffer[bufferlength];
  bool ret;
  while(ret=fgets(buffer,bufferlength,fp)){
    line.append(buffer);
    if(buffer[strlen(buffer)-1]=='\n') break;
  };
  if(*(line.end()-1)=='\n') line.erase(line.end()-1);
  return ret;
}

void Tools::trim(string & s){
  size_t n=s.find_last_not_of(" \t");
  s=s.substr(0,n+1);
}

void Tools::trimComments(string & s){
  size_t n=s.find_first_of("#");
  s=s.substr(0,n);
}

bool Tools::getKey(vector<string>& line,const string & key,string & s){
  s.clear();
  for(vector<string>::iterator p=line.begin();p!=line.end();++p){
    if((*p).length()==0) continue;
    string x=(*p).substr(0,key.length());
    if(x==key){
      if((*p).length()==key.length())return false;
      string tmp=(*p).substr(key.length(),(*p).length());
      line.erase(p);
      s=tmp;
      return true;
    }
  };
  return true;
}

void Tools::convert(int i,std::string & str){
        std::ostringstream ostr;
        ostr<<i;
        str=ostr.str();
}



