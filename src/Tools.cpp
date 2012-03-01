#include "Tools.h"
#include "AtomNumber.h"
#include "PlumedException.h"
#include <cstring>
#include <dirent.h>

using namespace PLMD;
using namespace std;

bool Tools::convert(const string & str,int & t){
        istringstream istr(str.c_str());
        bool ok=istr>>t;
        if(!ok) return false;
        string remaining;
        istr>>remaining;
        return remaining.length()==0;
}

bool Tools::convert(const string & str,unsigned & t){
        istringstream istr(str.c_str());
        bool ok=istr>>t;
        if(!ok) return false;
        string remaining;
        istr>>remaining;
        return remaining.length()==0;
}

bool Tools::convert(const string & str,AtomNumber &a){
  int i;
  bool r=convert(str,i);
  if(r) a.setSerial(i);
  return r;
}

bool Tools::convert(const string & str,double & t){
        if(str=="PI" || str=="+PI"){
          t=pi; return true;
        }else if(str=="2PI" || str=="+2PI"){
          t=2*pi; return true;
        }else if(str=="-PI"){
          t=-pi; return true;
        }else if(str=="-2PI"){
          t=-2*pi; return true;
        }
        istringstream istr(str.c_str());
        bool ok=istr>>t;
        if(!ok) return false;
        string remaining;
        istr>>remaining;
        return remaining.length()==0;
}


bool Tools::convert(const string & str,string & t){
        t=str;
        return true;
}

vector<string> Tools::getWords(const string & line,const char* separators){
  const string sep(separators);
  vector<string> words;
  string word;
  for(unsigned i=0;i<line.length();i++){
    bool found=false;
    for(unsigned j=0;j<sep.length();j++) if(line[i]==sep[j]) found=true;
    if(!found) word.push_back(line[i]);
    if(found && word.length()>0){
      words.push_back(word);
      word.clear();
    }
  }
  if(word.length()>0) words.push_back(word);
  return words;
}

bool Tools::getParsedLine(FILE* fp,vector<string> & words){
  string line("");
  words.clear();
  bool stat;
  bool inside=false;
  while((stat=getline(fp,line))){
    trimComments(line);
    trim(line);
    if(line.length()==0) continue;
    vector<string> w=getWords(line);
    if(w.size()==0) continue;
    if(inside && *(w.begin())=="..."){
      inside=false;
      if(w.size()==2) plumed_massert(w[1]==words[0],"second word in terminating \"...\" lines, if present, should be equal to first word of directive");
      plumed_massert(w.size()<=2,"terminating \"...\" lines cannot consist of more than two words");
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
  const int bufferlength=1024;
  char buffer[bufferlength];
  bool ret;
  for(int i=0;i<bufferlength;i++) buffer[i]='\0';
  while((ret=fgets(buffer,bufferlength,fp))){
    line.append(buffer);
    unsigned ss=strlen(buffer);
    if(ss>0) if(buffer[ss-1]=='\n') break;
  };
  if(line.length()>0) if(*(line.end()-1)=='\n') line.erase(line.end()-1);
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
  return false;
}

void Tools::interpretRanges(std::vector<std::string>&s){
  vector<string> news;
  for(vector<string>::iterator p=s.begin();p!=s.end();p++){
    vector<string> words;
    words=getWords(*p,"-");
    int a,b;
    if(words.size()==2 && convert(words[0],a) && convert(words[1],b)){
      plumed_massert(b>=a,"interpreting range \"" + words[0] + "-" + words[1] +"\", second atom should have higher index than first atom");
      for(int i=a;i<=b;i++){
        string ss;
        convert(i,ss);
        news.push_back(ss);
      }
    }else news.push_back(*p);
  }
  s=news;
}

void Tools::interpretLabel(vector<string>&s){
  if(s.size()<2)return;
  string s0=s[0];
  unsigned l=s0.length();
  if(l<1) return;
  if(s0[l-1]==':'){
    s[0]=s[1];
    s[1]="LABEL="+s0.substr(0,l-1);
  }
}

vector<string> Tools::ls(const string&d){
  DIR*dir;
  struct dirent *ent;
  vector<string> result;
  if ((dir=opendir(d.c_str()))){
    while ((ent = readdir (dir))) if(string(ent->d_name)!="." && string(ent->d_name)!="..") result.push_back(ent->d_name);
    closedir (dir);
  }
  return result;
}
