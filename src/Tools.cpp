#include "Tools.h"
#include "AtomNumber.h"
#include <cstring>
#include <cassert>
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

IMPLEMENT(int)

bool Tools::convert(const string & str,AtomNumber &a){
  int i;
  bool r=convert(str,i);
  a.setSerial(i);
  return r;
}

bool Tools::convert(const string & str,double & t){
        const double pi=3.141592653589793238462643383279502884197169399375105820974944592307;
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

void Tools::interpretRanges(std::vector<std::string>&s){
  vector<string> news;
  for(vector<string>::iterator p=s.begin();p!=s.end();p++){
    vector<string> words;
    words=getWords(*p,"-");
    int a,b;
    if(words.size()==2 && convert(words[0],a) && convert(words[1],b)){
      assert(b>=a);
      for(int i=a;i<=b;i++){
        string ss;
        convert(i,ss);
        news.push_back(ss);
      }
    }else news.push_back(*p);
  }
  s=news;
}

double Tools::switchingFunc(const double distance, const int nn, const int mm, const double r_0, const double d_0, double *dfunc){
  double result=0.;
  const double threshold=pow(0.00001,1./(nn-mm));
  const double rdist = (distance-d_0)/r_0;

  if(rdist<=0.){
     result+=1.;
     *dfunc=0.;
  }else if(rdist>(1.-epsilon) && rdist<(1+epsilon)){
     result+=nn/mm;
     *dfunc=0.5*nn*(nn-mm)/mm;
  }else if(rdist>threshold){
     *dfunc=0.;
  }else{
     double rNdist = pow(rdist, nn-1);
     double rMdist = pow(rdist, mm-1);
     double num = 1.-rNdist*rdist;
     double iden = 1./(1.-rMdist*rdist);
     double func = num*iden;
     result += func;
     *dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist))/(distance*r_0);
  }
  return result;
}
