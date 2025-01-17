/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "Tools.h"
#include "AtomNumber.h"
#include "Exception.h"
#include "IFile.h"
#include "lepton/Lepton.h"
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <map>
#if defined(__PLUMED_HAS_CHDIR) || defined(__PLUMED_HAS_GETCWD)
#include <unistd.h>
#endif

#include <iomanip>

namespace PLMD {

template<class T>
bool Tools::convertToAny(const std::string & str,T & t) {
  std::istringstream istr(str.c_str());
  bool ok=static_cast<bool>(istr>>t);
  if(!ok) return false;
  std::string remaining;
  istr>>remaining;
  return remaining.length()==0;
}

bool Tools::convertNoexcept(const std::string & str,int & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,long int & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,long long int & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,unsigned & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,long unsigned & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,long long unsigned & t) {
  return convertToInt(str,t);
}

bool Tools::convertNoexcept(const std::string & str,AtomNumber &a) {
  // Note: AtomNumber's are NOT converted as int, so as to
  // avoid using lepton conversions.
  unsigned i;
  bool r=convertToAny(str,i);
  if(r) a.setSerial(i);
  return r;
}

template<class T>
bool Tools::convertToInt(const std::string & str,T & t) {
  // First try standard conversion
  if(convertToAny(str,t)) return true;
  // Then use lepton
  try {
    double r=lepton::Parser::parse(str).evaluate(lepton::Constants());

    // now sanity checks on the resulting number

    // it should not overflow the requested int type:
    // (see https://stackoverflow.com/a/526092)
    if(r>std::nextafter(std::numeric_limits<T>::max(), 0)) return false;
    if(r<std::nextafter(std::numeric_limits<T>::min(), 0)) return false;

    // do the actual conversion
    auto tmp=static_cast<T>(std::round(r));

    // it should be *very close* to itself if converted back to double
    double diff= r-static_cast<double>(tmp);
    if(diff*diff > 1e-20) return false;
    // this is to accomodate small numerical errors and allow e.g. exp(log(7)) to be integer

    // it should be change if incremented or decremented by one (see https://stackoverflow.com/a/43656140)
    if(r == static_cast<double>(tmp-1)) return false;
    if(r == static_cast<double>(tmp+1)) return false;

    // everything is fine, then store in t
    t=tmp;
    return true;
  } catch(const PLMD::lepton::Exception& exc) {
  }
  return false;
}


template<class T>
bool Tools::convertToReal(const std::string & str,T & t) {
  if(convertToAny(str,t)) return true;
  if(str=="PI" || str=="+PI" || str=="+pi" || str=="pi") {
    t=pi; return true;
  } else if(str=="-PI" || str=="-pi") {
    t=-pi; return true;
  }
  try {
    t=lepton::Parser::parse(str).evaluate(lepton::Constants());
    return true;
  } catch(const PLMD::lepton::Exception& exc) {
  }
  if( str.find("PI")!=std::string::npos ) {
    std::size_t pi_start=str.find_first_of("PI");
    if(str.substr(pi_start)!="PI") return false;
    std::istringstream nstr(str.substr(0,pi_start));
    T ff=0.0; bool ok=static_cast<bool>(nstr>>ff);
    if(!ok) return false;
    t=ff*pi;
    std::string remains; nstr>>remains;
    return remains.length()==0;
  } else if( str.find("pi")!=std::string::npos ) {
    std::size_t pi_start=str.find_first_of("pi");
    if(str.substr(pi_start)!="pi") return false;
    std::istringstream nstr(str.substr(0,pi_start));
    T ff=0.0; bool ok=static_cast<bool>(nstr>>ff);
    if(!ok) return false;
    t=ff*pi;
    std::string remains; nstr>>remains;
    return remains.length()==0;
  } else if(str=="NAN") {
    t=std::numeric_limits<double>::quiet_NaN();
    return true;
  }
  return false;
}

bool Tools::convertNoexcept(const std::string & str,float & t) {
  return convertToReal(str,t);
}

bool Tools::convertNoexcept(const std::string & str,double & t) {
  return convertToReal(str,t);
}

bool Tools::convertNoexcept(const std::string & str,long double & t) {
  return convertToReal(str,t);
}

bool Tools::convertNoexcept(const std::string & str,std::string & t) {
  t=str;
  return true;
}

std::vector<std::string> Tools::getWords(const std::string & line,const char* separators,int * parlevel,const char* parenthesis, const bool& delete_parenthesis) {
  plumed_massert(std::strlen(parenthesis)==1,"multiple parenthesis type not available");
  plumed_massert(parenthesis[0]=='(' || parenthesis[0]=='[' || parenthesis[0]=='{',
                 "only ( [ { allowed as parenthesis");
  if(!separators) separators=" \t\n";
  const std::string sep(separators);
  char openpar=parenthesis[0];
  char closepar;
  if(openpar=='(') closepar=')';
  if(openpar=='[') closepar=']';
  if(openpar=='{') closepar='}';
  std::vector<std::string> words;
  std::string word;
  int parenthesisLevel=0;
  if(parlevel) parenthesisLevel=*parlevel;
  for(unsigned i=0; i<line.length(); i++) {
    bool found=false;
    bool onParenthesis=false;
    if( (line[i]==openpar || line[i]==closepar) && delete_parenthesis ) onParenthesis=true;
    if(line[i]==closepar) {
      parenthesisLevel--;
      plumed_massert(parenthesisLevel>=0,"Extra closed parenthesis in '" + line + "'");
    }
    if(parenthesisLevel==0) for(unsigned j=0; j<sep.length(); j++) if(line[i]==sep[j]) found=true;
// If at parenthesis level zero (outer)
    if(!(parenthesisLevel==0 && (found||onParenthesis))) word.push_back(line[i]);
    //if(onParenthesis) word.push_back(' ');
    if(line[i]==openpar) parenthesisLevel++;
    if(found && word.length()>0) {
      if(!parlevel) plumed_massert(parenthesisLevel==0,"Unmatching parenthesis in '" + line + "'");
      words.push_back(word);
      word.clear();
    }
  }
  if(word.length()>0) {
    if(!parlevel) plumed_massert(parenthesisLevel==0,"Unmatching parenthesis in '" + line + "'");
    words.push_back(word);
  }
  if(parlevel) *parlevel=parenthesisLevel;
  return words;
}

bool Tools::getParsedLine(IFile& ifile,std::vector<std::string> & words, bool trimcomments) {
  std::string line("");
  words.clear();
  bool stat;
  bool inside=false;
  int parlevel=0;
  bool mergenext=false;
  while((stat=ifile.getline(line))) {
    if(trimcomments) trimComments(line);
    trim(line);
    if(line.length()==0) continue;
    std::vector<std::string> w=getWords(line,NULL,&parlevel,"{",trimcomments);
    if(!w.empty()) {
      if(inside && *(w.begin())=="...") {
        inside=false;
        if(w.size()==2) plumed_massert(w[1]==words[0],"second word in terminating \"...\" "+w[1]+" line, if present, should be equal to first word of directive: "+words[0]);
        plumed_massert(w.size()<=2,"terminating \"...\" lines cannot consist of more than two words");
        w.clear(); if(!trimcomments) words.push_back("...");
      } else if(*(w.end()-1)=="...") {
        inside=true;
        w.erase(w.end()-1);
      };
      int i0=0;
      if(mergenext && words.size()>0 && w.size()>0) {
        words[words.size()-1]+=" "+w[0];
        i0=1;
      }
      for(unsigned i=i0; i<w.size(); ++i) words.push_back(w[i]);
    }
    mergenext=(parlevel>0);
    if(!inside)break;
    if(!trimcomments && parlevel==0) words.push_back("@newline");
    else if(!trimcomments) words[words.size()-1] += " @newline";
  }
  plumed_massert(parlevel==0,"non matching parenthesis");
  if(words.size()>0) return true;
  return stat;
}


bool Tools::getline(FILE* fp,std::string & line) {
  line="";
  const int bufferlength=1024;
  char buffer[bufferlength];
  bool ret;
  for(int i=0; i<bufferlength; i++) buffer[i]='\0';
  while((ret=fgets(buffer,bufferlength,fp))) {
    line.append(buffer);
    unsigned ss=std::strlen(buffer);
    if(ss>0) if(buffer[ss-1]=='\n') break;
  };
  if(line.length()>0) if(*(line.end()-1)=='\n') line.erase(line.end()-1);
  if(line.length()>0) if(*(line.end()-1)=='\r') line.erase(line.end()-1);
  return ret;
}

void Tools::trim(std::string & s) {
  auto n=s.find_last_not_of(" \t");
  if(n!=std::string::npos) s.resize(n+1);
}

void Tools::trimComments(std::string & s) {
  auto n=s.find_first_of("#");
  if(n!=std::string::npos) s.resize(n);
}

bool Tools::caseInSensStringCompare(const std::string & str1, const std::string &str2)
{
  return ((str1.size() == str2.size()) && std::equal(str1.begin(), str1.end(), str2.begin(), [](char c1, char c2) {
    return (c1 == c2 || std::toupper(c1) == std::toupper(c2));
  }));
}

bool Tools::getKey(std::vector<std::string>& line,const std::string & key,std::string & s,int rep) {
  s.clear();
  for(auto p=line.begin(); p!=line.end(); ++p) {
    if((*p).length()==0) continue;
    std::string x=(*p).substr(0,key.length());
    if(caseInSensStringCompare(x,key)) {
      if((*p).length()==key.length())return false;
      std::string tmp=(*p).substr(key.length(),(*p).length());
      line.erase(p);
      s=tmp;
      const std::string multi("@replicas:");
      if(rep>=0 && startWith(s,multi)) {
        s=s.substr(multi.length(),s.length());
        std::vector<std::string> words=getWords(s,"\t\n ,");
        plumed_massert(rep<static_cast<int>(words.size()),"Number of fields in " + s + " not consistent with number of replicas");
        s=words[rep];
      }
      return true;
    }
  };
  return false;
}

void Tools::interpretRanges(std::vector<std::string>&s) {
  std::vector<std::string> news;
  for(const auto & p :s) {
    news.push_back(p);
    size_t dash=p.find("-");
    if(dash==std::string::npos) continue;
    int first;
    if(!Tools::convertToAny(p.substr(0,dash),first)) continue;
    int stride=1;
    int second;
    size_t colon=p.substr(dash+1).find(":");
    if(colon!=std::string::npos) {
      if(!Tools::convertToAny(p.substr(dash+1).substr(0,colon),second) ||
          !Tools::convertToAny(p.substr(dash+1).substr(colon+1),stride)) continue;
    } else {
      if(!Tools::convertToAny(p.substr(dash+1),second)) continue;
    }
    news.resize(news.size()-1);
    if(first<=second) {
      plumed_massert(stride>0,"interpreting ranges "+ p + ", stride should be positive");
      for(int i=first; i<=second; i+=stride) {
        std::string ss;
        convert(i,ss);
        news.push_back(ss);
      }
    } else {
      plumed_massert(stride<0,"interpreting ranges "+ p + ", stride should be positive");
      for(int i=first; i>=second; i+=stride) {
        std::string ss;
        convert(i,ss);
        news.push_back(ss);
      }
    }
  }
  s=news;
}

void Tools::interpretLabel(std::vector<std::string>&s) {
  if(s.size()<2)return;
  std::string s0=s[0];
  unsigned l=s0.length();
  if(l<1) return;
  if(s0[l-1]==':') {
    s[0]=s[1];
    s[1]="LABEL="+s0.substr(0,l-1);
  }
  std::transform(s[0].begin(), s[0].end(), s[0].begin(), ::toupper);
}

std::vector<std::string> Tools::ls(const std::string&d) {
  DIR*dir;
  std::vector<std::string> result;
  if ((dir=opendir(d.c_str()))) {
#if defined(__PLUMED_HAS_READDIR_R)
    struct dirent ent;
#endif
    while(true) {
      struct dirent *res;
#if defined(__PLUMED_HAS_READDIR_R)
      readdir_r(dir,&ent,&res);
#else
      res=readdir(dir);
#endif
      if(!res) break;
      if(std::string(res->d_name)!="." && std::string(res->d_name)!="..") result.push_back(res->d_name);
    }
    closedir (dir);
  }
  return result;
}

void Tools::stripLeadingAndTrailingBlanks( std::string& str ) {
  std::size_t first=str.find_first_not_of(' ');
  std::size_t last=str.find_last_not_of(' ');
  if( first<=last && first!=std::string::npos) str=str.substr(first,last+1);
}

std::string Tools::extension(const std::string&s) {
  size_t n=s.find_last_of(".");
  std::string ext;
  if(n!=std::string::npos && n+1<s.length() && n+5>=s.length()) {
    ext=s.substr(n+1);
    if(ext.find("/")!=std::string::npos) ext="";
    std::string base=s.substr(0,n);
    if(base.length()==0) ext="";
    if(base.length()>0 && base[base.length()-1]=='/') ext="";
  }
  return ext;
}

double Tools::bessel0( const double& val ) {
  if (std::abs(val)<3.75) {
    double y = Tools::fastpow( val/3.75, 2 );
    return 1 + y*(3.5156229 +y*(3.0899424 + y*(1.2067492+y*(0.2659732+y*(0.0360768+y*0.0045813)))));
  }
  double ax=std::abs(val), y=3.75/ax, bx=std::exp(ax)/std::sqrt(ax);
  ax=0.39894228+y*(0.01328592+y*(0.00225319+y*(-0.00157565+y*(0.00916281+y*(-0.02057706+y*(0.02635537+y*(-0.01647633+y*0.00392377)))))));
  return ax*bx;
}

bool Tools::startWith(const std::string & full,const std::string &start) {
  return (full.substr(0,start.length())==start);
}

bool Tools::findKeyword(const std::vector<std::string>&line,const std::string&key) {
  const std::string search(key+"=");
  for(const auto & p : line) {
    if(startWith(p,search)) return true;
  }
  return false;
}

Tools::DirectoryChanger::DirectoryChanger(const char*path) {
  if(!path) return;
  if(std::strlen(path)==0) return;
#ifdef __PLUMED_HAS_GETCWD
  char* ret=getcwd(cwd,buffersize);
  plumed_assert(ret)<<"Name of current directory too long, increase buffer size";
#else
  plumed_error()<<"You are trying to use DirectoryChanger but your system does not support getcwd";
#endif
#ifdef __PLUMED_HAS_CHDIR
  int r=chdir(path);
  plumed_assert(r==0) <<"Cannot chdir to directory "<<path<<". The directory must exist!";
#else
  plumed_error()<<"You are trying to use DirectoryChanger but your system does not support chdir";
#endif
}

Tools::DirectoryChanger::~DirectoryChanger() {
#ifdef __PLUMED_HAS_CHDIR
  if(std::strlen(cwd)==0) return;
  int ret=chdir(cwd);
// we cannot put an assertion here (in a destructor) otherwise cppcheck complains
// we thus just report the problem
  if(ret!=0) std::fprintf(stderr,"+++ WARNING: cannot cd back to directory %s\n",cwd);
#endif
}

std::unique_ptr<std::lock_guard<std::mutex>> Tools::molfile_lock() {
  static std::mutex mtx;
  return Tools::make_unique<std::lock_guard<std::mutex>>(mtx);
}

/// Internal tool, I am keeping it private for now
namespace {

class process_one_exception {
  std::string & msg;
  bool first=true;
  void update() {
    if(!first) msg+="\n\nThe above exception was the direct cause of the following exception:\n";
    first=false;
  }
public:
  process_one_exception(std::string & msg):
    msg(msg)
  {}
  void operator()(const std::exception & e) {
    update();
    msg+=e.what();
  }
  void operator()(const std::string & e) {
    update();
    msg+=e;
  }
  void operator()(const char* e) {
    update();
    msg+=e;
  }
};

template<class T>
static void process_all_exceptions(T&& f) {
  try {
    // First throw the current exception
    throw;
  } catch(const std::nested_exception & e) {
    // If nested, we go recursive
    // notice that we apply function f only if exception is also a std::exception
    try {
      e.rethrow_nested();
    } catch(...) {
      process_all_exceptions(f);
    }
    auto d=dynamic_cast<const std::exception*>(&e);
    if(d) f(*d);
  } catch(const std::exception &e) {
    // If not nested, we end recursion
    f(e);
  } catch(const std::string &e) {
    // If not nested, we end recursion
    f(e);
  } catch(const char* e) {
    // If not nested, we end recursion
    f(e);
  } catch(...) {
    // If not nested and of unknown type, we stop the chain
  }
}

}

std::string Tools::concatenateExceptionMessages() {
  std::string msg;
  process_all_exceptions(process_one_exception(msg));
  return msg;
}

}
