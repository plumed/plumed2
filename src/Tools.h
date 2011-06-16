#ifndef __PLUMED_Tools_h
#define __PLUMED_Tools_h

#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <limits>

namespace PLMD{

/// Very small non-zero number
const double epsilon(std::numeric_limits<double>::epsilon());

/// Empty class which just contains several (static) tools
class Tools{
public:
/// Split the line in words using separators.
  static std::vector<std::string> getWords(const std::string & line,const char* sep=" \t\n");
  static std::vector<std::string> getWords(const std::string & line,const char* sep,const char* quote);
/// Get a line from the file pointer fp
  static bool getline(std::FILE* fp,std::string & line);
/// Get a parsed line from the file pointer fp.
/// This function already takes care of joining continued lines and splitting the
/// resulting line into an array of words
  static bool getParsedLine(std::FILE* fp,std::vector<std::string> & line);
/// Convert a string to a double, reading it
  static bool convert(const std::string & str,double & t);
  static bool convert(const std::string & str,int & t);
  static bool convert(const std::string & str,std::string & t);
  static void convert(int i,std::string & str);
  static void trim(std::string & s);
  static void trimComments(std::string & s);


  static double pbc(double);

  static bool getKey(std::vector<std::string>& line,const std::string & key,std::string & s);

  template <class T>
  static bool parse(std::vector<std::string>&line,const std::string&key,T&val);
  template <class T>
  static bool parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val);
  static bool parseFlag(std::vector<std::string>&line,const std::string&key,bool&val);

  template <class T>
  static double getArray(void*,int,int,double);
  template <class T>
  static double getArray(void*,void*,void*,int,int,double);
};

template <class T>
bool Tools::parse(std::vector<std::string>&line,const std::string&key,T&val){
  std::string s;
  if(!getKey(line,key+"=",s)) return false;
  if(s.length()>0 && !convert(s,val))return false;
  return true;
}

template <class T>
bool Tools::parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val){
  std::string s;
  if(!getKey(line,key+"=",s)) return false;
  if(s.length()==0) return true;
  val.clear();
  std::vector<std::string> words=getWords(s,",");
  for(unsigned i=0;i<words.size();++i){
    T v;
    if(!convert(words[i],v))return false;
    val.push_back(v);
  }
  return true;
}

inline
bool Tools::parseFlag(std::vector<std::string>&line,const std::string&key,bool&val){
  std::string s;
  for(std::vector<std::string>::iterator p=line.begin();p!=line.end();++p){
    if(key==*p){
      val=true;
      line.erase(p);
      return true;
    }
  }
  return true;
}

template <class T>
double Tools::getArray(void*p,int iat,int idim,double scale){
  return (static_cast<T*>(p))[3*iat+idim]*scale;
}

template <class T>
double Tools::getArray(void*px,void*py,void*pz,int iat,int idim,double scale){
  if(idim==0) return (static_cast<T*>(px))[iat]*scale;
  if(idim==1) return (static_cast<T*>(py))[iat]*scale;
              return (static_cast<T*>(pz))[iat]*scale;
}

inline
double Tools::pbc(double x){
  return x-floor(x+0.5);
}



}

#endif

