#ifndef __PLUMED_Tools_h
#define __PLUMED_Tools_h

#include "AtomNumber.h"
#include "Grid.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <limits>
#include <algorithm>

namespace PLMD{

/// Very small non-zero number
const double epsilon(std::numeric_limits<double>::epsilon());

/// Boltzman constant in kj/K
const double kBoltzmann(0.0083144621);

/// PI
const double pi(3.141592653589793238462643383279502884197169399375105820974944592307);

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
/// Convert a string to a int, reading it
  static bool convert(const std::string & str,int & t);
/// Convert a string to an unsigned int, reading it
  static bool convert(const std::string & str,unsigned & t);
/// Convert a string to a atom number, reading it
  static bool convert(const std::string & str,AtomNumber & t);
/// Convert a string to a string (i.e. copy)
  static bool convert(const std::string & str,std::string & t);
/// Convert an int to a string
  static void convert(int i,std::string & str);
/// Remove trailing blanks
  static void trim(std::string & s);
/// Remove trailing comments
  static void trimComments(std::string & s);
/// Apply pbc for a unitary cell
  static double pbc(double);
/// Retrieve a key from a vector of options.
/// It finds a key starting with "key=" or equal to "key" and copy the
/// part after the = on s. E.g.:
/// line.push_back("aa=xx");
/// getKey(line,"aa",s);
/// will set s="xx"
  static bool getKey(std::vector<std::string>& line,const std::string & key,std::string & s);
/// Find a keyword on the input line, eventually deleting it, and saving its value to val
  template <class T>
  static bool parse(std::vector<std::string>&line,const std::string&key,T&val);
/// Find a keyword on the input line, eventually deleting it, and saving its value to a vector
  template <class T>
  static bool parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val);
/// Find a keyword without arguments on the input line
  static bool parseFlag(std::vector<std::string>&line,const std::string&key,bool&val);
/// Interpret atom ranges
  static void interpretRanges(std::vector<std::string>&);
/// Remove duplicates from a vector of types <T> 
  template <typename T>
  static void removeDuplicates(std::vector<T>& vec);
/// interpret ":" syntax for labels
  static void interpretLabel(std::vector<std::string>&s);
/// read grid from file
  static Grid* readGridFromFile(FILE*,bool,bool,bool);
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

template<typename T>
void Tools::removeDuplicates(std::vector<T>& vec)
{
   std::sort(vec.begin(), vec.end());
   vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
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

inline
double Tools::pbc(double x){
  if(std::numeric_limits<int>::round_style == std::round_toward_zero) {
    const double offset=100.0;
    const double y=x+offset;
    if(y>=0) return y-int(y+0.5);
    else     return y-int(y-0.5);
  } else if(std::numeric_limits<int>::round_style == std::round_to_nearest) {
    return x-int(x);
  } else return x-floor(x+0.5);
}

}

#endif

