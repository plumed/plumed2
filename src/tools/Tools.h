/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_tools_Tools_h
#define __PLUMED_tools_Tools_h

#include "AtomNumber.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

namespace PLMD{

class IFile;

/// \ingroup TOOLBOX
/// Very small non-zero number
const double epsilon(std::numeric_limits<double>::epsilon());

/// \ingroup TOOLBOX
/// Boltzman constant in kj/K
const double kBoltzmann(0.0083144621);

/// \ingroup TOOLBOX
/// PI
const double pi(3.141592653589793238462643383279502884197169399375105820974944592307);

/// \ingroup TOOLBOX
/// Empty class which just contains several (static) tools
class Tools{
/// class to convert a string to a generic type T
  template<class T>
  static bool convertToAny(const std::string & str,T &t);
/// class to convert a string to a real type T.
/// T should be either float, double, or long double
  template<class T>
  static bool convertToReal(const std::string & str,T &t);
public:
/// Split the line in words using separators.
/// It also take into account parenthesis. Outer parenthesis found are removed from
/// output, and the text between them is considered as a single word. Only the
/// outer parenthesis are processed, to allow nesting them.
/// parlevel, if not NULL, is increased or decreased according to the number of opened/closed parenthesis
  static std::vector<std::string> getWords(const std::string & line,const char* sep=NULL,int* parlevel=NULL,const char* parenthesis="{");
/// Get a line from the file pointer ifile
  static bool getline(FILE*,std::string & line);
/// Get a parsed line from the file pointer ifile
/// This function already takes care of joining continued lines and splitting the
/// resulting line into an array of words
  static bool getParsedLine(IFile&ifile,std::vector<std::string> & line);
/// Convert a string to a double, reading it
  static bool convert(const std::string & str,double & t);
/// Convert a string to a long double, reading it
  static bool convert(const std::string & str,long double & t);
/// Convert a string to a float, reading it
  static bool convert(const std::string & str,float & t);
/// Convert a string to a int, reading it
  static bool convert(const std::string & str,int & t);
/// Convert a string to a long int, reading it
  static bool convert(const std::string & str,long int & t);
/// Convert a string to an unsigned int, reading it
  static bool convert(const std::string & str,unsigned & t);
/// Convert a string to a atom number, reading it
  static bool convert(const std::string & str,AtomNumber & t);
/// Convert a string to a string (i.e. copy)
  static bool convert(const std::string & str,std::string & t);
/// Convert anything into a string
  template<typename T>
  static void convert(T i,std::string & str);
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
/// Find a keyword on the input line, just reporting if it exists or not
  static bool findKeyword(const std::vector<std::string>&line,const std::string&key);
/// Interpret atom ranges
  static void interpretRanges(std::vector<std::string>&);
/// Remove duplicates from a vector of type T 
  template <typename T>
  static void removeDuplicates(std::vector<T>& vec);
/// interpret ":" syntax for labels
  static void interpretLabel(std::vector<std::string>&s);
/// list files in a directory
  static std::vector<std::string> ls(const std::string&);
/// removes leading and trailing blanks from a string
  static void stripLeadingAndTrailingBlanks( std::string& str );
/// Extract the extensions from a file name.
/// E.g.: extension("pippo.xyz")="xyz".
/// It only returns extensions with a length between 1 and 4
/// E.g.: extension("pippo.12345")="" whereas extenion("pippo.1234")="1234";
/// It is also smart enough to detect "/", so that 
/// extension("pippo/.t")="" whereas extension("pippo/a.t")="t"
  static std::string extension(const std::string&);
/// Fast int power
  static double fastpow(double base,int exp);
/// Check if a string full starts with string start.
/// Same as full.find(start)==0
  static bool startWith(const std::string & full,const std::string &start);
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
//  if(s.length()==0) return true;
  val.clear();
  std::vector<std::string> words=getWords(s,"\t\n ,");
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
  for(std::vector<std::string>::iterator p=line.begin();p!=line.end();++p){
    if(key==*p){
      val=true;
      line.erase(p);
      return true;
    }
  }
  return false;
}

/// beware: this brings any number into a pbc that ranges from -0.5 to 0.5
inline
double Tools::pbc(double x){
#ifdef __PLUMED_PBC_WHILE
  while (x>0.5) x-=1.0;
  while (x<-0.5) x+=1.0;
  return x;
#else
  if(std::numeric_limits<int>::round_style == std::round_toward_zero) {
    const double offset=100.0;
    const double y=x+offset;
    if(y>=0) return y-int(y+0.5);
    else     return y-int(y-0.5);
  } else if(std::numeric_limits<int>::round_style == std::round_to_nearest) {
    return x-int(x);
  } else return x-floor(x+0.5);
#endif
}

template<typename T>
void Tools::convert(T i,std::string & str){
        std::ostringstream ostr;
        ostr<<i;
        str=ostr.str();
}

inline
double Tools::fastpow(double base, int exp)
{
    if(exp<0){
      exp=-exp;
      base=1.0/base;
    }
    double result = 1.0;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}

}

#endif

