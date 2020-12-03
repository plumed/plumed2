/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include <cctype>
#include <cstdio>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>
#include <memory>
#include <cstddef>

namespace PLMD {

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
class Tools {
/// class to convert a string to a generic type T
  template<class T>
  static bool convertToAny(const std::string & str,T &t);
/// class to convert a string to a real type T.
/// T should be either float, double, or long double
  template<class T>
  static bool convertToReal(const std::string & str,T &t);
/// class to convert a string to a int type T
  template<class T>
  static bool convertToInt(const std::string & str,T &t);
public:
/// Split the line in words using separators.
/// It also take into account parenthesis. Outer parenthesis found are removed from
/// output, and the text between them is considered as a single word. Only the
/// outer parenthesis are processed, to allow nesting them.
/// parlevel, if not NULL, is increased or decreased according to the number of opened/closed parenthesis
  static std::vector<std::string> getWords(const std::string & line,const char* sep=NULL,int* parlevel=NULL,const char* parenthesis="{", const bool& delete_parenthesis=true);
/// Get a line from the file pointer ifile
  static bool getline(FILE*,std::string & line);
/// Get a parsed line from the file pointer ifile
/// This function already takes care of joining continued lines and splitting the
/// resulting line into an array of words
  static bool getParsedLine(IFile&ifile,std::vector<std::string> & line, const bool trimcomments=true);
/// compare two string in a case insensitive manner
  static bool caseInSensStringCompare(const std::string & str1, const std::string &str2);
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
/// Convert a string to a long unsigned int, reading it
  static bool convert(const std::string & str,long unsigned & t);
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
  static bool getKey(std::vector<std::string>& line,const std::string & key,std::string & s,int rep=-1);
/// Find a keyword on the input line, eventually deleting it, and saving its value to val
  template <class T>
  static bool parse(std::vector<std::string>&line,const std::string&key,T&val,int rep=-1);
/// Find a keyword on the input line, eventually deleting it, and saving its value to a vector
  template <class T>
  static bool parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val,int rep=-1);
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
/// Modified 0th-order Bessel function of the first kind
  static double bessel0(const double& val);
/// Check if a string full starts with string start.
/// Same as full.find(start)==0
  static bool startWith(const std::string & full,const std::string &start);
  /**
    Tool to create a vector of raw pointers from a vector of unique_pointers (const version).
  Returning a vector is fast in C++11. It can be used in order to feed a vector<unique_ptr<T>>
  to a function that takes a vector<T*>.
  \verbatim
  // some function that takes a vec
  void func(std::vector<Data*> & vec);
  std::vector<std::unique_ptr<Data>> vec;
  // func(vec); // does not compile
  func(Tools::unique2raw(vec)); // compiles
  \endverbatim
  Notice that the conversion is fast but takes
  some time to allocate the new vector and copy the pointers. In case the function
  acting on the vector<T*> is very fast and we do not want to add significant overhead,
  it might be convenient to store a separate set of raw pointers.
  \verbatim
  // some function that takes a vec
  void func(std::vector<Data*> & vec);
  std::vector<std::unique_ptr<Data>> vec;

  // conversion done only once:
  auto vec_ptr=Tools::unique2raw(vec);

  for(int i=0;i<1000;i++){
    func(vec_ptr);
  }
  \endverbatim
  */
  template <typename T>
  static std::vector<T*> unique2raw(const std::vector<std::unique_ptr<T>>&);
/// Tool to create a vector of raw pointers from a vector of unique_pointers.
/// See the non const version.
  template <typename T>
  static std::vector<const T*> unique2raw(const std::vector<std::unique_ptr<const T>>&);
/// Tiny class that changes directory and comes back when going out of scope.
/// In case system calls to change dir are not available it throws an exception.
/// \warning By construction, changing directory breaks thread safety! Use with care.
  class DirectoryChanger {
    static const std::size_t buffersize=4096;
    char cwd[buffersize]= {0};
  public:
    explicit DirectoryChanger(const char*path);
    ~DirectoryChanger();
  };
/// Mimic C++14 std::make_unique
  template<class T> struct _Unique_if {
    typedef std::unique_ptr<T> _Single_object;
  };
  template<class T> struct _Unique_if<T[]> {
    typedef std::unique_ptr<T[]> _Unknown_bound;
  };
  template<class T, std::size_t N> struct _Unique_if<T[N]> {
    typedef void _Known_bound;
  };
  template<class T, class... Args>
  static typename _Unique_if<T>::_Single_object
  make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
  template<class T>
  static typename _Unique_if<T>::_Unknown_bound
  make_unique(std::size_t n) {
    typedef typename std::remove_extent<T>::type U;
    return std::unique_ptr<T>(new U[n]());
  }
  template<class T, class... Args>
  static typename _Unique_if<T>::_Known_bound
  make_unique(Args&&...) = delete;
};

template <class T>
bool Tools::parse(std::vector<std::string>&line,const std::string&key,T&val,int rep) {
  std::string s;
  if(!getKey(line,key+"=",s,rep)) return false;
  if(s.length()>0 && !convert(s,val))return false;
  return true;
}

template <class T>
bool Tools::parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val,int rep) {
  std::string s;
  if(!getKey(line,key+"=",s,rep)) return false;
  val.clear();
  std::vector<std::string> words=getWords(s,"\t\n ,");
  for(unsigned i=0; i<words.size(); ++i) {
    T v;
    std::string s=words[i];
    const std::string multi("@replicas:");
    if(rep>=0 && startWith(s,multi)) {
      s=s.substr(multi.length(),s.length());
      std::vector<std::string> words=getWords(s,"\t\n ,");
      plumed_assert(rep<static_cast<int>(words.size()));
      s=words[rep];
    }
    if(!convert(s,v))return false;
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
bool Tools::parseFlag(std::vector<std::string>&line,const std::string&key,bool&val) {
  for(auto p=line.begin(); p!=line.end(); ++p) {
    if(key==*p) {
      val=true;
      line.erase(p);
      return true;
    }
  }
  return false;
}

/// beware: this brings any number into a pbc that ranges from -0.5 to 0.5
inline
double Tools::pbc(double x) {
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
void Tools::convert(T i,std::string & str) {
  std::ostringstream ostr;
  ostr<<i;
  str=ostr.str();
}

inline
double Tools::fastpow(double base, int exp)
{
  if(exp<0) {
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

template<typename T>
std::vector<T*> Tools::unique2raw(const std::vector<std::unique_ptr<T>> & x) {
  std::vector<T*> v(x.size());
  for(unsigned i=0; i<x.size(); i++) v[i]=x[i].get();
  return v;
}

template<typename T>
std::vector<const T*> Tools::unique2raw(const std::vector<std::unique_ptr<const T>> & x) {
  std::vector<const T*> v(x.size());
  for(unsigned i=0; i<x.size(); i++) v[i]=x[i].get();
  return v;
}

}

#endif

