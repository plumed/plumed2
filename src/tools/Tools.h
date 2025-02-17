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
#ifndef __PLUMED_tools_Tools_h
#define __PLUMED_tools_Tools_h

#include "AtomNumber.h"
#include "Vector.h"
#include "Tensor.h"
#include "small_vector/small_vector.h"
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
#include <queue>
#include <mutex>
#include <filesystem>
#include <utility>
#include <unordered_map>
#include <map>
#include <condition_variable>

namespace PLMD {

class IFile;

/// \ingroup TOOLBOX
/// Very small non-zero number
constexpr double epsilon(std::numeric_limits<double>::epsilon());

/// \ingroup TOOLBOX
/// Boltzman constant in kj/K
constexpr double kBoltzmann(0.0083144621);

/// \ingroup TOOLBOX
/// PI
constexpr double pi(3.141592653589793238462643383279502884197169399375105820974944592307);

constexpr double dp2cutoff(6.25);

constexpr double dp2cutoffA=1.00193418799744762399; // 1.0/(1-std::exp(-dp2cutoff));
constexpr double dp2cutoffB=-.00193418799744762399; // -std::exp(-dp2cutoff)/(1-std::exp(-dp2cutoff));

inline static bool dp2cutoffNoStretch() {
  static const auto* res=std::getenv("PLUMED_DP2CUTOFF_NOSTRETCH");
  return res;
}

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
/// @brief the  recursive part of the template fastpow implementation
  template <int exp, typename T=double, std::enable_if_t< (exp >=0), bool> = true>
  static inline /*consteval*/ T fastpow_rec(T base, T result);
public:
/// Split the line in words using separators.
/// It also take into account parenthesis. Outer parenthesis found are removed from
/// output, and the text between them is considered as a single word. Only the
/// outer parenthesis are processed, to allow nesting them.
/// parlevel, if not NULL, is increased or decreased according to the number of opened/closed parenthesis
  static std::vector<std::string> getWords(std::string_view line,const char* sep=NULL,int* parlevel=NULL,const char* parenthesis="{", const bool& delete_parenthesis=true);
/// Faster version
/// This version does not parse parenthesis and operates on a preallocated small_vector of string_view's
  static void getWordsSimple(gch::small_vector<std::string_view> & words,std::string_view line);
/// Get a line from the file pointer ifile
  static bool getline(FILE*,std::string & line);
/// Get a parsed line from the file pointer ifile
/// This function already takes care of joining continued lines and splitting the
/// resulting line into an array of words
  static bool getParsedLine(IFile&ifile,std::vector<std::string> & line, const bool trimcomments=true);
/// compare two string in a case insensitive manner
  static bool caseInSensStringCompare(const std::string & str1, const std::string &str2);
/// Convert a string to a double, reading it
  static bool convertNoexcept(const std::string & str,double & t);
/// Convert a string to a long double, reading it
  static bool convertNoexcept(const std::string & str,long double & t);
/// Convert a string to a float, reading it
  static bool convertNoexcept(const std::string & str,float & t);
/// Convert a string to a int, reading it
  static bool convertNoexcept(const std::string & str,int & t);
/// Convert a string to a long int, reading it
  static bool convertNoexcept(const std::string & str,long int & t);
/// Convert a string to a long long int, reading it
  static bool convertNoexcept(const std::string & str,long long int & t);
/// Convert a string to an unsigned int, reading it
  static bool convertNoexcept(const std::string & str,unsigned & t);
/// Convert a string to a long unsigned int, reading it
  static bool convertNoexcept(const std::string & str,long unsigned & t);
/// Convert a string to a long long unsigned int, reading it
  static bool convertNoexcept(const std::string & str,long long unsigned & t);
/// Convert a string to a atom number, reading it
  static bool convertNoexcept(const std::string & str,AtomNumber & t);
/// Convert a string to a string (i.e. copy)
  static bool convertNoexcept(const std::string & str,std::string & t);
/// Convert anything into a string
  template<typename T>
  static bool convertNoexcept(T i,std::string & str);
/// Convert anything into anything, throwing an exception in case there is an error
/// Remove trailing blanks
  static void trim(std::string & s);
/// Remove leading blanks
  static void ltrim(std::string & s);
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
  template <typename T,typename U>
  static void convert(const T & t,U & u) {
    plumed_assert(convertNoexcept(t,u)) <<"Error converting  "<<t;
  }
  template <typename T>
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
/// Fast int power for power known at compile time
  template <int exp, typename T=double>
  static inline /*consteval*/ T fastpow(T base);
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
    const std::filesystem::path path;
  public:
    explicit DirectoryChanger(const char*path);
    ~DirectoryChanger();
  };

  template<class T, class... Args>
  static auto make_unique(Args&&... args) {
    return std::make_unique<T>(std::forward<Args>(args)...);
  }

  static void set_to_zero(double*ptr,unsigned n) {
    for(unsigned i=0; i<n; i++) {
      ptr[i]=0.0;
    }
  }

  template<unsigned n>
  static void set_to_zero(std::vector<VectorGeneric<n>> & vec) {
    unsigned s=vec.size();
    if(s==0) {
      return;
    }
    set_to_zero(&vec[0][0],s*n);
  }

  template<unsigned n,unsigned m>
  static void set_to_zero(std::vector<TensorGeneric<n,m>> & vec) {
    unsigned s=vec.size();
    if(s==0) {
      return;
    }
    set_to_zero(&vec[0](0,0),s*n*m);
  }

  static std::unique_ptr<std::lock_guard<std::mutex>> molfile_lock();
  /// Build a concatenated exception message.
  /// Should be called with an in-flight exception.
  static std::string concatenateExceptionMessages();


  /// Tiny class implementing faster std::string_view access to an unordered_map
  /// It exposes a limited number of methods of std::unordered_map. Others could be added.
  /// Importantly, when it is accessed via a std::string_view, the access does not
  /// require constructing a std::string and is thus faster.
  /// Deletion would be slower instead. It's not even implemented yet.
  template<class T>
  class FastStringUnorderedMap {
    std::unordered_map<std::string_view,T> map;
    std::vector<std::unique_ptr<const char[]>> keys;

    // see https://stackoverflow.com/questions/34596768/stdunordered-mapfind-using-a-type-different-than-the-key-type
    std::unique_ptr<const char[]> conv(std::string_view str) {
      auto p=std::make_unique<char[]>(str.size()+1);
      std::memcpy(p.get(), str.data(), str.size()+1);
      return p;
    }

  public:

    FastStringUnorderedMap() = default;
    FastStringUnorderedMap(std::initializer_list<std::pair<const std::string_view,T>> init) {
      for(const auto & c : init) {
        (*this)[c.first]=c.second;
      }
    }

    T& operator[]( const std::string_view & key ) {
      auto f=map.find(key);
      if(f!=map.end()) {
        return f->second;
      }
      keys.push_back(conv(key));
      return map[keys.back().get()];
    }

    auto begin() {
      return map.begin();
    }
    auto end() {
      return map.end();
    }
    auto begin() const {
      return map.begin();
    }
    auto end() const {
      return map.end();
    }
    auto find(const std::string_view & key) {
      return map.find(key);
    }
    auto find(const std::string_view & key) const {
      return map.find(key);
    }
  };

  /// Utility to create named critical sections
  /// Key should be usable in a std::map
  template<class Key>
  class CriticalSectionWithKey {
    std::mutex mutex;
    std::condition_variable notify;
    std::map<Key, int> in_progress;
  public:
    void start(const Key & key) {
      std::unique_lock<std::mutex> lock(mutex);
      while (in_progress[key] > 0) {
        // Wait if this command is already in progress.
        notify.wait(lock);
      }
      // Mark this command as in progress.
      in_progress[key]++;
    }
    void stop(const Key & key) {
      std::unique_lock<std::mutex> lock(mutex);
      // Mark this command as completed.
      in_progress[key]--;
      // Notify other threads that may be waiting for this command to complete.
      notify.notify_all();
    }
    class Handler {
      CriticalSectionWithKey* section{nullptr};
      Key key;
      Handler(CriticalSectionWithKey* section,const Key& key):
        section(section),
        key(key) {
        section->start(key);
      }
      friend class CriticalSectionWithKey;
    public:
      /// Default constructor
      Handler() = default;
      /// Default copy constructor is deleted (not copyable)
      Handler(const Handler & handler) = delete;
      /// Default copy assignment is deleted (not copyable)
      Handler & operator=(const Handler & handler) = delete;
      /// Move constructor.
      Handler(Handler && handler) noexcept :
        section(handler.section),
        key(std::move(handler.key)) {
        handler.section=nullptr;
      };
      /// Move assignment.
      Handler & operator=(Handler && handler) noexcept {
        if(this!=&handler) {
          if(section) {
            section->stop(key);
          }
          section=handler.section;
          key=std::move(handler.key);
        }
        handler.section=nullptr;
        return *this;
      }
      /// Destructor
      ~Handler() {
        if(section) {
          section->stop(key);
        }
      }
    };

    Handler startStop(const Key & key) {
      return Handler(this,key);
    }

  };

};

template <class T>
bool Tools::parse(std::vector<std::string>&line,const std::string&key,T&val,int rep) {
  std::string s;
  if(!getKey(line,key+"=",s,rep)) {
    return false;
  }
  if(s.length()>0 && !convertNoexcept(s,val)) {
    return false;
  }
  return true;
}

template <class T>
bool Tools::parseVector(std::vector<std::string>&line,const std::string&key,std::vector<T>&val,int rep) {
  std::string s;
  if(!getKey(line,key+"=",s,rep)) {
    return false;
  }
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
    if(!convertNoexcept(s,v)) {
      return false;
    }
    val.push_back(v);
  }
  return true;
}

template<typename T>
void Tools::removeDuplicates(std::vector<T>& vec) {
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
  while (x>0.5) {
    x-=1.0;
  }
  while (x<-0.5) {
    x+=1.0;
  }
  return x;
#else
  if constexpr (std::numeric_limits<int>::round_style == std::round_toward_zero) {
    constexpr double offset=100.0;
    const double y=x+offset;
    if(y>=0) {
      return y-int(y+0.5);
    } else {
      return y-int(y-0.5);
    }
  } else if constexpr (std::numeric_limits<int>::round_style == std::round_to_nearest) {
    return x-int(x);
  } else {
    return x-floor(x+0.5);
  }
#endif
}

template<typename T>
bool Tools::convertNoexcept(T i,std::string & str) {
  std::ostringstream ostr;
  ostr<<i;
  str=ostr.str();
  return true;
}

inline
double Tools::fastpow(double base, int exp) {
  if(exp<0) {
    exp=-exp;
    base=1.0/base;
  }
  double result = 1.0;
  while (exp) {
    if (exp & 1) {
      result *= base;
    }
    exp >>= 1;
    base *= base;
  }

  return result;
}

template <int exp, typename T, std::enable_if_t< (exp >=0), bool>>
inline T Tools::fastpow_rec(T const base, T result) {
  if constexpr (exp == 0) {
    return result;
  }
  if constexpr (exp & 1) {
    result *= base;
  }
  return fastpow_rec<(exp>>1),T> (base*base, result);
}

template <int exp, typename T>
inline T Tools::fastpow(T const base) {
  if constexpr (exp<0) {
    return  fastpow_rec<-exp,T>(1.0/base,1.0);
  } else {
    return fastpow_rec<exp,T>(base, 1.0);
  }
}

template<typename T>
std::vector<T*> Tools::unique2raw(const std::vector<std::unique_ptr<T>> & x) {
  std::vector<T*> v(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    v[i]=x[i].get();
  }
  return v;
}

template<typename T>
std::vector<const T*> Tools::unique2raw(const std::vector<std::unique_ptr<const T>> & x) {
  std::vector<const T*> v(x.size());
  for(unsigned i=0; i<x.size(); i++) {
    v[i]=x[i].get();
  }
  return v;
}

}

#endif
