#ifndef __PLUMED_Keywords_h
#define __PLUMED_Keywords_h
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include "Tools.h"
#include "Log.h"

namespace PLMD{

/// This class lets me pass keyword types easily
class KeyType{
friend class Keyword;
private:
  enum {compulsory,flag,optional,atoms,numbered,nohtml,hidden} style;
public:
  KeyType( const std::string& type );
  bool isCompulsory() const { return (style==compulsory); }
  bool isFlag() const { return (style==flag); }
  bool isOptional() const { return (style==optional); }
  bool isAtomList() const { return (style==atoms); }
  bool isNumbered() const { return (style==numbered); }
  bool isNoHTML() const { return (style==nohtml); }
  bool isHidden() const { return (style==hidden); }
};

/// This class holds the keywords and their documentation
class Keywords{
friend class Action;
private:
/// Whether the keyword is compulsory, optional...
  std::vector<KeyType> types;
/// The names of the keyword
  std::vector<std::string> keys;
/// The documentation for the keywords
  std::vector<std::string> documentation;
/// The default values for the flags (are they on or of)
  std::map<std::string,bool> booldefs; 
/// The default values (if there are default values) for compulsory keywords
  std::map<std::string,std::string> numdefs;
/// Print the documentation for the jth keyword in html
  void print_html_item( const unsigned& j ) const;
/// Print the documentation to the log file (used by PLMD::Action::error)
  void print( Log& log ) const ;
/// find out whether flag key is on or off by default.
  bool getLogicalDefault( std::string key, bool& def ) const ;
/// Get the value of the default for the keyword named key
  bool getDefaultValue( std::string key, std::string& def ) const ;
/// Clear all the keywords
  void clear();
/// Return the number of defined keywords 
  unsigned size() const;
public:
/// Add a new keyword of type t with name k and description d
  void add( const std::string t, const std::string k, const std::string d );
/// Add a new compulsory keyword (t must equal compulsory) with name k, default value def and description d
  void add( const std::string t, const std::string k, const std::string def, const std::string d );
/// Add a falg with name k that is by default on if def is true and off if def is false.  d should provide a description of the flag
  void addFlag( const std::string k, const bool def, const std::string d );
/// Remove the keyword with name k
  void remove( const std::string k );
// /// Clear all the keywords
//   void clear();
// /// Return the number of defined keywords 
//   unsigned size() const;
/// Check if there is a keyword with name k
  bool exists( const std::string k ) const ;
/// Check if the keyword with name k has style t
  bool style( const std::string k, const std::string t ) const ;
// /// Print the documentation to the log file (used by PLMD::Action::error)
//   void print( Log& log ) const ;
/// Print an html version of the documentation
  void print_html() const ;
// /// find out whether flag key is on or off by default.
//   bool getLogicalDefault( std::string key, bool& def ) const ;
// /// Get the value of the default for the keyword named key
//   bool getDefaultValue( std::string key, std::string& def ) const ;
};

}

#endif
