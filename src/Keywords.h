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
  std::vector<KeyType> types;
  std::vector<std::string> keys;
  std::vector<std::string> documentation;
  std::map<std::string,bool> booldefs; 
  std::map<std::string,std::string> numdefs;
public:
  void add( const std::string, const std::string, const std::string );
  void add( const std::string t, const std::string k, const std::string def, const std::string d );
  void addFlag( const std::string k, const bool def, const std::string d );
  void remove( const std::string );
  void clear();
  unsigned size() const;
  bool exists( const std::string k ) const ;
  bool style( const std::string k, const std::string t ) const ;
  void print( Log& log ) const ;
  void print_html() const ;
  void print_html_item( const unsigned& j ) const;
  bool getLogicalDefault( std::string key, bool& def ) const ;
  bool getDefaultValue( std::string key, std::string& def ) const ;
};

}

#endif
