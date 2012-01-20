#ifndef __PLUMED_Keywords_h
#define __PLUMED_Keywords_h
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include "Tools.h"
#include "Log.h"

namespace PLMD{

/// This class lets me pass keyword types easily
class KeyType{
friend class Keyword;
private:
  enum {compulsory,flag,optional,input,numbered,modifier} style;
public:
  KeyType( const std::string& type );
  bool isCompulsory() const { return (style==compulsory); }
  bool isFlag() const { return (style==flag); }
  bool isOptional() const { return (style==optional); }
  bool isInput() const { return (style==input); }
  bool isNumbered() const { return (style==numbered); }
  bool isModifier() const { return (style==modifier); }
};

/// This function can be used to print out the keyword and its documentation in html
//void print_html( Log& log, const std::string& key, const std::string& documentation );

/// This class holds the keywords and their documentation
class Keywords{
friend class Action;
private:
  std::vector<KeyType> types;
  std::vector<std::string> keys;
  std::vector<std::string> documentation;
  std::vector<bool> defaults;
public:
  void add( const std::string, const std::string, const std::string );
  void addFlag( const std::string k, const bool def, const std::string d );
  void remove( const std::string );
  void clear();
  bool exists( const std::string k ) const ;
  KeyType style( const std::string ) const ;
  void print( Log& log ) const ;
  void print_html() const ;
  void print_html_item( const unsigned& j ) const;
};

}

#endif
