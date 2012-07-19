/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_Keywords_h
#define __PLUMED_Keywords_h
#include <vector>
#include <string>
#include <set>
#include <map>

namespace PLMD{

class Log;

/// This class lets me pass keyword types easily
class KeyType{
friend class Keyword;
private:
  enum {compulsory,flag,optional,atoms,nohtml,hidden} style;
public:
  KeyType( const std::string& type );
  bool isCompulsory() const { return (style==compulsory); }
  bool isFlag() const { return (style==flag); }
  bool isOptional() const { return (style==optional); }
  bool isAtomList() const { return (style==atoms); }
  bool isNoHTML() const { return (style==nohtml); }
  bool isHidden() const { return (style==hidden); }
};

/// This class holds the keywords and their documentation
class Keywords{
friend class Action;
private:
/// Whether the keyword is compulsory, optional...
  std::vector<KeyType> reserved_types;
/// The names of the keyword
  std::vector<std::string> reserved_keys;
/// The allowance of stuff like key1, key2 etc
  std::vector<bool> reserved_allowmultiple;
/// The documentation for the keywords
  std::vector<std::string> reserved_documentation;
/// Whether the keyword is compulsory, optional...
  std::vector<KeyType> types;
/// The names of the keyword
  std::vector<std::string> keys;
/// Do we allow stuff like key1, key2 etc
  std::vector<bool> allowmultiple;
/// The documentation for the keywords
  std::vector<std::string> documentation;
/// The default values for the flags (are they on or of)
  std::map<std::string,bool> booldefs; 
/// The default values (if there are default values) for compulsory keywords
  std::map<std::string,std::string> numdefs;
/// Print the documentation for the jth keyword in html
  void print_html_item( const unsigned& j ) const;
/// Print a particular keyword
  void printKeyword( const unsigned& j, Log& log ) const ;
/// find out whether flag key is on or off by default.
  bool getLogicalDefault( std::string key, bool& def ) const ;
/// Get the value of the default for the keyword named key
  bool getDefaultValue( std::string key, std::string& def ) const ;
/// Clear all the keywords
  void clear();
/// Return the number of defined keywords 
  unsigned size() const;
/// Check if numbered keywords are allowed for this action
  bool numbered( const std::string & k ) const ;
public:
/// Print the documentation to the log file (used by PLMD::Action::error)
  void print( Log& log ) const ;
/// Reserve a keyword 
  void reserve( const std::string & t, const std::string & k, const std::string & d );
/// Use one of the reserved keywords
  void use( const std::string  k );
/// Add a new keyword of type t with name k and description d
  void add( const std::string & t, const std::string & k, const std::string & d );
/// Add a new compulsory keyword (t must equal compulsory) with name k, default value def and description d
  void add( const std::string & t, const std::string & k, const std::string & def, const std::string & d );
/// Add a falg with name k that is by default on if def is true and off if def is false.  d should provide a description of the flag
  void addFlag( const std::string & k, const bool def, const std::string & d );
/// Remove the keyword with name k
  void remove( const std::string & k );
/// Check if there is a keyword with name k
  bool exists( const std::string & k ) const ;
/// Check the keyword k has been reserved
  bool reserved( const std::string & k ) const ;
/// Check if the keyword with name k has style t
  bool style( const std::string & k, const std::string & t ) const ;
/// Print an html version of the documentation
  void print_html() const ;
/// Change the style of a keyword
  void reset_style( const std::string & k, const std::string & style );
};

}

#endif
