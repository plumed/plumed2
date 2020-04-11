/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_tools_Keywords_h
#define __PLUMED_tools_Keywords_h
#include <vector>
#include <string>
#include <set>
#include <map>

#include "Exception.h"

namespace PLMD {

class Log;

/// This class holds the keywords and their documentation
class Keywords {
/// This class lets me pass keyword types easily
  class KeyType {
  public:
    enum {hidden,compulsory,flag,optional,atoms,vessel} style;
    explicit KeyType( const std::string& type );
    void setStyle( const std::string& type );
    bool isCompulsory() const { return (style==compulsory); }
    bool isFlag() const { return (style==flag); }
    bool isOptional() const { return (style==optional); }
    bool isAtomList() const { return (style==atoms); }
    bool isVessel() const { return (style==vessel); }
    std::string toString() const {
      if(style==compulsory) return "compulsory";
      else if(style==optional) return "optional";
      else if(style==atoms) return "atoms";
      else if(style==flag) return "flag";
      else if(style==hidden) return "hidden";
      else if(style==vessel) return "vessel";
      else plumed_assert(0);
      return "";
    }
  };
  friend class Action;
private:
/// Is this an action or driver (this bool affects what style==atoms does in print)
  bool isaction;
/// This allows us to overwrite the behavior of the atoms type in analysis actions
  bool isatoms;
/// The names of the allowed keywords
  std::vector<std::string> keys;
/// The names of the reserved keywords
  std::vector<std::string> reserved_keys;
/// Whether the keyword is compulsory, optional...
  std::map<std::string,KeyType> types;
/// Do we allow stuff like key1, key2 etc
  std::map<std::string,bool> allowmultiple;
/// The documentation for the keywords
  std::map<std::string,std::string> documentation;
/// The default values for the flags (are they on or of)
  std::map<std::string,bool> booldefs;
/// The default values (if there are default values) for compulsory keywords
  std::map<std::string,std::string> numdefs;
/// The tags for atoms - we use this so the manual can differentiate between different ways of specifying atoms
  std::map<std::string,std::string> atomtags;
/// The string that should be printed out to describe how the components work for this particular action
  std::string cstring;
/// The names of all the possible components for an action
  std::vector<std::string> cnames;
/// The keyword that turns on a particular component
  std::map<std::string,std::string> ckey;
/// The documentation for a particular component
  std::map<std::string,std::string> cdocs;
/// Print the documentation for the jth keyword in html
  void print_html_item( const std::string& ) const;
/// Print a particular keyword
  void printKeyword( const std::string& j, Log& log ) const ;
/// Print a particular keyword (copy of the above that works with files)
  void printKeyword( const std::string& j, FILE* out ) const ;
public:
/// Constructor
  Keywords() : isaction(true), isatoms(true) {}
///
  void isDriver() { isaction=false; }
///
  void isAnalysis() { isatoms=false; }
/// find out whether flag key is on or off by default.
  bool getLogicalDefault( std::string key, bool& def ) const ;
/// Get the value of the default for the keyword named key
  bool getDefaultValue( std::string key, std::string& def ) const ;
/// Return the number of defined keywords
  unsigned size() const;
/// Check if numbered keywords are allowed for this action
  bool numbered( const std::string & k ) const ;
/// Return the ith keyword
  std::string getKeyword( const unsigned i ) const ;
/// Print the documentation to the log file (used by PLMD::Action::error)
  void print( Log& log ) const ;
/// Print the documentation to a file (use by PLUMED::CLTool::readCommandLineArgs)
  void print( FILE* out ) const ;
/// Print a file containing the list of keywords for a particular action (used for spell checking)
  void print_spelling() const ;
/// Reserve a keyword
  void reserve( const std::string & t, const std::string & k, const std::string & d );
/// Reserve a flag
  void reserveFlag( const std::string & k, const bool def, const std::string & d );
/// Use one of the reserved keywords
  void use( const std::string  & k );
/// Get the ith keyword
  std::string get( const unsigned k ) const ;
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
/// Print keywords in form readable by vim
  void print_vim() const ;
/// Print the template version for the documenation
  void print_template( const std::string& actionname, bool include_optional) const ;
/// Change the style of a keyword
  void reset_style( const std::string & k, const std::string & style );
/// Add keywords from one keyword object to another
  void add( const Keywords& keys );
/// Copy the keywords data
  void copyData( std::vector<std::string>& kk, std::vector<std::string>& rk, std::map<std::string,KeyType>& tt, std::map<std::string,bool>& am,
                 std::map<std::string,std::string>& docs, std::map<std::string,bool>& bools, std::map<std::string,std::string>& nums,
                 std::map<std::string,std::string>& atags, std::vector<std::string>& cnam, std::map<std::string,std::string>& ck,
                 std::map<std::string,std::string>& cd ) const ;
/// Clear everything from the keywords object.
/// Not actually needed if your Keywords object is going out of scope.
  void destroyData();
/// Set the text that introduces how the components for this action are introduced
  void setComponentsIntroduction( const std::string& instr );
/// Add a potential component which can be output by this particular action
  void addOutputComponent( const std::string& name, const std::string& key, const std::string& descr );
/// Has a component with this name been added?
  bool outputComponentExists( const std::string& name, const bool& custom ) const ;
/// Remove a component with a particular name from the keywords
  void removeComponent( const std::string& name );
/// Reference to keys
  std::vector<std::string> getKeys() const { return keys; }
};

}

#endif
