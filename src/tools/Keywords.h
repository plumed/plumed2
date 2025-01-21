/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include <string_view>
#include <map>
#include <variant>

#include "Exception.h"
#include "BitmaskEnum.h"

namespace PLMD {

class Log;

/// This class holds the keywords and their documentation
class Keywords {
/// This class lets me pass keyword types easily
  struct KeyType {
    enum class keyStyle {hidden,compulsory,flag,optional,atoms,vessel,unknown} style;
    static keyStyle keyStyleFromString(std::string_view type );
    explicit KeyType( keyStyle type );
    explicit KeyType( std::string_view type );
    void setStyle( std::string_view type );
    bool isCompulsory() const {
      return (style==keyStyle::compulsory);
    }
    bool isFlag() const {
      return (style==keyStyle::flag);
    }
    bool isOptional() const {
      return (style==keyStyle::optional);
    }
    bool isAtomList() const {
      return (style==keyStyle::atoms);
    }
    bool isVessel() const {
      return (style==keyStyle::vessel);
    }
    bool isHidden() const {
      return (style==keyStyle::hidden);
    }
    std::string toString() const {
      //if you add a style and you forget to update this function the compiler will refuse to compile
      switch(style) {
      case keyStyle::compulsory:
        return "compulsory";
      case keyStyle::optional:
        return "optional";
      case keyStyle::atoms:
        return "atoms";
      case keyStyle::flag:
        return "flag";
      case keyStyle::hidden:
        return "hidden";
      case keyStyle::vessel:
        return "vessel";
      default:
        plumed_massert(false,"unknown keyword type");
      }
      return "unknown";
    }
  };

public:
  enum class argType {scalar=1,grid=1<<2,vector=1<<3,matrix=1<<4};
  enum class componentType {scalar=1,grid=1<<2,vector=1<<3,matrix=1<<4,atoms=1<<5,atom=1<<6};
private:
/// Is this an action or driver (this bool affects what style==atoms does in print)
  bool isaction;
/// This allows us to overwrite the behavior of the atoms type in analysis actions
  bool isatoms;
/// The name of the action that has this set of keywords
  std::string thisactname;

  struct keyInfo {
    /// Whether the keyword is compulsory, optional...
    KeyType type;
    /// The documentation for the keyword
    std::string docstring;
    /// The default values (if there are default values) for compulsory keywords or flags
    std::variant<std::monostate,std::string,bool> defaultValue;
    /// Do we allow stuff like key1, key2 etc
    bool allowmultiple;
    keyInfo();
    //these functions are not neeeded (this is a struct), but are useful in constructing the key
    keyInfo& setType(KeyType t);
    keyInfo& setDocString(std::string_view d);
    keyInfo& setDefaultValue(std::string_view d);
    keyInfo& setAllowMultiple(bool a);
    keyInfo& setDefaultFlag(bool a);
  };
  //std::less<void> make some magic and makes find and [] work with string_view
/// Stores the keywords along with their settings
  std::map<std::string,keyInfo,std::less<void>> keywords;
/// The names of the allowed keywords, in order of declaration
  std::vector<std::string> keys;
/// The names of the reserved keywords, in order of declaration
  std::vector<std::string> reserved_keys;
/// The type for the arguments in this action
  std::map<std::string,argType,std::less<void>> argument_types;
/// The tags for atoms - we use this so the manual can differentiate between different ways of specifying atoms
  std::map<std::string,std::string,std::less<void>> atomtags;
  struct component {
    /// The keyword that turns on this component
    std::string key;
    /// The documentation for the component
    std::string docstring;
    /// The type of the component
    componentType type;
    component();
    //these functions are not neeeded (this is a struct), but are useful in constructing the component
    component& setKey(std::string k);
    component& setDocstring(std::string d);
    component& setType(componentType t);
  };
  //the "exists component" is stored in the map keys
  std::map<std::string,component,std::less<void>> components;
/// The string that should be printed out to describe how the components work for this particular action
  std::string cstring;
/// The names of all the possible components for an action, in order of their (first) declaration
  std::vector<std::string> cnames;
/// The list of actions that are needed by this action
  std::vector<std::string> neededActions;
/// List of suffixes that can be used with this action
  std::vector<std::string> actionNameSuffixes;
/// Print the documentation for the named keyword in html
  void print_html_item( const std::string& ) const;
public:
/// Constructor
  Keywords() : isaction(true), isatoms(true) {}
///
  void isDriver() {
    isaction=false;
  }
///
  void isAnalysis() {
    isatoms=false;
  }
/// find out whether flag key is on or off by default.
  bool getLogicalDefault(const std::string & key, bool& def ) const ;
/// Get the value of the default for the keyword named key
  bool getDefaultValue(const std::string & key, std::string& def ) const ;
/// Return the number of defined keywords
  unsigned size() const;
/// Check if numbered keywords are allowed for this action
  bool numbered( const std::string & k ) const;
  /// Get the ordered list of active keywords (not the reserved ones)
  const std::vector<std::string>& getKeys() const {
    return keys;
  }
/// Return the ith keyword
  std::string getKeyword( const unsigned i ) const ;
/// Get the documentation for a particular keyword
  std::string getKeywordDocs( const std::string& key ) const ;
/// Print the documentation to the log file (used by PLMD::Action::error)
  void print( Log& log ) const ;
/// Print the documentation to a file (use by PLUMED::CLTool::readCommandLineArgs)
  void print( FILE* out ) const ;
/// Get the help string
  std::string getHelpString() const ;
/// Print a file containing the list of keywords for a particular action (used for spell checking)
  void print_spelling() const ;
/// Reserve a keyword
  void reserve( const std::string & keytype, const std::string & key, const std::string & docstring );
/// Reserve a flag
  void reserveFlag( const std::string & k, const bool def, const std::string & d );
/// Use one of the reserved keywords
  void use( std::string_view k );
  /// append the data from another Keywords object (**only** keywords, reserved keywords and components)
  void add( const Keywords& keys );
/// Add a new keyword of type t with name k and description d
  void add( std::string_view keytype, std::string_view key, std::string_view docstring );
/// Add a new compulsory keyword (t must equal compulsory) with name k, default value def and description d
  void add( std::string_view keytype, std::string_view key, std::string_view defaultValue, std::string_view docstring );
/// Add a falg with name k that is by default on if def is true and off if def is false.  d should provide a description of the flag
  void addFlag(std::string_view key, bool defaultValue, std::string_view docstring);
/// Remove the keyword with name k
  void remove( const std::string & k );
/// Check if there is a keyword with name k
  bool exists( std::string_view k ) const ;
/// Check the keyword k has been reserved
  bool reserved( std::string_view k ) const ;
/// Get the type for the keyword with string k
  std::string getStyle( const std::string & k ) const ;
/// Check if the keyword with name k has style t
  bool style( const std::string & k, const std::string & t ) const ;
/// Print an html version of the documentation
  void print_html() const ;
/// Print keywords in form readable by vim
  void print_vim() const ;
/// Print the template version for the documentation
  void print_template( const std::string& actionname, bool include_optional) const ;
/// Change the style of a keyword
  void reset_style( const std::string & k, const std::string & style );
/// Clear everything from the keywords object.
/// Not actually needed if your Keywords object is going out of scope.
  void destroyData();
/// Set the text that introduces how the components for this action are introduced
  void setComponentsIntroduction( const std::string& instr );
/// Add the description of the value
  void setValueDescription( const std::string& type, const std::string& descr );
/// Add a potential component which can be output by this particular action
  [[deprecated("Use addOutputComponent with four argument and specify valid types for value from scalar/vector/matrix/grid")]]
  void addOutputComponent( const std::string& name, const std::string& key, const std::string& descr );
/// Add a potential component which can be output by this particular action
  void addOutputComponent( const std::string& name, const std::string& key, const std::string& type, const std::string& descr );
/// Remove a component that can be output by this particular action
  void removeOutputComponent( const std::string& name );
/// Has a component with this name been added?
  bool outputComponentExists( const std::string& name ) const ;
/// Check that type for component has been documented correctly
  bool componentHasCorrectType( const std::string& name, const std::size_t& rank, const bool& hasderiv ) const ;
/// Create the documentation for a keyword that reads arguments
/// @todo prepend [[deprecated("Please specify the data type for the argument from scalar/vector/matrix/grid with the Keywords::argType enum")]]
  void addInputKeyword( const std::string & keyType, const std::string & key, const std::string & dataType, const std::string & docstring );
  /// Create the documentation for a keyword that reads arguments
  void addInputKeyword( const std::string & keyType, const std::string & key, argType dataType, const std::string & docstring );
  /// Create the documentation for a keyword that reads arguments
  /// @todo prepend[[deprecated("Please specify the data type for the argument from scalar/vector/matrix/grid with the Keywords::argType enum")]]
  void addInputKeyword( const std::string & keyType, const std::string & key, const std::string & dataType, const std::string& defaultValue, const std::string & docstring );
  /// Create the documentation for a keyword that reads arguments
  void addInputKeyword( const std::string & keyType, const std::string & key, argType dataType, const std::string& defaultValue, const std::string & docstring );
/// Check the documentation of the argument types
  bool checkArgumentType( const std::size_t& rank, const bool& hasderiv ) const ;
/// Get the valid types that can be used as argument for this keyword
  std::string getArgumentType( const std::string& name ) const ;
/// Get the flag that forces thie named component to be calculated
  std::string getOutputComponentFlag( const std::string& name ) const ;
/// Get the type for the named output component
  std::string getOutputComponentType( const std::string& name ) const ;
/// Get the description of the named component
  std::string getOutputComponentDescription( const std::string& name ) const ;
/// Get the full ordered list of output components
  const std::vector<std::string>& getOutputComponents() const {
    return cnames;
  }
/// Get the description of a particular keyword
  std::string getKeywordDescription( const std::string& name ) const ;
/// Remove a component with a particular name from the keywords
  void removeComponent( const std::string& name );
/// Get the description of a particular keyword
  std::string getTooltip( const std::string& name ) const ;
/// Note that another actions is required to create this shortcut
  void needsAction( const std::string& name );
/// Check if the requested action is in the list of the needed actions
  bool isActionNeeded( std::string_view name ) const ;
/// Add a suffix to the list of possible action name suffixes
  void addActionNameSuffix( const std::string& suffix );
  /** @brief checks that name is is a composition of basename and one of the possible suffixes

  Cycles on the registered suffixed and return true if the combination
  `basename+suffix` equals to the passed name
  */
  bool isActionSuffixed( std::string_view name, std::string_view basename) const ;
/// Get the list of keywords that are needed by this action
  const std::vector<std::string>& getNeededKeywords() const ;
/// Return the name of the action that has this set of keywords
  std::string getDisplayName() const ;
/// Set the display name
  void setDisplayName( const std::string& name );
};

//the following templates specializations make the bitmask enum work with the
// bitwise operators `|`, `&` and the "valid" function (valid converts to bool the result of a "mask operation")
template<>
struct enum_traits::BitmaskEnum< Keywords::componentType > {
  static constexpr bool has_valid = true;
  static constexpr bool has_bit_or = true;
  static constexpr bool has_bit_and = true;
};

template<>
struct enum_traits::BitmaskEnum< Keywords::argType > {
  static constexpr bool has_valid = true;
  static constexpr bool has_bit_or = true;
  static constexpr bool has_bit_and = true;
};

std::string toString(Keywords::argType at);
/**
 * Converts a string to the corresponding Keywords::argType.
 *
 * @param str The string to convert.
 * @return The Keywords::argType corresponding to the string.
 * @throws std::invalid_argument If the string does not match any enum value.
 */
Keywords::argType stoat(std::string_view str);
std::string toString(Keywords::componentType at);

/**
 * Converts a string to the corresponding Keywords::componentType.
 * @param str The string to convert.
 * @return The Keywords::componentType corresponding to the string.
 * @throws std::invalid_argument if the string does not match any enum value.
 */
Keywords::componentType stoct(std::string_view str) ;

} // namespace PLMD

#endif
