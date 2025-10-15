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
#include "Keywords.h"
#include "Log.h"
#include "Tools.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

template <typename T>
void erase_remove(std::vector<T>& vec, const T& value) {
  vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
}

void erase_remove(std::string& vec, const char value) {
  vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
}

//few definition to avoid rewriting the too many times the same docstring
#define NUMBERED_DOCSTRING(key) ". You can use multiple instances of this keyword i.e. " \
        + std::string(key) +"1, " + std::string(key) + "2, " + std::string(key) + "3..."


namespace PLMD {

std::string toString(Keywords::argType at) {
  //the simple cases
  switch (at) {
  case Keywords::argType::scalar:
    return  "scalar";

  case Keywords::argType::vector:
    return  "vector";

  case Keywords::argType::matrix:
    return  "matrix";

  case Keywords::argType::grid:
    return  "grid";
  }
  //the not simple cases
  {
    std::string ret="";
    std::string next="";
    if(valid(at & Keywords::argType::scalar)) {
      ret+="scalar";
      next="/";
    }
    if(valid(at & Keywords::argType::vector)) {
      ret+=next+"vector";
      next="/";
    }
    if(valid(at & Keywords::argType::matrix)) {
      ret+=next+"matrix";
      next="/";
    }
    if(valid(at & Keywords::argType::grid)) {
      ret+=next+"grid";
    }
    return ret;
  }
  //the return is outsids so the compile should block the compilation
  //when expanding the enum without updating the toString
  return "";
}

Keywords::argType stoat(std::string_view str) {
  if(auto pos = str.find("/"); pos!=str.npos) {
    //here we can express that we do not want certain combinations
    auto val=stoat(str.substr(0,pos));
    return val | stoat(str.substr(pos+1));
  }
  if (str == "scalar") {
    return Keywords::argType::scalar;
  }
  if (str == "vector") {
    return Keywords::argType::vector;
  }
  if (str == "matrix") {
    return Keywords::argType::matrix;
  }
  if (str == "grid") {
    return Keywords::argType::grid;
  }
  // Handle the case where the string does not match any enum value.
  plumed_massert(false,"invalid argType specifier " + std::string(str));
}

std::string toString(Keywords::componentType at) {
  switch (at) {
  case Keywords::componentType::scalar:
    return  "scalar";

  case Keywords::componentType::vector:
    return  "vector";

  case Keywords::componentType::matrix:
    return  "matrix";

  case Keywords::componentType::grid:
    return  "grid";

  case Keywords::componentType::atom:
    return  "atom";

  case Keywords::componentType::atoms:
    return  "atoms";
  }
  //the not simple cases
  {
    std::string ret="";
    std::string next="";
    if(valid(at & Keywords::componentType::scalar)) {
      ret+="scalar";
      next="/";
    }
    if(valid(at & Keywords::componentType::vector)) {
      ret+=next+"vector";
      next="/";
    }
    if(valid(at & Keywords::componentType::matrix)) {
      ret+=next+"matrix";
      next="/";
    }
    if(valid(at & Keywords::componentType::grid)) {
      ret+=next+"grid";
      next="/";
    }
    //I do not think these two are necessary
    if(valid(at & Keywords::componentType::atom)) {
      ret+=next+"atom";
      next="/";
    }
    if(valid(at & Keywords::componentType::atoms)) {
      ret+=next+"atoms";
    }
    return ret;
  }
  //the return is outsids so the compile should block the compilation
  //when expanding the enum without updating the toString
  return "";
}

inline Keywords::componentType stoct(std::string_view str) {
  if(auto pos = str.find("/"); pos!=str.npos) {
    //here we can express that we do not want certain combinations
    auto val=stoct(str.substr(0,pos));
    return val | stoct(str.substr(pos+1));
  }
  if (str == "scalar") {
    return Keywords::componentType::scalar;
  }
  if (str == "grid") {
    return Keywords::componentType::grid;
  }
  if (str == "vector") {
    return Keywords::componentType::vector;
  }
  if (str == "matrix") {
    return Keywords::componentType::matrix;
  }
  if (str == "atom") {
    return Keywords::componentType::atom;
  }
  if (str == "atoms") {
    return Keywords::componentType::atoms;
  }

  plumed_massert(false,"invalid componentType specifier " + std::string(str));
}

Keywords::KeyType::keyStyle Keywords::KeyType::keyStyleFromString(std::string_view type ) {
  if( type=="compulsory" ) {
    return keyStyle::compulsory;
  } else if( type=="flag" ) {
    return keyStyle::flag;
  } else if( type=="optional" ) {
    return keyStyle::optional;
    //this is special: some atoms keywords have extra characters usually a "-" followed by a number
  } else if( type.find("atoms")!=type.npos || type.find("residues")!=type.npos) {
    return keyStyle::atoms;
  } else if( type=="hidden" ) {
    return keyStyle::hidden;
  } else if( type=="deprecated" ) {
    return keyStyle::deprecated;
  } else {
    plumed_massert(false,"invalid keyword specifier " + std::string(type));
  }
}

Keywords::KeyType::KeyType( std::string_view type )
  : style(keyStyleFromString(type)) {}

Keywords::KeyType::KeyType( Keywords::KeyType::keyStyle type )
  : style(type) {}

void Keywords::KeyType::setStyle( std::string_view type ) {
  style=keyStyleFromString(type);
}

Keywords::keyInfo::keyInfo()
  : type(Keywords::KeyType::keyStyle::unknown),
    docstring(""),
    defaultValue(std::monostate()),
    allowmultiple(false)
{}

Keywords::keyInfo& Keywords::keyInfo::setType(Keywords::KeyType t) {
  type=t;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setDocString(std::string_view d) {
  docstring=d;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setDefaultValue(std::string_view d) {
  defaultValue=std::string(d);
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setDefaultFlag(bool a) {
  defaultValue=a;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setArgumentType(argType a) {
  argument_type=a;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setAllowMultiple(bool a) {
  allowmultiple=a;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setLinkedAction(std::string_view a) {
  linkaction=a;
  return *this;
}
Keywords::keyInfo& Keywords::keyInfo::setLinkedPage(std::string_view a) {
  linkpage=a;
  return *this;
}
bool Keywords::keyInfo::isArgument() const {
  return std::holds_alternative<argType>(argument_type);
}

Keywords::component::component():
//the 0 ensures something that always fails unles explicitly set
  type(static_cast<componentType>(0)) {}
Keywords::component& Keywords::component::setKey(std::string_view k) {
  key=k;
  return *this;
}
Keywords::component& Keywords::component::setDocstring(std::string_view d) {
  docstring=d;
  return *this;
}
Keywords::component& Keywords::component::setType(componentType t) {
  type=t;
  return *this;
}

std::string Keywords::getStyle( const std::string&  k ) const {
  plumed_massert( exists(k)||reserved(k), "Did not find keyword " + k );
  return (keywords.at(k).type).toString();
}

void Keywords::add( const Keywords& newkeys ) {
  //copies data from
  //loop on the declared keys
  for(const auto& thiskey:newkeys.keys) {
    plumed_massert( !exists(thiskey), "keyword " + thiskey + " is in twice" );
    plumed_massert( !reserved(thiskey), "keyword " + thiskey + " is in twice" );
    keywords[thiskey] = newkeys.keywords.at(thiskey);
    keys.emplace_back( thiskey );
  }
  //loop on the reserved keys
  for (const auto&thiskey : newkeys.reserved_keys) {
    plumed_massert( !exists(thiskey), "keyword " + thiskey + " is in twice" );
    plumed_massert( !reserved(thiskey), "keyword " + thiskey + " is in twice" );

    keywords[thiskey] = newkeys.keywords.at(thiskey);
    reserved_keys.emplace_back( thiskey );
  }
  for (const auto& thisnam : newkeys.cnames) {
    plumed_massert( components.find(thisnam)!=components.end(), "keyword for component" + thisnam + " is in twice" );
    cnames.push_back( thisnam );
    components[thisnam]=newkeys.components.at(thisnam);
  }
}

void Keywords::addOrReserve( std::string_view keytype,
                             std::string_view key,
                             std::string_view docstring,
                             const bool reserve ) {
  plumed_massert(!exists(key),  "keyword " + std::string(key) + " has already been registered");
  plumed_massert(!reserved(key),"keyword " + std::string(key) + " has already been reserved");
  std::string t_type{keytype};
  bool isNumbered = keytype=="numbered";
  if( isNumbered ) {
    t_type="optional";
  }
  //let's fail asap in case of typo
  auto type = KeyType(t_type);
  if (!reserve) {
    plumed_massert( !type.isFlag(),   "use addFlag() to register a flag keyword (" + std::string(key) + ")");
  }

  std::string fd{docstring};
  bool allowMultiple= false;
  if( isNumbered ) {
    fd += NUMBERED_DOCSTRING(key);
    allowMultiple = true;
  }

  keywords[std::string(key)] = keyInfo()
                               .setType(type)
                               .setDocString(fd)
                               .setAllowMultiple(allowMultiple)
                               .setLinkedAction("none");
  if( type.isAtomList() ) {
    //keytype may be "residues" or something like "atoms-3"
    keywords.find(key)->second.atomtag=keytype;
    if (isaction && keytype == "atoms") { //this narrow the doctrstring ONLY to atoms
      keywords.find(key)->second.docstring+= ".  For more information on how to specify lists of atoms see \\ref Group";
    }
  }
  if (reserve) {
    reserved_keys.emplace_back(key);
  } else {
    keys.emplace_back(key);
  }
}

void Keywords::reserve( std::string_view keytype,
                        std::string_view key,
                        std::string_view docstring ) {
  //If you modify this function, please update also the add with three arguments
  addOrReserve(keytype,key,docstring,true);
}

void Keywords::reserveFlag(const std::string & key, const bool defaultValue, const std::string & docstring ) {
  plumed_assert( !exists(key) && !reserved(key) );
  std::string defstr;
  if( defaultValue ) {
    defstr="( default=on ) ";
  } else {
    defstr="( default=off ) ";
  }

  keywords[key] = keyInfo()
                  .setType(KeyType{KeyType::keyStyle::flag})
                  .setDocString(defstr + docstring)
                  .setAllowMultiple(false)
                  .setDefaultFlag(defaultValue)
                  .setLinkedAction("none");
  reserved_keys.emplace_back(key);
}

///this "copies" a reserved key into the keylist so it can be used
void Keywords::use(std::string_view  k ) {
  plumed_massert( reserved(k), "the " + std::string(k) + " keyword is not reserved");
  keys.emplace_back(k);
}

void Keywords::reset_style( const std::string & k, const std::string & style ) {
  plumed_massert( exists(k) || reserved(k), "no " + k + " keyword" );
  //Adding this two feels correct, but breaks some actions where a numbered keyword is changed to compulsory
  //So also the atomtag is removed
  //keywords.at(k).atomtag="";
  //keywords.at(k).allowmultiple=false;
  if( style=="numbered" ) {
    keywords.at(k).allowmultiple=true;
    return;
  }
  keywords.at(k).type.setStyle(style);
  if( (keywords.at(k).type).isAtomList() ) {
    keywords.at(k).atomtag=style;
  }
}

void Keywords::add(std::string_view keytype,
                   std::string_view key,
                   std::string_view docstring ) {
  //the 'false' deactivates the "reserve mode"
  addOrReserve(keytype,key,docstring,false);
}

void Keywords::addDeprecatedFlag( const std::string& key,
                                  const std::string& replacement ) {
  addDeprecatedKeyword(key, replacement);
  keywords.at(key).setDefaultFlag(false);
}

void Keywords::addDeprecatedKeyword( std::string_view key,
                                     const std::string& replacement ) {
  if( exists(replacement) ) {
    std::string docs = "You should use " + replacement + " instead of this keyword which was used in older versions of PLUMED and is provided for back compatibility only.";
    add("deprecated", key, docs);
  } else {
    add("deprecated", key, "Including this keyword in the input to this action makes no difference to the calculation performed it was used in older versions of PLUMED and is provided here for back compatibility only.");
  }
}

void Keywords::addInputKeyword( const std::string & keyType,
                                const std::string & key,
                                const std::string & datatype,
                                const std::string & docstring ) {
  addInputKeyword(keyType,key,stoat(datatype),docstring);
}

void Keywords::addInputKeyword( const std::string & keyType,
                                const std::string & key,
                                argType datatype,
                                const std::string & docstring ) {
  if( exists(key) ) {
    remove(key);
  }
  //insert({k,datatype}) Inserts element(s) into the container, if the container doesn't already contain an element with an equivalent key.[cit.]
  //operator[] inserts if the key doesn't exist, or overwrites if it does
  add( keyType, key, docstring );
  keywords.at(key).setArgumentType(datatype);
}

void Keywords::addInputKeyword( const std::string & keyType,
                                const std::string & key,
                                const std::string & datatype,
                                const std::string & defaultV,
                                const std::string & docstring ) {
  addInputKeyword(keyType,key,stoat(datatype),defaultV,docstring);
}

void Keywords::addInputKeyword( const std::string & keyType,
                                const std::string & key,
                                argType datatype,
                                const std::string & defaultV,
                                const std::string & docstring ) {
  if( exists(key) ) {
    remove(key);
  }
  add( keyType, key, defaultV, docstring );
  keywords[key].setArgumentType(datatype);
}

void Keywords::add( std::string_view keytype,
                    std::string_view key,
                    std::string_view defaultValue,
                    std::string_view docstring ) {
  //let's fail asap in case of typo
  auto type = KeyType(keytype);

  plumed_massert( !exists(key) && !reserved(key), "failing on keyword " + std::string(key) );
  // An optional keyword can't have a default
  plumed_massert(type.isCompulsory() || type.isHidden(), "You can't set a default value for an optional keyword, failing on " + std::string(key));
  keywords[std::string(key)] = keyInfo()
                               .setType(type)
                               .setDefaultValue(defaultValue)
                               .setDocString("( default=" + std::string(defaultValue) + " ) " + std::string(docstring) )
                               .setAllowMultiple(false)
                               .setLinkedAction("none");

  keys.emplace_back(key);
}

void Keywords::addFlag(std::string_view key, bool defaultValue, std::string_view docstring) {
  plumed_massert( !exists(key) && !reserved(key), "keyword " + std::string(key) + " has already been registered");
  plumed_massert( !defaultValue, "the second argument to addFlag must be false " + std::string(key) );
  std::string defstr="( default=off ) ";
  keywords[std::string(key)] = keyInfo()
                               .setType(KeyType("flag"))
                               .setDefaultFlag(false)
                               .setDocString(std::string(defstr) + std::string(docstring))
                               .setAllowMultiple(false)
                               .setLinkedAction("none");

  keys.emplace_back(key);
}

void Keywords::remove( const std::string & k ) {
  bool found=false;
  if(exists(k)) {
    erase_remove(keys,k);
    found=true;
  } else if(reserved(k)) {
    erase_remove(reserved_keys,k);
    found=true;
  }
  plumed_massert(found,"You are trying to forbid " + k + " a keyword that isn't there");
  // Delete documentation, type and so on from the description
  keywords.erase(k);

  // Remove any output components that this keyword creates
  //we need the double loop because we should not remove and iterate on the map at the same time
  std::vector<std::string> markForRemoval{};
  for(const auto& dkey : components ) {
    if( dkey.second.key==k ) {
      markForRemoval.push_back(dkey.first);
    }
  }
  for(const auto& toremove : markForRemoval ) {
    removeOutputComponent( toremove );
  }
}

bool Keywords::numbered( const std::string & k ) const {
  if( style( k,"atoms") ) {
    return true;
  }
  //should I add also the "reserved(k)" to the test?
  plumed_massert( exists(k), "Did not find keyword " + k );
  return keywords.at(k).allowmultiple;
}

bool Keywords::style( const std::string & k, const std::string & t ) const {
  if( getStyle(k)==t ) {
    return true;
  }
  return false;
}

unsigned Keywords::size() const {
  return keys.size();
}

std::string Keywords::getKeyword( const unsigned i ) const {
  plumed_assert( i<size() );
  return keys[i];
}

bool Keywords::exists( std::string_view k ) const {
  return std::find(keys.begin(), keys.end(), k) != keys.end();
}

bool Keywords::reserved( std::string_view k ) const {
  return std::find(reserved_keys.begin(), reserved_keys.end(), k) != reserved_keys.end();
}

void Keywords::print_template(const std::string& actionname, bool include_optional) const {
  std::printf("%s",actionname.c_str());
  {
    std::string prevtag="start";
    for(const auto& key : keys) {
      if( keywords.at(key).type.isAtomList() ) {
        plumed_massert( keywords.at(key).atomtag!="", "keyword " + key + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        const auto & currentTag=keywords.at(key).atomtag;
        if( prevtag!="start" && prevtag!=currentTag ) {
          break;
        }
        if( currentTag.find("residues")!=std::string::npos) {
          std::printf(" %s=<residue selection>", key.c_str() );
        } else {
          std::printf(" %s=<atom selection>", key.c_str() );
        }
        prevtag=currentTag;
      }
    }
  }

  for(const auto& key : keys) {
    if ( keywords.at(key).type.isCompulsory() ) {
      std::string def;
      if( getDefaultValue( key, def) ) {
        std::printf(" %s=%s ", key.c_str(), def.c_str() );
      } else {
        std::printf(" %s=    ", key.c_str() );
      }
    } else if (include_optional) {
      // TG no defaults for optional keywords?
      std::printf(" [%s]", key.c_str() );

    }
  }
  std::printf("\n");
  std::flush(std::cout);
}

void Keywords::print_vim() const {
  for(const auto& key : keys) {
    if( keywords.at(key).type.isFlag() ) {
      std::printf( ",flag:%s", key.c_str() );
    } else {
      if( keywords.at(key).allowmultiple ) {
        std::printf(",numbered:%s",key.c_str() );
      } else {
        std::printf(",option:%s",key.c_str() );
      }
    }
  }
  std::fprintf(stdout, "\n%s", getHelpString().c_str() );
}

void Keywords::print_html() const {
// This is the part that outputs the details of the components
  if( cnames.size()>0 ) {
    unsigned ndef=0;
    //running on the order of insertion
    for(const auto& cname : cnames) {
      if(components.at(cname).key=="default") {
        ndef++;
      }
    }

    if( ndef>0 ) {
      std::cout<<"\\par Description of components\n\n";
      std::cout<<cstring<<"\n\n";
      std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
      std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Description </b> </td> </tr>\n");
      unsigned nndef=0;
      for(const auto& cname : cnames) {
        //plumed_assert( components.at(cname).key=="default" );
        if( components.at(cname).key!="default" ) {
          nndef++;
          continue;
        }
        std::printf("<tr>\n");
        std::printf("<td width=15%%> <b> %s </b></td>\n",cname.c_str() );
        std::printf("<td> %s </td>\n",(components.at(cname).docstring).c_str() );
        std::printf("</tr>\n");
      }
      std::cout<<"</table>\n\n";
      if( nndef>0 ) {
        std::cout<<"In addition the following quantities can be calculated by employing the keywords listed below"<<std::endl;
        std::cout<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(const auto& cname : cnames) {
          if( components.at(cname).key!="default") {
            std::printf("<tr>\n");
            std::printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                        cname.c_str(),(components.at(cname).key).c_str() );
            std::printf("<td> %s </td>\n",(components.at(cname).docstring).c_str() );
            std::printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    } else {
      unsigned nregs=0;
      for(const auto& cname : cnames) {
        if( exists(components.at(cname).key) ) {
          nregs++;
        }
      }
      if( nregs>0 ) {
        std::cout<<"\\par Description of components\n\n";
        std::cout<<cstring<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(const auto& cname : cnames) {
          if( exists(components.at(cname).key) ) {
            std::printf("<tr>\n");
            std::printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                        cname.c_str(),(components.at(cname).key).c_str() );
            std::printf("<td> %s </td>\n",(components.at(cname).docstring).c_str() );
            std::printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    }
  }

  unsigned nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isAtomList() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    if(isaction && isatoms) {
      std::cout<<"\\par The atoms involved can be specified using\n\n";
    } else if(isaction) {
      std::cout<<"\\par The data to analyze can be the output from another analysis algorithm\n\n";
    } else {
      std::cout<<"\\par The input trajectory is specified using one of the following\n\n";
    }
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    std::string prevtag="start";
    unsigned counter=0;
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isAtomList() ) {
        const auto& currentTag = keywords.at(key).atomtag;
        plumed_massert( currentTag!="", "keyword " + key + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        if( prevtag!="start" && prevtag!=currentTag && isaction ) {
          std::cout<<"</table>\n\n";
          if( isatoms ) {
            std::cout<<"\\par Or alternatively by using\n\n";
          } else if( counter==0 ) {
            std::cout<<"\\par Alternatively data can be collected from the trajectory using \n\n";
            counter++;
          } else {
            std::cout<<"\\par Lastly data collected in a previous analysis action can be reanalyzed by using the keyword \n\n";
          }
          std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        }
        print_html_item( key );
        prevtag=currentTag;
      }
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isCompulsory() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    if(isaction) {
      std::cout<< "\\par Compulsory keywords\n\n";
    } else {
      std::cout<<"\\par The following must be present\n\n";
    }
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isCompulsory() ) {
        print_html_item( key );
      }
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isFlag() || keywords.at(key).type.isOptional() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    if(isaction) {
      std::cout<<"\\par Options\n\n";
    } else {
      std::cout<<"\\par The following options are available\n\n";
    }
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isFlag() ) {
        print_html_item( key );
      }
    }
    std::cout<<"\n";
  }
  nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isOptional() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isOptional() ) {
        print_html_item( key );
      }
    }
  }
  std::cout<<"</table>\n\n";
}

void Keywords::print_spelling() const {
  for(const auto& key : keys) {
    std::printf("%s\n", key.c_str() );
  }
  for(const auto& cname : cnames) {
    std::printf("%s\n",cname.c_str() );
  }
}

std::string Keywords::getKeywordDocs( const std::string& key ) const {
  bool killdot=( keywords.at(key).docstring.find("\\f$")!=std::string::npos ); // Check for latex
  std::vector<std::string> w=Tools::getWords( keywords.at(key).docstring );
  std::stringstream sstr;
  sstr<<std::setw(23)<<key<<" - ";
  unsigned nl=0;
  std::string blank=" ";
  for(unsigned i=0; i<w.size(); ++i) {
    nl+=w[i].length() + 1;
    if( nl>60 ) {
      sstr<<"\n"<<std::setw(23)<<blank<<"   "<<w[i]<<" ";
      nl=0;
    } else {
      sstr<<w[i]<<" ";
    }
    if( killdot && w[i].find(".")!=std::string::npos ) {
      break;  // If there is latex only write up to first dot
    }
  }
  sstr<<"\n";
  return sstr.str();
}

std::string Keywords::getHelpString() const {
  std::string helpstr;
  unsigned nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isAtomList() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    helpstr += "The input trajectory can be in any of the following formats: \n\n";
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isAtomList() ) {
        helpstr += getKeywordDocs( key );
      }
    }
  }
  unsigned ncompulsory=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isCompulsory() ) {
      ncompulsory++;
    }
  }
  if( ncompulsory>0 ) {
    helpstr += "\nThe following arguments are compulsory: \n\n";
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isCompulsory() ) {
        helpstr += getKeywordDocs( key );
      }
    }
  }
  nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isFlag() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    if(ncompulsory>0) {
      helpstr += "\nIn addition you may use the following options: \n\n";
    } else {
      helpstr += "\nThe following options are available\n\n";
    }
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isFlag() ) {
        helpstr += getKeywordDocs( key ).c_str();
      }
    }
  }
  nkeys=0;
  for(const auto& key : keys) {
    if ( keywords.at(key).type.isOptional() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    for(const auto& key : keys) {
      if ( keywords.at(key).type.isOptional() ) {
        helpstr += getKeywordDocs( key );
      }
    }
    helpstr += "\n";
  }
  return helpstr;
}

void Keywords::print( Log& log ) const {
  log.printf("%s", getHelpString().c_str() );
}

void Keywords::print( FILE* out ) const {
  fprintf( out,"%s", getHelpString().c_str() );
}

std::string Keywords::getTooltip( const std::string& name ) const {
  std::size_t dd=name.find_first_of("0123456789");
  std::string kname=name.substr(0,dd);
  if( !exists(kname) ) {
    return "<b> could not find this keyword </b>";
  }
  std::string mystring;
  std::string docstr = keywords.at(kname).docstring;
  if( keywords.at(kname).type.isCompulsory() ) {
    mystring += "<b>compulsory keyword ";
    if( docstr.find("default")!=std::string::npos ) {
      std::size_t bra = docstr.find_first_of(")");
      mystring += docstr.substr(0,bra+1);
      docstr = docstr.substr(bra+1);
    }
    mystring += "</b>\n";
  }
  std::vector<std::string> w=Tools::getWords( docstr );
  unsigned nl=0;
  for(unsigned i=0; i<w.size(); ++i) {
    nl+=w[i].length() + 1;
    if( nl>80 ) {
      mystring += w[i] + "\n";
      nl=0;
    } else {
      mystring += w[i] + " ";
    }
    if( w[i].find(".")!=std::string::npos ) {
      break;  // Only write up the the first dot
    }
  }
  return mystring;
}

void Keywords::print_html_item( const std::string& key ) const {
  std::printf("<tr>\n");
  std::printf("<td width=15%%> <b> %s </b></td>\n",key.c_str() );
  std::printf("<td> %s </td>\n",(keywords.at(key).docstring).c_str() );
  std::printf("</tr>\n");
}

bool Keywords::getLogicalDefault(const std::string & key, bool& def ) const {
  // plumed_massert(exists(key)||reserved(key),"You can't ask for the default value of a keyword that doesn't exist("+key+")");
  if (std::holds_alternative<bool>(keywords.at(key).defaultValue)) {
    def = std::get<bool>(keywords.at(key).defaultValue);
    return true;
  } else {
    return false;
  }
}

bool Keywords::getDefaultValue(const std::string & key, std::string& def ) const {
  plumed_massert( style(key,"compulsory") || style(key,"hidden"),"You can't ask for the default value of a keyword that doesn't have one ("+key+")" );
  if (std::holds_alternative<std::string>(keywords.at(key).defaultValue)) {
    def = std::get<std::string>(keywords.at(key).defaultValue);
    return true;
  } else {
    return false;
  }
}

void Keywords::destroyData() {
  keys.clear();
  reserved_keys.clear();
  keywords.clear();
  components.clear();
  //cname was missing before, is it wanted or not?
  cnames.clear();
}

void Keywords::setComponentsIntroduction( const std::string& instr ) {
  cstring = instr;
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, const std::string& descr ) {
  addOutputComponent( name, key, "scalar", descr );
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, const std::string& type, const std::string& descr ) {
  addOutputComponent(name,key,stoct(type),descr);
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, componentType type, const std::string& descr ) {
  plumed_assert( !outputComponentExists(name) );
  plumed_massert( name!=".#!value", name + " is reserved for storing description of value" );
  plumed_massert( name.find("-")==std::string::npos,"dash is reseved character in component names" );

  std::size_t num2=name.find_first_of("_");
  if( num2!=std::string::npos ) {
    char uu = '_';
    plumed_massert( std::count(name.begin(),name.end(), uu)==1, "underscore is reserved character in component names and there should only be one");
    plumed_massert( num2==0, "underscore is reserved character in component names that has special meaning");
  }
  if( key=="default" ) {
    cstring = "By default this Action calculates the following quantities. These quantities can "
              "be referenced elsewhere in the input by using this Action's label followed by a "
              "dot and the name of the quantity required from the list below.";
  }
  components[name] = component()
                     .setKey(key)
                     .setDocstring(descr)
                     .setType(type);
  cnames.emplace_back(name);
}

void Keywords::setValueDescription( const std::string& type, const std::string& descr ) {
  setValueDescription(stoct (type),descr);
}

void Keywords::setValueDescription( componentType type, const std::string& descr ) {
  if( !outputComponentExists(".#!value") ) {
    components[".#!value"] =component()
                            .setKey("default")
                            .setDocstring(descr)
                            .setType(type);
    cnames.emplace_back(".#!value");
  } else {
    components[".#!value"].docstring = descr;
    components[".#!value"].type = type;
  }
}

bool Keywords::outputComponentExists( const std::string& name ) const {
  if( cstring.find("customize")!=std::string::npos ) {
    return true;
  }

  std::string sname;
  std::size_t num=name.find_first_of("-");
  std::size_t num2=name.find_last_of("_");

  if( num2!=std::string::npos ) {
    sname=name.substr(num2);
  } else if( num!=std::string::npos ) {
    sname=name.substr(0,num);
  } else {
    sname=name;
  }

  return components.find(sname)!=components.end();
}

bool Keywords::componentHasCorrectType( const std::string& name, const std::size_t& rank, const bool& hasderiv ) const {
  if( cstring.find("customize")!=std::string::npos ) {
    return true;
  }

  std::string sname;
  std::size_t num=name.find_first_of("-");
  std::size_t num2=name.find_last_of("_");
  if( num2!=std::string::npos ) {
    sname=name.substr(num2);
  } else if( num!=std::string::npos ) {
    sname=name.substr(0,num);
  } else {
    sname=name;
  }

  // using valid(components.at(sname).type & (componentType::atom | componentType::atoms) will have a sligthly different flavour
  // == means "is exactly", the valid(&) construct instead measn "can be different, but must contain the asked flag"
  if( thisactname=="CENTER" && (components.at(sname).type == componentType::atom || components.at(sname).type == componentType::atoms)) {
    return true;
  }

  if( rank==0 ) {
    return (valid(components.at(sname).type | componentType::scalar));
  } else if( hasderiv ) {
    return (valid(components.at(sname).type | componentType::grid));
  } else if( rank==1 ) {
    return (valid(components.at(sname).type | componentType::vector));
  } else if( rank==2 ) {
    return (valid(components.at(sname).type | componentType::matrix ));
  }
  return false;
}

std::vector<std::string> Keywords::getArgumentKeys() const {
  std::vector<std::string> arguments;
  std::copy_if(keys.begin(), keys.end(),std::back_inserter(arguments),
  [this](auto const& kw) {
    return keywords.at(kw).isArgument();
  });
  return arguments;
}

bool Keywords::checkArgumentType( const std::size_t& rank, const bool& hasderiv ) const {
  bool allArgumentsAreCorrect = true;
  for(auto const& kw : getArgumentKeys() ) {
    const auto & at = std::get<argType>(keywords.at(kw).argument_type);
    bool kwIsCorrect = false;
    if( rank==0  && valid(at | argType::scalar)) {
      kwIsCorrect = true;
    }
    if( hasderiv && valid(at | argType::grid)) {
      kwIsCorrect = true;
    }
    if( rank==1  && valid(at | argType::vector)) {
      kwIsCorrect = true;
    }
    if( rank==2  && valid(at | argType::matrix)) {
      kwIsCorrect = true;
    }
    allArgumentsAreCorrect &= kwIsCorrect;
  }
  return allArgumentsAreCorrect;
}

std::string Keywords::getArgumentType( const std::string& name ) const {
  auto argument_keys = getArgumentKeys();
  if( find(argument_keys.begin(),argument_keys.end(),name)==argument_keys.end() ) {
    return "";
  }
  return toString(std::get<argType>(keywords.at(name).argument_type));
}

std::string Keywords::getOutputComponentFlag( const std::string& name ) const {
  return components.at(name).key;
}

std::string Keywords::getOutputComponentType( const std::string& name ) const {
  //return toString( components.find(name)->second.type); brings to segfault in case name is ot present
  //at at least throws
  return toString( components.at(name).type);
}

std::string Keywords::getOutputComponentDescription( const std::string& name ) const {
  std::string checkname = name;
  std::size_t hyp=name.find_first_of("-");
  if( hyp!=std::string::npos ) {
    checkname = name.substr(0,hyp);
  }

  bool found=components.find(checkname)!=components.end();

  if( !found ) {
    if( name==".#!value" ) {
      return "the value calculated by this action";
    }
    if( outputComponentExists( name ) ) {
      plumed_merror("cannot find description for component " + name + " that allegedly exists. Gareth Tribello might know what the fuck that is about.");
    }
    plumed_merror("could not find output component named " + name );
  }
  return components.at(checkname).docstring;
}

///////////DUPLICATED??????????///////
void Keywords::removeOutputComponent( const std::string& name ) {
  components.erase(name);
  erase_remove(cnames,name);
}

std::string Keywords::getKeywordDescription( const std::string& key ) const {
  plumed_assert( exists( key ) );
  return keywords.at(key).docstring;
}

void Keywords::needsAction( const std::string& name ) {
  if( std::find(neededActions.begin(), neededActions.end(), name )!=neededActions.end() ) {
    return;
  }
  neededActions.push_back( name );
}

bool Keywords::isActionNeeded( std::string_view name ) const {
  return std::find(neededActions.begin(), neededActions.end(), name )!=neededActions.end();
}

const std::vector<std::string>& Keywords::getNeededKeywords() const {
  return neededActions;
}

void Keywords::addActionNameSuffix( const std::string& suffix ) {
  if( std::find(actionNameSuffixes.begin(), actionNameSuffixes.end(), suffix )!=actionNameSuffixes.end() ) {
    return;
  }
  actionNameSuffixes.push_back( suffix );
}

bool Keywords::isActionSuffixed( std::string_view name, std::string_view basename) const {
  std::string bname{basename};
  return std::any_of(actionNameSuffixes.begin(),
                     actionNameSuffixes.end(),
  [name,&bname](const std::string& suffix)->bool{
    return (bname + suffix)==name ;
  }
                    );
}

void Keywords::setDisplayName( const std::string& name ) {
  thisactname = name;
}

std::string Keywords::getDisplayName() const {
  return thisactname;
}

void Keywords::setDeprecated( const std::string& name ) {
  replaceaction = name;
}

std::string Keywords::getReplacementAction() const {
  return replaceaction;
}

void Keywords::addDOI( const std::string& doi ) {
  doilist.push_back( doi );
}

const std::vector<std::string>& Keywords::getDOIList() const {
  return doilist;
}

void Keywords::linkActionInDocs( const std::string& k, const std::string& action ) {
  plumed_massert( exists(k), "no " + k + " keyword" );
  keywords.at(k).setLinkedAction(action);
}

void Keywords::addLinkInDocForFlag( const std::string& k, const std::string& page ) {
  plumed_massert( exists(k), "no " + k + " keyword" );
  plumed_massert( style(k,"flag"), k + " is not a flag" );
  keywords.at(k).setLinkedPage(page);
}

std::string Keywords::getLinkedActions( const std::string& key ) const {
  plumed_assert( exists( key ) );
  return keywords.at(key).linkaction;
}

std::string Keywords::getLinkedPages( const std::string& key ) const {
  plumed_assert( exists( key ) );
  return keywords.at(key).linkpage;
}

}// namespace PLMD
