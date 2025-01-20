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
namespace PLMD {

std::string toString(Keywords::argType at) {
  //the simple cases
  switch (at) {
  case Keywords::argType::scalar:
    return  "scalar";

  case Keywords::argType::grid:
    return  "grid";

  case Keywords::argType::vector:
    return  "vector";

  case Keywords::argType::matrix:
    return  "matrix";
  }
  //the not simple cases
  {
    std::string ret="";
    std::string next="";
    if(valid(at & Keywords::argType::scalar)) {
      ret+="scalar";
      next="/";
    }
    if(valid(at & Keywords::argType::grid)) {
      ret+=next+"grid";
      next="/";
    }
    if(valid(at & Keywords::argType::vector)) {
      ret+=next+"vector";
      next="/";
    }
    if(valid(at & Keywords::argType::matrix)) {
      ret+=next+"matrix";
    }
    return ret;
  }
  //the return is outsids so the compile should block the compilation
  //when expanding the enum without updating the toString
  return "";
}

Keywords::argType stoat(std::string_view str) {
  using namespace std::literals;
  if(auto pos = str.find("/"sv); pos!=str.npos) {
    //here we can express that we do not want certain combinations
    auto val=stoat(str.substr(0,pos));
    return val | stoat(str.substr(pos+1));
  }
  if (str == "scalar") {
    return Keywords::argType::scalar;
  }
  if (str == "grid") {
    return Keywords::argType::grid;
  }
  if (str == "vector") {
    return Keywords::argType::vector;
  }
  if (str == "matrix") {
    return Keywords::argType::matrix;
  }
  // Handle the case where the string does not match any enum value.
  plumed_massert(false,"invalid argType specifier " + std::string(str));
}

std::string toString(Keywords::componentType at) {
  switch (at) {
  case Keywords::componentType::scalar:
    return  "scalar";

  case Keywords::componentType::grid:
    return  "grid";

  case Keywords::componentType::vector:
    return  "vector";

  case Keywords::componentType::matrix:
    return  "matrix";

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
    if(valid(at & Keywords::componentType::grid)) {
      ret+=next+"grid";
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
    //I do not think these two are necessary
    if(valid(at & Keywords::componentType::atom)) {
      ret+=next+"atom";
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
  using namespace std::literals;
  if(auto pos = str.find("/"sv); pos!=str.npos) {
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
  } else if( type=="vessel" ) {
    return keyStyle::vessel;
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

std::string Keywords::getStyle( const std::string & k ) const {
  plumed_massert( types.count(k), "Did not find keyword " + k );
  return (types.find(k)->second).toString();
}

void Keywords::add( const Keywords& newkeys ) {
  newkeys.copyData( keys, reserved_keys, types, allowmultiple, documentation, booldefs, numdefs, atomtags, cnames, ckey, cdocs  );
}

void Keywords::copyData( std::vector<std::string>& kk,
                         std::vector<std::string>& rk,
                         std::map<std::string,KeyType,std::less<void>>& tt, std::map<std::string,bool>& am,
                         std::map<std::string,std::string>& docs, std::map<std::string,bool>& bools, std::map<std::string,std::string>& nums,
                         std::map<std::string,std::string>& atags, std::vector<std::string>& cnam, std::map<std::string,std::string>& ck,
                         std::map<std::string,std::string>& cd ) const {
  for(unsigned i=0; i<keys.size(); ++i) {
    std::string thiskey=keys[i];
    for(unsigned j=0; j<kk.size(); ++j) {
      plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
    }
    for(unsigned j=0; j<rk.size(); ++j) {
      plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
    }
    kk.push_back( thiskey );
    plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
    tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
    if( (types.find(thiskey)->second).isAtomList() ) {
      atags.insert( std::pair<std::string,std::string>( thiskey,atomtags.find(thiskey)->second) );
    }
    plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" );
    am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
    plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" );
    docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
    if( booldefs.count( thiskey ) ) {
      bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) );
    }
    if( numdefs.count( thiskey ) ) {
      nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
    }
  }
  for(unsigned i=0; i<reserved_keys.size(); ++i) {
    std::string thiskey=reserved_keys[i];
    for(unsigned j=0; j<kk.size(); ++j) {
      plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
    }
    for(unsigned j=0; j<rk.size(); ++j) {
      plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
    }
    rk.push_back( thiskey );
    plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
    tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
    if( (types.find(thiskey)->second).isAtomList() ) {
      atags.insert( std::pair<std::string,std::string>( thiskey,atomtags.find(thiskey)->second) );
    }
    plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" );
    am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
    plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" );
    docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
    if( booldefs.count( thiskey ) ) {
      bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) );
    }
    if( numdefs.count( thiskey ) ) {
      nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
    }
  }
  for(unsigned i=0; i<cnames.size(); ++i) {
    std::string thisnam=cnames[i];
    for(unsigned j=0; j<cnam.size(); ++j) {
      plumed_massert( thisnam!=cnam[j], "component " + thisnam + " is in twice" );
    }
    cnam.push_back( thisnam );
    plumed_massert( ckey.count( thisnam ), "no keyword data on component " + thisnam + " to copy" );
    ck.insert( std::pair<std::string,std::string>( thisnam, ckey.find(thisnam)->second) );
    plumed_massert( cdocs.count( thisnam ), "no documentation on component " + thisnam + " to copy" );
    cd.insert( std::pair<std::string,std::string>( thisnam, cdocs.find(thisnam)->second) );
  }
}

#define NUMBERED_DOCSTRING(key) ". You can use multiple instances of this keyword i.e. " + std::string(key) +"1, " + std::string(key) + "2, " + std::string(key) + "3..."
void Keywords::reserve( const std::string & keytype,
                        const std::string & key,
                        const std::string & docstring ) {
  plumed_assert( !exists(key) && !reserved(key) );
  std::string t_type{keytype};
  bool isNumbered = keytype=="numbered";
  if( isNumbered ) {
    t_type="optional";
  }
  //let's fail asap in case of typo
  auto type = KeyType(t_type);
  plumed_assert( !exists(key) && !reserved(key) );

  std::string fd{docstring};
  bool allowMultiple= false;
  if( type.isVessel() ) {
    // Convert to lower case
    std::string lowkey{key};
    std::transform(lowkey.begin(),lowkey.end(),lowkey.begin(),[](unsigned char c) {
      return std::tolower(c);
    });
    // Remove any underscore characters
    lowkey.erase(std::remove(lowkey.begin(), lowkey.end(), '_'), lowkey.end());

    fd += " The final value can be referenced using <em>label</em>." + lowkey;
    if(docstring.find("flag")==std::string::npos) {
      fd += NUMBERED_DOCSTRING(key) "  The corresponding values are then "
            "referenced using <em>label</em>."+ lowkey +"-1,  <em>label</em>." + lowkey +
            "-2,  <em>label</em>." + lowkey + "-3...";
    }
    allowMultiple = true;
  } else if( isNumbered ) {
    fd += NUMBERED_DOCSTRING(key);
    allowMultiple = true;
  } else {
    //if( type.isAtomList() && isaction ) {//<- why not this? atoms could also be "residues" or "atoms-n"
    if( keytype=="atoms" && isaction ) {
      fd += ".  For more information on how to specify lists of atoms see \\ref Group";
    }
    if( type.isAtomList() ) {
      atomtags.insert( std::pair<std::string,std::string>(key,keytype) );
    }
  }
  types.insert( std::pair<std::string,KeyType>(key,type) );
  allowmultiple.insert( std::pair<std::string,bool>(key,allowMultiple) );
  documentation.insert( std::pair<std::string,std::string>(key,fd) );
  reserved_keys.emplace_back(key);
}

void Keywords::reserveFlag(const std::string & k, const bool def, const std::string & d ) {
  plumed_assert( !exists(k) && !reserved(k) );
  std::string defstr;
  if( def ) {
    defstr="( default=on ) ";
  } else {
    defstr="( default=off ) ";
  }
  types.insert( std::pair<std::string,KeyType>(k,KeyType::keyStyle::flag) );
  std::string fd,lowkey=k;
  std::transform(lowkey.begin(),lowkey.end(),lowkey.begin(),[](unsigned char c) {
    return std::tolower(c);
  });
  fd=defstr + d;
  documentation.insert( std::pair<std::string,std::string>(k,fd) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) );
  reserved_keys.push_back(k);
}

void Keywords::use(std::string_view  k ) {
  plumed_massert( reserved(k), "the " + std::string(k) + " keyword is not reserved");
  keys.emplace_back(k);
  //reserved(k) verifies reserved_keys[i]==k, so no need to traverse again the reserved keys?
  //since, from reserve(), the reserved keys are unique (shall we use a set?)
  // for(unsigned i=0; i<reserved_keys.size(); ++i) {
  //   if(reserved_keys[i]==k) keys.push_back( reserved_keys[i] );
  // }
}

void Keywords::reset_style( const std::string & k, const std::string & style ) {
  plumed_massert( exists(k) || reserved(k), "no " + k + " keyword" );
  if( style=="numbered" ) {
    allowmultiple[k]=true;
    return;
  }
  (types.find(k)->second).setStyle(style);
  if( (types.find(k)->second).isVessel() ) {
    allowmultiple[k]=true;
  }
  if( (types.find(k)->second).isAtomList() ) {
    atomtags.insert( std::pair<std::string,std::string>(k,style) );
  }
}

void Keywords::add(std::string_view keytype,
                   std::string_view key,
                   std::string_view docstring ) {
  std::string t_type{keytype};
  bool isNumbered = keytype=="numbered";
  if( isNumbered ) {
    t_type="optional";
  }
  //let's fail asap in case of typo
  auto type = KeyType(t_type);
  plumed_massert( !exists(key) && (!type.isFlag()) && !reserved(key) && (!type.isVessel()),
                  "keyword " + std::string(key) + " has already been registered");
  std::string fd;
  fd=docstring;
  if( isNumbered ) {
    fd += NUMBERED_DOCSTRING(key);
  } else {
    if( (types.find(key)->second).isAtomList() ) {
      //keytype may be "residues" or something like "atoms-3"
      atomtags.insert( std::pair<std::string,std::string>(key,keytype) );
    }
  }
  allowmultiple.insert( std::pair<std::string,bool>(key,isNumbered) );
  types.insert( std::pair<std::string,KeyType>(key, type) );
  if( type.isAtomList() && isaction ) {
    fd += ".  For more information on how to specify lists of atoms see \\ref Group";
  }
  documentation.insert( std::pair<std::string,std::string>(key,fd) );
  keys.emplace_back(key);
}

void Keywords::addInputKeyword( const std::string & typekey,
                                const std::string & key,
                                const std::string & datatype,
                                const std::string & docstring ) {
  addInputKeyword(typekey,key,stoat(datatype),docstring);
}

void Keywords::addInputKeyword( const std::string & typekey,
                                const std::string & key,
                                argType datatype,
                                const std::string & docstring ) {
  if( exists(key) ) {
    remove(key);
  }
  //insert({k,datatype}) Inserts element(s) into the container, if the container doesn't already contain an element with an equivalent key.[cit.]
  //operator[] inserts if the key doesn't exist, or overwrites if it does
  argument_types[key] = datatype;
  add( typekey, key, docstring );
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
  argument_types[key] = datatype;
  add( keyType, key, defaultV, docstring );
}

void Keywords::add( std::string_view keytype,
                    std::string_view key,
                    std::string_view  defaultValue,
                    std::string_view docstring ) {
  //let's fail asap in case of typo
  auto type = KeyType(keytype);
  // An optional keyword can't have a default
  plumed_massert( !exists(key) && !reserved(key) &&
                  (type.isCompulsory() || type.isHidden() ), "failing on keyword " + std::string(key) );
  types.insert(  std::pair<std::string,KeyType> {key, type} );
  documentation.insert( std::pair<std::string,std::string>(key,"( default=" + std::string(defaultValue) + " ) " + std::string(docstring) ));
  allowmultiple.insert( std::pair<std::string,bool>(key,false) );
  numdefs.insert( std::pair<std::string,std::string>(key,defaultValue) );
  keys.emplace_back(key);
}

void Keywords::addFlag( const std::string & k, const bool def, const std::string & d ) {
  plumed_massert( !exists(k) && !reserved(k), "keyword " + k + " has already been registered");
  std::string defstr;
  plumed_massert( !def, "the second argument to addFlag must be false " + k );
  defstr="( default=off ) ";
  types.insert( std::pair<std::string,KeyType>(k,KeyType("flag")) );
  documentation.insert( std::pair<std::string,std::string>(k,defstr + d) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) );
  keys.push_back(k);
}

void Keywords::remove( const std::string & k ) {
  bool found=false;
  unsigned j=0, n=0;

  while(true) {
    for(j=0; j<keys.size(); j++) {
      if(keys[j]==k) {
        break;
      }
    }
    for(n=0; n<reserved_keys.size(); n++) {
      if(reserved_keys[n]==k) {
        break;
      }
    }
    if(j<keys.size()) {
      keys.erase(keys.begin()+j);
      found=true;
    } else if(n<reserved_keys.size()) {
      reserved_keys.erase(reserved_keys.begin()+n);
      found=true;
    } else {
      break;
    }
  }
  // Delete documentation, type and so on from the description
  types.erase(k);
  documentation.erase(k);
  allowmultiple.erase(k);
  booldefs.erase(k);
  numdefs.erase(k);
  // Remove any output comonents that this keyword creates
  for(const auto& dkey : ckey ) {
    if( dkey.second==k ) {
      removeOutputComponent( dkey.first );
    }
  }
  plumed_massert(found,"You are trying to forbid " + k + " a keyword that isn't there"); // You have tried to forbid a keyword that isn't there
}

bool Keywords::numbered( const std::string & k ) const {
  if( style( k,"atoms") ) {
    return true;
  }
  plumed_massert( allowmultiple.count(k), "Did not find keyword " + k );
  return allowmultiple.find(k)->second;
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
  unsigned nkeys=0;
  std::printf("%s",actionname.c_str());
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    std::string prevtag="start";
    for(unsigned i=0; i<keys.size(); ++i) {
      if( (types.find(keys[i])->second).isAtomList() ) {
        plumed_massert( atomtags.count(keys[i]), "keyword " + keys[i] + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        if( prevtag!="start" && prevtag!=atomtags.find(keys[i])->second ) {
          break;
        }
        if( (atomtags.find(keys[i])->second).find("residues")!=std::string::npos) {
          std::printf(" %s=<residue selection>", keys[i].c_str() );
        } else {
          std::printf(" %s=<atom selection>", keys[i].c_str() );
        }
        prevtag=atomtags.find(keys[i])->second;
      }
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( include_optional || \
         (types.find(keys[i])->second).isCompulsory() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) {
        std::string def;
        if( getDefaultValue( keys[i], def) ) {
          std::printf(" %s=%s ", keys[i].c_str(), def.c_str() );
        } else {
          std::printf(" %s=    ", keys[i].c_str() );
        }
      } else if (include_optional) {
        // TG no defaults for optional keywords?
        std::printf(" [%s]", keys[i].c_str() );
      }
    }
  }
  std::printf("\n");
  std::flush(std::cout);
}

void Keywords::print_vim() const {
  for(unsigned i=0; i<keys.size(); ++i) {
    if( (types.find(keys[i])->second).isFlag() ) {
      std::printf( ",flag:%s", keys[i].c_str() );
    } else {
      if( allowmultiple.find(keys[i])->second ) {
        std::printf(",numbered:%s",keys[i].c_str() );
      } else {
        std::printf(",option:%s",keys[i].c_str() );
      }
    }
  }
  std::fprintf(stdout, "\n%s", getHelpString().c_str() );
}

void Keywords::print_html() const {

// This is the part that outputs the details of the components
  if( cnames.size()>0 ) {
    unsigned ndef=0;
    for(unsigned i=0; i<cnames.size(); ++i) {
      if(ckey.find(cnames[i])->second=="default") {
        ndef++;
      }
    }

    if( ndef>0 ) {
      std::cout<<"\\par Description of components\n\n";
      std::cout<<cstring<<"\n\n";
      std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
      std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Description </b> </td> </tr>\n");
      unsigned nndef=0;
      for(unsigned i=0; i<cnames.size(); ++i) {
        //plumed_assert( ckey.find(cnames[i])->second=="default" );
        if( ckey.find(cnames[i])->second!="default" ) {
          nndef++;
          continue;
        }
        std::printf("<tr>\n");
        std::printf("<td width=15%%> <b> %s </b></td>\n",cnames[i].c_str() );
        std::printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
        std::printf("</tr>\n");
      }
      std::cout<<"</table>\n\n";
      if( nndef>0 ) {
        std::cout<<"In addition the following quantities can be calculated by employing the keywords listed below"<<std::endl;
        std::cout<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(unsigned i=0; i<cnames.size(); ++i) {
          if( ckey.find(cnames[i])->second!="default") {
            std::printf("<tr>\n");
            std::printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                        cnames[i].c_str(),(ckey.find(cnames[i])->second).c_str() );
            std::printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
            std::printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    } else {
      unsigned nregs=0;
      for(unsigned i=0; i<cnames.size(); ++i) {
        if( exists(ckey.find(cnames[i])->second) ) {
          nregs++;
        }
      }
      if( nregs>0 ) {
        std::cout<<"\\par Description of components\n\n";
        std::cout<<cstring<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        std::printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(unsigned i=0; i<cnames.size(); ++i) {
          if( exists(ckey.find(cnames[i])->second) ) {
            std::printf("<tr>\n");
            std::printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                        cnames[i].c_str(),(ckey.find(cnames[i])->second).c_str() );
            std::printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
            std::printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    }
  }

  unsigned nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) {
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
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isAtomList() ) {
        plumed_massert( atomtags.count(keys[i]), "keyword " + keys[i] + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        if( prevtag!="start" && prevtag!=atomtags.find(keys[i])->second && isaction ) {
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
        print_html_item( keys[i] );
        prevtag=atomtags.find(keys[i])->second;
      }
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isCompulsory() ) {
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
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) {
        print_html_item( keys[i] );
      }
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isFlag() || (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) {
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
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isFlag() ) {
        print_html_item( keys[i] );
      }
    }
    std::cout<<"\n";
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) {
        print_html_item( keys[i] );
      }
    }
  }
  std::cout<<"</table>\n\n";
}

void Keywords::print_spelling() const {
  for(unsigned i=0; i<keys.size(); ++i) {
    std::printf("%s\n", keys[i].c_str() );
  }
  for(unsigned i=0; i<cnames.size(); ++i) {
    std::printf("%s\n",cnames[i].c_str() );
  }
}

std::string Keywords::getKeywordDocs( const std::string& key ) const {
  bool killdot=( (documentation.find(key)->second).find("\\f$")!=std::string::npos ); // Check for latex
  std::vector<std::string> w=Tools::getWords( documentation.find(key)->second );
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
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    helpstr += "The input trajectory can be in any of the following formats: \n\n";
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isAtomList() ) {
        helpstr += getKeywordDocs( keys[i] );
      }
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isCompulsory() ) {
      nkeys++;
    }
  }
  unsigned ncompulsory=nkeys;
  if( nkeys>0 ) {
    helpstr += "\nThe following arguments are compulsory: \n\n";
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) {
        helpstr += getKeywordDocs( keys[i] );
      }
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isFlag() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    if(ncompulsory>0) {
      helpstr += "\nIn addition you may use the following options: \n\n";
    } else {
      helpstr += "\nThe following options are available\n\n";
    }
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isFlag() ) {
        helpstr += getKeywordDocs( keys[i] ).c_str();
      }
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) {
      nkeys++;
    }
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) {
        helpstr += getKeywordDocs( keys[i] );
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
  std::string mystring, docstr = documentation.find(kname)->second;
  if( types.find(kname)->second.isCompulsory() ) {
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
  std::printf("<td> %s </td>\n",(documentation.find(key)->second).c_str() );
  std::printf("</tr>\n");
}

std::string Keywords::get( const unsigned k ) const {
  plumed_assert( k<size() );
  return keys[k];
}

bool Keywords::getLogicalDefault(const std::string & key, bool& def ) const {
  if( booldefs.find(key)!=booldefs.end() ) {
    def=booldefs.find(key)->second;
    return true;
  } else {
    return false;
  }
}

bool Keywords::getDefaultValue(const std::string & key, std::string& def ) const {
  plumed_assert( style(key,"compulsory") || style(key,"hidden") );

  if( numdefs.find(key)!=numdefs.end() ) {
    def=numdefs.find(key)->second;
    return true;
  } else {
    return false;
  }
}

void Keywords::destroyData() {
  keys.clear();
  reserved_keys.clear();
  types.clear();
  allowmultiple.clear();
  documentation.clear();
  booldefs.clear();
  numdefs.clear();
  atomtags.clear();
  ckey.clear();
  cdocs.clear();
  ckey.clear();
}

void Keywords::setComponentsIntroduction( const std::string& instr ) {
  cstring = instr;
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, const std::string& descr ) {
  addOutputComponent( name, key, "scalar", descr );
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, const std::string& type, const std::string& descr ) {
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

  ckey.insert( std::pair<std::string,std::string>(name,key) );
  cdocs.insert( std::pair<std::string,std::string>(name,descr) );
  ctypes.insert( std::pair<std::string,componentType>(name,stoct(type)) );
  cnames.push_back(name);
}

void Keywords::removeOutputComponent( const std::string& name ) {
  unsigned j=0;
  while(true) {
    for(j=0; j<cnames.size(); j++)
      if(cnames[j]==name) {
        break;
      }
    if(j<cnames.size()) {
      cnames.erase(cnames.begin()+j);
    } else {
      break;
    }
  }
  cdocs.erase(name);
}

void Keywords::setValueDescription( const std::string& type, const std::string& descr ) {
  if( !outputComponentExists(".#!value") ) {
    ckey.insert( std::pair<std::string,std::string>(".#!value","default") );
    cdocs.insert( std::pair<std::string,std::string>(".#!value",descr) );
    ctypes.insert( std::pair<std::string,componentType>(".#!value",stoct (type)) );
    cnames.push_back(".#!value");
  } else {
    cdocs[".#!value"] = descr;
    ctypes[".#!value"] = stoct(type);
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

  for(unsigned i=0; i<cnames.size(); ++i) {
    if( sname==cnames[i] ) {
      return true;
    }
  }
  return false;
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

  if( thisactname=="CENTER" && (components.at(sname).type== componentType::atom || components.at(sname).type== componentType::atoms) ) {
    return true;
  }

  if( rank==0 ) {
    return (valid(ctypes.find(sname)->second | componentType::scalar));
  } else if( hasderiv ) {
    return (valid(ctypes.find(sname)->second | componentType::grid));
  } else if( rank==1 ) {
    return (valid(ctypes.find(sname)->second | componentType::vector));
  } else if( rank==2 ) {
    return (valid(ctypes.find(sname)->second | componentType::matrix ));
  }
  return false;
}

bool Keywords::checkArgumentType( const std::size_t& rank, const bool& hasderiv ) const {
  for(auto const& x : argument_types ) {
    if( rank==0  && valid(x.second | argType::scalar)) {
      return true;
    }
    if( hasderiv && valid(x.second | argType::grid)) {
      return true;
    }
    if( rank==1  && valid(x.second | argType::vector)) {
      return true;
    }
    if( rank==2  && valid(x.second | argType::matrix)) {
      return true;
    }
  }
  plumed_merror("WARNING: type for input argument has not been specified");
  return false;
}

std::string Keywords::getArgumentType( const std::string& name ) const {
  if( argument_types.find(name)==argument_types.end() ) {
    return "";
  }
  return toString(argument_types.find(name)->second);
}

std::string Keywords::getOutputComponentFlag( const std::string& name ) const {
  return ckey.find(name)->second;
}

std::string Keywords::getOutputComponentType( const std::string& name ) const {
  return toString( components.find(name)->second.type);
}

std::string Keywords::getOutputComponentDescription( const std::string& name ) const {
  std::string checkname = name;
  std::size_t hyp=name.find_first_of("-");
  if( hyp!=std::string::npos ) {
    checkname = name.substr(0,hyp);
  }

  bool found=false;
  for(unsigned i=0; i<cnames.size(); ++i) {
    if( checkname==cnames[i] ) {
      found=true;
    }
  }
  if( !found ) {
    if( name==".#!value" ) {
      return "the value calculated by this action";
    }
    if( outputComponentExists( name ) ) {
      plumed_merror("cannot find description for component " + name + " that allegedly exists. Gareth Tribello might know what the fuck that is about.");
    }
    plumed_merror("could not find output component named " + name );
  }
  return cdocs.find(checkname)->second;
}

void Keywords::removeComponent( const std::string& name ) {
  bool found=false;

  while(true) {
    unsigned j;
    for(j=0; j<cnames.size(); j++)
      if(cnames[j]==name) {
        break;
      }
    if(j<cnames.size()) {
      cnames.erase(cnames.begin()+j);
      found=true;
    } else {
      break;
    }
  }
  // Delete documentation, type and so on from the description
  cdocs.erase(name);
  ckey.erase(name);
  plumed_massert(found,"You are trying to remove " + name + " a component that isn't there");
}

std::vector<std::string> Keywords::getOutputComponents() const {
  return cnames;
}

std::string Keywords::getKeywordDescription( const std::string& key ) const {
  plumed_assert( exists( key ) );
  return documentation.find(key)->second;
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
  ///@todo convert into a find, or asses that this is more readable
  for(auto const& suffix : actionNameSuffixes ) {
    if( (bname + suffix)==name ) {
      return true;
    }
  }
  return false;
}

void Keywords::setDisplayName( const std::string& name ) {
  thisactname = name;
}

std::string Keywords::getDisplayName() const {
  return thisactname;
}

}
