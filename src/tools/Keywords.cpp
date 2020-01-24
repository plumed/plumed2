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
#include "Keywords.h"
#include "Log.h"
#include "Tools.h"
#include <iostream>

namespace PLMD {

Keywords::KeyType::KeyType( const std::string& type ) {
  if( type=="compulsory" ) {
    style=compulsory;
  } else if( type=="flag" ) {
    style=flag;
  } else if( type=="optional" ) {
    style=optional;
  } else if( type.find("atoms")!=std::string::npos || type.find("residues")!=std::string::npos ) {
    style=atoms;
  } else if( type=="hidden" ) {
    style=hidden;
  } else if( type=="vessel" ) {
    style=vessel;
  } else {
    plumed_massert(false,"invalid keyword specifier " + type);
  }
}

void Keywords::KeyType::setStyle( const std::string& type ) {
  if( type=="compulsory" ) {
    style=compulsory;
  } else if( type=="flag" ) {
    style=flag;
  } else if( type=="optional" ) {
    style=optional;
  } else if( type.find("atoms")!=std::string::npos || type.find("residues")!=std::string::npos ) {
    style=atoms;
  } else if( type=="hidden" ) {
    style=hidden;
  } else if( type=="vessel" ) {
    style=vessel;
  } else {
    plumed_massert(false,"invalid keyword specifier " + type);
  }
}

void Keywords::add( const Keywords& newkeys ) {
  newkeys.copyData( keys, reserved_keys, types, allowmultiple, documentation, booldefs, numdefs, atomtags, cnames, ckey, cdocs  );
}

void Keywords::copyData( std::vector<std::string>& kk, std::vector<std::string>& rk, std::map<std::string,KeyType>& tt, std::map<std::string,bool>& am,
                         std::map<std::string,std::string>& docs, std::map<std::string,bool>& bools, std::map<std::string,std::string>& nums,
                         std::map<std::string,std::string>& atags, std::vector<std::string>& cnam, std::map<std::string,std::string>& ck,
                         std::map<std::string,std::string>& cd ) const {
  for(unsigned i=0; i<keys.size(); ++i) {
    std::string thiskey=keys[i];
    for(unsigned j=0; j<kk.size(); ++j) plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
    for(unsigned j=0; j<rk.size(); ++j) plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
    kk.push_back( thiskey );
    plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
    tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
    if( (types.find(thiskey)->second).isAtomList() ) atags.insert( std::pair<std::string,std::string>( thiskey,atomtags.find(thiskey)->second) );
    plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" );
    am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
    plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" );
    docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
    if( booldefs.count( thiskey ) ) bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) );
    if( numdefs.count( thiskey ) ) nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
  }
  for(unsigned i=0; i<reserved_keys.size(); ++i) {
    std::string thiskey=reserved_keys[i];
    for(unsigned j=0; j<kk.size(); ++j) plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
    for(unsigned j=0; j<rk.size(); ++j) plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
    rk.push_back( thiskey );
    plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
    tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
    if( (types.find(thiskey)->second).isAtomList() ) atags.insert( std::pair<std::string,std::string>( thiskey,atomtags.find(thiskey)->second) );
    plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" );
    am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
    plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" );
    docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
    if( booldefs.count( thiskey ) ) bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) );
    if( numdefs.count( thiskey ) ) nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
  }
  for(unsigned i=0; i<cnames.size(); ++i) {
    std::string thisnam=cnames[i];
    for(unsigned j=0; j<cnam.size(); ++j) plumed_massert( thisnam!=cnam[j], "component " + thisnam + " is in twice" );
    cnam.push_back( thisnam );
    plumed_massert( ckey.count( thisnam ), "no keyword data on component " + thisnam + " to copy" );
    ck.insert( std::pair<std::string,std::string>( thisnam, ckey.find(thisnam)->second) );
    plumed_massert( cdocs.count( thisnam ), "no documentation on component " + thisnam + " to copy" );
    cd.insert( std::pair<std::string,std::string>( thisnam, cdocs.find(thisnam)->second) );
  }
}

void Keywords::reserve( const std::string & t, const std::string & k, const std::string & d ) {
  plumed_assert( !exists(k) && !reserved(k) );
  std::string fd, lowkey=k;
  // Convert to lower case
  std::transform(lowkey.begin(),lowkey.end(),lowkey.begin(),tolower);
// Remove any underscore characters
  for(unsigned i=0;; ++i) {
    std::size_t num=lowkey.find_first_of("_");
    if( num==std::string::npos ) break;
    lowkey.erase( lowkey.begin() + num, lowkey.begin() + num + 1 );
  }
  if( t=="vessel" ) {
    fd = d + " The final value can be referenced using <em>label</em>." + lowkey;
    if(d.find("flag")==std::string::npos) fd += ".  You can use multiple instances of this keyword i.e. " +
          k +"1, " + k + "2, " + k + "3...  The corresponding values are then "
          "referenced using <em>label</em>."+ lowkey +"-1,  <em>label</em>." + lowkey +
          "-2,  <em>label</em>." + lowkey + "-3...";
    allowmultiple.insert( std::pair<std::string,bool>(k,true) );
    types.insert( std::pair<std::string,KeyType>(k,KeyType("vessel")) );
  } else if( t=="numbered" ) {
    fd = d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
    allowmultiple.insert( std::pair<std::string,bool>(k,true) );
    types.insert( std::pair<std::string,KeyType>(k,KeyType("optional")) );
  } else {
    fd = d;
    if( t=="atoms" && isaction ) fd = d + ".  For more information on how to specify lists of atoms see \\ref Group";
    allowmultiple.insert( std::pair<std::string,bool>(k,false) );
    types.insert( std::pair<std::string,KeyType>(k,KeyType(t)) );
    if( (types.find(k)->second).isAtomList() ) atomtags.insert( std::pair<std::string,std::string>(k,t) );
  }
  documentation.insert( std::pair<std::string,std::string>(k,fd) );
  reserved_keys.push_back(k);
}

void Keywords::reserveFlag( const std::string & k, const bool def, const std::string & d ) {
  plumed_assert( !exists(k) && !reserved(k) );
  std::string defstr;
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.insert( std::pair<std::string,KeyType>(k,KeyType("flag")) );
  std::string fd,lowkey=k; std::transform(lowkey.begin(),lowkey.end(),lowkey.begin(),tolower);
  fd=defstr + d;
  documentation.insert( std::pair<std::string,std::string>(k,fd) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) );
  reserved_keys.push_back(k);
}

void Keywords::use( const std::string & k ) {
  plumed_massert( reserved(k), "the " + k + " keyword is not reserved");
  for(unsigned i=0; i<reserved_keys.size(); ++i) {
    if(reserved_keys[i]==k) keys.push_back( reserved_keys[i] );
  }
}

void Keywords::reset_style( const std::string & k, const std::string & style ) {
  plumed_massert( exists(k) || reserved(k), "no " + k + " keyword" );
  (types.find(k)->second).setStyle(style);
  if( (types.find(k)->second).isVessel() ) allowmultiple[k]=true;
  if( (types.find(k)->second).isAtomList() ) atomtags.insert( std::pair<std::string,std::string>(k,style) );
}

void Keywords::add( const std::string & t, const std::string & k, const std::string & d ) {
  plumed_massert( !exists(k) && t!="flag" && !reserved(k) && t!="vessel", "keyword " + k + " has already been registered");
  std::string fd;
  if( t=="numbered" ) {
    fd=d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
    allowmultiple.insert( std::pair<std::string,bool>(k,true) );
    types.insert( std::pair<std::string,KeyType>(k, KeyType("optional")) );
  } else {
    fd=d;
    allowmultiple.insert( std::pair<std::string,bool>(k,false) );
    types.insert( std::pair<std::string,KeyType>(k,KeyType(t)) );
    if( (types.find(k)->second).isAtomList() ) atomtags.insert( std::pair<std::string,std::string>(k,t) );
  }
  if( t=="atoms" && isaction ) fd = d + ".  For more information on how to specify lists of atoms see \\ref Group";
  documentation.insert( std::pair<std::string,std::string>(k,fd) );
  keys.push_back(k);
}

void Keywords::add( const std::string & t, const std::string & k, const std::string &  def, const std::string & d ) {
  plumed_assert( !exists(k) && !reserved(k) &&  (t=="compulsory" || t=="hidden" )); // An optional keyword can't have a default
  types.insert(  std::pair<std::string,KeyType>(k, KeyType(t)) );
  documentation.insert( std::pair<std::string,std::string>(k,"( default=" + def + " ) " + d) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  numdefs.insert( std::pair<std::string,std::string>(k,def) );
  keys.push_back(k);
}

void Keywords::addFlag( const std::string & k, const bool def, const std::string & d ) {
  plumed_massert( !exists(k) && !reserved(k), "keyword " + k + " has already been registered");
  std::string defstr;
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.insert( std::pair<std::string,KeyType>(k,KeyType("flag")) );
  documentation.insert( std::pair<std::string,std::string>(k,defstr + d) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) );
  keys.push_back(k);
}

void Keywords::remove( const std::string & k ) {
  bool found=false; unsigned j=0, n=0;

  while(true) {
    for(j=0; j<keys.size(); j++) if(keys[j]==k)break;
    for(n=0; n<reserved_keys.size(); n++) if(reserved_keys[n]==k)break;
    if(j<keys.size()) {
      keys.erase(keys.begin()+j);
      found=true;
    } else if(n<reserved_keys.size()) {
      reserved_keys.erase(reserved_keys.begin()+n);
      found=true;
    } else break;
  }
  // Delete documentation, type and so on from the description
  types.erase(k); documentation.erase(k); allowmultiple.erase(k); booldefs.erase(k); numdefs.erase(k);
  plumed_massert(found,"You are trying to forbid " + k + " a keyword that isn't there"); // You have tried to forbid a keyword that isn't there
}

bool Keywords::numbered( const std::string & k ) const {
  if( style( k,"atoms") ) return true;
  plumed_massert( allowmultiple.count(k), "Did not find keyword " + k );
  return allowmultiple.find(k)->second;
}

bool Keywords::style( const std::string & k, const std::string & t ) const {
  plumed_massert( types.count(k), "Did not find keyword " + k );

  if( (types.find(k)->second).toString()==t ) return true;
  return false;
}

unsigned Keywords::size() const {
  return keys.size();
}

std::string Keywords::getKeyword( const unsigned i ) const {
  plumed_assert( i<size() );
  return keys[i];
}

bool Keywords::exists( const std::string & k ) const {
  for(unsigned i=0; i<keys.size(); ++i) {
    if( keys[i]==k ) return true;
  }
  return false;
}

bool Keywords::reserved( const std::string & k ) const {
  for(unsigned i=0; i<reserved_keys.size(); ++i) {
    if( reserved_keys[i]==k ) return true;
  }
  return false;
}

void Keywords::print_template(const std::string& actionname, bool include_optional) const {
  unsigned nkeys=0;
  printf("%s",actionname.c_str());
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if( nkeys>0 ) {
    std::string prevtag="start";
    for(unsigned i=0; i<keys.size(); ++i) {
      if( (types.find(keys[i])->second).isAtomList() ) {
        plumed_massert( atomtags.count(keys[i]), "keyword " + keys[i] + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        if( prevtag!="start" && prevtag!=atomtags.find(keys[i])->second ) break;
        if( (atomtags.find(keys[i])->second).find("residues")!=std::string::npos) printf(" %s=<residue selection>", keys[i].c_str() );
        else printf(" %s=<atom selection>", keys[i].c_str() );
        prevtag=atomtags.find(keys[i])->second;
      }
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( include_optional || \
         (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) {
        std::string def;
        if( getDefaultValue( keys[i], def) ) {
          printf(" %s=%s ", keys[i].c_str(), def.c_str() );
        } else {
          printf(" %s=    ", keys[i].c_str() );
        }
      } else if (include_optional) {
        // TG no defaults for optional keywords?
        printf(" [%s]", keys[i].c_str() );
      }
    }
  }
  printf("\n");
}

void Keywords::print_vim() const {
  for(unsigned i=0; i<keys.size(); ++i) {
    if( (types.find(keys[i])->second).isFlag() ) {
      printf( ",flag:%s", keys[i].c_str() );
    } else {
      if( allowmultiple.find(keys[i])->second ) printf(",numbered:%s",keys[i].c_str() );
      else printf(",option:%s",keys[i].c_str() );
    }
  }
  fprintf(stdout,"\n");
  print(stdout);
}

void Keywords::print_html() const {

// This is the part that outputs the details of the components
  if( cnames.size()>0 ) {
    unsigned ndef=0;
    for(unsigned i=0; i<cnames.size(); ++i) {
      if(ckey.find(cnames[i])->second=="default") ndef++;
    }

    if( ndef>0 ) {
      std::cout<<"\\par Description of components\n\n";
      std::cout<<cstring<<"\n\n";
      std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
      printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Description </b> </td> </tr>\n");
      unsigned nndef=0;
      for(unsigned i=0; i<cnames.size(); ++i) {
        //plumed_assert( ckey.find(cnames[i])->second=="default" );
        if( ckey.find(cnames[i])->second!="default" ) { nndef++; continue; }
        printf("<tr>\n");
        printf("<td width=15%%> <b> %s </b></td>\n",cnames[i].c_str() );
        printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
        printf("</tr>\n");
      }
      std::cout<<"</table>\n\n";
      if( nndef>0 ) {
        std::cout<<"In addition the following quantities can be calculated by employing the keywords listed below"<<std::endl;
        std::cout<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(unsigned i=0; i<cnames.size(); ++i) {
          if( ckey.find(cnames[i])->second!="default") {
            printf("<tr>\n");
            printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                   cnames[i].c_str(),(ckey.find(cnames[i])->second).c_str() );
            printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
            printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    } else {
      unsigned nregs=0;
      for(unsigned i=0; i<cnames.size(); ++i) {
        if( exists(ckey.find(cnames[i])->second) ) nregs++;
      }
      if( nregs>0 ) {
        std::cout<<"\\par Description of components\n\n";
        std::cout<<cstring<<"\n\n";
        std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
        printf("<tr> <td width=5%%> <b> Quantity </b> </td> <td> <b> Keyword </b> </td> <td> <b> Description </b> </td> </tr>\n");
        for(unsigned i=0; i<cnames.size(); ++i) {
          if( exists(ckey.find(cnames[i])->second) ) {
            printf("<tr>\n");
            printf("<td width=5%%> <b> %s </b></td> <td width=10%%> <b> %s </b> </td> \n",
                   cnames[i].c_str(),(ckey.find(cnames[i])->second).c_str() );
            printf("<td> %s </td>\n",(cdocs.find(cnames[i])->second).c_str() );
            printf("</tr>\n");
          }
        }
        std::cout<<"</table>\n\n";
      }
    }
  }

  unsigned nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if( nkeys>0 ) {
    if(isaction && isatoms) std::cout<<"\\par The atoms involved can be specified using\n\n";
    else if(isaction) std::cout<<"\\par The data to analyze can be the output from another analysis algorithm\n\n";
    else std::cout<<"\\par The input trajectory is specified using one of the following\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    std::string prevtag="start"; unsigned counter=0;
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isAtomList() ) {
        plumed_massert( atomtags.count(keys[i]), "keyword " + keys[i] + " allegedly specifies atoms but no tag has been specified. Please email Gareth Tribello");
        if( prevtag!="start" && prevtag!=atomtags.find(keys[i])->second && isaction ) {
          std::cout<<"</table>\n\n";
          if( isatoms ) std::cout<<"\\par Or alternatively by using\n\n";
          else if( counter==0 ) { std::cout<<"\\par Alternatively data can be collected from the trajectory using \n\n"; counter++; }
          else std::cout<<"\\par Lastly data collected in a previous analysis action can be reanalyzed by using the keyword \n\n";
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
    if ( (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ) {
    if(isaction) std::cout<< "\\par Compulsory keywords\n\n";
    else std::cout<<"\\par The following must be present\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) print_html_item( keys[i] );
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isFlag() || (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) nkeys++;
  }
  if( nkeys>0 ) {
    if(isaction) std::cout<<"\\par Options\n\n";
    else std::cout<<"\\par The following options are available\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isFlag() ) print_html_item( keys[i] );
    }
    std::cout<<"\n";
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) nkeys++;
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) print_html_item( keys[i] );
    }
  }
  std::cout<<"</table>\n\n";
}

void Keywords::print_spelling() const {
  for(unsigned i=0; i<keys.size(); ++i) printf("%s\n", keys[i].c_str() );
  for(unsigned i=0; i<cnames.size(); ++i) printf("%s\n",cnames[i].c_str() );
}

void Keywords::print( FILE* out ) const {
  unsigned nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if( nkeys>0 ) {
    fprintf(out,"The input trajectory can be in any of the following formats: \n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isAtomList() ) printKeyword( keys[i], out );
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  unsigned ncompulsory=nkeys;
  if( nkeys>0 ) {
    fprintf(out,"\nThe following arguments are compulsory: \n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) printKeyword( keys[i], out );   //log.printKeyword( keys[i], documentation[i] );
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isFlag() ) nkeys++;
  }
  if( nkeys>0 ) {
    if(ncompulsory>0) fprintf( out,"\nIn addition you may use the following options: \n\n");
    else fprintf( out,"\nThe following options are available\n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isFlag() ) printKeyword( keys[i], out );   //log.printKeyword( keys[i], documentation[i] );
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) nkeys++;
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isOptional() || (types.find(keys[i])->second).isVessel() ) printKeyword( keys[i], out );   //log.printKeyword( keys[i], documentation[i] );
    }
    fprintf(out,"\n");
  }
}

void Keywords::printKeyword( const std::string& key, FILE* out ) const {
  bool killdot=( (documentation.find(key)->second).find("\\f$")!=std::string::npos ); // Check for latex
  std::vector<std::string> w=Tools::getWords( documentation.find(key)->second );
  fprintf(out,"%23s - ", key.c_str() );
  unsigned nl=0; std::string blank=" ";
  for(unsigned i=0; i<w.size(); ++i) {
    nl+=w[i].length() + 1;
    if( nl>60 ) {
      fprintf(out,"\n%23s   %s ", blank.c_str(), w[i].c_str() ); nl=0;
    } else {
      fprintf(out,"%s ", w[i].c_str() );
    }
    if( killdot && w[i].find(".")!=std::string::npos ) break; // If there is latex only write up to first dot
  }
  fprintf(out,"\n");
}

void Keywords::print( Log& log ) const {
  unsigned nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if (nkeys>0 ) {
    log.printf( "The input for this keyword can be specified using one of the following \n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isAtomList() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ) {
    log.printf( "\n The compulsory keywords for this action are: \n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isCompulsory() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
    }
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isFlag() ) nkeys++;
  }
  if( nkeys>0 ) {
    log.printf( "\n The following options are available: \n\n");
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isFlag() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
    }
    log.printf("\n");
  }
  nkeys=0;
  for(unsigned i=0; i<keys.size(); ++i) {
    if ( (types.find(keys[i])->second).isOptional() ) nkeys++;
  }
  if( nkeys>0 ) {
    for(unsigned i=0; i<keys.size(); ++i) {
      if ( (types.find(keys[i])->second).isOptional() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
    }
    log.printf("\n");
  }
}

void Keywords::printKeyword( const std::string& key, Log& log ) const {
  bool killdot=( (documentation.find(key)->second).find("\\f$")!=std::string::npos ); // Check for latex
  std::vector<std::string> w=Tools::getWords( documentation.find(key)->second );
  log.printf("%23s - ", key.c_str() );
  unsigned nl=0; std::string blank=" ";
  for(unsigned i=0; i<w.size(); ++i) {
    nl+=w[i].length() + 1;
    if( nl>60 ) {
      log.printf("\n%23s   %s ", blank.c_str(), w[i].c_str() ); nl=0;
    } else {
      log.printf("%s ", w[i].c_str() );
    }
    if( killdot && w[i].find(".")!=std::string::npos ) break; // If there is latex only write up to first dot
  }
  log.printf("\n");
}

void Keywords::print_html_item( const std::string& key ) const {
  printf("<tr>\n");
  printf("<td width=15%%> <b> %s </b></td>\n",key.c_str() );
  printf("<td> %s </td>\n",(documentation.find(key)->second).c_str() );
  printf("</tr>\n");
}

std::string Keywords::get( const unsigned k ) const {
  plumed_assert( k<size() );
  return keys[k];
}

bool Keywords::getLogicalDefault( std::string key, bool& def ) const {
  if( booldefs.find(key)!=booldefs.end() ) {
    def=booldefs.find(key)->second;
    return true;
  } else {
    return false;
  }
}

bool Keywords::getDefaultValue( std::string key, std::string& def ) const {
  plumed_assert( style(key,"compulsory") || style(key,"hidden") );

  if( numdefs.find(key)!=numdefs.end() ) {
    def=numdefs.find(key)->second;
    return true;
  } else {
    return false;
  }
}

void Keywords::destroyData() {
  keys.clear(); reserved_keys.clear(); types.clear();
  allowmultiple.clear(); documentation.clear();
  booldefs.clear(); numdefs.clear(); atomtags.clear();
  ckey.clear(); cdocs.clear(); ckey.clear();
}

void Keywords::setComponentsIntroduction( const std::string& instr ) {
  cstring = instr;
}

void Keywords::addOutputComponent( const std::string& name, const std::string& key, const std::string& descr ) {
  plumed_assert( !outputComponentExists( name, false ) );
  plumed_massert( name.find("-")==std::string::npos,"dash is reseved character in component names" );

  std::size_t num2=name.find_first_of("_");
  if( num2!=std::string::npos ) plumed_massert( num2==0, "underscore is reserved character in component names that has special meaning");

  ckey.insert( std::pair<std::string,std::string>(name,key) );
  cdocs.insert( std::pair<std::string,std::string>(name,descr) );
  cnames.push_back(name);
}

bool Keywords::outputComponentExists( const std::string& name, const bool& custom ) const {
  if( custom && cstring.find("customize")!=std::string::npos ) return true;

  std::string sname;
  std::size_t num=name.find_first_of("-");
  std::size_t num2=name.find_last_of("_");

  if( num2!=std::string::npos ) sname=name.substr(num2);
  else if( num!=std::string::npos ) sname=name.substr(0,num);
  else sname=name;

  for(unsigned i=0; i<cnames.size(); ++i) {
    if( sname==cnames[i] ) return true;
  }
  return false;
}

void Keywords::removeComponent( const std::string& name ) {
  bool found=false; unsigned j=0, n=0;

  while(true) {
    for(j=0; j<cnames.size(); j++) if(cnames[j]==name)break;
    if(j<cnames.size()) {
      cnames.erase(cnames.begin()+j);
      found=true;
    } else break;
  }
  // Delete documentation, type and so on from the description
  cdocs.erase(name); ckey.erase(name);
  plumed_massert(found,"You are trying to remove " + name + " a component that isn't there");
}

}
