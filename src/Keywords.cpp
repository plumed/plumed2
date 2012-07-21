/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include <iostream>

using namespace PLMD;

KeyType::KeyType( const std::string& type ){
  if( type=="compulsory" ){
      style=compulsory;
  } else if( type=="flag" ){
      style=flag;
  } else if( type=="optional" ){
      style=optional;
  } else if( type=="atoms" ){
      style=atoms;
  } else if( type=="hidden" ){
      style=hidden;
  } else {
      plumed_massert(false,"invalid keyword specifier " + type);   
  }
}

void KeyType::setStyle( const std::string& type ){
  if( type=="compulsory" ){
      style=compulsory;
  } else if( type=="flag" ){
      style=flag;
  } else if( type=="optional" ){
      style=optional;
  } else if( type=="atoms" ){
      style=atoms;
  } else if( type=="hidden" ){
      style=hidden;
  } else {
      plumed_massert(false,"invalid keyword specifier " + type);    
  }  
}

void Keywords::add( const Keywords& newkeys ){
  newkeys.copyData( keys, reserved_keys, types, allowmultiple, documentation, booldefs, numdefs  ); 
}

void Keywords::copyData( std::vector<std::string>& kk, std::vector<std::string>& rk, std::map<std::string,KeyType>& tt, std::map<std::string,bool>& am, 
                         std::map<std::string,std::string>& docs, std::map<std::string,bool>& bools, std::map<std::string,std::string>& nums ) const {
  for(unsigned i=0;i<keys.size();++i){
     std::string thiskey=keys[i];
     for(unsigned j=0;j<kk.size();++j) plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
     for(unsigned j=0;j<rk.size();++j) plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
     kk.push_back( thiskey ); 
     plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
     tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
     plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" ); 
     am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
     plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" ); 
     docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
     if( booldefs.count( thiskey ) ) bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) ); 
     if( numdefs.count( thiskey ) ) nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
  }
  for(unsigned i=0;i<reserved_keys.size();++i){
     std::string thiskey=reserved_keys[i];
     for(unsigned j=0;j<kk.size();++j) plumed_massert( thiskey!=kk[j], "keyword " + thiskey + " is in twice" );
     for(unsigned j=0;j<rk.size();++j) plumed_massert( thiskey!=rk[j], "keyword " + thiskey + " is in twice" );
     rk.push_back( thiskey ); 
     plumed_massert( types.count( thiskey ), "no type data on keyword " + thiskey + " to copy" );
     tt.insert( std::pair<std::string,KeyType>( thiskey,types.find(thiskey)->second) );
     plumed_massert( allowmultiple.count( thiskey ), "no numbered data on keyword " + thiskey + " to copy" ); 
     am.insert( std::pair<std::string,bool>(thiskey,allowmultiple.find(thiskey)->second) );
     plumed_massert( documentation.count( thiskey ), "no documentation for keyword " + thiskey + " to copy" ); 
     docs.insert( std::pair<std::string,std::string>(thiskey,documentation.find(thiskey)->second) );
     if( booldefs.count( thiskey ) ) bools.insert( std::pair<std::string,bool>( thiskey,booldefs.find(thiskey)->second) ); 
     if( numdefs.count( thiskey ) ) nums.insert( std::pair<std::string,std::string>( thiskey,numdefs.find(thiskey)->second) );
  }
}   

void Keywords::reserve( const std::string & t, const std::string & k, const std::string & d ){
  plumed_assert( !exists(k) && !reserved(k) );
  std::string fd;
  if( t=="numbered" ){
     fd=d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
     allowmultiple.insert( std::pair<std::string,bool>(k,true) );
     types.insert( std::pair<std::string,KeyType>(k,KeyType("optional")) );
  } else {
     fd=d;
     allowmultiple.insert( std::pair<std::string,bool>(k,false) );
     types.insert( std::pair<std::string,KeyType>(k,KeyType(t)) );
  }
  documentation.insert( std::pair<std::string,std::string>(k,fd) ); 
  reserved_keys.push_back(k); 
}

void Keywords::reserveFlag( const std::string & k, const bool def, const std::string & d ){
  plumed_assert( !exists(k) && !reserved(k) ); 
  std::string defstr, flag="flag";
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.insert( std::pair<std::string,KeyType>(k,KeyType("flag")) );
  documentation.insert( std::pair<std::string,std::string>(k,defstr + d) ); 
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) ); 
  reserved_keys.push_back(k); 
}

void Keywords::use( const std::string k ){
  plumed_massert( reserved(k), "the " + k + " keyword is not reserved");
  for(unsigned i=0;i<reserved_keys.size();++i){
     if(reserved_keys[i]==k) keys.push_back( reserved_keys[i] ); 
  }
}

void Keywords::reset_style( const std::string & k, const std::string & style ){
  plumed_assert( exists(k) || reserved(k) );
  (types.find(k)->second).setStyle(style); 
}

void Keywords::add( const std::string & t, const std::string & k, const std::string & d ){
  plumed_assert( !exists(k) && t!="flag" && !reserved(k) );  
  std::string fd;
  if( t=="numbered" ){
     fd=d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
     allowmultiple.insert( std::pair<std::string,bool>(k,true) );
     types.insert( std::pair<std::string,KeyType>(k, KeyType("optional")) );
  } else { 
     fd=d;
     allowmultiple.insert( std::pair<std::string,bool>(k,false) );
     types.insert( std::pair<std::string,KeyType>(k,KeyType(t)) );
  }
  documentation.insert( std::pair<std::string,std::string>(k,fd) );
  keys.push_back(k);  
}

void Keywords::add( const std::string & t, const std::string & k, const std::string &  def, const std::string & d ){
  plumed_assert( !exists(k) && !reserved(k) &&  t=="compulsory" ); // An optional keyword can't have a default
  types.insert(  std::pair<std::string,KeyType>(k, KeyType(t)) ); 
  documentation.insert( std::pair<std::string,std::string>(k,"( default=" + def + " ) " + d) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  numdefs.insert( std::pair<std::string,std::string>(k,def) );
  keys.push_back(k); 
} 

void Keywords::addFlag( const std::string & k, const bool def, const std::string & d ){
  plumed_assert( !exists(k) && !reserved(k) ); std::string defstr, flag="flag";
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.insert( std::pair<std::string,KeyType>(k,KeyType("flag")) );
  documentation.insert( std::pair<std::string,std::string>(k,defstr + d) );
  allowmultiple.insert( std::pair<std::string,bool>(k,false) );
  booldefs.insert( std::pair<std::string,bool>(k,def) ); 
  keys.push_back(k);  
} 

void Keywords::remove( const std::string & k ){
  bool found=false; unsigned j=0;

  while(true){
    for(j=0;j<keys.size();j++) if(keys[j]==k)break;
    if(j<keys.size()){
      keys.erase(keys.begin()+j);
      found=true;
    } else break;
  }
  plumed_massert(found,"You are trying to forbid " + k + " a keyword that isn't there"); // You have tried to forbid a keyword that isn't there
}

bool Keywords::numbered( const std::string & k ) const {
  unsigned j=0; 
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
  for(unsigned i=0;i<keys.size();++i){
     if( keys[i]==k ) return true;
  }
  return false;
}

bool Keywords::reserved( const std::string & k ) const {
  for(unsigned i=0;i<reserved_keys.size();++i){
     if( reserved_keys[i]==k ) return true;
  }
  return false;
}

void Keywords::print_html() const {
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if( nkeys>0 ){
    std::cout<<"\\par Specifying the atoms involved\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isAtomList() ) print_html_item( keys[i] );
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<< "\\par Compulsory keywords\n\n";
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isCompulsory() ) print_html_item( keys[i] );
     }
     std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isFlag() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<<"\\par Options\n\n";
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isFlag() ) print_html_item( keys[i] );
     }
     std::cout<<"\n";
  }
  std::cout<<"</table>\n\n";
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isOptional() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isOptional() ) print_html_item( keys[i] );
     }
     std::cout<<"\n";
  }
  std::cout<<"</table>\n\n";
}

void Keywords::print( Log& log ) const {
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isAtomList() ) nkeys++;
  }
  if (nkeys>0 ){
     log.printf( "The input for this keyword can be specified using one of the following \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isAtomList() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
     }
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The compulsory keywords for this action are: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isCompulsory() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
     }
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isFlag() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The following options are available: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( (types.find(keys[i])->second).isFlag() ) printKeyword( keys[i], log );   //log.printKeyword( keys[i], documentation[i] );
     }
     log.printf("\n");
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( (types.find(keys[i])->second).isOptional() ) nkeys++;
  }
  if( nkeys>0 ){
     for(unsigned i=0;i<keys.size();++i){
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
  for(unsigned i=0;i<w.size();++i){
      nl+=w[i].length() + 1;
      if( nl>60 ){
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

bool Keywords::getLogicalDefault( std::string key, bool& def ) const {
   if( booldefs.find(key)!=booldefs.end() ){ 
     def=booldefs.find(key)->second;
     return true;
   } else {
     return false;
   }
}

bool Keywords::getDefaultValue( std::string key, std::string& def ) const {
   plumed_assert( style(key,"compulsory") );

   if( numdefs.find(key)!=numdefs.end() ){
      def=numdefs.find(key)->second;
      return true;
   } else {
      return false;
   }
}

void Keywords::destroyData(){
   keys.clear(); reserved_keys.clear(); types.clear();
   allowmultiple.clear(); documentation.clear(); 
   booldefs.clear(); numdefs.clear();
}

