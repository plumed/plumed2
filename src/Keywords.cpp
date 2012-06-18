#include "Keywords.h"
#include "DistributionFunctions.h"
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
  } else if( type=="nohtml" ){
      style=nohtml;
  } else if( type=="hidden" ){
      style=hidden;
  } else {
      plumed_assert(false);    // Invalid keyword specifier
  }
}

void Keywords::reserve( const std::string & t, const std::string & k, const std::string & d ){
  plumed_assert( !exists(k) );  plumed_assert( t!="flag");  // Cannot reserve flags at the moment - if you need it let me know
  plumed_assert( !reserved(k) );
  std::string fd;
  if( t=="numbered" ){
     fd=d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
     reserved_allowmultiple.push_back(true);
     reserved_types.push_back( KeyType("optional") );
  } else {
     fd=d;
     reserved_allowmultiple.push_back(false);
     reserved_types.push_back( KeyType(t) );
  }
  reserved_keys.push_back(k); reserved_documentation.push_back(fd);
}

void Keywords::use( const std::string & k ){
  plumed_assert( reserved(k) );
  for(unsigned i=0;i<reserved_keys.size();++i){
     if(reserved_keys[i]==k){ 
       types.push_back( reserved_types[i] ); keys.push_back( reserved_keys[i] ); 
       documentation.push_back( reserved_documentation[i] );
       allowmultiple.push_back( reserved_allowmultiple[i] );
     }
  }
}

void Keywords::reset_style( const std::string & k, const std::string & style ){
  if( exists(k) ){ 
    for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ){ types[i]=KeyType(style); }
    }
  } else if ( reserved(k) ){
    for(unsigned i=0;i<reserved_keys.size();++i){
       if( reserved_keys[i]==k ){ reserved_types[i]=KeyType(style); }
    }
  } else {
    plumed_assert(0);
  }
}

void Keywords::add( const std::string & t, const std::string & k, const std::string & d ){
  plumed_assert( !exists(k) ); plumed_assert( t!="flag");  // Use addFlag to add flags
  plumed_assert( !reserved(k) );
  std::string fd;
  if( t=="numbered" ){
     fd=d + " You can use multiple instances of this keyword i.e. " + k +"1, " + k + "2, " + k + "3...";
     allowmultiple.push_back(true);
     types.push_back( KeyType("optional") );
  } else { 
     fd=d;
     allowmultiple.push_back(false);
     types.push_back( KeyType(t) );
  }
  keys.push_back(k); documentation.push_back(fd); 
}

void Keywords::add( const std::string & t, const std::string & k, const std::string &  def, const std::string & d ){
  plumed_assert( !exists(k) ); plumed_assert( t=="compulsory" ); // An optional keyword can't have a default
  plumed_assert( !reserved(k) );
  types.push_back( KeyType(t) ); keys.push_back(k); 
  documentation.push_back( "( default=" + def + " ) " + d ); 
  allowmultiple.push_back(false);
  numdefs.insert( std::pair<std::string,std::string>(k,def) );
} 

void Keywords::addFlag( const std::string & k, const bool def, const std::string & d ){
  plumed_assert( !exists(k) ); std::string defstr, flag="flag";
  plumed_assert( !reserved(k) );
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.push_back( KeyType(flag) ); 
  keys.push_back(k); 
  documentation.push_back( defstr + d ); //defaults.push_back(def);
  allowmultiple.push_back(false); 
  booldefs.insert( std::pair<std::string,bool>(k,def) );
} 

void Keywords::remove( const std::string & k ){
  bool found=false; unsigned j=0;

  while(true){
    for(j=0;j<keys.size();j++) if(keys[j]==k)break;
    if(j<keys.size()){
      types.erase(types.begin()+j);
      keys.erase(keys.begin()+j);
      allowmultiple.erase(allowmultiple.begin()+j); 
      documentation.erase(documentation.begin()+j);
      found=true;
    } else break;
  }
  plumed_massert(found,"You are trying to forbid " + k + " a keyword that isn't there"); // You have tried to forbid a keyword that isn't there
}

bool Keywords::numbered( const std::string & k ) const {
  unsigned j=0; 
  while(true){
    for(j=0;j<keys.size();j++) if(keys[j]==k) break;
    if( j<keys.size() ){
        return allowmultiple[j];
    } else break;
  }
  plumed_massert(0,"Did not find keyword " + k );
  return false;
}

void Keywords::clear() {
  types.clear(); keys.clear(); documentation.clear(); allowmultiple.clear();   //defaults.clear();
}

bool Keywords::style( const std::string & k, const std::string & t ) const {
  plumed_assert( exists(k) );

  if( t=="compulsory" ){
     for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isCompulsory();
     }
     plumed_assert(false);
  } else if( t=="flag" ){
     for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isFlag();
     }
     plumed_assert(false);
  } else if( t=="optional" ){
     for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isOptional();
     }
     plumed_assert(false);
  } else if( t=="atoms" ){
     for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isAtomList();
     }
     plumed_assert(false);
  } else if( t=="nohtml" ){
     for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isNoHTML();
     }
     plumed_assert(false);
  } else if( t=="hidden" ){
      for(unsigned i=0;i<keys.size();++i){
        if( keys[i]==k ) return types[i].isHidden();
     }
     plumed_assert(false);
  } else {
     plumed_assert(false);
  }
// this is to avoid warnings:
  return false;
}

unsigned Keywords::size() const {
  plumed_assert( keys.size()==documentation.size() );
  plumed_assert( keys.size()==types.size() );
  plumed_assert(  keys.size()==allowmultiple.size() );
  return keys.size();
}

bool Keywords::exists( const std::string & k ) const {
  plumed_massert( keys.size()==documentation.size(), "documentation doesn't match keys" ); 
  plumed_massert( keys.size()==types.size(), "types doesn't match keys" );
  plumed_assert(  keys.size()==allowmultiple.size() );

  for(unsigned i=0;i<keys.size();++i){
     if( keys[i]==k ) return true;
  }
  return false;
}

bool Keywords::reserved( const std::string & k ) const {
  plumed_massert( reserved_keys.size()==reserved_documentation.size(), "documentation doesn't match keys" );
  plumed_massert( reserved_keys.size()==reserved_types.size(), "types doesn't match keys" );
  plumed_massert( reserved_keys.size()==reserved_allowmultiple.size(),"allowmultiple doesn't match keys" );

  for(unsigned i=0;i<reserved_keys.size();++i){
     if( reserved_keys[i]==k ) return true;
  }
  return false;
}

void Keywords::print_html() const {
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isAtomList() ) nkeys++;
  }
  if( nkeys>0 ){
    std::cout<<"\\par Specifying the atoms involved\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isAtomList() ) print_html_item( i );
    }
    std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<< "\\par Compulsory keywords\n\n";
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isCompulsory() ) print_html_item( i );
     }
     std::cout<<"</table>\n\n";
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isFlag() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<<"\\par Options\n\n";
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isFlag() ) print_html_item( i );
     }
     std::cout<<"\n";
  }
  std::cout<<"</table>\n\n";
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isOptional() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isOptional() ) print_html_item( i );
     }
     std::cout<<"\n";
  }
  std::cout<<"</table>\n\n";
}

void Keywords::print( Log& log ) const {
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isAtomList() ) nkeys++;
  }
  if (nkeys>0 ){
     log.printf( "The input for this keyword can be specified using one of the following \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isAtomList() ) printKeyword( i, log );   //log.printKeyword( keys[i], documentation[i] );
     }
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The compulsory keywords for this action are: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isCompulsory() ) printKeyword( i, log );   //log.printKeyword( keys[i], documentation[i] );
     }
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isFlag() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The following options are available: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isFlag() ) printKeyword( i, log );   //log.printKeyword( keys[i], documentation[i] );
     }
     log.printf("\n");
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isOptional() || types[i].isNoHTML() ) nkeys++;
  }
  if( nkeys>0 ){
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isOptional() || types[i].isNoHTML() ) printKeyword( i, log );   //log.printKeyword( keys[i], documentation[i] );
     }
     log.printf("\n");
  }
}

void Keywords::printKeyword( const unsigned& j, Log& log ) const {
  bool killdot=( documentation[j].find("\\f$")!=std::string::npos ); // Check for latex
  std::vector<std::string> w=Tools::getWords( documentation[j] );
  log.printf("%23s - ", keys[j].c_str() );
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

void Keywords::print_html_item( const unsigned& j ) const {
  printf("<tr>\n");
  printf("<td width=15%%> <b> %s </b></td>\n",keys[j].c_str() );
  printf("<td> %s </td>\n",documentation[j].c_str() );
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

