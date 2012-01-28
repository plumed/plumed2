#include "Keywords.h"

using namespace PLMD;

KeyType::KeyType( const std::string& type ){
  if( type=="compulsory" ){
      style=compulsory;
  } else if( type=="flag" ){
      style=flag;
  } else if( type=="optional" ){
      style=optional;
  } else if( type=="input" ){
      style=input;
  } else if( type=="numbered" ){
      style=numbered;
  } else if( type=="modifier" ){
      style=modifier;
  } else if( type=="hidden" ){
      style=hidden;
  } else {
      plumed_assert(false);    // Invalid keyword specifier
  }
}

void Keywords::add( const std::string t, const std::string k, const std::string d ){
  plumed_assert( !exists(k) ); plumed_assert( t!="flag");  // Use addFlag to add flags
  types.push_back( KeyType(t) ); keys.push_back(k); documentation.push_back(d); 
  //defaults.push_back(false);
}

void Keywords::add( const std::string t, const std::string k, const std::string  def, const std::string d ){
  plumed_assert( !exists(k) ); plumed_assert( t=="compulsory" ); // An optional keyword can't have a default
  types.push_back( KeyType(t) ); keys.push_back(k); documentation.push_back( "( default=" + def + " ) " + d ); 
  numdefs.insert( std::pair<std::string,std::string>(k,def) );
} 

void Keywords::addFlag( const std::string k, const bool def, const std::string d ){
  plumed_assert( !exists(k) ); std::string defstr, flag="flag";
  if( def ) { defstr="( default=on ) "; } else { defstr="( default=off ) "; }
  types.push_back( KeyType(flag) ); keys.push_back(k); documentation.push_back( defstr + d ); //defaults.push_back(def);
  booldefs.insert( std::pair<std::string,bool>(k,def) );
} 

void Keywords::remove( const std::string k ){
  bool found=false; unsigned j=0;

  while(true){
    for(j=0;j<keys.size();j++) if(keys[j]==k)break;
    if(j<keys.size()){
      types.erase(types.begin()+j);
      keys.erase(keys.begin()+j);
      documentation.erase(documentation.begin()+j);
      //defaults.erase(defaults.begin()+j);
      found=true;
    } else break;
  }
  plumed_assert(found); // You have tried to forbid a keyword that isn't there
}

void Keywords::clear() {
  types.clear(); keys.clear(); documentation.clear(); //defaults.clear();
}

KeyType Keywords::style( const std::string k ) const {
  plumed_assert( exists(k) );
  for(unsigned i=0;i<keys.size();++i){
     if( keys[i]==k ) return types[i];
  }
  plumed_assert(false);
}

unsigned Keywords::size() const {
  plumed_assert( keys.size()==documentation.size() );
  plumed_assert( keys.size()==types.size() );
//  plumed_assert( keys.size()==defaults.size() );
  return keys.size();
}

bool Keywords::exists( const std::string k ) const {
  plumed_assert( keys.size()==documentation.size() ); 
  plumed_assert( keys.size()==types.size() );
//  plumed_assert( keys.size()==defaults.size() );

  for(unsigned i=0;i<keys.size();++i){
     if( keys[i]==k ) return true;
  }
  return false;
}

void Keywords::print_html() const {
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isInput() ) nkeys++;
  }
  if( nkeys>0 ){
    std::cout<<"\\par Specifying the atoms involved\n\n";
    std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
    for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isInput() ) print_html_item( i );
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
     if ( types[i].isNumbered() ) nkeys++;
  }
  if( nkeys>0 ){
     std::cout<<" <table align=center frame=void width=95%% cellpadding=5%%> \n";
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isNumbered() ) print_html_item( i );
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
  log.printf( "The input for this keyword can be specified using one of the following \n\n");
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isInput() ) log.printKeyword( keys[i], documentation[i] );
  }
  unsigned nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isCompulsory() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The compulsory keywords for this action are: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isCompulsory() ) log.printKeyword( keys[i], documentation[i] );
     }
  }
  // This is a special option specifically for steered MD
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isNumbered() ) nkeys++;
  }
  if( nkeys>0 ){
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isNumbered() ) log.printKeyword( keys[i], documentation[i] );
     }
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isFlag() ) nkeys++;
  }
  if( nkeys>0 ){
     log.printf( "\n The following options are available: \n\n");
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isFlag() ) log.printKeyword( keys[i], documentation[i] );
     }
     log.printf("\n");
  }
  nkeys=0;
  for(unsigned i=0;i<keys.size();++i){
     if ( types[i].isOptional() ) nkeys++;
  }
  if( nkeys>0 ){
     for(unsigned i=0;i<keys.size();++i){
        if ( types[i].isOptional() ) log.printKeyword( keys[i], documentation[i] );
     }
     log.printf("\n");
  }
}

void Keywords::print_html_item( const unsigned& j ) const {
  printf("<tr>\n");
  printf("<td width=15%%> <b> %s </b></td>\n",keys[j].c_str() );
  printf("<td> %s </td>\n",documentation[j].c_str() );
  printf("</tr>\n");
}

