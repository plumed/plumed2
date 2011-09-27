#include "Action.h"
#include "PlumedMain.h"
#include "Log.h"
#include <cassert>

using namespace PLMD;

ActionOptions::ActionOptions(PlumedMain&p,const std::vector<std::string>&l):
plumed(p),
line(l)
{
}

Keyword::Keyword( const unsigned& c, const std::string& k, const std::string& d ):
compulsory(c),
key(k),
documentation(d),
forbidden(false)
{
  assert( key.length()<23 );  // If you find yourself here then you have created a keyword that is too long for the formatting
  if( compulsory==2 ) forbidden=true;
}

void Keyword::print( Log& log, const unsigned len ) const {
  unsigned nlines; nlines=floor( documentation.length() / 60 );
  if ( nlines==0 ){
     log.printf("%23s - %-60s \n", key.c_str(), documentation.c_str() );
  } else {
     std::vector<std::string> w=Tools::getWords( documentation );
     std::vector<unsigned> lens( nlines + 1 ); 
     unsigned ll=1, nl=0;
     for(unsigned i=0;i<w.size();++i){
        nl+=w[i].length() + 1;
        if( nl>=ll*60 ){ lens[ll]=nl-1-w[i].length(); ll++; } 
     } 

     log.printf("%23s - %-60s \n", key.c_str(), documentation.substr(0,lens[1]).c_str() );
     std::string blank=" ";
     for(unsigned i=1;i<nlines;++i){
        log.printf("%23s   %-60s \n", blank.c_str(), documentation.substr(lens[i],lens[i+1]-lens[i]).c_str() );
     }
     log.printf("%23s   %-60s  \n", blank.c_str(), documentation.substr(lens[nlines]).c_str() );
  }
} 

Action::Action(const ActionOptions&ao):
  name(ao.line[0]),
  line(ao.line),
  active(false),
  stride(0),
  plumed(ao.plumed),
  log(plumed.getLog()),
  comm(plumed.comm)
{
  registerKeyword( 0, "LABEL", "a label for the action so that it can be referenced" );
  registerKeyword( 1, "STRIDE", "the frequency with which this particular action should be performed"); 
  line.erase(line.begin());
  log.printf("Action %s\n",name.c_str());
}

void Action::readAction(){
  parse("LABEL",label);
  if(label.length()==0){
    std::string s; Tools::convert(static_cast<int>( plumed.getActionSet().size() ),s);
    label="@"+s;
  } 
  if (plumed.getActionSet().selectWithLabel<Action*>(label)) error("there is already an action with label " + label);
  log.printf("  with label %s\n",label.c_str());

  parse("STRIDE",stride);
  if(stride!=0) log.printf("  with stride %d\n",stride);
}

FILE* Action::fopen(const char *path, const char *mode){
  FILE*fp=std::fopen(const_cast<char*>(path),const_cast<char*>(mode));
  files.insert(fp);
  return fp;
}

int Action::fclose(FILE*fp){
  files.erase(fp);
  return std::fclose(fp);
}

void Action::addCitation( const std::string& citation ){
  cites.push_back(citation);
  plumed.addCitation( citation);
}

void Action::registerKeyword( const unsigned& compulsory, const std::string& key, const std::string& description ){
  for(unsigned i=0;i<keys.size();++i){
     if( keys[i].compulsory>1 && key.find(keys[i].key)!=std::string::npos ){
         assert(false);  // A keyword is already in use
     } else if ( key==keys[i].key ){
         assert(false);  // A keyword is already in use
     }
  }
  Keyword tmpkey( compulsory, key, description ); keys.push_back( tmpkey );
}

void Action::forbidKeyword( const std::string& key ){
  bool removed=false;
  for(unsigned i=0;i<keys.size();++i){
     if( key==keys[i].key ){ keys[i].forbidden=true; removed=true; break; }
  }
  assert(removed); // You can't forbid something that isn't there
}

void Action::allowKeyword( const std::string& key ){
  bool allowed=false;
  for(unsigned i=0;i<keys.size();++i){
     if( key==keys[i].key ){ keys[i].forbidden=false; allowed=true; break; }
  }
  assert(allowed); // You can't allow something that isn't there
}

void Action::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    std::fflush((*p));
  }
}

void Action::parseFlag(const std::string&key,bool & t){
  bool forbid, found=false;
  for(unsigned i=0;i<keys.size();++i){
     if ( key==keys[i].key ){
          if(keys[i].compulsory!=0 ) error("a flag cannot be a compulsary keyword register with 0");
          forbid=keys[i].forbidden; found=true;
     }
  }
  assert(found);
  Tools::parseFlag(line,key,t);
  if( t && forbid ) error("keyword " + key + " is forbidden for this action");
}

bool Action::testForKey( const std::string& key) const {
   bool forbid, found=false;
   for( unsigned i=0;i<keys.size();++i){
       // Test for key does a search for the pattern and should only be used with 
       // keywords that are written like ATOMS0 ATOMS1 ATOMS2...
       if( key==keys[i].key ){ assert( keys[i].compulsory>1 ); forbid=keys[i].forbidden; found=true; }
       // If you want to do a parse for an optinal keyword just use parse. 
       // If you have registered the keyword properly parse is clever enough to
       // work out what is optional and what is compulsory.
   }
   assert(found);
   
   if( forbid ){
      for(unsigned i=0;i<line.size();++i){
          if( line[i].find(key)!=std::string::npos ) error("keyword " + key + " is forbidden for this action");
      }
   } else {
      for(unsigned i=0;i<line.size();++i){
          if( line[i].find(key)!=std::string::npos ) return true;
      }
   }
   return false;
}  
   
bool Action::testForNumberedKeys( const std::string& key) const {
   bool forbid, found=false;
   for( unsigned i=0;i<keys.size();++i){
       // Test for key does a search for the pattern and should only be used with 
       // keywords that are written like ATOMS0 ATOMS1 ATOMS2...
       if( key==keys[i].key ){ assert( keys[i].compulsory>1 ); forbid=keys[i].forbidden; found=true; }
       // If you want to do a parse for an optinal keyword just use parse. 
       // If you have registered the keyword properly parse is clever enough to
       // work out what is optional and what is compulsory.
   }   
   assert(found);
   
   if( forbid ){
      for(unsigned i=0;i<line.size();++i){
          if( line[i].find(key)!=std::string::npos ) error("keyword " + key + " is forbidden for this action");
      }
   } else {
      for(unsigned i=0;i<line.size();++i){
          if( line[i]==(key+"1") ) return true;
      } 
   }     
   return false;
}  

void Action::addDependency(Action*action){
  assert(this->after.count(action)==action->before.count(this));
  this->after.insert(action);
  action->before.insert(this);
}

void Action::activate(){
// preparation step is called only the first time an Action is activated.
// since it could change its dependences (e.g. in an ActionAtomistic which is
// accessing to a virtual atom), this is done just before dependencies are
// activated
  if(!active){
    this->unlockRequests();
    prepare();
    this->lockRequests();
  }
  for(Dependencies::iterator p=after.begin();p!=after.end();p++) (*p)->activate();
  active=true;
}

void Action::clearDependencies(){
  for(Dependencies::iterator p=after.begin();p!=after.end();p++){
    (*p)->before.erase(this);
  };
  this->after.clear();
}

void Action::checkRead(){
  if(line.size()>0){
    std::string shit; for(unsigned i=0;i<line.size();++i) shit=shit+line[i];
    error("found the following unregistered keywords on the input line : " + shit);
  }
  keys.clear(); cites.clear();
}

int Action::getStep()const{
  return plumed.getStep();
}

double Action::getTime()const{
  return plumed.getAtoms().getTimeStep()*getStep();
}

double Action::getTimeStep()const{
  return plumed.getAtoms().getTimeStep();
}

void Action::calculateNumericalDerivatives(){
  assert(0);
}

void Action::prepare(){
  return;
}

void Action::error( const std::string& stream ) const {
  log.printf("ERROR in action %s with label %s : %s\n\n",name.c_str(),label.c_str(),stream.c_str() );
  log.printf( "The input for this keyword can be specified using one of the following \n\n");
  for(unsigned i=0;i<keys.size();++i){
     if ( keys[i].compulsory==2 && !keys[i].forbidden ) keys[i].print( log, 60 ); 
  }
  log.printf( "\n The compulsory keywords for this action are: \n\n");
  for(unsigned i=0;i<keys.size();++i){
     if ( keys[i].compulsory==1 && !keys[i].forbidden ) keys[i].print( log, 60 );  
  }
  log.printf( "\n The optional keywords for this action are: \n\n");
  for(unsigned i=0;i<keys.size();++i){
     if ( keys[i].compulsory==0 && !keys[i].forbidden ) keys[i].print( log, 60 );  
  }
  log.printf("\n");
  log.printf("\n The following optional keywords may be use multiple times by using <keyword>1, <keyword>2, ... \n or just once by using <keyword>\n\n");
  for(unsigned i=0;i<keys.size();++i){
     if ( keys[i].compulsory==3 && !keys[i].forbidden ) keys[i].print( log, 60 );   
  }
  plumed.exit(1);
}

void Action::warning( const std::string& stream ) const {
  log.printf("WARNING in action %s with label %s : %s\n",name.c_str(),label.c_str(),stream.c_str() );
}   





