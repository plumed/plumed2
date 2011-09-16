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
documentation(d)
{
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

  registerKeyword( 1, "LABEL", "a label for the action so that it can be referenced" );
  registerKeyword( 0, "STRIDE", "the frequency with which this particular action should be performed in the absence of dependencies"); 
  line.erase(line.begin());
  log.printf("Action %s\n",name.c_str());
}


void Action::strideKeywordIsCompulsory(){
  for(unsigned i=0;i<keys.size();++i){
      if( keys[i].key=="STRIDE" ){ keys[i].compulsory=1; break; }
  }
}

void Action::readAction(){
  parse("LABEL",label);
  if(label.length()==0){
    std::string s; Tools::convert(plumed.getActionSet().size(),s);
    label="@"+s;
  } 
  if (plumed.getActionSet().selectWithLabel<Action*>(label)) error("there is already an action with label " + label);
  log.printf("  with label %s\n",label.c_str());

  // GAT - eventually move something like this into parse
  bool needstride=false;
  for(unsigned i=0;i<keys.size();++i){
      if( keys[i].key=="STRIDE" && keys[i].compulsory!=0 ){ needstride=true; break; }
  }

  if( needstride || Tools::testForKey(line, "STRIDE=") ){ 
    parse("STRIDE",stride); 
    log.printf("  with stride %d\n",stride);
  }
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

void Action::registerKeyword( const unsigned& compulsory, const std::string& key, const std::string& description ){
  Keyword tmpkey( compulsory, key, description ); keys.push_back( tmpkey );
}

void Action::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    std::fflush((*p));
  }
}

void Action::parseFlag(const std::string&key,bool & t){
  if(!Tools::parseFlag(line,key,t)) error("problem parsing keyword " + key );
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

std::string Action::getDocumentation()const{
  return std::string("UNDOCUMENTED ACTION");
}

void Action::checkRead(){
  if(line.size()>0){
    log.printf("ERROR READING INPUT FILE\n");
    log.printf("I CANNOT UNDERSTAND THE FOLLOWING WORDS:\n");
    for(unsigned i=0;i<line.size();i++) log.printf("  %s\n",line[i].c_str());
    log.printf("STOP!!\n");
    plumed.exit(1);
  }
  keys.resize(0);
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

void Action::error( const std::string& stream ){
	log.printf("ERROR in action %s with label %s : %s\n\n",name.c_str(),label.c_str(),stream.c_str() );
	log.printf( "%s\n\n",getDocumentation().c_str() );
	log.printf( "Compulsory keywords for this action are: \n\n");
        for(unsigned i=0;i<keys.size();++i){
	   if ( keys[i].compulsory==1 ) log.printf(" %s - %s \n",keys[i].key.c_str(), keys[i].documentation.c_str() );
	}
        log.printf( "\n Optional keywords for this action are: \n\n");
        for(unsigned i=0;i<keys.size();++i){
           if ( keys[i].compulsory==0 ) log.printf(" %s - %s \n",keys[i].key.c_str(), keys[i].documentation.c_str() );
        }
	log.printf("\n");
	plumed.exit(1);
}

void Action::warning( const std::string& stream ){
	log.printf("WARNING in action %s with label %s : %s\n",name.c_str(),label.c_str(),stream.c_str() );
}   





