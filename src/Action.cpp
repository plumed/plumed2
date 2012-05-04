#include "Action.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "Log.h"
#include "PlumedException.h"
#include "Atoms.h"
#include "ActionSet.h"
#include <iostream>

using namespace PLMD;

Keywords ActionOptions::emptyKeys;

ActionOptions::ActionOptions(PlumedMain&p,const std::vector<std::string>&l):
plumed(p),
line(l),
keys(emptyKeys)
{
}

ActionOptions::ActionOptions(const ActionOptions&ao,const Keywords&keys):
plumed(ao.plumed),
line(ao.line),
keys(keys)
{
}

void Action::registerKeywords( Keywords& keys ){
  plumed_assert( keys.size()==0 );
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual componets are then refered to using label.component" );
}

Action::Action(const ActionOptions&ao):
  name(ao.line[0]),
  line(ao.line),
  active(false),
  plumed(ao.plumed),
  log(plumed.getLog()),
  comm(plumed.comm),
  keywords(ao.keys)
{
  line.erase(line.begin());
  log.printf("Action %s\n",name.c_str());

  if ( keywords.exists("LABEL") ){ parse("LABEL",label); }

  if(label.length()==0){
    std::string s; Tools::convert(plumed.getActionSet().size(),s);
    label="@"+s;
  }
  plumed_massert(!plumed.getActionSet().selectWithLabel<Action*>(label),
                "label " + label + " has been already used");
  log.printf("  with label %s\n",label.c_str());
}

Action::~Action(){
  plumed_massert(files.size()==0,"some files open in action "+getLabel()+" where not properly closed.");
}

FILE* Action::fopen(const char *path, const char *mode){
  bool write(false);
  for(const char*p=mode;*p;p++) if(*p=='w' || *p=='a' || *p=='+') write=true;
  FILE* fp;
  if(write && comm.Get_rank()!=0) fp=plumed.fopen("/dev/null",mode);
  else      fp=plumed.fopen(path,mode);
  files.insert(fp);
  return fp;
}

int Action::fclose(FILE*fp){
  files.erase(fp);
  return plumed.fclose(fp);
}

void Action::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    std::fflush((*p));
  }
}

void Action::parseFlag(const std::string&key,bool & t){
  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  // Check keyword is a flag
  if(!keywords.style(key,"nohtml")){
     plumed_massert(keywords.style(key,"flag"), "keyword " + key + " is not a flag");
  }

  // Read in the flag otherwise get the default value from the keywords object
  if(!Tools::parseFlag(line,key,t)){
     if( keywords.style(key,"nohtml") ){ 
        t=false; 
     } else if ( !keywords.getLogicalDefault(key,t) ){
        log.printf("ERROR in action %s with label %s : flag %s has no default",name.c_str(),label.c_str(),key.c_str() );
        plumed_assert(0);
     } 
  }
}

void Action::addDependency(Action*action){
  plumed_massert(this->after.count(action)==action->before.count(this),"internal consistency of dependencies");
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
  for(Dependencies::iterator p=after.begin();p!=after.end();++p) (*p)->activate();
  active=true;
}

void Action::setOption(const std::string &s){
// This overloads the action and activate some options  
  options.insert(s);
  for(Dependencies::iterator p=after.begin();p!=after.end();++p) (*p)->setOption(s);
}

void Action::clearOptions(){
// This overloads the action and activate some options  
  options.clear();
}


void Action::clearDependencies(){
  for(Dependencies::iterator p=after.begin();p!=after.end();++p){
    (*p)->before.erase(this);
  };
  this->after.clear();
}

std::string Action::getDocumentation()const{
  return std::string("UNDOCUMENTED ACTION");
}

void Action::checkRead(){
  if(!line.empty()){
    std::string msg="cannot understand the following words from the input line : ";
    for(unsigned i=0;i<line.size();i++) msg = msg + line[i] + ", ";
    error(msg);
  }
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



void Action::exit(int c){
  plumed.exit(c);
}

void Action::calculateNumericalDerivatives( ActionWithValue* a ){
  plumed_merror("if you get here it means that you are trying to use numerical derivatives for a class that does not implement them");
}

void Action::prepare(){
  return;
}

void Action::error( const std::string & msg ){
  log.printf("ERROR in input to action %s with label %s : %s \n \n", name.c_str(), label.c_str(), msg.c_str() );
  keywords.print( log );
  plumed_merror("ERROR in input to action " + name + " with label " + label + " : " + msg );
  this->exit(1);
}

void Action::warning( const std::string & msg ){
  log.printf("WARNING for action %s with label %s : %s \n", name.c_str(), label.c_str(), msg.c_str() ); 
}





