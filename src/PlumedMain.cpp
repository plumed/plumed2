#include "PlumedMain.h"
#include "Tools.h"
#include <cstring>
#include <cassert>
#include "ActionWithValue.h"
#include "ActionWithExternalArguments.h"
#include "Atoms.h"
#include <set>
#include "PlumedConfig.h"
#include "ColvarEnergy.h"

#include <cstdlib>

#include "ActionRegister.h"

using namespace PLMD;

// !!!!!!!!!!!!!!!!!!!!!!    DANGER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// THE FOLLOWING ARE UTILITIES WHICH ARE NECESSARY FOR DYNAMIC LOADING OF THE PLUMED KERNEL:
// This section should be consistent with the Plumed.h file.
// Since the Plumed.h file may be included in host MD codes, **NEVER** MODIFY THE CODE DOWN HERE

/* Holder for plumedmain function pointers */
typedef struct {
  void*(*create)();
  void(*cmd)(void*,const char*,const void*);
  void(*finalize)(void*);
} plumed_plumedmain_function_holder;

extern "C" void*plumedmain_create();
extern "C" void plumedmain_cmd(void*plumed,const char*key,const void*val);
extern "C" void plumedmain_finalize(void*plumed);

void*plumedmain_create(){
  return new PlumedMain;
}

void plumedmain_cmd(void*plumed,const char*key,const void*val){
  assert(plumed);
  static_cast<PlumedMain*>(plumed)->cmd(key,val);
}

void plumedmain_finalize(void*plumed){
  assert(plumed);
  delete static_cast<PlumedMain*>(plumed);
}

extern "C" plumed_plumedmain_function_holder* plumed_kernel_register(const plumed_plumedmain_function_holder*);
extern "C" void* plumed_dlopen(const char*);
extern "C" const char* plumed_dlerror(void);

namespace PLMD{

/// Static object which registers Plumed.
/// This is a static object which, during its construction at startup,
/// registers the pointers to plumedmain_create, plumedmain_cmd and plumedmain_finalize
/// to the plumed_kernel_register function
static class PlumedMainInitializer{
  public:
  PlumedMainInitializer(){
    plumed_plumedmain_function_holder fh={plumedmain_create,plumedmain_cmd,plumedmain_finalize};
    plumed_kernel_register(&fh);
  };
} RegisterMe;

}

// END OF DANGER
////////////////////////////////////////////////////////////


PlumedMain::PlumedMain():
  initialized(false),
  log(comm),
  step(0),
  active(false),
  atoms(*this),
  actionSet((*this)),
  novirial(false)
{
}

/////////////////////////////////////////////////////////////
//  MAIN INTERPRETER

void PlumedMain::cmd(const std::string & word,void*val){

  if(false){
// for efficiency, words frequently used are checked first

// words used at every MD steps:
  } else if(word=="setBox") {
       assert(initialized);
       assert(val);
       atoms.setBox(val);
  } else if(word=="setPositions") {
       assert(initialized);
       assert(val);
       atoms.setPositions(val);
  } else if(word=="setMasses") {
       assert(initialized);
       assert(val);
       atoms.setMasses(val);
  } else if(word=="setCharges") {
       assert(initialized);
       assert(val);
       atoms.setCharges(val);
  } else if(word=="setPositionsX") {
       assert(initialized);
       assert(val);
       atoms.setPositions(val,0);
  } else if(word=="setPositionsY") {
       assert(initialized);
       assert(val);
       atoms.setPositions(val,1);
  } else if(word=="setPositionsZ") {
       assert(initialized);
       assert(val);
       atoms.setPositions(val,2);
  } else if(word=="setVirial") {
       assert(initialized);
       assert(val);
       atoms.setVirial(val);
  } else if(word=="setEnergy") {
       assert(initialized);
       assert(val);
       atoms.setEnergy(val);
  } else if(word=="setForces") {
       assert(initialized);
       assert(val);
       atoms.setForces(val);
  } else if(word=="setForcesX") {
       assert(initialized);
       assert(val);
       atoms.setForces(val,0);
  } else if(word=="setForcesY") {
       assert(initialized);
       assert(val);
       atoms.setForces(val,1);
  } else if(word=="setForcesZ") {
       assert(initialized);
       assert(val);
       atoms.setForces(val,2);
  } else if(word=="calc") {
       assert(initialized);
       calc();
  } else if(word=="prepareDependencies") {
       assert(initialized);
       prepareDependencies();
  } else if(word=="shareData") {
       assert(initialized);
       shareData();
  } else if(word=="prepareCalc") {
       assert(initialized);
       prepareCalc();
  } else if(word=="performCalc") {
       assert(initialized);
       performCalc();
  } else if(word=="setStep") {
       assert(initialized);
       step=(*static_cast<int*>(val));
// words used less frequently:
  } else if(word=="setAtomsNlocal"){
       assert(initialized);
       assert(val);
       atoms.setAtomsNlocal(*static_cast<int*>(val));
  } else if(word=="setAtomsGatindex"){
       assert(initialized);
       assert(val);
       atoms.setAtomsGatindex(static_cast<int*>(val));
  } else if(word=="setAtomsContiguous"){
       assert(initialized);
       assert(val);
       atoms.setAtomsContiguous(*static_cast<int*>(val));
  } else if(word=="createFullList"){
       assert(initialized);
       atoms.createFullList(static_cast<int*>(val));
  } else if(word=="getFullList"){
       assert(initialized);
       atoms.getFullList(static_cast<int**>(val));
  } else if(word=="clearFullList"){
       assert(initialized);
       atoms.clearFullList();
  } else if(word=="read"){
       assert(initialized);
       if(val)readInputFile(static_cast<char*>(val));
       else   readInputFile("plumed.dat");
  } else if(word=="clear"){
       assert(initialized);
       actionSet.clearDelete();
  } else if(word=="getApiVersion"){
       *(static_cast<int*>(val))=1;
// commands which can be used only before initialization:
  } else if(word=="init"){
       assert(!initialized);
       init();
  } else if(word=="setRealPrecision"){
       assert(!initialized);
       assert(val);
       atoms.setRealPrecision(*static_cast<int*>(val));
  } else if(word=="setMDLengthUnits"){
       assert(!initialized);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDLengthUnits(d);
  } else if(word=="setMDEnergyUnits"){
       assert(!initialized);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDEnergyUnits(d);
  } else if(word=="setMDTimeUnits"){
       assert(!initialized);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDTimeUnits(d);
  } else if(word=="setPlumedDat"){
       assert(!initialized);
       plumedDat=static_cast<char*>(val);
  } else if(word=="setMPIComm"){
       assert(!initialized);
       comm.Set_comm(val);
       atoms.setDomainDecomposition(comm);
  } else if(word=="setMPIFComm"){
       assert(!initialized);
       comm.Set_fcomm(val);
       atoms.setDomainDecomposition(comm);
  } else if(word=="setNatoms"){
       assert(!initialized);
       assert(val);
       atoms.setNatoms(*static_cast<int*>(val));
  } else if(word=="setTimestep"){
       assert(!initialized);
       assert(val);
       atoms.setTimeStep(val);
  } else if(word=="setMDEngine"){
       assert(!initialized);
       assert(val);
       assert(MDEngine.length()==0);
       MDEngine=static_cast<char*>(val);
  } else if(word=="setLog"){
       assert(!initialized);
       log.set(static_cast<FILE*>(val));
  } else if(word=="setLogFile"){
       assert(!initialized);
       log.setFile(static_cast<char*>(val));
  } else if(word=="setKBoltzman"){
       assert(val);
//
  } else {
// multi word commands

     std::vector<std::string> words=Tools::getWords(word);
     int nw=words.size();
   
     if(false){
     } else if(nw==2 && words[0]=="checkAction"){
       int check=0;
       if(actionRegister().check(words[1])) check=1;
       *(static_cast<int*>(val))=check;
     } else{
   // error
       fprintf(stderr,"+++ PLUMED ERROR\n");
       fprintf(stderr,"+++ CANNOT INTERPRET CALL TO cmd() ROUTINE WITH ARG %s\n",word.c_str());
       fprintf(stderr,"+++ There might be a mistake in the MD code\n");
       fprintf(stderr,"+++ or you may be using an out-dated plumed version\n");
       exit(1);
     };
  };
}

/////
////////////////////////////////////////////////////////////////////////

void PlumedMain::init(){
// check that initialization just happens once
  initialized=true;
  atoms.init();
  log<<"PLUMED is starting\n";
  log<<"****  THIS IS AN EXPERIMENTAL VERSION ****\n";
  log<<"PLUMED compiled on " __DATE__ " at " __TIME__ "\n";
  log.printf("Please read and cite:\n");
  log.printf("  M. Bonomi, D. Branduardi, G. Bussi, C. Camilloni, D. Provasi, P. Raiteri,\n");
  log.printf("  D. Donadio, F. Marinelli, F. Pietrucci, R. A. Broglia and M. Parrinello\n");
  log.printf("  PLUMED: a portable plugin for free-energy calculations with molecular dynamics\n");
  log.printf("  Comp. Phys. Comm. 180, 1961 (2009)\n");
  log.printf("For further information see the PLUMED web page at www.plumed-code.org\n");
  log.printf("List of registered actions:\n");
  log<<actionRegister();
  log.printf("Molecular dynamics engine: %s\n",MDEngine.c_str());
  log.printf("Precision of reals: %d\n",atoms.getRealPrecision());
  log.printf("Running over %d %s\n",comm.Get_size(),(comm.Get_size()>1?"nodes":"node"));
  log.printf("Number of atoms: %d\n",atoms.getNatoms());
  if(plumedDat.length()>0){
    readInputFile(plumedDat);
    plumedDat="";
  }
  atoms.updateUnits();
  log.printf("Timestep: %f\n",atoms.getTimeStep());
}

void PlumedMain::readInputFile(std::string str){
  assert(initialized);
  log.printf("FILE: %s\n",str.c_str());
  FILE*fp=fopen(str.c_str(),"r");
  std::vector<std::string> words;
  while(Tools::getParsedLine(fp,words)){
    if(words.size()==0)continue;
    else if(words[0]=="ENDPLUMED") break;
    else if(words[0]=="LOAD"){
      std::string s=words[1];
      assert(words.size()==2);
      void *p=plumed_dlopen(s.c_str());
      if(!p){
// try with different extension
        size_t n=s.find_last_of(".");
        if(n==std::string::npos) s+=".";
        else s=s.substr(0,n+1);
        s+=soext;
        p=plumed_dlopen(s.c_str());
      }
      if(!p){
        log<<"ERROR\n";
        log<<"I cannot load library "<<words[1].c_str()<<"\n";
        log<<plumed_dlerror();
        log<<"\n";
        this->exit(1);
      }
      log<<"Loading shared library "<<s.c_str()<<"\n";
      log<<"Here is the new list of available actions\n";
      log<<actionRegister();
    } else if(words[0]=="INCLUDE"){
      assert(words.size()==2);
      readInputFile(words[1]);
      continue;
    } else {
      Tools::interpretLabel(words);
      Action* action=actionRegister().create(ActionOptions(*this,words));
      if(!action){
        log<<"ERROR\n";
        log<<"I cannot understand line:";
        for(unsigned i=0;i<words.size();++i) log<<" "<<words[i];
        log<<"\n";
        exit(1);
      };
      actionSet.push_back(action);
    };
  };
  fclose(fp);
  log.printf("END FILE: %s\n",str.c_str());
}

////////////////////////////////////////////////////////////////////////

void PlumedMain::exit(int c){
  comm.Abort(c);
}

Log& PlumedMain::getLog(){
  return log;
}





void PlumedMain::calc(){
  prepareCalc();
  performCalc();
}

void PlumedMain::prepareCalc(){
  prepareDependencies();
  shareData();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// here we have the main steps in "calc()"
// they can be called individually, but the standard thing is to
// traverse them in this order:
void PlumedMain::prepareDependencies(){

// activate all the actions which are on step
// activation is recursive and enables also the dependencies
// before doing that, the prepare() method is called to see if there is some
// new/changed dependency (up to now, only useful for dependences on virtual atoms,
// which can be dynamically changed).
//
// for optimization, an "active" flag remains false if no action at all is active
  active=false;
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
     if( (*p)->onStep() ){ (*p)->activate(); active=true; }
  };

// also, if one of them is the total energy, tell to atoms that energy should be collected
  bool collectEnergy=false;
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();p++){
    if((*p)->isActive()){
      if(ColvarEnergy* c=dynamic_cast<ColvarEnergy*>(*p)) collectEnergy=true; 
    }
  }
  atoms.setCollectEnergy(collectEnergy);

}

void PlumedMain::shareData(){
// atom positions are shared (but only if there is something to do)
  if(!active)return;
  atoms.share();
}


void PlumedMain::performCalc(){

  if(!active)return;
  atoms.wait();

// calculate the active actions in order (assuming *backward* dependence)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    {
      ActionWithValue* a=dynamic_cast<ActionWithValue*>(*p);
      if(a) a->clearInputForces();
      if(a) a->clearDerivatives();
    }
    {
      ActionWithExternalArguments* a=dynamic_cast<ActionWithExternalArguments*>(*p);
      if(a) a->clearOutputForces();
      if(a) if(a->isActive()) a->retrieveData();
    }
    if((*p)->isActive()){
      if( (*p)->checkNumericalDerivatives() ) (*p)->calculateNumericalDerivatives();
      else (*p)->calculate();
    }
  }
  
// apply them in reverse order
  for(ActionSet::reverse_iterator p=actionSet.rbegin();p!=actionSet.rend();++p){
    if((*p)->isActive()) (*p)->apply();
//    ActionAtomistic*a=dynamic_cast<ActionAtomistic*>(*p);
// still ActionAtomistic has a special treatment, since they may need to add forces on atoms
//    if(a) if(a->isActive()) a->applyForces();
  }

// this is updating the MD copy of the forces
  atoms.updateForces();

// Finally switch off all actions (they will be switched on at next step
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();p++) (*p)->deactivate();

}

// GAT perhaps need to think a bit more how to do this better
void PlumedMain::addCitation( const std::string& citation ){
  cites.push_back( citation );
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



