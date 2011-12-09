#include "PlumedMain.h"
#include "Tools.h"
#include <cstring>
#include <cassert>
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "Atoms.h"
#include <set>
#include "PlumedConfig.h"
#include "Colvar.h"
#include <cstdlib>
#include "ActionRegister.h"
#include "GREX.h"

using namespace PLMD;
using namespace std;

PlumedMain::PlumedMain():
  initialized(false),
  grex(NULL),
  log(comm),
  step(0),
  active(false),
  atoms(*this),
  actionSet((*this)),
  bias(0.0),
  novirial(false)
{}

PlumedMain::~PlumedMain(){
  if(grex) delete grex;
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
  } else if(word=="setNaturalUnits"){
// set the boltzman constant for MD in natural units (kb=1)
// only needed in LJ codes if the MD is passing temperatures to plumed (so, not yet...)
// use as cmd("setNaturalUnits")
       assert(!initialized);
       assert(!val);
       atoms.setMDNaturalUnits(true);
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
  } else {
// multi word commands

     std::vector<std::string> words=Tools::getWords(word);
     int nw=words.size();
   
     if(false){
     } else if(nw==2 && words[0]=="checkAction"){
       int check=0;
       if(actionRegister().check(words[1])) check=1;
       *(static_cast<int*>(val))=check;
     } else if(nw>1 && words[0]=="GREX"){
       if(!grex) grex=new GREX(*this);
       assert(grex);
       std::string kk=words[1];
       for(int i=2;i<words.size();i++) kk+=" "+words[i];
       grex->cmd(kk.c_str(),val);
     } else{
   // error
       fprintf(stderr,"+++ PLUMED ERROR\n");
       fprintf(stderr,"+++ CANNOT INTERPRET CALL TO cmd() ROUTINE WITH ARG '%s'\n",word.c_str());
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
  if(grex) log.printf("GROMACS-like replica exchange is on\n");
  log.printf("File suffix: %s\n",getSuffix().c_str());
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
    else if(words[0]=="LOAD") load(words);
    else if(words[0]=="_SET_SUFFIX"){
      assert(words.size()==2);
      setSuffix(words[1]);
    }
    else if(words[0]=="INCLUDE"){
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

  pilots=actionSet.select<ActionPilot*>();
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

// First switch off all actions
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();p++) (*p)->deactivate();

// for optimization, an "active" flag remains false if no action at all is active
  active=false;
  for(unsigned i=0;i<pilots.size();++i){
    if(pilots[i]->onStep()){
      pilots[i]->activate();
      active=true;
     }
  };

// also, if one of them is the total energy, tell to atoms that energy should be collected
  bool collectEnergy=false;
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();p++){
    if((*p)->isActive()){
      if(Colvar *c=dynamic_cast<Colvar*>(*p)) {
        if(c->checkIsEnergy()) collectEnergy=true;
      }
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
  waitData();
  justCalculate();
  justApply();
}

void PlumedMain::waitData(){
  if(!active)return;
  atoms.wait();
}


void PlumedMain::justCalculate(){

  bias=0.0;

// calculate the active actions in order (assuming *backward* dependence)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    ActionWithValue*av=dynamic_cast<ActionWithValue*>(*p);
    ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(*p);
    {
      if(av) av->clearInputForces();
      if(av) av->clearDerivatives();
    }
    {
      if(aa) aa->clearOutputForces();
      if(aa) if(aa->isActive()) aa->retrieveAtoms();
    }
    if((*p)->isActive()){
      if((*p)->checkNumericalDerivatives()) (*p)->calculateNumericalDerivatives();
      else (*p)->calculate();
      if(av)for(int i=0;i<av->getNumberOfValues();++i){
        if(av->getValue(i)->getName()=="bias") bias+=av->getValue(i)->get();
      }
    }
  }
}

void PlumedMain::justApply(){
  
// apply them in reverse order
  for(ActionSet::reverse_iterator p=actionSet.rbegin();p!=actionSet.rend();++p){
    if((*p)->isActive()) (*p)->apply();
    ActionAtomistic*a=dynamic_cast<ActionAtomistic*>(*p);
// still ActionAtomistic has a special treatment, since they may need to add forces on atoms
    if(a) if(a->isActive()) a->applyForces();
  }

// this is updating the MD copy of the forces
  atoms.updateForces();

// update step (for statistics, etc)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive()) (*p)->update();
  }
}

void PlumedMain::load(std::vector<std::string> & words){
  if(DLLoader::installed()){
     string s=words[1];
     assert(words.size()==2);
     size_t n=s.find_last_of(".");
     string extension="";
     string base=s;
     if(n!=std::string::npos && n<s.length()-1) extension=s.substr(n+1);
     if(n!=std::string::npos && n<s.length())   base=s.substr(0,n);
     if(extension=="cpp"){
       string cmd="plumed mklib "+s;
       log<<"Executing: "<<cmd;
       if(comm.Get_size()>0) log<<" (only on master node)";
       log<<"\n";
       if(comm.Get_rank()==0) system(cmd.c_str());
       comm.Barrier();
       base="./"+base;
     }
     s=base+"."+soext;
     void *p=dlloader.load(s);
     if(!p){
       log<<"ERROR\n";
       log<<"I cannot load library "<<words[1].c_str()<<"\n";
       log<<dlloader.error();
       log<<"\n";
       this->exit(1);
     }
     log<<"Loading shared library "<<s.c_str()<<"\n";
     log<<"Here is the new list of available actions\n";
     log<<actionRegister();
  } else assert(0); // Loading not enabled; please recompile with -D__PLUMED_HAS_DLOPEN
}

double PlumedMain::getBias() const{
  return bias;
}

FILE* PlumedMain::fopen(const char *path, const char *mode){
  std::string mmode(mode);
  std::string ppath(path);
  std::string suffix(getSuffix());
  std::string ppathsuf=ppath+suffix;
  FILE*fp=std::fopen(const_cast<char*>(ppathsuf.c_str()),const_cast<char*>(mmode.c_str()));
  if(!fp) fp=std::fopen(const_cast<char*>(ppath.c_str()),const_cast<char*>(mmode.c_str()));
  assert(fp);
  return fp;
}

int PlumedMain::fclose(FILE*fp){
  return std::fclose(fp);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



