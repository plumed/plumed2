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
#include "PlumedMain.h"
#include "Tools.h"
#include <cstring>
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include <set>
#include "PlumedConfig.h"
#include "Colvar.h"
#include <cstdlib>
#include "ActionRegister.h"
#include "GREX.h"
#include "PlumedException.h"
#include "Atoms.h"
#include "ActionSet.h"
#include "Log.h"
#include "DLLoader.h"
#include "PlumedCommunicator.h"
#include "CLToolMain.h"
#include "Stopwatch.h"
#include "Citations.h"

using namespace PLMD;
using namespace std;

PlumedMain::PlumedMain():
  comm(*new PlumedCommunicator),
  dlloader(*new DLLoader),
  cltool(NULL),
  stopwatch(*new Stopwatch),
  grex(NULL),
  initialized(false),
  log(*new Log),
  citations(*new Citations),
  step(0),
  active(false),
  atoms(*new Atoms(*this)),
  actionSet(*new ActionSet(*this)),
  bias(0.0),
  novirial(false)
{
  log.link(comm);
  log.setLinePrefix("PLUMED: ");
  log.link(stdout);
  stopwatch.start();
  stopwatch.pause();
}

PlumedMain::~PlumedMain(){
  stopwatch.start();
  stopwatch.stop();
  if(initialized) log<<stopwatch;
  delete &actionSet;
  delete &citations;
  delete &atoms;
  delete &log;
  if(grex)  delete grex;
  delete &stopwatch;
  if(cltool) delete cltool;
  delete &dlloader;
  delete &comm;
}

/////////////////////////////////////////////////////////////
//  MAIN INTERPRETER

#define CHECK_INIT(ini,word) plumed_massert(ini,"cmd(\"" + word +"\") should be only used after plumed initialization")
#define CHECK_NOTINIT(ini,word) plumed_massert(!(ini),"cmd(\"" + word +"\") should be only used before plumed initialization")
#define CHECK_NULL(val,word) plumed_massert(val,"NULL pointer received in cmd(\"" + word + "\")");

void PlumedMain::cmd(const std::string & word,void*val){

  stopwatch.start();

  if(false){
// for efficiency, words frequently used are checked first

// words used at every MD steps:
  } else if(word=="setBox") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setBox(val);
  } else if(word=="setPositions") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setPositions(val);
  } else if(word=="setMasses") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setMasses(val);
  } else if(word=="setCharges") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setCharges(val);
  } else if(word=="setPositionsX") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setPositions(val,0);
  } else if(word=="setPositionsY") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setPositions(val,1);
  } else if(word=="setPositionsZ") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setPositions(val,2);
  } else if(word=="setVirial") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setVirial(val);
  } else if(word=="setEnergy") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setEnergy(val);
  } else if(word=="setForces") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setForces(val);
  } else if(word=="setForcesX") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setForces(val,0);
  } else if(word=="setForcesY") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setForces(val,1);
  } else if(word=="setForcesZ") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setForces(val,2);
  } else if(word=="calc") {
       CHECK_INIT(initialized,word);
       calc();
  } else if(word=="prepareDependencies") {
       CHECK_INIT(initialized,word);
       prepareDependencies();
  } else if(word=="shareData") {
       CHECK_INIT(initialized,word);
       shareData();
  } else if(word=="prepareCalc") {
       CHECK_INIT(initialized,word);
       prepareCalc();
  } else if(word=="performCalc") {
       CHECK_INIT(initialized,word);
       performCalc();
  } else if(word=="setStep") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       step=(*static_cast<int*>(val));
// words used less frequently:
  } else if(word=="setAtomsNlocal"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setAtomsNlocal(*static_cast<int*>(val));
  } else if(word=="setAtomsGatindex"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setAtomsGatindex(static_cast<int*>(val));
  } else if(word=="setAtomsContiguous"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setAtomsContiguous(*static_cast<int*>(val));
  } else if(word=="createFullList"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.createFullList(static_cast<int*>(val));
  } else if(word=="getFullList"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.getFullList(static_cast<int**>(val));
  } else if(word=="clearFullList"){
       CHECK_INIT(initialized,word);
       atoms.clearFullList();
  } else if(word=="read"){
       CHECK_INIT(initialized,word);
       if(val)readInputFile(static_cast<char*>(val));
       else   readInputFile("plumed.dat");
  } else if(word=="clear"){
       CHECK_INIT(initialized,word);
       actionSet.clearDelete();
  } else if(word=="getApiVersion"){
       CHECK_NULL(val,word);
       *(static_cast<int*>(val))=1;
// commands which can be used only before initialization:
  } else if(word=="init"){
       CHECK_NOTINIT(initialized,word);
       init();
  } else if(word=="setRealPrecision"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setRealPrecision(*static_cast<int*>(val));
  } else if(word=="setMDLengthUnits"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDLengthUnits(d);
  } else if(word=="setMDEnergyUnits"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDEnergyUnits(d);
  } else if(word=="setMDTimeUnits"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       double d;
       atoms.MD2double(val,d);
       atoms.setMDTimeUnits(d);
  } else if(word=="setNaturalUnits"){
// set the boltzman constant for MD in natural units (kb=1)
// only needed in LJ codes if the MD is passing temperatures to plumed (so, not yet...)
// use as cmd("setNaturalUnits")
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setMDNaturalUnits(true);
  } else if(word=="setPlumedDat"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       plumedDat=static_cast<char*>(val);
  } else if(word=="setMPIComm"){
       CHECK_NOTINIT(initialized,word);
       comm.Set_comm(val);
       atoms.setDomainDecomposition(comm);
  } else if(word=="setMPIFComm"){
       CHECK_NOTINIT(initialized,word);
       comm.Set_fcomm(val);
       atoms.setDomainDecomposition(comm);
  } else if(word=="setNatoms"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setNatoms(*static_cast<int*>(val));
  } else if(word=="setTimestep"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setTimeStep(val);
  } else if(word=="setMDEngine"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       MDEngine=static_cast<char*>(val);
  } else if(word=="setLog"){
       CHECK_NOTINIT(initialized,word);
       log.link(static_cast<FILE*>(val));
  } else if(word=="setLogFile"){
       CHECK_NOTINIT(initialized,word);
       CHECK_NULL(val,word);
       log.open(static_cast<char*>(val),"w");
  } else if(word=="getExchangesFlag"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangepatterns.getFlag((*static_cast<int*>(val)));
  } else if(word=="setExchangesSeed"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangepatterns.setSeed((*static_cast<int*>(val)));
  } else if(word=="setNumberOfReplicas"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangepatterns.setNofR((*static_cast<int*>(val)));
  } else if(word=="getExchangesList"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangepatterns.getList((static_cast<int*>(val)));
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
       plumed_massert(grex,"error allocating grex");
       std::string kk=words[1];
       for(unsigned i=2;i<words.size();i++) kk+=" "+words[i];
       grex->cmd(kk.c_str(),val);
     } else if(nw>1 && words[0]=="CLTool"){
       CHECK_NOTINIT(initialized,word);
       if(!cltool) cltool=new CLToolMain;
       std::string kk=words[1];
       for(unsigned i=2;i<words.size();i++) kk+=" "+words[i];
       cltool->cmd(kk.c_str(),val);
     } else{
       plumed_merror("cannot interpret cmd(\"" + word + "\"). check plumed developers manual to see the available commands.");
     };
  };

 stopwatch.pause();
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
  log<<"Relevant bibliography:\n";
  log<<citations;
  log<<"Please read and cite where appropriate!\n";
  log<<"Finished setup\n";
}

void PlumedMain::readInputFile(std::string str){
  plumed_assert(initialized);
  log.printf("FILE: %s\n",str.c_str());
  PlumedIFile ifile;
  ifile.link(*this);
  ifile.open(str,"r");
  std::vector<std::string> words;
  exchangepatterns.setFlag(exchangepatterns.NONE);
  while(Tools::getParsedLine(ifile,words)){
    if(words.empty())continue;
    else if(words[0]=="ENDPLUMED") break;
    else if(words[0]=="LOAD") load(words);
    else if(words[0]=="_SET_SUFFIX"){
      plumed_assert(words.size()==2);
      setSuffix(words[1]);
    }
    else if(words[0]=="RANDOM_EXCHANGES"){
      exchangepatterns.setFlag(exchangepatterns.RANDOM);
      // I convert the seed to -seed because I think it is more general to use a positive seed in input
      if(words.size()>2&&words[1]=="SEED") {int seed; Tools::convert(words[2],seed); exchangepatterns.setSeed(-seed); }
    }
    else if(words[0]=="INCLUDE"){
      plumed_assert(words.size()==2);
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
      action->checkRead();
      actionSet.push_back(action);
    };
  };
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

  stopwatch.start("1 Prepare dependencies");

// activate all the actions which are on step
// activation is recursive and enables also the dependencies
// before doing that, the prepare() method is called to see if there is some
// new/changed dependency (up to now, only useful for dependences on virtual atoms,
// which can be dynamically changed).
//

// First switch off all actions
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
     (*p)->deactivate();
     (*p)->clearOptions();
  }

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
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive()){
      if(Colvar *c=dynamic_cast<Colvar*>(*p)) {
        if(c->checkIsEnergy()) collectEnergy=true;
      }
      if((*p)->checkNeedsGradients()) (*p)->setOption("GRADIENTS");
    }
  }
  atoms.setCollectEnergy(collectEnergy);

  stopwatch.stop("1 Prepare dependencies");

}

void PlumedMain::shareData(){
// atom positions are shared (but only if there is something to do)
  if(!active)return;
  stopwatch.start("2 Sharing data");
  atoms.share();
  stopwatch.stop("2 Sharing data");
}

void PlumedMain::performCalc(){
  waitData();
  justCalculate();
  justApply();
}

void PlumedMain::waitData(){
  if(!active)return;
  stopwatch.start("3 Waiting for data");
  atoms.wait();
  stopwatch.stop("3 Waiting for data");
}


void PlumedMain::justCalculate(){

  stopwatch.start("4 Calculating (forward loop)");
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
      // This retrieves components called bias 
      if(av) bias+=av->getOutputQuantity("bias");
      if(av)av->setGradientsIfNeeded();	
      ActionWithVirtualAtom*avv=dynamic_cast<ActionWithVirtualAtom*>(*p);
      if(avv)avv->setGradientsIfNeeded();	
    }
  }
  stopwatch.stop("4 Calculating (forward loop)");
}

void PlumedMain::justApply(){
  
  stopwatch.start("5 Applying (backward loop)");
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
  stopwatch.stop("5 Applying (backward loop)");
}

void PlumedMain::load(std::vector<std::string> & words){
  if(DLLoader::installed()){
     string s=words[1];
     plumed_assert(words.size()==2);
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
  } else plumed_merror("loading not enabled, please recompile with -D__PLUMED_HAS_DLOPEN");
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
  plumed_massert(fp,"file " + ppath + "cannot be found");
  return fp;
}

int PlumedMain::fclose(FILE*fp){
  return std::fclose(fp);
}

std::string PlumedMain::cite(const std::string&item){
  return citations.cite(item);
}

void PlumedMain::fflush(){
  for(files_iterator p=files.begin();p!=files.end();++p){
    (*p)->flush();
  }
}

void PlumedMain::insertFile(PlumedFileBase&f){
  files.insert(&f);
}

void PlumedMain::eraseFile(PlumedFileBase&f){
  files.erase(&f);
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



