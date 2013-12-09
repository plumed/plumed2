/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "PlumedMain.h"
#include "tools/Tools.h"
#include <cstring>
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithVirtualAtom.h"
#include "Atoms.h"
#include <set>
#include "config/Config.h"
#include <cstdlib>
#include "ActionRegister.h"
#include "GREX.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "ActionSet.h"
#include "tools/Log.h"
#include "tools/DLLoader.h"
#include "tools/Communicator.h"
#include "CLToolMain.h"
#include "tools/Stopwatch.h"
#include "tools/Citations.h"
#include "ExchangePatterns.h"
#include "tools/IFile.h"

using namespace std;

namespace PLMD{

PlumedMain::PlumedMain():
  comm(*new Communicator),
  multi_sim_comm(*new Communicator),
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
  exchangePatterns(*new(ExchangePatterns)),
  exchangeStep(false),
  restart(false),
  stopFlag(NULL),
  stopNow(false),
  novirial(false),
  detailedTimers(false)
{
  log.link(comm);
  log.setLinePrefix("PLUMED: ");
  stopwatch.start();
  stopwatch.pause();
}

PlumedMain::~PlumedMain(){
  stopwatch.start();
  stopwatch.stop();
  if(initialized) log<<stopwatch;
  delete &exchangePatterns;
  delete &actionSet;
  delete &citations;
  delete &atoms;
  delete &log;
  if(grex)  delete grex;
  delete &stopwatch;
  if(cltool) delete cltool;
  delete &dlloader;
  delete &comm;
  delete &multi_sim_comm;
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
       atoms.setPositions(val);
  } else if(word=="setMasses") {
       CHECK_INIT(initialized,word);
       atoms.setMasses(val);
  } else if(word=="setCharges") {
       CHECK_INIT(initialized,word);
       atoms.setCharges(val);
  } else if(word=="setPositionsX") {
       CHECK_INIT(initialized,word);
       atoms.setPositions(val,0);
  } else if(word=="setPositionsY") {
       CHECK_INIT(initialized,word);
       atoms.setPositions(val,1);
  } else if(word=="setPositionsZ") {
       CHECK_INIT(initialized,word);
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
       atoms.setForces(val);
  } else if(word=="setForcesX") {
       CHECK_INIT(initialized,word);
       atoms.setForces(val,0);
  } else if(word=="setForcesY") {
       CHECK_INIT(initialized,word);
       atoms.setForces(val,1);
  } else if(word=="setForcesZ") {
       CHECK_INIT(initialized,word);
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
       atoms.startStep();
 } else if(word=="setStepLong") {
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       step=(*static_cast<long int*>(val));
       atoms.startStep();
// words used less frequently:
  } else if(word=="setAtomsNlocal"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       atoms.setAtomsNlocal(*static_cast<int*>(val));
  } else if(word=="setAtomsGatindex"){
       CHECK_INIT(initialized,word);
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
       atoms.setMDNaturalUnits(true);
  } else if(word=="setNoVirial"){
       CHECK_NOTINIT(initialized,word);
       novirial=true;
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
  } else if(word=="setMPImultiSimComm"){
       CHECK_NOTINIT(initialized,word);
       multi_sim_comm.Set_comm(val);
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
// other commands that should be used after initialization:
  } else if(word=="setStopFlag"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       stopFlag=static_cast<int*>(val);
  } else if(word=="getExchangesFlag"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangePatterns.getFlag((*static_cast<int*>(val)));
  } else if(word=="setExchangesSeed"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangePatterns.setSeed((*static_cast<int*>(val)));
  } else if(word=="setNumberOfReplicas"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangePatterns.setNofR((*static_cast<int*>(val)));
  } else if(word=="getExchangesList"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       exchangePatterns.getList((static_cast<int*>(val)));
  } else if(word=="runFinalJobs"){  
       CHECK_INIT(initialized,word);
       runJobsAtEndOfCalculation();
  } else if(word=="isEnergyNeeded"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       if(atoms.isEnergyNeeded()) *(static_cast<int*>(val))=1;
       else                       *(static_cast<int*>(val))=0;
  } else if(word=="getBias"){
       CHECK_INIT(initialized,word);
       CHECK_NULL(val,word);
       *(static_cast<double*>(val))=getBias(); 
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
  if(!log.isOpen()) log.link(stdout);
  log<<"PLUMED is starting\n";
  log<<"PLUMED compiled on " __DATE__ " at " __TIME__ "\n";
  log<<"Please cite this paper when using PLUMED ";
  log<<cite("Tribello, Bonomi, Branduardi, Camilloni, and Bussi, Comput. Phys. Commun. 185, 604 (2014)");
  log<<"\n";
  log<<"For further information see the PLUMED web page at http://www.plumed-code.org\n";
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
  IFile ifile;
  ifile.link(*this);
  ifile.open(str);
  std::vector<std::string> words;
  exchangePatterns.setFlag(exchangePatterns.NONE);
  while(Tools::getParsedLine(ifile,words) && words[0]!="ENDPLUMED") readInputWords(words);
  log.printf("END FILE: %s\n",str.c_str());
  log.flush();	

  pilots=actionSet.select<ActionPilot*>();
}

void PlumedMain::readInputWords(const std::vector<std::string> & words){
  plumed_assert(initialized);
  if(words.empty())return;
  else if(words[0]=="ENDPLUMED") return;
  else if(words[0]=="_SET_SUFFIX"){
    plumed_assert(words.size()==2);
    setSuffix(words[1]);
  } else {
    std::vector<std::string> interpreted(words);
    Tools::interpretLabel(interpreted);
    Action* action=actionRegister().create(ActionOptions(*this,interpreted));
    if(!action){
      log<<"ERROR\n";
      log<<"I cannot understand line:";
      for(unsigned i=0;i<interpreted.size();++i) log<<" "<<interpreted[i];
      log<<"\n";
      exit(1);
    };
    action->checkRead();
    actionSet.push_back(action);
  };

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
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive()){
      if((*p)->checkNeedsGradients()) (*p)->setOption("GRADIENTS");
    }
  }

  stopwatch.stop("1 Prepare dependencies");

}

void PlumedMain::shareData(){
// atom positions are shared (but only if there is something to do)
  if(!active)return;
  stopwatch.start("2 Sharing data");
  if(atoms.getNatoms()>0) atoms.share();
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
  if(atoms.getNatoms()>0) atoms.wait();
  stopwatch.stop("3 Waiting for data");
}


void PlumedMain::justCalculate(){
  if(!active)return;
  stopwatch.start("4 Calculating (forward loop)");
  bias=0.0;

  int iaction=0;
// calculate the active actions in order (assuming *backward* dependence)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    std::string actionNumberLabel;
    if(detailedTimers){
      Tools::convert(iaction,actionNumberLabel);
      actionNumberLabel="4A "+actionNumberLabel+" "+(*p)->getLabel();
      stopwatch.start(actionNumberLabel);
    }
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

    if(detailedTimers) stopwatch.stop(actionNumberLabel);
    iaction++;
  }
  stopwatch.stop("4 Calculating (forward loop)");
}

void PlumedMain::justApply(){
  
  if(!active)return;
  int iaction=0;
  stopwatch.start("5 Applying (backward loop)");
// apply them in reverse order
  for(ActionSet::reverse_iterator p=actionSet.rbegin();p!=actionSet.rend();++p){
    if((*p)->isActive()){

      std::string actionNumberLabel;
      if(detailedTimers){
        Tools::convert(iaction,actionNumberLabel);
        actionNumberLabel="5A "+actionNumberLabel+" "+(*p)->getLabel();
        stopwatch.start(actionNumberLabel);
      }

      (*p)->apply();
      ActionAtomistic*a=dynamic_cast<ActionAtomistic*>(*p);
// still ActionAtomistic has a special treatment, since they may need to add forces on atoms
      if(a) a->applyForces();

      if(detailedTimers) stopwatch.stop(actionNumberLabel);
    }
    iaction++;
  }

// this is updating the MD copy of the forces
  if(detailedTimers) stopwatch.start("5B Update forces");
  if(atoms.getNatoms()>0) atoms.updateForces();
  if(detailedTimers) stopwatch.stop("5B Update forces");

  if(detailedTimers) stopwatch.start("5C Update");
// update step (for statistics, etc)
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
    if((*p)->isActive()) (*p)->update();
  }
  if(detailedTimers) stopwatch.stop("5C Update");
// Check that no action has told the calculation to stop
  if(stopNow){
     if(stopFlag) (*stopFlag)=1;
     else plumed_merror("your md code cannot handle plumed stop events - add a call to plumed.comm(stopFlag,stopCondition)");
  }  
  stopwatch.stop("5 Applying (backward loop)");

// flush by default every 10000 steps
// hopefully will not affect performance
  if(step%10000==0){
    fflush();
    log.flush();
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p) (*p)->fflush();
  }
}

void PlumedMain::load(const std::string& ss){
  if(DLLoader::installed()){
     string s=ss;
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
     s=base+"."+config::getSoExt();
     void *p=dlloader.load(s);
     if(!p){
       log<<"ERROR\n";
       log<<"I cannot load library "<<ss<<"\n";
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
  plumed_massert(fp,"file " + ppath + " cannot be found");
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

void PlumedMain::insertFile(FileBase&f){
  files.insert(&f);
}

void PlumedMain::eraseFile(FileBase&f){
  files.erase(&f);
}

void PlumedMain::stop(){ 
  stopNow=true;
}

void PlumedMain::runJobsAtEndOfCalculation(){
  for(ActionSet::iterator p=actionSet.begin();p!=actionSet.end();++p){
      (*p)->runFinalJobs();
  }
} 

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



