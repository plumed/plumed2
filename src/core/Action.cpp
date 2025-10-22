/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "Action.h"
#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithArguments.h"
#include "ActionWithVirtualAtom.h"
#include "ActionForInterface.h"
#include "DomainDecomposition.h"
#include "PbcAction.h"
#include "ActionToPutData.h"
#include "ActionToGetData.h"
#include "PlumedMain.h"
#include "tools/Log.h"
#include "tools/Exception.h"
#include "tools/Communicator.h"
#include "ActionSet.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
namespace PLMD {

Keywords ActionOptions::emptyKeys;

ActionOptions::ActionOptions(PlumedMain&p,const std::vector<std::string>&l):
  plumed(p),
  line(l),
  keys(emptyKeys) {
}

ActionOptions::ActionOptions(const ActionOptions&ao,const Keywords&mykeys):
  plumed(ao.plumed),
  line(ao.line),
  keys(mykeys) {
}

void Action::registerKeywords( Keywords& keys ) {
  plumed_assert( keys.size()==0 );
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual components are then referred to using label.component" );
  keys.reserve("optional","UPDATE_FROM","Only update this action from this time");
  keys.reserve("optional","UPDATE_UNTIL","Only update this action until this time");
  keys.reserve("optional","RESTART","allows per-action setting of restart (YES/NO/AUTO)");
}

Action::Action(const ActionOptions&ao):
  actionName(ao.line[0]),
  //skipping the name with the +1
  linemap(ao.line.begin()+1,ao.line.end()),
  update_from(std::numeric_limits<double>::max()),
  update_until(std::numeric_limits<double>::max()),
  timestep(0),
  active(false),
  restart(ao.plumed.getRestart()),
  doCheckPoint(ao.plumed.getCPT()),
  never_activate(actionName=="CONSTANT"),
  plumed(ao.plumed),
  log(plumed.getLog()),
  comm(plumed.comm),
  multi_sim_comm(plumed.multi_sim_comm),
  keywords(ao.keys) {
  // Retrieve the timestep and save it
  resetStoredTimestep();

  //line.erase(line.begin());
  //making the line map
  if( !keywords.exists("NO_ACTION_LOG") ) {
    log.printf("Action %s\n",actionName.c_str());
    if(ao.fullPath.length()>0) {
      log<<"  from library: "<<ao.fullPath<<"\n";
    }
  }

  if(comm.Get_rank()==0) {
    replica_index=multi_sim_comm.Get_rank();
  }
  comm.Bcast(replica_index,0);

  if ( keywords.exists("LABEL") ) {
    parse("LABEL",actionLabel);
  }
  if(actionLabel.length()==0) {
    std::string s;
    Tools::convert(plumed.getActionSet().size()-plumed.getActionSet().select<ActionForInterface*>().size(),s);
    actionLabel="@"+s;
  } else if ( actionLabel.find(".")!=std::string::npos ) {
    warning("using full stop in an action label should be avaoided as . has a special meaning in PLUMED action labels");
  }
  if( plumed.getActionSet().selectWithLabel<Action*>(actionLabel) ) {
    error("label " + actionLabel + " has been already used");
  }
  if( !keywords.exists("NO_ACTION_LOG") ) {
    log.printf("  with label %s\n",actionLabel.c_str());
  }
  if ( keywords.exists("UPDATE_FROM") ) {
    parse("UPDATE_FROM",update_from);
  }
  if( !keywords.exists("NO_ACTION_LOG") && update_from!=std::numeric_limits<double>::max()) {
    log.printf("  only update from time %f\n",update_from);
  }
  if ( keywords.exists("UPDATE_UNTIL") ) {
    parse("UPDATE_UNTIL",update_until);
  }
  if( !keywords.exists("NO_ACTION_LOG") && update_until!=std::numeric_limits<double>::max()) {
    log.printf("  only update until time %f\n",update_until);
  }
  if ( keywords.exists("RESTART") ) {
    std::string srestart="AUTO";
    parse("RESTART",srestart);
    if( plumed.parseOnlyMode() ) {
      restart=false;
    } else if(srestart=="YES") {
      restart=true;
    } else if(srestart=="NO") {
      restart=false;
    } else if(srestart=="AUTO") {
      // do nothing, this is the default
    } else {
      error("RESTART should be either YES, NO, or AUTO");
    }
  }
}

void Action::resetStoredTimestep() {
  ActionWithValue* ts = plumed.getActionSet().selectWithLabel<ActionWithValue*>("timestep");
  if( ts ) {
    timestep = (ts->copyOutput(0))->get();
  }
}

Action::~Action() {
  if(files.size()!=0) {
    std::cerr<<"WARNING: some files open in action "+getLabel()+" where not properly closed. This could lead to data loss!!\n";
  }
}

FILE* Action::fopen(const char *path, const char *mode) {
  bool write(false);
  for(const char*p=mode; *p; p++)
    if(*p=='w' || *p=='a' || *p=='+') {
      write=true;
    }
  FILE* fp;
  if(write && comm.Get_rank()!=0) {
    fp=plumed.fopen("/dev/null",mode);
  } else {
    fp=plumed.fopen(path,mode);
  }
  files.insert(fp);
  return fp;
}

int Action::fclose(FILE*fp) {
  files.erase(fp);
  return plumed.fclose(fp);
}

void Action::fflush() {
  for(const auto & p : files) {
    std::fflush(p);
  }
}

std::string Action::getKeyword(const std::string& key) {
  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");

  std::string outkey=linemap.getKeyword(key);
  if(!outkey.empty()) {
    return outkey;
  }
  if( keywords.style(key,"compulsory") ) {
    if( keywords.getDefaultValue(key,outkey) ) {
      if( outkey.length()==0 ) {
        error("keyword " + key + " has weird default value");
      }
      return key + "=" +  outkey;
    } else {
      error("keyword " + key + " is compulsory for this action");
    }
  }
  return "";
}

void Action::parseFlag(const std::string&key,bool & t) {
  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  // Check keyword is a flag
  plumed_massert( keywords.style(key,"deprecated") || keywords.style(key,"flag") || keywords.style(key,"hidden"), "keyword " + key + " is not a flag");

  // Read in the flag otherwise get the default value from the keywords object
  t = linemap.readAndRemoveFlag(key);
  if(!t) {
    if( keywords.style(key,"nohtml") ) {
      t=false;
    } else if ( !keywords.getLogicalDefault(key,t) ) {
      log.printf("ERROR in action %s with label %s : flag %s has no default",actionName.c_str(),actionLabel.c_str(),key.c_str() );
      plumed_error();
    }
  }
}

void Action::addDependency(Action*action) {
  bool found=false;
  for(const auto & d : after ) {
    if( action==d ) {
      found=true;
      break;
    }
  }
  if( !found ) {
    after.push_back(action);
  }
}

bool Action::checkForDependency( Action* action ) {
  for(const auto & d : after) {
    if( action==d ) {
      return true;
    }
    if( d->checkForDependency(action) ) {
      return true;
    }
  }
  return false;
}

void Action::activate() {
// This is set to true if actions are only need to be computed in setup (during checkRead)
  if( never_activate ) {
    return;
  }
// preparation step is called only the first time an Action is activated.
// since it could change its dependences (e.g. in an ActionAtomistic which is
// accessing to a virtual atom), this is done just before dependencies are
// activated
  if(!active) {
    this->unlockRequests();
    prepare();
    this->lockRequests();
  } else {
    return;
  }
  for(const auto & p : after) {
    p->activate();
  }
  active=true;
}

void Action::setOption(const std::string &s) {
// This overloads the action and activate some options
  actionOptions.insert(s);
  for(const auto & p : after) {
    p->setOption(s);
  }
}

void Action::clearOptions() {
// This overloads the action and activate some options
  actionOptions.clear();
}


void Action::clearDependencies() {
  after.clear();
}

void Action::checkRead() {
  if (!linemap.empty()) {
    std::string msg="cannot understand the following words from the input line : ";
    msg += linemap.keyList(", ");
    error(msg);
  }

  setupConstantValues(false);
}

void Action::setupConstantValues( const bool& have_atoms ) {
  if( have_atoms ) {
    // This ensures that we switch off actions that only depend on constant when passed from the
    // MD code on the first step
    ActionAtomistic* at = castToActionAtomistic();
    ActionWithValue* av = castToActionWithValue();
    if( at && av ) {
      never_activate=av->getNumberOfComponents()>0;
      for(unsigned i=0; i<av->getNumberOfComponents(); ++i) {
        if( !av->copyOutput(i)->isConstant() ) {
          never_activate=false;
          break;
        }
      }
    }
  }
  ActionWithArguments* aa = castToActionWithArguments();
  if( aa && aa->getNumberOfArguments()>0 && getName()!="BIASVALUE" ) {
    never_activate = aa->calculateConstantValues( have_atoms );
  }
}

long long int Action::getStep()const {
  return plumed.getStep();
}

double Action::getTime()const {
  return timestep*getStep();
}

double Action::getTimeStep()const {
  return timestep;
}

double Action::getkBT() {
  double temp=-1.0;
  if( keywords.exists("TEMP") ) {
    parse("TEMP",temp);
  }
  if(temp>=0.0 && keywords.style("TEMP","optional") ) {
    return getKBoltzmann()*temp;
  }
  ActionForInterface* kb=plumed.getActionSet().selectWithLabel<ActionForInterface*>("kBT");
  double kbt=0;
  if(kb) {
    kbt=(kb->copyOutput(0))->get();
  }
  if( temp>=0 && keywords.style("TEMP","compulsory") ) {
    double kB=getKBoltzmann();
    if( kbt>0 && std::abs(kbt-kB*temp)>1e-4) {
      std::string strt1, strt2;
      Tools::convert( temp, strt1 );
      Tools::convert( kbt/kB, strt2 );
      warning("using TEMP=" + strt1 + " while MD engine uses " + strt2 + "\n");
    }
    kbt = kB*temp;
    plumed_massert(kbt>0,"your MD engine does not pass the temperature to plumed, you must specify it using TEMP");
    return kbt;
  }
  return kbt;
}

void Action::exit(int c) {
  plumed.exit(c);
}

void Action::calculateNumericalDerivatives( ActionWithValue* a ) {
  plumed_merror("if you get here it means that you are trying to use numerical derivatives for a class that does not implement them");
}

void Action::prepare() {
  return;
}

[[noreturn]] void Action::error( const std::string & msg ) const {
  log.printf("ERROR in input to action %s with label %s : %s \n \n", actionName.c_str(), actionLabel.c_str(), msg.c_str() );
  plumed_merror("ERROR in input to action " + actionName + " with label " + actionLabel + " : " + msg );
}

void Action::warning( const std::string & msg ) {
  log.printf("WARNING for action %s with label %s : %s \n", actionName.c_str(), actionLabel.c_str(), msg.c_str() );
}

void Action::calculateFromPDB( const PDB& pdb ) {
  activate();
  for(const auto & p : after) {
    ActionWithValue*av=castToActionWithValue();
    if(av) {
      av->clearInputForces();
      av->clearDerivatives();
    }
    p->readAtomsFromPDB( pdb );
    p->calculate();
  }
  readAtomsFromPDB( pdb );
  calculate();
}

bool Action::getExchangeStep()const {
  return plumed.getExchangeStep();
}

std::string Action::cite(const std::string&s) {
  return plumed.cite(s);
}

/// Check if action should be updated.
bool Action::checkUpdate()const {
  double t=getTime();
  if(t<update_until && (update_from==std::numeric_limits<double>::max() || t>=update_from)) {
    return true;
  } else {
    return false;
  }
}

bool Action::getCPT() const {
  return plumed.getCPT();
}

const Units& Action::getUnits() const {
  return plumed.getUnits();
}

bool Action::usingNaturalUnits() const {
  return plumed.usingNaturalUnits();
}

double Action::getKBoltzmann() const {
  if( usingNaturalUnits() ) {
    return 1.0;
  } else {
    return kBoltzmann/getUnits().getEnergy();
  }
}

std::string Action::writeInGraph() const {
  std::string nam=getName();
  std::size_t u=nam.find_last_of("_");
  std::string sub=nam.substr(u+1);
  if( sub=="SCALAR" || sub=="VECTOR" || sub=="GRID" ) {
    return nam.substr(0,u);
  }
  return nam;
}

}

