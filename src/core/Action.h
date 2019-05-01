/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#ifndef __PLUMED_core_Action_h
#define __PLUMED_core_Action_h
#include <vector>
#include <string>
#include <set>
#include "tools/Keywords.h"
#include "tools/Tools.h"
#include "tools/Log.h"

namespace PLMD {

class PDB;
class PlumedMain;
class Communicator;
class ActionWithValue;

/// This class is used to bring the relevant information to the Action constructor.
/// Only Action and ActionRegister class can access to its content, which is
/// kept private to other classes, and may change in the future.
class ActionOptions {
  friend class Action;
  friend class ActionRegister;
/// Reference to main PlumedMain object
  PlumedMain& plumed;
/// Input line which sets up the action
  std::vector<std::string> line;
/// The documentation for this action
  const Keywords& keys;
  static Keywords emptyKeys;
public:
/// Constructor
  ActionOptions(PlumedMain&p,const std::vector<std::string>&);
  ActionOptions(const ActionOptions&,const Keywords& keys);
};

/// Base class for all the input Actions.
/// The input Actions are more or less corresponding to the directives
/// in the plumed.dat file and are applied in order at each time-step.
class Action {
  friend class ActionShortcut;

/// Name of the directive in the plumed.dat file.
  const std::string name;

/// Label of the Action, as set with LABEL= in the plumed.dat file.
  std::string label;

/// Directive line.
/// This line is progressively erased during Action construction
/// so as to check if all the present keywords are correct.
  std::vector<std::string> line;

/// Update only after this time.
  double update_from;

/// Update only until this time.
  double update_until;

public:

/// Check if action should be updated.
  bool checkUpdate()const;

public:
  typedef std::vector<Action*> Dependencies;

private:
/// Actions on which this Action depends.
  Dependencies after;

/// Switch to activate Action on this step.
  bool active;

/// Option that you might have enabled
  std::set<std::string> options;

  bool restart;

  bool doCheckPoint;

public:

/// Reference to main plumed object
  PlumedMain& plumed;

/// Reference to the log stream
  Log& log;

/// Specify that this Action depends on another one
  void addDependency(Action*);

/// Clear the dependence list for this Action
  void clearDependencies();

/// Return the present timestep
  long int getStep()const;

/// Return the present time
  double getTime()const;

/// Return the timestep
  double getTimeStep()const;

/// Return true if we are doing a restart
  bool getRestart()const;

/// Return true if we are doing at a checkpoint step
  bool getCPT()const;

/// Just read one of the keywords and return the whole thing as a string
  std::string getKeyword(const std::string& key);

/// Parse one keyword as generic type
  template<class T>
  void parse(const std::string&key,T&t);

/// Parse one numbered keyword as generic type
  template<class T>
  bool parseNumbered(const std::string&key, const int no, T&t);

/// Parse one keyword as std::vector
  template<class T>
  void parseVector(const std::string&key,std::vector<T>&t);

/// Parse a vector with a number
  template<class T>
  bool parseNumberedVector(const std::string& key, const int no, std::vector<T>&t);

/// Parse one keyword as boolean flag
  void parseFlag(const std::string&key,bool&t);

/// Crash calculation and print documentation
  void error( const std::string & msg ) const;

/// Issue a warning
  void warning( const std::string & msg );

/// Exit with error code c
  void exit(int c=0);

///
  std::set<FILE*> files;

public:
/// Standard constructor from ActionOptions
  explicit Action(const ActionOptions&);
/// Destructor
  virtual ~Action();
private:
/// Copy constructor is deleted
  Action(const Action&a) = delete;
/// Assignment operator is deleted
  Action& operator=(const Action&a) = delete;
  int replica_index;
public:
/// Check if Action was properly read.
/// This checks if Action::line is empty. It must be called after
/// a final Action has been initialized
  void checkRead();

  Communicator& comm;
  Communicator& multi_sim_comm;

  const Keywords& keywords;
/// Prepare an Action for calculation
/// This can be used by Action if they need some special preparation
/// before calculation. Typical case is for collective variables
/// which would like to change their list of requested atoms.
/// By default (if not overridden) does nothing.
  virtual void prepare();

/// Register all the relevant keywords for the action
  static void registerKeywords( Keywords& keys );

  virtual void lockRequests() {}
  virtual void unlockRequests() {}

/// Calculate an Action.
/// This method is called one or more times per step.
/// The set of all Actions is calculated in forward order.
  virtual void calculate()=0;

/// Apply an Action.
/// This method is called one time per step.
/// The set of all Actions is applied in backward order.
  virtual void apply()=0;

/// Before Update.
/// This is a special method that is called just
/// before the update() method. It can be used by
/// actions that want to do something irrespectively
/// of the fact that update() is active or not.
/// In other words, this is *always* called, even when action
/// is not active.
  virtual void beforeUpdate() {}

/// Update.
/// This method is called one time per step.
/// The set of all Actions is updated in forward order.
  virtual void update() {}

/// RunFinalJobs
/// This method is called once at the very end of the calculation.
/// The set of all Actions in run for the final time in forward order.
  virtual void runFinalJobs() {}

/// Tell to the Action to flush open files
  void fflush();

  virtual std::string getDocumentation()const;

/// Returns the label
  const std::string & getLabel()const;

/// Returns the name
  const std::string & getName()const;

/// Set action to active
  virtual void activate();

///
  virtual void setOption(const std::string &s);

  virtual void clearOptions();

/// Set action to inactive
  virtual void deactivate();

/// Check if action is active
  bool isActive()const;

/// Check if an option is on
  bool isOptionOn(const std::string &s)const;

/// Return dependencies
  const Dependencies & getDependencies()const {return after;}

/// Check if numerical derivatives should be performed
  virtual bool checkNumericalDerivatives()const {return false;}

/// Check if the action needs gradient
  virtual bool checkNeedsGradients()const {return false;}

/// Perform calculation using numerical derivatives
/// N.B. only pass an ActionWithValue to this routine if you know exactly what you
/// are doing.
  virtual void calculateNumericalDerivatives( ActionWithValue* a=NULL );

/// Opens a file.
/// This is similar to plain fopen, but with some extra functionality.
/// * When opened for writing, processors other than the one with rank 0 just open /dev/null
/// * PlumedMain::fopen is used, so that other tricks may appear (see \ref PlumedMain::fopen)
  FILE *fopen(const char *path, const char *mode);
/// Closes a file opened with Action::fclose().
  int   fclose(FILE*fp);

/// Calculate the action given a pdb file as input.  This is used to initialize
/// things like distance from a point in CV map space given a pdb as an input file
  void calculateFromPDB( const PDB&  );
/// This is overwritten in ActionAtomistic so that we can read
/// the atoms from the pdb input file rather than taking them from the
/// MD code
  virtual void readAtomsFromPDB( const PDB&  ) {}
/// Check if we are on an exchange step
  bool getExchangeStep()const;

/// Cite a paper see PlumedMain::cite
  std::string cite(const std::string&s);
};

/////////////////////
// FAST INLINE METHODS

inline
const std::string & Action::getLabel()const {
  return label;
}

inline
const std::string & Action::getName()const {
  return name;
}

template<class T>
void Action::parse(const std::string&key,T&t) {
//  if(!Tools::parse(line,key,t)){
//    log.printf("ERROR parsing keyword %s\n",key.c_str());
//    log.printf("%s\n",getDocumentation().c_str());
//    this->exit(1);
//  }
  // Check keyword has been registered
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");

  // Now try to read the keyword
  std::string def;
  bool present=Tools::findKeyword(line,key);
  bool found=Tools::parse(line,key,t,replica_index);
  if(present && !found) error("keyword " + key +" could not be read correctly");

  // If it isn't read and it is compulsory see if a default value was specified
  if ( !found && (keywords.style(key,"compulsory") || keywords.style(key,"hidden")) ) {
    if( keywords.getDefaultValue(key,def) ) {
      if( def.length()==0 || !Tools::convert(def,t) ) {
        log.printf("ERROR in action %s with label %s : keyword %s has weird default value",name.c_str(),label.c_str(),key.c_str() );
        this->exit(1);
      }
    } else if( keywords.style(key,"compulsory") ) {
      error("keyword " + key + " is compulsory for this action");
    }
  }
}

template<class T>
bool Action::parseNumbered(const std::string&key, const int no, T&t) {
  // Check keyword has been registered
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");
  if( !keywords.numbered(key) ) error("numbered keywords are not allowed for " + key );

  // Now try to read the keyword
  std::string num; Tools::convert(no,num);
  return Tools::parse(line,key+num,t,replica_index);
}

template<class T>
void Action::parseVector(const std::string&key,std::vector<T>&t) {
//  if(!Tools::parseVector(line,key,t)){
//    log.printf("ERROR parsing keyword %s\n",key.c_str());
//    log.printf("%s\n",getDocumentation().c_str());
//    this->exit(1);
//  }

  // Check keyword has been registered
  plumed_massert(keywords.exists(key), "keyword " + key + " has not been registered");
  unsigned size=t.size(); bool skipcheck=false;
  if(size==0) skipcheck=true;

  // Now try to read the keyword
  std::string def; T val;
  bool present=Tools::findKeyword(line,key);
  bool found=Tools::parseVector(line,key,t,replica_index);
  if(present && !found) error("keyword " + key +" could not be read correctly");

  // Check vectors size is correct (not if this is atoms or ARG)
  if( !keywords.style(key,"atoms") && found ) {
//     bool skipcheck=false;
//     if( keywords.style(key,"compulsory") ){ keywords.getDefaultValue(key,def); skipcheck=(def=="nosize"); }
    if( !skipcheck && t.size()!=size ) error("vector read in for keyword " + key + " has the wrong size");
  }

  // If it isn't read and it is compulsory see if a default value was specified
  if ( !found && (keywords.style(key,"compulsory") || keywords.style(key,"hidden")) ) {
    if( keywords.getDefaultValue(key,def) ) {
      if( def.length()==0 || !Tools::convert(def,val) ) {
        log.printf("ERROR in action %s with label %s : keyword %s has weird default value",name.c_str(),label.c_str(),key.c_str() );
        this->exit(1);
      } else {
        if(t.size()>0) for(unsigned i=0; i<t.size(); ++i) t[i]=val;
        else t.push_back(val);
      }
    } else if( keywords.style(key,"compulsory") ) {
      error("keyword " + key + " is compulsory for this action");
    }
  } else if ( !found ) {
    t.resize(0);
  }
}

template<class T>
bool Action::parseNumberedVector(const std::string&key, const int no, std::vector<T>&t) {
  plumed_massert(keywords.exists(key),"keyword " + key + " has not been registered");
  if( !keywords.numbered(key) ) error("numbered keywords are not allowed for " + key );

  unsigned size=t.size(); bool skipcheck=false;
  if(size==0) skipcheck=true;
  std::string num; Tools::convert(no,num);
  bool present=Tools::findKeyword(line,key);
  bool found=Tools::parseVector(line,key+num,t,replica_index);
  if(present && !found) error("keyword " + key +" could not be read correctly");

  if(  keywords.style(key,"compulsory") ) {
    if (!skipcheck && found && t.size()!=size ) error("vector read in for keyword " + key + num + " has the wrong size");
  } else if ( !found ) {
    t.resize(0);
  }
  return found;
}

inline
void Action::deactivate() {
  options.clear();
  active=false;
}

inline
bool Action::isActive()const {
  return active;
}

inline
bool Action::isOptionOn(const std::string &s)const {
  return options.count(s);
}

inline
bool Action::getRestart()const {
  return restart;
}

}
#endif

