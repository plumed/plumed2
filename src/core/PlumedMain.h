/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_core_PlumedMain_h
#define __PLUMED_core_PlumedMain_h

#include "WithCmd.h"
#include "tools/ForwardDecl.h"
#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <stack>
#include <memory>
#include <map>

// !!!!!!!!!!!!!!!!!!!!!!    DANGER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// THE FOLLOWING ARE DEFINITIONS WHICH ARE NECESSARY FOR DYNAMIC LOADING OF THE PLUMED KERNEL:
// This section should be consistent with the Plumed.h file.
// Since the Plumed.h file may be included in host MD codes, **NEVER** MODIFY THE CODE DOWN HERE

/* Generic function pointer */
typedef void (*plumed_function_pointer)(void);

/* Holder for function pointer */
typedef struct {
  plumed_function_pointer p;
} plumed_function_holder;

// END OF DANGER
////////////////////////////////////////////////////////////

namespace PLMD {



class ActionAtomistic;
class ActionPilot;
class Log;
class Atoms;
class ActionSet;
class DLLoader;
class Communicator;
class Stopwatch;
class Citations;
class ExchangePatterns;
class FileBase;
class DataFetchingObject;

/**
Main plumed object.
In MD engines this object is not manipulated directly but it is wrapped in
plumed or PLMD::Plumed objects. Its main method is cmd(),
which defines completely the external plumed interface.
It does not contain any static data.
*/
class PlumedMain:
  public WithCmd
{
/// Pointers to files opened in actions associated to this object.
/// Notice that with the current implementation this should be at the top of this
/// structure. Indeed, this should be destroyed *after* all the actions allocated
/// in this PlumedMain object have been destroyed.
  std::set<FileBase*> files;
/// Forward declaration.
  ForwardDecl<Communicator> comm_fwd;
public:
/// Communicator for plumed.
/// Includes all the processors used by plumed.
  Communicator&comm=*comm_fwd;

private:
/// Forward declaration.
  ForwardDecl<Communicator> multi_sim_comm_fwd;
public:
  Communicator&multi_sim_comm=*multi_sim_comm_fwd;

private:
/// Error handler.
/// Pointer to a function that is called an exception thrown within
/// the library is about to leave the library.
/// Can be used to remap exceptions in case the plumed wrapper was compiled
/// with a different version of the C++ standard library.
/// Should only be called from \ref plumed_plumedmain_cmd().
  typedef struct {
    void* ptr;
    void(*handler)(void* ptr,int code,const char*);
  } plumed_error_handler;

  plumed_error_handler error_handler= {NULL,NULL};

/// Forward declaration.
  ForwardDecl<DLLoader> dlloader_fwd;
  DLLoader& dlloader=*dlloader_fwd;

  std::unique_ptr<WithCmd> cltool;

  std::unique_ptr<WithCmd> grex;
/// Flag to avoid double initialization
  bool  initialized;
/// Name of MD engine
  std::string MDEngine;

/// Forward declaration.
  ForwardDecl<Log> log_fwd;
/// Log stream
  Log& log=*log_fwd;

/// Forward declaration.
/// Should be placed after log since its constructor takes a log reference as an argument.
  ForwardDecl<Stopwatch> stopwatch_fwd;
  Stopwatch& stopwatch=*stopwatch_fwd;

/// Forward declaration.
  ForwardDecl<Citations> citations_fwd;
/// tools/Citations.holder
  Citations& citations=*citations_fwd;

/// Present step number.
  long int step;

/// Condition for plumed to be active.
/// At every step, PlumedMain is checking if there are Action's requiring some work.
/// If at least one Action requires some work, this variable is set to true.
  bool active;

/// Name of the input file
  std::string plumedDat;

/// Object containing data we would like to grab and pass back
  std::unique_ptr<DataFetchingObject> mydatafetcher;

/// End of input file.
/// Set to true to terminate reading
  bool endPlumed;

/// Forward declaration.
  ForwardDecl<Atoms> atoms_fwd;
/// Object containing information about atoms (such as positions,...).
  Atoms&    atoms=*atoms_fwd;           // atomic coordinates

/// Forward declaration.
  ForwardDecl<ActionSet> actionSet_fwd;
/// Set of actions found in plumed.dat file
  ActionSet& actionSet=*actionSet_fwd;

/// Set of Pilot actions.
/// These are the action the, if they are Pilot::onStep(), can trigger execution
  std::vector<ActionPilot*> pilots;

/// Suffix string for file opening, useful for multiple simulations in the same directory
  std::string suffix;

/// The total bias (=total energy of the restraints)
  double bias;

/// The total work.
/// This computed by accumulating the change in external potentials.
  double work;

/// Forward declaration.
  ForwardDecl<ExchangePatterns> exchangePatterns_fwd;
/// Class of possible exchange patterns, used for BIASEXCHANGE but also for future parallel tempering
  ExchangePatterns& exchangePatterns=*exchangePatterns_fwd;

/// Set to true if on an exchange step
  bool exchangeStep;

/// Flag for restart
  bool restart;

/// Flag for checkpointig
  bool doCheckPoint;


/// Stuff to make plumed stop the MD code cleanly
  int* stopFlag;
  bool stopNow;

/// Stack for update flags.
/// Store information used in class \ref generic::UpdateIf
  std::stack<bool> updateFlags;

public:
/// Flag to switch off virial calculation (for debug and MD codes with no barostat)
  bool novirial;

/// Flag to switch on detailed timers
  bool detailedTimers;

/// Generic map string -> double
/// intended to pass information across Actions
  std::map<std::string,double> passMap;

/// Add a citation, returning a string containing the reference number, something like "[10]"
  std::string cite(const std::string&);

/// Get number of threads that can be used by openmp
  unsigned getNumThreads()const;

/// Get a reasonable number of threads so as to access to an array of size s located at x
  template<typename T>
  unsigned getGoodNumThreads(const T*x,unsigned s)const;

/// Get a reasonable number of threads so as to access to vector v;
  template<typename T>
  unsigned getGoodNumThreads(const std::vector<T> & v)const;

public:
  PlumedMain();
// this is to access to WithCmd versions of cmd (allowing overloading of a virtual method)
  using WithCmd::cmd;
  /**
   cmd method, accessible with standard Plumed.h interface.
   \param key The name of the command to be executed.
   \param val The argument of the command to be executed.
   It is called as plumed_cmd() or as PLMD::Plumed::cmd()
   It is the interpreter for plumed commands. It basically contains the definition of the plumed interface.
   If you want to add a new functionality to the interface between plumed
   and an MD engine, this is the right place
   Notice that this interface should always keep retro-compatibility
  */
  void cmd(const std::string&key,void*val=NULL) override;
  ~PlumedMain();
  /**
    Read an input file.
    \param str name of the file
  */
  void readInputFile(std::string str);
  /**
    Read an input string.
    \param str name of the string
  */
  void readInputWords(const std::vector<std::string> &  str);

  /**
    Read an input string.
    \param str name of the string
    At variance with readInputWords(), this is splitting the string into words
  */
  void readInputLine(const std::string & str);

  /**
    Read an input buffer.
    \param str name of the string
    Same as readInputFile, but first write str on a temporary file and then read
    that files. At variance with readInputLine, it can take care of comments and
    continuation lines.
  */
  void readInputLines(const std::string & str);

  /**
    Initialize the object.
    Should be called once.
  */
  void init();
  /**
    Prepare the calculation.
    Here it is checked which are the active Actions and communication of the relevant atoms is initiated.
    Shortcut for prepareDependencies() + shareData()
  */
  void prepareCalc();
  /**
    Prepare the list of active Actions and needed atoms.
    Scan the Actions to see which are active and which are not, so as to prepare a list of
    the atoms needed at this step.
  */
  void prepareDependencies();
  /**
    Share the needed atoms.
    In asynchronous implementations, this method sends the required atoms to all the plumed processes,
    without waiting for the communication to complete.
  */
  void shareData();
  /**
    Perform the calculation.
    Shortcut for waitData() + justCalculate() + justApply().
    Equivalently: waitData() + justCalculate() + backwardPropagate() + update().
  */
  void performCalc();
  /**
    Perform the calculation without update()
    Shortcut for: waitData() + justCalculate() + backwardPropagate()
  */
  void performCalcNoUpdate();
  /**
    Complete PLUMED calculation.
    Shortcut for prepareCalc() + performCalc()
  */
  void calc();
  /**
    Scatters the needed atoms.
    In asynchronous implementations, this method waits for the communications started in shareData()
    to be completed. Otherwise, just send around needed atoms.
  */
  void waitData();
  /**
    Perform the forward loop on active actions.
  */
  void justCalculate();
  /**
    Backward propagate and update.
    Shortcut for backwardPropagate() + update()
    I leave it here for backward compatibility
  */
  void justApply();
  /**
    Perform the backward loop on active actions.
    Needed to apply the forces back.
  */
  void backwardPropagate();
  /**
    Call the update() method.
  */
  void update();
  /**
    If there are calculations that need to be done at the very end of the calculations this
    makes sures they are done
  */
  void runJobsAtEndOfCalculation();
/// Reference to atoms object
  Atoms& getAtoms();
/// Reference to the list of Action's
  const ActionSet & getActionSet()const;
/// Referenge to the log stream
  Log & getLog();
/// Return the number of the step
  long int getStep()const {return step;}
/// Stop the run
  void exit(int c=0);
/// Load a shared library
  void load(const std::string&);
/// Get the suffix string
  const std::string & getSuffix()const;
/// Set the suffix string
  void setSuffix(const std::string&);
/// get the value of the bias
  double getBias()const;
/// get the value of the work
  double getWork()const;
/// Opens a file.
/// Similar to plain fopen, but, if it finds an error in opening the file, it also tries with
/// path+suffix.  This trick is useful for multiple replica simulations.
  FILE* fopen(const char *path, const char *mode);
/// Closes a file opened with PlumedMain::fopen()
  int fclose(FILE*fp);
/// Insert a file
  void insertFile(FileBase&);
/// Erase a file
  void eraseFile(FileBase&);
/// Flush all files
  void fflush();
/// Check if restarting
  bool getRestart()const;
/// Set restart flag
  void setRestart(bool f) {restart=f;}
/// Check if checkpointing
  bool getCPT()const;
/// Set exchangeStep flag
  void setExchangeStep(bool f);
/// Get exchangeStep flag
  bool getExchangeStep()const;
/// Stop the calculation cleanly (both the MD code and plumed)
  void stop();
/// Enforce active flag.
/// This is a (bit dirty) hack to solve a bug. When there is no active ActionPilot,
/// several shortcuts are used. However, these shortcuts can block GREX module.
/// This function allows to enforce active plumed when doing exchanges,
/// thus fixing the bug.
  void resetActive(bool active);

/// Access to exchange patterns
  ExchangePatterns& getExchangePatterns() {return exchangePatterns;}

/// Push a state to update flags
  void updateFlagsPush(bool);
/// Pop a state from update flags
  void updateFlagsPop();
/// Get top of update flags
  bool updateFlagsTop();
/// Set end of input file
  void setEndPlumed();
/// Call error handler.
/// Should only be called from \ref plumed_plumedmain_cmd().
/// If the error handler was not set, returns false.
  bool callErrorHandler(int code,const char* msg)const;
};

/////
// FAST INLINE METHODS:

inline
const ActionSet & PlumedMain::getActionSet()const {
  return actionSet;
}

inline
Atoms& PlumedMain::getAtoms() {
  return atoms;
}

inline
const std::string & PlumedMain::getSuffix()const {
  return suffix;
}

inline
void PlumedMain::setSuffix(const std::string&s) {
  suffix=s;
}

inline
bool PlumedMain::getRestart()const {
  return restart;
}

inline
bool PlumedMain::getCPT()const {
  return doCheckPoint;
}

inline
void PlumedMain::setExchangeStep(bool s) {
  exchangeStep=s;
}

inline
bool PlumedMain::getExchangeStep()const {
  return exchangeStep;
}

inline
void PlumedMain::resetActive(bool active) {
  this->active=active;
}

inline
void PlumedMain::updateFlagsPush(bool on) {
  updateFlags.push(on);
}

inline
void PlumedMain::updateFlagsPop() {
  updateFlags.pop();
}

inline
bool PlumedMain::updateFlagsTop() {
  return updateFlags.top();
}

inline
void PlumedMain::setEndPlumed() {
  endPlumed=true;
}

inline
bool PlumedMain::callErrorHandler(int code,const char* msg)const {
  if(error_handler.handler) {
    error_handler.handler(error_handler.ptr,code,msg);
    return true;
  } else return false;
}


}

#endif

