#ifndef __PLUMED_PlumedMain_h
#define __PLUMED_PlumedMain_h

#include "WithCmd.h"
#include <cstdio>
#include <string>
#include <vector>


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
class PlumedCommunicator;

/// Main plumed object.
/// In MD engines this object is not manipulated directly but it is wrapped in
/// plumed or PLMD::Plumed objects. Its main method is cmd(),
/// which defines completely the external plumed interface.
/// It does not contain any static data.
class PlumedMain:
  public WithCmd
{
public:
/// Communicator for plumed.
/// Includes all the processors used by plumed.
  PlumedCommunicator&comm;

private:
  DLLoader& dlloader;

  WithCmd* grex;
/// Flag to avoid double initialization
  bool  initialized;
/// Name of MD engine
  std::string MDEngine;
/// Log stream
  Log& log;


/// Present step number.
  int step;

/// Condition for plumed to be active.
/// At every step, PlumedMain is checking if there are Action's requiring some work.
/// If at least one Action requires some work, this variable is set to true.
  bool active;

/// Name of the input file
  std::string plumedDat;

/// Object containing information about atoms (such as positions,...).
  Atoms&    atoms;           // atomic coordinates

/// Set of actions found in plumed.dat file
  ActionSet& actionSet;

/// Set of Pilot actions.
/// These are the action the, if they are Pilot::onStep(), can trigger execution
  std::vector<ActionPilot*> pilots;

/// Suffix string for file opening, useful for multiple simulations in the same directory
  std::string suffix;

/// The total bias (=total energy of the restraints)
  double bias;

public:
/// Flag to switch off virial calculation (for debug)
  bool novirial;

public:
  PlumedMain();
// this is to access to WithCmd versions of cmd (allowing overloading of a virtual method)
  using WithCmd::cmd;
/// cmd method, accessible with standard Plumed.h interface.
/// It is called as plumed_cmd() or as PLMD::Plumed::cmd()
/// It is the interpreter for plumed commands. It basically contains the definition of the plumed interface.
/// If you want to add a new functionality to the interface between plumed
/// and an MD engine, this is the right place
/// Notice that this interface should always keep retro-compatibility
  void cmd(const std::string&key,void*val=NULL);
  ~PlumedMain();
/// Read an input file.
/// \param str name of the file
  void readInputFile(std::string str);
/// Initialize the object
  void init();
/// Prepare the calculation.
/// Here it is checked which are the active Actions and communication of the relevant atoms is initiated
  void prepareCalc();
  void prepareDependencies();
  void shareData();
/// Perform the calculation.
/// Here atoms coordinates are received and all Action's are applied in order
  void performCalc();
/// Shortcut for prepareCalc() + performCalc()
  void calc();

  void waitData();
  void justCalculate();
  void justApply();
/// Reference to atoms object
  Atoms& getAtoms();
/// Reference to the list of Action's
  const ActionSet & getActionSet()const;
/// Referenge to the log stream
  Log & getLog();
/// Return the number of the step
  const int & getStep()const{return step;};
/// Stop the run
  void exit(int c=0);
/// Load a shared library
  void load(std::vector<std::string> & words);
/// Get the suffix string
  const std::string & getSuffix()const;
/// Set the suffix string
  void setSuffix(const std::string&);
/// get the value of the bias
  double getBias()const;
  FILE* fopen(const char *path, const char *mode);
  int fclose(FILE*fp);
};

/////
// FAST INLINE METHODS:

inline
const ActionSet & PlumedMain::getActionSet()const{
  return actionSet;
}

inline
Atoms& PlumedMain::getAtoms(){
  return atoms;
}

inline
const std::string & PlumedMain::getSuffix()const{
  return suffix;
}

inline
void PlumedMain::setSuffix(const std::string&s){
  suffix=s;
}


}

#endif

