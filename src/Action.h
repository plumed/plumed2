#ifndef __PLUMED_Action_h
#define __PLUMED_Action_h
#include <vector>
#include <string>
#include <set>
#include "Tools.h"
#include "Log.h"

namespace PLMD{

class PlumedMain;
class PlumedCommunicator;
class Log;

/// This class is used to bring the relevant information to the Action constructor.
/// Only Action and ActionRegister class can access to its content, which is 
/// kept private to other classes, and may change in the future.
class ActionOptions{
  friend class Action;
  friend class ActionRegister;
/// Reference to main PlumedMain object
  PlumedMain& plumed;
/// Input line which sets up the action
  std::vector<std::string> line;
public:
/// Constructor
  ActionOptions(PlumedMain&p,const std::vector<std::string>&);
};

/// This class holds the keywords and their documentation
class Keyword{
friend class Action;
private:
  unsigned compulsory;
  std::string key;
  std::string documentation;
public:
  Keyword(){}
  Keyword( const unsigned& , const std::string&, const std::string& );
};

/// Base class for all the input Actions.
/// The input Actions are more or less corresponding to the directives
/// in the plumed.dat file and are applied in order at each time-step.
class Action
{

/// Name of the directive in the plumed.dat file.
  const std::string name;

/// Label of the Action, as set with LABEL= in the plumed.dat file.
  std::string label;

/// Directive line.
/// This line is progressively erased during Action construction
/// so as to check if all the present keywords are correct.
  std::vector<std::string> line;

public:
  typedef std::set<Action*> Dependencies;

private:
/// Actions on which this Action depends.
  Dependencies after;
/// Actions depending on this Action.
  Dependencies before;

/// The keywords that are available to this particular action 
/// This list is created using registerKeyword.
  std::vector<Keyword> keys;
 
/// Switch to activate Action on this step.
  bool active;

/// The frequency with which the action is performed (in the absence of any dependencies)
  int stride;

protected:

/// Reference to main plumed object
  PlumedMain& plumed;

/// Reference to the log stream
  Log& log;

/// Specify that this Action depends on another one
  void addDependency(Action*);

/// Clear the dependence list for this Action
  void clearDependencies();

/// Return the present timestep
  int getStep()const; 

/// Return the present time
  double getTime()const;

/// Return the timestep
  double getTimeStep()const;

/// Return the value of the stride parameter
  int getStride() const;

/// Make stride a compulsory keyword (for actions like BIAS, PRINT and so on)
  void strideKeywordIsCompulsory();

/// Read in the keywords relevant to action (LABEL and STRIDE)
  void readAction();

/// Check if Action was properly read.
/// This checks if Action::line is empty. It must be called after
/// a final Action has been initialized
  void checkRead();

/// Add a keyword to the list of keywords for this action
  void registerKeyword( const unsigned& , const std::string& , const std::string& ); 

/// Parse one keyword as generic type
  template<class T>
  void parse(const std::string&key,T&t);

/// Parse one keyword as std::vector
  template<class T>
  void parseVector(const std::string&key,std::vector<T>&t);

/// Parse one keyword as boolean flag
  void parseFlag(const std::string&key,bool&t);
	
/// Print errors and die 
  void error( const std::string& stream );
  
/// Print a warning but don't die	
  void warning( const std::string& stream );	

///
  std::set<FILE*> files;
  typedef std::set<FILE*>::iterator files_iterator;

public:
  Action(const ActionOptions&);
  virtual ~Action(){};

  PlumedCommunicator& comm;

/// Prepare an Action for calculation
/// This can be used by Action if they need some special preparation
/// before calculation. Typical case is for collective variables
/// which would like to change their list of requested atoms.
/// By default (if not overridden) does nothing.
  virtual void prepare();

  virtual void lockRequests(){};
  virtual void unlockRequests(){};

/// Calculate an Action.
/// This method is called one or more times per step.
/// The set of all Actions is calculated in forward order.
  virtual void calculate()=0;

/// Apply an Action.
/// This method is called one time per step.
/// The set of all Actions is applied in backward order.
  virtual void apply()=0;

/// Tell to the Action to flush open files
  void fflush();

  virtual std::string getDocumentation()const;

/// Returns the label
  const std::string & getLabel()const;

/// Returns the name
  const std::string & getName()const;

/// Is it time to peform this action
  virtual bool onStep() const;

/// Set action to active
  virtual void activate();

/// Set action to inactive
  virtual void deactivate();

/// Check if action is active
  bool isActive()const;

/// Return dependencies
  const Dependencies & getDependencies()const{return after;}

/// Check if numerical derivatives should be performed
  virtual bool checkNumericalDerivatives()const{return false;}

/// Perform calculation using numerical derivatives
  virtual void calculateNumericalDerivatives();

  FILE *fopen(const char *path, const char *mode);
  int   fclose(FILE*fp);
};

/////////////////////
// FAST INLINE METHODS

inline
const std::string & Action::getLabel()const{
  return label;
}

inline
const std::string & Action::getName()const{
  return name;
}

inline
bool Action::onStep() const {
  return ( stride>0 && getStep()%stride==0 );
}

inline
int Action::getStride() const {
  assert(stride>0); return stride;
}

template<class T>
void Action::parse(const std::string&key,T&t){
  if(!Tools::parse(line,key,t)) error("missing " + key + " keyword");
}

template<class T>
void Action::parseVector(const std::string&key,std::vector<T>&t){
  if(!Tools::parseVector(line,key,t)) error("missing " + key + " keyword");
}

inline
void Action::deactivate(){
  active=false;
}

inline
bool Action::isActive()const{
  return active;
}

}
#endif

