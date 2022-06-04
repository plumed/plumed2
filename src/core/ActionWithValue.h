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
#ifndef __PLUMED_core_ActionWithValue_h
#define __PLUMED_core_ActionWithValue_h

#include "Action.h"
#include "Value.h"
#include "tools/Exception.h"
#include "tools/MultiValue.h"
#include "tools/ForwardDecl.h"
#include <vector>
#include <memory>

namespace PLMD {
class Stopwatch;

/**
\ingroup MULTIINHERIT
Used to create a PLMD::Action that has some scalar or vectorial output that may or may not have some derivatives.
This is used for PLMD::Bias, PLMD::Colvar and PLMD::Function

The vast majority of the PLMD::Action objects that are implemented in
plumed calculate some quantity or a set of quantities.  This could be
the value of a CV, the value of a function or the potential due to a bias.
PLMD::ActionWithValue provides the functionality for storing these quantities
and (in tandem with PLMD::ActionWithArguments) the functionality for passing
quantities between PLMD::Actions.  When you are deciding what quantities
your new PLMD::Action will need to store using PLMD::ActionWithValue you must
ask yourself the following two questions:

- Do I need to differentiate my output quantities
- Is my PLMD::Action calculating a single thing or does the output have multiple components

If the answer to the first of these questions is yes then you must setup your values
you using either PLMD::ActionWithValue::addValueWithDerivatives() or
PLMD::ActionWithValue::addComponentWithDerivatives.  If the answer is no you
can set up values using PLMD::ActionWithValue::addValue() or PLMD::ActionWithValue::addComponent().
The precise routine you use to setup your values will depend on your answer to the
second question.  As you are probably aware if the output of your PLMD::Action is a
single quantity you can reference that quantity in the input file using the label of the
PLMD::Action it was calculated in.  If your action <b> outputs only one quantity </b>
we call that quantity the <b> value </b> of the Action.  To set the <b> value </b> and get pointers to it
you should <b> use the set of routines that have the word value in the name </b>.  If, by contrast,
your PLMD::Action calculates multiple quantities then these quantities are referenced in input using the
label.component syntax.  We refer to these <b> multiple quantities </b> the <b> components </b>
of the PLMD::Action.  Perhaps unsurprisingly, when you manipulate the <b> components </b> of an
PLMD::Action you should use <b> the routines with the word component in the name. </b>
*/

class ActionAtomistic;
class ActionWithArguments;

class ActionWithValue :
  public virtual Action
{
  friend class AverageBase;
  friend class CollectFrames;
  friend class ActionWithArguments;
private:
/// An array containing the values for this action
  std::vector<std::unique_ptr<Value>> values;
/// Are we skipping the calculation of the derivatives
  bool noderiv;
/// Are we using numerical derivatives to differentiate
  bool numericalDerivatives;
/// Are we no using open mp for some reason
  bool no_openmp;
/// Are we running calculations in serial
  bool serial;
/// Are we using timers
  bool timers;
/// Can the input be a mixture of history dependent quantities and non-history dependent quantities
  bool allow_mixed_history_input;
/// The stopwatch that times the different parts of the calculation
  ForwardDecl<Stopwatch> stopwatch_fwd;
  Stopwatch& stopwatch=*stopwatch_fwd;
/// The current number of active tasks
  unsigned nactive_tasks;
/// Stores the labels of all the actions that are in this chain.  Used when setting up what to calculate
  std::vector<std::string> actionsLabelsInChain;
/// The indices of the tasks in the full list of tasks
  std::vector<unsigned> indexOfTaskInFullList;
/// The list of currently active tasks  
  std::vector<unsigned> partialTaskList;
/// Ths full list of tasks we have to perform
  std::vector<unsigned> fullTaskList;
/// This list is used to update the active tasks in the list
  std::vector<unsigned> taskFlags;
/// The buffer that we use (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
/// A pointer to this but with actionWithArgument so we can avoid lots of dynamic_cast
  const ActionWithArguments* thisAsActionWithArguments;
/// Action that must be done before this one
  ActionWithValue* action_to_do_before;
/// Actions that must be done after this one
  ActionWithValue* action_to_do_after;
  ActionAtomistic* atom_action_to_do_after;
/// Check if the input to this action is a time series
  bool inputIsTimeSeries() const ;
/// Return the index for the component named name
  int getComponent( const std::string& name ) const;
////
  void selectActiveTasks( const std::vector<std::string>& actionLabelsInChain, bool& forceAllTasks,
                          std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
///
  void prepareForTaskLoop( const unsigned& nactive, const std::vector<unsigned>& pTaskList );
///  Run the task
  void runTask( const unsigned& index, const unsigned& taskno, MultiValue& myvals ) const ;
/// Retrieve the forces for a particualr task
  void gatherAccumulators( const unsigned& index, const MultiValue& myvals, std::vector<double>& buf ) const ;
  void clearAllForcesInChain();
  bool checkForGrids() const ;
  void getNumberOfStreamedDerivatives( unsigned& nderivatives ) const ;
  void getNumberOfStreamedQuantities( unsigned& nquants, unsigned& ncols, unsigned& nmat ) const ;
  unsigned getGridArgumentIndex( const ActionWithArguments* aa ) const ;
public:

// -------- The action has one value only  ---------------- //

/// Add a value with the name label
  void addValue( const std::vector<unsigned>& shape=std::vector<unsigned>() );
/// Add a value with the name label that has derivatives
  void addValueWithDerivatives( const std::vector<unsigned>& shape=std::vector<unsigned>() );
/// Set your default value to have no periodicity
  void setNotPeriodic();
/// Set the value to be periodic with a particular domain
  void setPeriodic( const std::string& min, const std::string& max );
  virtual void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize );
protected:
/// Get a pointer to the default value
  Value* getPntrToValue();
///
  void getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const ;
/// Set the default value (the one without name)
  void setValue(const double& d);
// -------- The action has multiple components ---------- //

///
  unsigned getTaskCode( const unsigned& ii ) const ;
///
  bool checkUsedOutsideOfChain( const std::vector<std::string>& actionLabelsInChain, const std::string& parent, 
                                std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
///
  unsigned setTaskFlags( std::vector<unsigned>& tflags, std::vector<unsigned>&  pTaskList, std::vector<unsigned>& pIndexList );
/// Run all the tasks in the list
  void runAllTasks();
/// Run all calculations in serial
  bool runInSerial() const ;
///  Ruan all calculations without open MP
  bool runWithoutOpenMP() const ;
/// Gather all the data in one row of the matrix
  void gatherMatrixRow( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                        const unsigned& bufstart, std::vector<double>& buffer ) const ;
public:
/// Get the action that does the calculation
  ActionWithValue* getActionThatCalculates();
/// This allows us to have components and values in special circumstances
  virtual bool allowComponentsAndValue() const { return false; }
/// Get the list of actions that are required to calculate this particular action
  void getAllActionsRequired( std::vector<const ActionWithValue*>& allvals ) const ;
/// Add a value with a name like label.name
  void addComponent( const std::string& name, const std::vector<unsigned>& shape=std::vector<unsigned>() );
/// Add a value with a name like label.name that has derivatives
  void addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape=std::vector<unsigned>() );
/// Set your value component to have no periodicity
  void componentIsNotPeriodic( const std::string& name );
/// Set the value to be periodic with a particular domain
  void componentIsPeriodic( const std::string& name, const std::string& min, const std::string& max );
///
  void addTaskToList( const unsigned& taskCode );
/// Retrieve all the scalar values calculated in the loop
  void retrieveAllScalarValuesInLoop( const std::string& ulab, unsigned& nargs, std::vector<Value*>& myvals );
///
  const std::vector<unsigned>& getCurrentTasks() const ;
protected:
/// Get a pointer to the output value
  Value* getPntrToOutput( const unsigned& ind ) const ;
/// Return a pointer to the component by index
  Value* getPntrToComponent(int i);
/// Return a pointer to the value by name
  Value* getPntrToComponent(const std::string& name);
public:
  explicit ActionWithValue(const ActionOptions&ao);
  ~ActionWithValue();

/// Register all the relevant keywords for the action
  static void registerKeywords( Keywords& keys );
/// Insist that numerical derivatives should always be used for an action and make this fact appear in the manual
  static void noAnalyticalDerivatives(Keywords& keys);
/// Puts a message into the manual that the components always output
  static void componentsAreNotOptional(Keywords& keys);
/// The components in the action will depend on the user
  static void useCustomisableComponents(Keywords& keys);
/// Are we not calculating derivatives
  virtual bool doNotCalculateDerivatives() const ;
/// Get the value of one of the components of the PLMD::Action
  double getOutputQuantity( const unsigned j ) const ;
/// Get the value with a specific name (N.B. if there is no such value this returns zero)
  double getOutputQuantity( const std::string& name ) const ;

//  --- Routines for passing stuff to ActionWithArguments -- //

/// Check if a value with a particular name is present.  This is only used in PLMD::ActionWithArguments.
/// You should not use it when manipulating components.
  bool exists( const std::string& name ) const;
/// Return a pointer to the value with name (this is used to retrieve values in other PLMD::Actions)
/// You should NEVER use this routine to refer to the components of your PLMD::Action.  Use
/// getPntrToComponent instead.
  Value* copyOutput( const std::string&name ) const;
/// Return a pointer to the value with this number (this is used to retrieve values in other PLMD::Actions)
/// You should NEVER use this routine to refer to the components of your PLMD::Action.  Use
/// getPntrToComponent instead.
  Value* copyOutput( const unsigned& n ) const;
/// get a string that contains all the available components
  std::string getComponentsList( ) const ;
/// get a vector that contains the label for all the components
  std::vector<std::string> getComponentsVector( ) const ;


// -- Routines for everything else -- //

/// Returns the number of values defined
  int getNumberOfComponents() const ;
///
  bool actionInChain() const ;
///
  virtual bool canChainFromThisAction() const { return true; }
/// Clear the forces on the values
  void clearInputForces();
/// Clear the derivatives of values wrt parameters
  virtual void clearDerivatives( const bool& force=false );
/// Calculate the gradients and store them for all the values (need for projections)
  virtual void setGradientsIfNeeded();
/// Set the value
  void setValue(Value*,double);
/// Check if numerical derivatives should be used
  bool checkNumericalDerivatives() const override;
/// This forces the class to use numerical derivatives
  void useNumericalDerivatives();
// These are things for using vectors of values as fields
  virtual void checkFieldsAllowed() { error("cannot use this action as a field"); }
  virtual unsigned getNumberOfDerivatives() const = 0;
/// Activate the calculation of derivatives
  virtual void turnOnDerivatives();
///
  unsigned getFullNumberOfTasks() const ;
/// Reperform one of the tasks
  void rerunTask( const unsigned& task_index, MultiValue& myvals ) const ;
/// Get the number of columns for the matrix
  virtual unsigned getNumberOfColumns() const { plumed_merror("in " + getName() + " method for number of columns is not defined"); }
/// Run a task for a matrix element
  void runTask( const std::string& controller, const unsigned& task_index, const unsigned& current, const unsigned colno, MultiValue& myvals ) const ;
  void clearMatrixElements( MultiValue& myvals ) const ;
///
  virtual void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {}
///
  virtual void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                     std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                     std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
///
  virtual void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  virtual void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ; 
/// Make sure all tasks required for loop are done before loop starts
  virtual void prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {}
///
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const { plumed_error(); }
///
  virtual void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                  const unsigned& bufstart, std::vector<double>& buffer ) const ; 
/// This one calculates matrix elements
  virtual bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const { return true; }
///
  virtual void gatherForces( const unsigned& task_index, const MultiValue& myvals, std::vector<double>& forces ) const ;
///
  virtual void finishComputations( const std::vector<double>& buf );
///
  bool addActionToChain( const std::vector<std::string>& alabels, ActionWithValue* act );
///
  virtual bool canBeAfterInChain( ActionWithValue* av ) { return true; }
///
  virtual bool valuesComputedInChain() const { return true; }
///
  virtual void transformFinalValueAndDerivatives( const std::vector<double>& buf  ) {};
/// Retrieve the forces acting on all values
  bool getForcesFromValues( std::vector<double>& forces );
///
  virtual std::string writeInGraph() const { return getName(); }
///
  void generateGraphNodes( OFile& ofile, std::vector<std::string>& graph_actions ) const ;
  static std::string getCleanGraphLabel( const std::string& glab );
};

inline
double ActionWithValue::getOutputQuantity(const unsigned j) const {
  plumed_massert(j<values.size(),"index requested is out of bounds");
  return values[j]->get();
}

inline
double ActionWithValue::getOutputQuantity( const std::string& name ) const {
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->name==thename ) return values[i]->data[0];
  }
  return 0.0;
}

inline
void ActionWithValue::setValue(const double& d) {
  plumed_massert(values.size()==1, "cannot use setValue in multi-component actions");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->set(d);
}

inline
int ActionWithValue::getNumberOfComponents() const {
  return values.size();
}

inline
void ActionWithValue::useNumericalDerivatives() {
  plumed_massert( keywords.exists("NUMERICAL_DERIVATIVES"), "numerical derivatives are not permitted for this action" );
  numericalDerivatives=true;
}

inline
bool ActionWithValue::checkNumericalDerivatives() const {
  return numericalDerivatives;
}

inline
bool ActionWithValue::doNotCalculateDerivatives() const {
  return noderiv;
}

inline
unsigned ActionWithValue::getFullNumberOfTasks() const {
  return fullTaskList.size();
}

inline
unsigned ActionWithValue::getTaskCode( const unsigned& ii ) const {
  plumed_dbg_assert( ii<fullTaskList.size() );
  return fullTaskList[ii];
}

inline
bool ActionWithValue::runInSerial() const {
  return serial;
}

inline
bool ActionWithValue::runWithoutOpenMP() const {
  return no_openmp;
}

inline
bool ActionWithValue::actionInChain() const {
  return (action_to_do_before!=NULL);
}

inline
const std::vector<unsigned>& ActionWithValue::getCurrentTasks() const {
  return indexOfTaskInFullList;
}



}

#endif
