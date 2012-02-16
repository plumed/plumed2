#ifndef __PLUMED_ActionWithValue_h
#define __PLUMED_ActionWithValue_h

#include "Action.h"
#include "Value.h"
#include "PlumedException.h"
#include <vector>

namespace PLMD{

//+DEVELDOC MULTI-INHERIT ActionWithValue
/**
Used to create a PLMD::Action that has some scalar or vectorial output that may or may not have some derivatives.  
This is used for PLMD::Bias, PLMD::Colvar and PLMD::Function 
*/
//+ENDDEVELDOC

/// The vast majority of the PLMD::Action objects that are implemented in
/// plumed calculate some quantity or a set of quantities.  This could be
/// the value of a CV, the value of a function or the potential due to a bias.
/// PLMD::ActionWithValue provides the functionality for storing these quantities
/// and (in tandem with PLMD::ActionWithArguments) the functionality for passing
/// quantities between PLMD::Actions.  When you are deciding what quantities
/// your new PLMD::Action will need to store using PLMD::ActionWithValue you must
/// ask yourself the following two questions:

/// - Do I need to differentiate my output quantities 
/// - Is my PLMD::Action calculating a single thing or does the output have multiple components  

/// If the answer to the first of these questions is yes then you must setup your values
/// you using either PLMD::ActionWithValue::addValueWithDerivatives() or 
/// PLMD::ActionWithValue::addComponentWithDerivatives.  If the answer is no you 
/// can set up values using PLMD::ActionWithValue::addValue() or PLMD::ActionWithValue::addComponent().
/// The precise routine you use to setup your values will depend on your answer to the
/// second question.  As you are probably aware if the output of your PLMD::Action is a 
/// single quantity you can reference that quantity in the input file using the label of the 
/// PLMD::Action it was calculated in.  If your action <b> outputs only one quantity </b> 
/// we call that quantity the <b> value </b> of the Action.  To set the <b> value </b> and get pointers to it 
/// you should <b> use the set of routines that have the word value in the name </b>.  If, by contrast, 
/// your PLMD::Action calculates multiple quantities then these quantities are referenced in input using the
/// label.component syntax.  We refer to these <b> multiple quantities </b> the <b> components </b> 
/// of the PLMD::Action.  Perhaps unsurprisingly, when you manipulate the <b> components </b> of an
/// PLMD::Action you should use <b> the routines with the word component in the name. </b>

class ActionWithValue : 
  public virtual Action
{
friend class ActionWithDistribution;
private:
/// An array containing the values for this action
  std::vector<Value*> values;
/// Are we using numerical derivatives to differentiate
  bool numericalDerivatives;
/// Return the index for the component named name
  int getComponent( const std::string& name ) const;
protected:

// -------- The action has one value only  ---------------- //

/// Add a value with the name <label>
  void addValue();
/// Add a value with the name <label> that has derivatives
  void addValueWithDerivatives();
/// Set your default value to have no periodicity
  void setNotPeriodic();
/// Set the value to be periodic with a particular domain
  void setPeriodic( const double min, const double max );
/// Get a pointer to the default value
  Value* getPntrToValue();
/// Set the default value (the one without name)
  void setValue(const double& d);

// -------- The action has multiple components ---------- //

/// Add a value with a name like <label>.name
  void addComponent( const std::string& name );
/// Add a value with a name like <label>.name that has derivatives
  void addComponentWithDerivatives( const std::string& name );
/// Set your value component to have no periodicity
  void componentIsNotPeriodic( const std::string& name );
/// Set the value to be periodic with a particular domain
  void componentIsPeriodic( const std::string& name, const double min, const double max );
/// Return a pointer to the component by index
  Value* getPntrToComponent(int i);
/// Return a pointer to the value by name
  Value* getPntrToComponent(const std::string& name);
public:
  ActionWithValue(const ActionOptions&ao);
  ~ActionWithValue();

/// Register all the relevant keywords for the action  
  static void registerKeywords( Keywords& keys );
/// Insist that numerical derivatives should always be used for an action and make this fact appear in the manual
  static void noAnalyticalDerivatives(Keywords& keys);
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

// -- Routines for everything else -- //

/// Returns the number of values defined
  int getNumberOfComponents() const ;
/// Clear the forces on the values
  void clearInputForces();
/// Clear the derivatives of values wrt parameters
  void clearDerivatives();
/// Calculate the gradients and store them for all the values (need for projections) 
  void setGradientsIfNeeded();
/// Set the value
  void setValue(Value*,double);
/// Check if numerical derivatives should be used
  bool checkNumericalDerivatives() const ;
};

inline
double ActionWithValue::getOutputQuantity(const unsigned j) const {
  plumed_massert(j<values.size(),"index requested is out of bounds");
  return values[j]->get();
}

inline
double ActionWithValue::getOutputQuantity( const std::string& name ) const {
  std::string thename; thename=getLabel() + "." + name;
  if( exists(thename) ){
      for(unsigned i=0;i<values.size();++i){
         if( values[i]->name==thename ) return values[i]->value;
      }
  }
  return 0.0;
}

inline
void ActionWithValue::setValue(const double& d){
  plumed_massert(values.size()==1, "cannot use setValue in multi-component actions");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->set(d);
}

inline
int ActionWithValue::getNumberOfComponents() const {
  return values.size();
}

inline
bool ActionWithValue::checkNumericalDerivatives() const {
  return numericalDerivatives;
}


}

#endif
