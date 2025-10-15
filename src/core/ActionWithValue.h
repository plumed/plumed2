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
#ifndef __PLUMED_core_ActionWithValue_h
#define __PLUMED_core_ActionWithValue_h

#include "Action.h"
#include "Value.h"
#include "tools/Exception.h"
#include <vector>
#include <memory>
#include <cstring>

namespace PLMD {

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

class ActionWithValue :
  public virtual Action {
  friend class ActionWithVector;
  friend class ActionWithArguments;
private:
/// This finishes setup on first step to check if actions are calculated during update
  bool firststep;
/// An array containing the values for this action
  std::vector<std::unique_ptr<Value>> values;
/// A vector that is used to hold the forces that we will apply on the input quantities
  std::vector<double> forcesForApply;
  std::vector<unsigned> valsToForce;
/// Are we skipping the calculation of the derivatives
  bool noderiv;
/// Are we using numerical derivatives to differentiate
  bool numericalDerivatives;
/// Return the index for the component named name
  int getComponent( const std::string& valname ) const;
public:

// -------- The action has one value only  ---------------- //

/// Add a value with the name label
  void addValue( const std::vector<std::size_t>& shape=std::vector<std::size_t>() );
/// Add a value with the name label that has derivatives
  virtual void addValueWithDerivatives( const std::vector<std::size_t>& shape=std::vector<std::size_t>() );
/// Set your default value to have no periodicity
  void setNotPeriodic();
/// Set the value to be periodic with a particular domain
  void setPeriodic( const std::string& min, const std::string& max );
protected:
/// Get a pointer to the default value
  Value* getPntrToValue();
/// Set the default value (the one without name)
  void setValue(const double& d);

// -------- The action has multiple components ---------- //

public:
/// Add a value with a name like label.name
  void addComponent( const std::string& valname, const std::vector<std::size_t>& shape=std::vector<std::size_t>() );
/// Add a value with a name like label.name that has derivatives
  virtual void addComponentWithDerivatives( const std::string& valname, const std::vector<std::size_t>& shape=std::vector<std::size_t>() );
/// Set your value component to have no periodicity
  void componentIsNotPeriodic( const std::string& valname );
/// Set the value to be periodic with a particular domain
  void componentIsPeriodic( const std::string& valname, const std::string& min, const std::string& max );
/// Get the description of this component
  virtual std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const ;
/// Get a const pointer to the ith component
  const Value* getConstPntrToComponent(unsigned i) const;
protected:
/// Return a pointer to the component by index
  Value* getPntrToComponent(unsigned i);
/// Return a pointer to the value by name
  Value* getPntrToComponent(const std::string& valname);
/// Accumulate the forces from the Values
  bool checkForForces();
/// Get the forces to apply
  const std::vector<double>& getForcesToApply() const;
public:
  explicit ActionWithValue(const ActionOptions&ao);
  ~ActionWithValue();

/// Register all the relevant keywords for the action
  static void registerKeywords( Keywords& keys );
/// Insist that numerical derivatives should always be used for an action and make this fact appear in the manual
  static void noAnalyticalDerivatives(Keywords& keys);
/// The components in the action will depend on the user
  static void useCustomisableComponents(Keywords& keys);
/// Are we not calculating derivatives
  virtual bool doNotCalculateDerivatives() const ;
/// Get the value of one of the components of the PLMD::Action
  double getOutputQuantity( const unsigned j ) const ;
/// Get the value with a specific name (N.B. if there is no such value this returns zero)
  double getOutputQuantity( const std::string& valname ) const ;

//  --- Routines for passing stuff to ActionWithArguments -- //

/// Check if a value with a particular name is present.  This is only used in PLMD::ActionWithArguments.
/// You should not use it when manipulating components.
  bool exists( const std::string& valname ) const;
/// Return a pointer to the value with name (this is used to retrieve values in other PLMD::Actions)
/// You should NEVER use this routine to refer to the components of your PLMD::Action.  Use
/// getPntrToComponent instead.
  Value* copyOutput( const std::string& valname ) const;
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
  std::size_t getNumberOfComponents() const ;
/// Clear the forces on the values
  virtual void clearInputForces( const bool& force=false );
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
  virtual void checkFieldsAllowed() {
    error("cannot use this action as a field");
  }
  virtual unsigned getNumberOfDerivatives()=0;
/// Activate the calculation of derivatives
  virtual void turnOnDerivatives();
/// Get the titles to use for the columns of the matrix
  virtual void getMatrixColumnTitles( std::vector<std::string>& argnames ) const ;
/// This is used to check if we run calculate during the update step
  virtual bool calculateOnUpdate();
  ActionWithValue* castToActionWithValue() noexcept final {
    return this;
  }
};

inline
double ActionWithValue::getOutputQuantity(const unsigned j) const {
  plumed_massert(j<values.size(),"index requested is out of bounds");
  return values[j]->get();
}

inline
double ActionWithValue::getOutputQuantity( const std::string& quantityName ) const {
  auto offset=getLabel().size();
  for(unsigned i=0; i<values.size(); ++i) {
    const std::string & valname=values[i]->name;
    if(valname.size()>offset+1 && valname[offset]=='.' ) {
      plumed_dbg_assert(Tools::startWith(valname,getLabel()));
      if(!std::strcmp(valname.c_str()+offset+1,quantityName.c_str())) {
        return values[i]->get();
      }
    }
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
std::size_t ActionWithValue::getNumberOfComponents() const {
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
const std::vector<double>& ActionWithValue::getForcesToApply() const {
  return forcesForApply;
}

inline
Value* ActionWithValue::getPntrToValue() {
  plumed_dbg_massert(values.size()==1,"The number of components is not equal to one");
  plumed_dbg_massert(values[0]->name==getLabel(), "The value you are trying to retrieve is not the default");
  return values[0].get();
}

}

#endif
