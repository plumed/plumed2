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
#include "ActionWithValue.h"
#include "tools/PlumedException.h"

using namespace std;
using namespace PLMD;

void ActionWithValue::registerKeywords(Keywords& keys){
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
}

void ActionWithValue::noAnalyticalDerivatives(Keywords& keys){
   keys.remove("NUMERICAL_DERIVATIVES");
   keys.addFlag("NUMERICAL_DERIVATIVES",true,"analytical derivatives are not implemented for this keyword so numerical derivatives are always used");
}

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  numericalDerivatives(false)
{
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if(numericalDerivatives) log.printf("  using numerical derivatives\n");
}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i)delete values[i];
}

void ActionWithValue::clearInputForces(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearInputForce();
}
void ActionWithValue::clearDerivatives(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearDerivatives();
} 

// -- These are the routine for copying the value pointers to other classes -- //

bool ActionWithValue::exists( const std::string& name ) const {
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==name) return true;
  }
  return false;
}

Value* ActionWithValue::copyOutput( const std::string& name ) const {
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==name) return values[i];
  }
  plumed_merror("there is no pointer with name " + name);
  return NULL;
}

Value* ActionWithValue::copyOutput( const unsigned& n ) const {
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n];
}

// -- HERE WE HAVE THE STUFF FOR THE DEFAULT VALUE -- //

void ActionWithValue::addValue(){
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.push_back(new Value(this,getLabel(), false ) );
}

void ActionWithValue::addValueWithDerivatives(){
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.push_back(new Value(this,getLabel(), true ) );
}

void ActionWithValue::setNotPeriodic(){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=0; values[0]->max=0;
  values[0]->setupPeriodicity();
}

void ActionWithValue::setPeriodic( const std::string& min, const std::string& max ){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->setDomain( min, max );
}

Value* ActionWithValue::getPntrToValue(){
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to retrieve is not the default");
  return values[0];
}

// -- HERE WE HAVE THE STUFF FOR NAMED VALUES / COMPONENTS -- //

void ActionWithValue::addComponent( const std::string& name ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
     plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.push_back(new Value(this,thename, false ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

void ActionWithValue::addComponentWithDerivatives( const std::string& name ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
     plumed_massert(values[i]->name!=thename,"there is already a value with this name");
  }
  values.push_back(new Value(this,thename, true ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

int ActionWithValue::getComponent( const std::string& name ) const {
  plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  plumed_massert(name!=getLabel(),"You should never be calling this routine to retrieve the value");
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     if (values[i]->name==thename) return i;
  }
  plumed_massert(0,"there is no component with name" + name);
  return -1;
}

std::string ActionWithValue::getComponentsList( ) const {
  std::string complist;
  plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  for(unsigned i=0;i<values.size();++i){
     complist+=values[i]->name+" ";
  }
  return complist;
}

void ActionWithValue::componentIsNotPeriodic( const std::string& name ){
  int kk=getComponent(name);
  values[kk]->min=0; values[kk]->max=0;
  values[kk]->setupPeriodicity();
}

void ActionWithValue::setGradientsIfNeeded(){
  if(isOptionOn("GRADIENTS")) {
     for(unsigned i=0;i<values.size();i++) values[i]->setGradients();
  }
}

void ActionWithValue::componentIsPeriodic( const std::string& name, const std::string& min, const std::string& max ){
  int kk=getComponent(name);
  values[kk]->setDomain(min,max);
}

Value* ActionWithValue::getPntrToComponent( const std::string& name ){
  int kk=getComponent(name);
  return values[kk]; 
}

Value* ActionWithValue::getPntrToComponent( int n ){
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n];
}

void ActionWithValue::mergeFieldDerivatives( const std::vector<double>& derivatives, Value* val_out ){
  plumed_assert( derivatives.size()==getNumberOfComponents() );
  for(unsigned i=0;i<derivatives.size();++i){
      for(unsigned j=0;j<getPntrToComponent(i)->getNumberOfDerivatives();++j){
          val_out->addDerivative( j, derivatives[i]*getPntrToComponent(i)->getDerivative(j) );
      }
  }    
}
