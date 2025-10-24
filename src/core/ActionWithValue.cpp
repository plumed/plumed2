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
#include "ActionWithValue.h"
#include "ActionWithArguments.h"
#include "ActionAtomistic.h"
#include "tools/Exception.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"
#include "blas/blas.h"

namespace PLMD {

void ActionWithValue::registerKeywords(Keywords& keys) {
  keys.setComponentsIntroduction("By default the value of the calculated quantity can be referenced elsewhere in the "
                                 "input file by using the label of the action.  Alternatively this Action can be used "
                                 "to calculate the following quantities by employing the keywords listed "
                                 "below.  These quantities can be referenced elsewhere in the input by using this Action's "
                                 "label followed by a dot and the name of the quantity required from the list below.");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.add("hidden","HAS_VALUES","this is used in json output to determine those actions that have values");
}

void ActionWithValue::noAnalyticalDerivatives(Keywords& keys) {
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addFlag("NUMERICAL_DERIVATIVES",false,"analytical derivatives are not implemented for this keyword so numerical derivatives are always used");
}

void ActionWithValue::useCustomisableComponents(Keywords& keys) {
  if( !keys.outputComponentExists(".#!custom") ) {
    keys.addOutputComponent(".#!custom","default","scalar","the names of the output components for this action depend on the actions input file see the example inputs below for details");
  }
  keys.setComponentsIntroduction("The names of the components in this action can be customized by the user in the "
                                 "actions input file.  However, in addition to the components that can be customized the "
                                 "following quantities will always be output");
}

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  firststep(true),
  noderiv(true),
  numericalDerivatives(false) {
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) {
    parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  }
  if(!keywords.exists("NO_ACTION_LOG") && numericalDerivatives) {
    log.printf("  using numerical derivatives\n");
  }
}

ActionWithValue::~ActionWithValue() {
// empty destructor to delete unique_ptr
}

void ActionWithValue::clearInputForces( const bool& force ) {
  for(unsigned i=0; i<values.size(); i++) {
    values[i]->clearInputForce();
  }
}

void ActionWithValue::clearDerivatives( const bool& force ) {
#ifdef _OPENMP
  //nt is unused if openmp is not declared
  const unsigned nt=OpenMP::getNumThreads();
#endif //_OPENMP
  #pragma omp parallel num_threads(nt)
  {
    #pragma omp for
    for(unsigned i=0; i<values.size(); i++) {
      values[i]->clearDerivatives();
    }
  }
}

// -- These are the routine for copying the value pointers to other classes -- //

bool ActionWithValue::exists( const std::string& valname ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==valname) {
      return true;
    }
  }
  return false;
}

void ActionWithValue::getMatrixColumnTitles( std::vector<std::string>& argnames ) const {
  plumed_assert( getNumberOfComponents()==1 && getConstPntrToComponent(0)->getRank()==2 );
  unsigned nargs = getConstPntrToComponent(0)->getShape()[1];
  std::string aname = getConstPntrToComponent(0)->getName();
  for(unsigned j=0; j<nargs; ++j) {
    std::string nn;
    Tools::convert( j+1, nn );
    argnames.push_back( aname + "." + nn );
  }
}

Value* ActionWithValue::copyOutput( const std::string& valname ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==valname) {
      return values[i].get();
    }
  }
  plumed_merror("there is no pointer with name " + valname);
}

Value* ActionWithValue::copyOutput( const unsigned& n ) const {
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

// -- HERE WE HAVE THE STUFF FOR THE DEFAULT VALUE -- //

void ActionWithValue::addValue( const std::vector<std::size_t>& shape ) {
  if( !keywords.outputComponentExists(".#!value") ) {
    warning("documentation for the value calculated by this action has not been included");
  } else {
    plumed_massert( keywords.componentHasCorrectType(".#!value",shape.size(),false), "documentation for type of value is incorrect");
  }
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(Tools::make_unique<Value>(this,getLabel(), false, shape ) );
}

void ActionWithValue::addValueWithDerivatives( const std::vector<std::size_t>& shape ) {
  if( !keywords.outputComponentExists(".#!value") ) {
    warning("documentation for the value calculated by this action has not been included");
  } else {
    plumed_massert( keywords.componentHasCorrectType(".#!value",shape.size(),true), "documentation for type of value is incorrect");
  }
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(Tools::make_unique<Value>(this,getLabel(), true, shape ) );
}

void ActionWithValue::setNotPeriodic() {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=0;
  values[0]->max=0;
  values[0]->setupPeriodicity();
}

void ActionWithValue::setPeriodic( const std::string& min, const std::string& max ) {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->setDomain( min, max );
}

// -- HERE WE HAVE THE STUFF FOR NAMED VALUES / COMPONENTS -- //

void ActionWithValue::addComponent( const std::string& valname, const std::vector<std::size_t>& shape ) {
  if( !keywords.outputComponentExists(valname) ) {
    plumed_merror("a description of component " + valname + " has not been added to the manual. Components should be registered like keywords in "
                  "registerKeywords as described in the developer docs.");
  }
  plumed_massert( keywords.componentHasCorrectType(valname,shape.size(),false), "documentation for type of component " + valname + " is incorrect");
  std::string thename;
  thename=getLabel() + "." + valname;
  for(unsigned i=0; i<values.size(); ++i) {
    plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name: "+thename);
    plumed_massert(values[i]->name!=thename&&valname!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponent(\"bias\") from your bias.");
  }
  values.emplace_back(Tools::make_unique<Value>(this,thename, false, shape ) );
  log.printf("  added component to this action: %s \n", thename.c_str() );
}

void ActionWithValue::addComponentWithDerivatives( const std::string& valname, const std::vector<std::size_t>& shape ) {
  if( !keywords.outputComponentExists(valname) ) {
    plumed_merror("a description of component " + valname + " has not been added to the manual. Components should be registered like keywords in "
                  "registerKeywords as described in the developer doc.");
  }
  plumed_massert( keywords.componentHasCorrectType(valname,shape.size(),true), "documentation for type of component " + valname + " is incorrect");
  std::string thename;
  thename=getLabel() + "." + valname;
  for(unsigned i=0; i<values.size(); ++i) {
    plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name: "+thename);
    plumed_massert(values[i]->name!=thename&&valname!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponentWithDerivatives(\"bias\") from your bias.");
  }
  values.emplace_back(Tools::make_unique<Value>(this,thename, true, shape ) );
  log.printf("  added component to this action: %s \n", thename.c_str() );
}

std::string ActionWithValue::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( keys.outputComponentExists(".#!custom") ) {
    return "a quantity calculated by the action " + getName() + " with label " + getLabel();
  }
  std::size_t und=cname.find_last_of("_");
  std::size_t hyph=cname.find_first_of("-");
  if( und!=std::string::npos ) {
    return keys.getOutputComponentDescription(cname.substr(und)) + " This particular component measures this quantity for the input CV named " + cname.substr(0,und);
  }
  if( hyph!=std::string::npos ) {
    return keys.getOutputComponentDescription(cname.substr(0,hyph)) + "  This is the " + cname.substr(hyph+1) + "th of these quantities";
  }
  plumed_massert( keys.outputComponentExists(cname), "component " + cname + " does not exist in " + keys.getDisplayName() + " if the component names are customizable then you should override this function" );
  return keys.getOutputComponentDescription( cname );
}

int ActionWithValue::getComponent( const std::string& valname ) const {
  plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  std::string thename;
  thename=getLabel() + "." + valname;
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==thename) {
      return i;
    }
  }
  plumed_merror("there is no component with name " + valname);
}

std::string ActionWithValue::getComponentsList( ) const {
  std::string complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist+=values[i]->name+" ";
  }
  return complist;
}

std::vector<std::string> ActionWithValue::getComponentsVector( ) const {
  std::vector<std::string> complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist.push_back(values[i]->name);
  }
  return complist;
}

void ActionWithValue::componentIsNotPeriodic( const std::string& valname ) {
  int kk=getComponent(valname);
  values[kk]->min=0;
  values[kk]->max=0;
  values[kk]->setupPeriodicity();
}

void ActionWithValue::componentIsPeriodic( const std::string& valname, const std::string& min, const std::string& max ) {
  int kk=getComponent(valname);
  values[kk]->setDomain(min,max);
}

void ActionWithValue::setGradientsIfNeeded() {
  if(isOptionOn("GRADIENTS")) {
    ActionAtomistic* aa=castToActionAtomistic();
    if(aa) {
      for(unsigned i=0; i<values.size(); i++) {
        unsigned start=0;
        values[i]->gradients.clear();
        values[i]->setGradients( aa, start );
      }
    } else {
      ActionWithArguments* aarg = castToActionWithArguments();
      if( !aarg ) {
        plumed_merror( "failing in " + getLabel() );
      }
      for(unsigned i=0; i<values.size(); i++) {
        unsigned start=0;
        values[i]->gradients.clear();
        aarg->setGradients( values[i].get(), start );
      }
    }
  }
}

void ActionWithValue::turnOnDerivatives() {
  // Turn on the derivatives
  noderiv=false;
  // Resize the derivatives
  for(unsigned i=0; i<values.size(); ++i) {
    values[i]->resizeDerivatives( getNumberOfDerivatives() );
  }
  // And turn on the derivatives in all actions on which we are dependent
  for(unsigned i=0; i<getDependencies().size(); ++i) {
    ActionWithValue* vv=getDependencies()[i]->castToActionWithValue();
    if(vv) {
      vv->turnOnDerivatives();
    }
  }
}

Value* ActionWithValue::getPntrToComponent( const std::string& valname ) {
  int kk=getComponent(valname);
  return values[kk].get();
}

const Value* ActionWithValue::getConstPntrToComponent(unsigned n) const {
  plumed_dbg_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

Value* ActionWithValue::getPntrToComponent( unsigned n ) {
  plumed_dbg_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

bool ActionWithValue::calculateOnUpdate() {
  if( firststep ) {
    ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>(this);
    if(aa) {
      const std::vector<Value*> & args(aa->getArguments());
      for(const auto & p : args ) {
        if( p->calculateOnUpdate() ) {
          for(unsigned i=0; i<values.size(); ++i) {
            values[i]->setValType("calcFromAverage");
          }
          break;
        }
      }
    }
    firststep=false;
  }
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->calculateOnUpdate() ) {
      return true;
    }
  }
  return false;
}

bool ActionWithValue::checkForForces() {
  const unsigned    ncp=getNumberOfComponents();
  unsigned    nder=getNumberOfDerivatives();
  if( ncp==0 || nder==0 ) {
    return false;
  }

  unsigned nvalsWithForce=0;
  valsToForce.resize(ncp);
  for(unsigned i=0; i<ncp; ++i) {
    if( values[i]->hasForce && !values[i]->isConstant() ) {
      valsToForce[nvalsWithForce]=i;
      nvalsWithForce++;
    }
  }
  if( nvalsWithForce==0 ) {
    return false;
  }

  // Make sure forces to apply is empty of forces
  if( forcesForApply.size()!=nder ) {
    forcesForApply.resize( nder );
  }
  std::fill(forcesForApply.begin(),forcesForApply.end(),0);

  unsigned stride=1;
  unsigned rank=0;
  if(ncp>static_cast<unsigned>(4*comm.Get_size())) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();
  if(nt>ncp/(4*stride)) {
    nt=1;
  }

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_f;
    if( nt>1 ) {
      omp_f.resize(nder,0);
    }
    #pragma omp for
    for(unsigned i=rank; i<nvalsWithForce; i+=stride) {
      double ff=values[valsToForce[i]]->inputForce[0];
      std::vector<double> & thisderiv( values[valsToForce[i]]->data );
      int nn=nder;
      int one1=1;
      int one2=1;
      if( nt>1 ) {
        plumed_blas_daxpy(&nn,&ff,thisderiv.data()+1,&one1,omp_f.data(),&one2);
      } else {
        plumed_blas_daxpy(&nn,&ff,thisderiv.data()+1,&one1,forcesForApply.data(),&one2);
      }
      // if( nt>1 ) for(unsigned j=0; j<nder; ++j) omp_f[j] += ff*thisderiv[1+j];
      //else for(unsigned j=0; j<nder; ++j) forcesForApply[j] += ff*thisderiv[1+j];
    }
    #pragma omp critical
    {
      if( nt>1 ) {
        int nn=forcesForApply.size();
        double one0=1.0;
        int one1=1;
        int one2=1;
        plumed_blas_daxpy(&nn,&one0,omp_f.data(),&one1,forcesForApply.data(),&one2);
      }
      // for(unsigned j=0; j<forcesForApply.size(); ++j) {
      // forcesForApply[j]+=omp_f[j];
      // }
    }
  }

  if(ncp>static_cast<unsigned>(4*comm.Get_size())) {
    comm.Sum(&forcesForApply[0],nder);
  }
  return true;
}

}
