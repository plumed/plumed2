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
#include "Value.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"
#include "ActionWithVirtualAtom.h"
#include "tools/Exception.h"
#include "Atoms.h"
#include "PlumedMain.h"

namespace PLMD {

Value::Value():
  action(NULL),
  value_set(false),
  value(0.0),
  inputForce(0.0),
  hasForce(false),
  hasDeriv(true),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
}

Value::Value(ActionWithValue* av, const std::string& name, const bool withderiv):
  action(av),
  value_set(false),
  value(0.0),
  inputForce(0.0),
  hasForce(false),
  name(name),
  hasDeriv(withderiv),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{
}

void Value::setupPeriodicity() {
  if( min==0 && max==0 ) {
    periodicity=notperiodic;
  } else {
    periodicity=periodic;
    max_minus_min=max-min;
    plumed_massert(max_minus_min>0, "your function has a very strange domain?");
    inv_max_minus_min=1.0/max_minus_min;
  }
}

bool Value::isPeriodic()const {
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

bool Value::applyForce(std::vector<double>& forces ) const {
  if( !hasForce ) return false;
  plumed_dbg_massert( derivatives.size()==forces.size()," forces array has wrong size" );
  const unsigned N=derivatives.size();
  for(unsigned i=0; i<N; ++i) forces[i]=inputForce*derivatives[i];
  return true;
}

void Value::setNotPeriodic() {
  min=0; max=0; periodicity=notperiodic;
}

void Value::setDomain(const std::string& pmin,const std::string& pmax) {
  str_min=pmin;
  if( !Tools::convert(str_min,min) ) action->error("could not convert period string " + str_min + " to real");
  str_max=pmax;
  if( !Tools::convert(str_max,max) ) action->error("could not convert period string " + str_max + " to read");
  setupPeriodicity();
}

void Value::getDomain(std::string&minout,std::string&maxout) const {
  plumed_massert(periodicity==periodic,"function should be periodic");
  minout=str_min;
  maxout=str_max;
}

void Value::getDomain(double&minout,double&maxout) const {
  plumed_massert(periodicity==periodic,"function should be periodic");
  minout=min;
  maxout=max;
}

void Value::setGradients() {
  // Can't do gradients if we don't have derivatives
  if( !hasDeriv ) return;
  gradients.clear();
  ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(action);
  ActionWithArguments*aw=dynamic_cast<ActionWithArguments*>(action);
  if(aa) {
    Atoms&atoms((aa->plumed).getAtoms());
    for(unsigned j=0; j<aa->getNumberOfAtoms(); ++j) {
      AtomNumber an=aa->getAbsoluteIndex(j);
      if(atoms.isVirtualAtom(an)) {
        const ActionWithVirtualAtom* a=atoms.getVirtualAtomsAction(an);
        for(const auto & p : a->getGradients()) {
// controllare l'ordine del matmul:
          gradients[p.first]+=matmul(Vector(derivatives[3*j],derivatives[3*j+1],derivatives[3*j+2]),p.second);
        }
      } else {
        for(unsigned i=0; i<3; i++) gradients[an][i]+=derivatives[3*j+i];
      }
    }
  } else if(aw) {
    std::vector<Value*> values=aw->getArguments();
    for(unsigned j=0; j<derivatives.size(); j++) {
      for(const auto & p : values[j]->gradients) {
        AtomNumber iatom=p.first;
        gradients[iatom]+=p.second*derivatives[j];
      }
    }
  } else plumed_error();
}

double Value::projection(const Value& v1,const Value&v2) {
  double proj=0.0;
  const std::map<AtomNumber,Vector> & grad1(v1.gradients);
  const std::map<AtomNumber,Vector> & grad2(v2.gradients);
  for(const auto & p1 : grad1) {
    AtomNumber a=p1.first;
    const auto p2=grad2.find(a);
    if(p2!=grad2.end()) {
      proj+=dotProduct(p1.second,(*p2).second);
    }
  }
  return proj;
}

ActionWithValue* Value::getPntrToAction() {
  plumed_assert( action!=NULL );
  return action;
}

void copy( const Value& val1, Value& val2 ) {
  unsigned nder=val1.getNumberOfDerivatives();
  if( nder!=val2.getNumberOfDerivatives() ) { val2.resizeDerivatives( nder ); }
  val2.clearDerivatives();
  for(unsigned i=0; i<val1.getNumberOfDerivatives(); ++i) val2.addDerivative( i, val1.getDerivative(i) );
  val2.set( val1.get() );
}

void copy( const Value& val1, Value* val2 ) {
  unsigned nder=val1.getNumberOfDerivatives();
  if( nder!=val2->getNumberOfDerivatives() ) { val2->resizeDerivatives( nder ); }
  val2->clearDerivatives();
  for(unsigned i=0; i<val1.getNumberOfDerivatives(); ++i) val2->addDerivative( i, val1.getDerivative(i) );
  val2->set( val1.get() );
}

void add( const Value& val1, Value* val2 ) {
  plumed_assert( val1.getNumberOfDerivatives()==val2->getNumberOfDerivatives() );
  for(unsigned i=0; i<val1.getNumberOfDerivatives(); ++i) val2->addDerivative( i, val1.getDerivative(i) );
  val2->set( val1.get() + val2->get() );
}

}
