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
#include "ValueFromMDCode.h"

namespace PLMD {

ValueFromMDCode::ValueFromMDCode(const std::string& name, const std::vector<unsigned>& shape):
fixed(false),
collect(false),
set(false),
domain_decomposed(false)
{
  value=std::unique_ptr<Value>(new Value(NULL, name, false, shape));
  value->created_in_plumedmain=true;
  value->setShape( shape ); 
}

void ValueFromMDCode::setupPeriodicity( const bool& isperiodic, const std::string& min, const std::string& max ) {
  if( isperiodic ) { value->setDomain( min, max ); } 
  else { value->min=0; value->max=0; value->setupPeriodicity(); } 
}

template <class T>
class ValueFromMDCodeTyped : public ValueFromMDCode {
private:
/// A pointer to the value in the MD code that we get the value from 
  T* v;
/// A pointer to the value in the MD code that we add the force to
  T* f;
public:
  explicit ValueFromMDCodeTyped(const std::string& name, const std::vector<unsigned>& shape);
  void setv(void*p) override; 
  void setf(void*f) override;
  void gather() override;
  void updateForces() override;
};

template <class T>
ValueFromMDCodeTyped<T>::ValueFromMDCodeTyped(const std::string& name, const std::vector<unsigned>& shape):
ValueFromMDCode(name,shape)
{
}

template <class T>
void ValueFromMDCodeTyped<T>::setv(void*p) {
  // Clear the forces on this guy
  value->clearInputForce();
  // Note that it has been set
  set=true;
  // Set the pointer so we have it
  this->v=static_cast<T*>(p); 
}

template <class T>
void ValueFromMDCodeTyped<T>::gather() {
  if( !collect || domain_decomposed ) return;
  unsigned nvals = value->getSize(); 
  for(unsigned i=0;i<nvals;++i) value->set( i, this->v[i] );
} 

template <class T>
void ValueFromMDCodeTyped<T>::setf(void*f) { 
  // Set the pointer so we have it
  this->f=static_cast<T*>(f); 
}

template <class T>
void ValueFromMDCodeTyped<T>::updateForces(){
  if( fixed || !value->forcesWereAdded() ) return;
  // This adds forces if we more than one value
  unsigned nvals=value->getNumberOfValues( value->getName() );
  for(unsigned i=0;i<nvals;++i) this->f[i] += T(value->getForce(i));
}

std::unique_ptr<ValueFromMDCode> ValueFromMDCode::create(unsigned p, const std::string& name, const std::vector<unsigned>& shape ) {
  if(p==sizeof(double)) {
    return std::unique_ptr<ValueFromMDCodeTyped<double>>(new ValueFromMDCodeTyped<double>(name,shape));
  } else if (p==sizeof(float)) {
    return std::unique_ptr<ValueFromMDCodeTyped<float>>(new ValueFromMDCodeTyped<float>(name,shape));
  }
  std::string pp; Tools::convert(p,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

}
