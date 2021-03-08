/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "DataPassingObject.h"
#include "tools/Tools.h"

namespace PLMD {

template <class T>
class DataPassingObjectTyped : public DataPassingObject {
private:
/// A pointer to the value
  T* v;
/// A pointer to the force 
  T* f;
public:
  static std::unique_ptr<DataPassingObject> create(unsigned n);
/// Set the pointer to the value
  void setValuePointer( void* p ) override;
/// Set the pointer to the force
  void setForcePointer( void* p ) override;
/// Share the data and put it in the value
  void share_data( Value* vv ) override;
/// Pass the force from the value to the output value
  void set_force( Value* vv ) override;
/// This transfers everything to the output
  void setData( const std::vector<double>& data ) override;
};

std::unique_ptr<DataPassingObject> DataPassingObject::create(unsigned n) {
  if(n==sizeof(double)) {
    return std::unique_ptr<DataPassingObject>(new DataPassingObjectTyped<double>());
  } else  if(n==sizeof(float)) {
    return std::unique_ptr<DataPassingObject>(new DataPassingObjectTyped<float>());
  }
  std::string pp; Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

template <class T>
void DataPassingObjectTyped<T>::setValuePointer( void* p ) {
   v=static_cast<T*>(p); 
}

template <class T>
void DataPassingObjectTyped<T>::setForcePointer( void* p ) {
   f=static_cast<T*>(p);
}

template <class T>
void DataPassingObjectTyped<T>::setData( const std::vector<double>& data ) {
   for(unsigned i=0;i<data.size();++i) v[i] = static_cast<T>( data[i] );
}

template <class T>
void DataPassingObjectTyped<T>::share_data( Value* value ) {
   unsigned nvals = value->getSize();
   for(unsigned i=0;i<nvals;++i) value->set( i, unit*this->v[i] );
}

template <class T>
void DataPassingObjectTyped<T>::set_force( Value* value ) {
   unsigned nvals=value->getNumberOfValues( value->getName() );
   for(unsigned i=0;i<nvals;++i) this->f[i] = T(value->getForce(i));
}
 
}
