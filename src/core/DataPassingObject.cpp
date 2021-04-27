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
#include "tools/OpenMP.h"
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
/// Set the pointer to the value
  void setValuePointer( void* p ) override;
/// Set the pointer to the force
  void setForcePointer( void* p ) override;
/// Share the data and put it in the value from sequential data
  void share_data( const unsigned& j, const unsigned& k, Value* value ) override;
/// Share the data and put it in the value from a scattered data
  void share_data( const std::set<AtomNumber>&index, const std::vector<unsigned>& i, Value* value ) override;
/// Pass the force from the value to the output value
  void add_force( Value* vv ) override;
  void add_force( const std::vector<int>& index, Value* value ) override;
  void add_force( const std::set<AtomNumber>& index, const std::vector<unsigned>& i, Value* value ) override;
/// Rescale the force on the output value
  void rescale_force( const unsigned& n, const double& factor, Value* value ) override;
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
void DataPassingObjectTyped<T>::share_data( const unsigned& j, const unsigned& k, Value* value ) {
  #pragma omp parallel for num_threads(value->getGoodNumThreads(j,k))
  for(unsigned i=j; i<k; ++i) { value->set( i, unit*this->v[i*stride] ); }
}

template <class T>
void DataPassingObjectTyped<T>::share_data( const std::set<AtomNumber>&index, const std::vector<unsigned>& i, Value* value ) {
  // cannot be parallelized with omp because access to data is not ordered
  unsigned k=0; for(const auto & p : index) { value->set( p.index(), unit*this->v[i[k]*stride] ); k++; }
}

template <class T>
void DataPassingObjectTyped<T>::add_force( Value* value ) {
   unsigned nvals=value->getNumberOfValues( value->getName() );
   #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(this->f,nvals))
   for(unsigned i=0;i<nvals;++i) this->f[i*stride] += T(funit*value->getForce(i));
}

template <class T>
void DataPassingObjectTyped<T>::add_force( const std::vector<int>& index, Value* value ) {
   #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(this->f,index.size()))
   for(unsigned i=0; i<index.size(); ++i) this->f[i*stride] += T(funit*value->getForce(index[i])); 
}

template <class T>
void DataPassingObjectTyped<T>::add_force( const std::set<AtomNumber>& index, const std::vector<unsigned>& i, Value* value ) {
   unsigned k=0; for(const auto & p : index) { this->f[stride*i[k]] += T(funit*value->getForce(p.index())); k++; }
}

template <class T>
void DataPassingObjectTyped<T>::rescale_force( const unsigned& n, const double& factor, Value* value ) {
   #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(this->f,n))
   for(unsigned i=0;i<n;++i) this->f[i*stride] *= T(factor);
}
 
}
