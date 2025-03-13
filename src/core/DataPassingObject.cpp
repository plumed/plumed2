/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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

template<typename T>
static void getPointer(const TypesafePtr & p, const std::vector<std::size_t>& shape, const unsigned& start, const unsigned& stride, T*&pp ) {
  if( shape.size()==1 && stride==1 ) {
    pp=p.get<T*>( {shape[0]} );
  } else if( shape.size()==1 ) {
    auto p_=p.get<T*>( {shape[0],stride} );
    pp = p_+start;
  } else if( shape.size()==2 ) {
    pp=p.get<T*>( {shape[0],shape[1]} );
  }
}

template <class T>
class DataPassingObjectTyped : public DataPassingObject {
private:
/// A pointer to the value
  TypesafePtr v;
/// A pointer to the force
  TypesafePtr f;
public:
/// This convers a number from the MD code into a double
  double MD2double(const TypesafePtr &) const override ;
/// This is used when you want to save the passed object to a double variable in PLUMED rather than the pointer
/// this can be used even when you don't pass a pointer from the MD code
  void saveValueAsDouble( const TypesafePtr & val ) override;
/// Set the pointer to the value
  void setValuePointer( const TypesafePtr & p, const std::vector<std::size_t>& shape, const bool& isconst ) override;
/// Set the pointer to the force
  void setForcePointer( const TypesafePtr & p, const std::vector<std::size_t>& shape ) override;
/// This gets the data in the pointer and passes it to the output value
  void share_data( std::vector<double>& values ) const override ;
/// Share the data and put it in the value from sequential data
  void share_data( const unsigned& j, const unsigned& k, Value* value ) override;
/// Share the data and put it in the value from a scattered data
  void share_data( const std::vector<AtomNumber>&index, const std::vector<unsigned>& i, Value* value ) override;
/// Pass the force from the value to the output value
  void add_force( Value* vv ) override;
  void add_force( const std::vector<int>& index, Value* value ) override;
  void add_force( const std::vector<AtomNumber>& index, const std::vector<unsigned>& i, Value* value ) override;
/// Rescale the force on the output value
  void rescale_force( const unsigned& n, const double& factor, Value* value ) override;
/// This transfers everything to the output
  void setData( Value* value ) override;
};

std::unique_ptr<DataPassingObject> DataPassingObject::create(unsigned n) {
  if(n==sizeof(double)) {
    return std::make_unique<DataPassingObjectTyped<double>>();
  } else  if(n==sizeof(float)) {
    return std::make_unique<DataPassingObjectTyped<float>>();
  }
  std::string pp;
  Tools::convert(n,pp);
  plumed_merror("cannot create an MD interface with sizeof(real)=="+ pp);
  return NULL;
}

template <class T>
double DataPassingObjectTyped<T>::MD2double(const TypesafePtr & m) const {
  double d=double(m.template get<T>());
  return d;
}

template <class T>
void DataPassingObjectTyped<T>::saveValueAsDouble( const TypesafePtr & val ) {
  hasbackup=true;
  bvalue=double(val.template get<T>());
}

template <class T>
void DataPassingObjectTyped<T>::setValuePointer( const TypesafePtr & val, const std::vector<std::size_t>& shape, const bool& isconst ) {
  if( shape.size()==0 ) {
    if( isconst ) {
      val.get<const T*>();
    } else {
      val.get<T*>();  // just check type and discard pointer
    }
  } else if( shape.size()==1 ) {
    if( isconst )
      val.get<const T*>({shape[0]});
    else
      val.get<T*>({shape[0]});  // just check type and discard pointer
  } else if( shape.size()==2 ) {
    if( isconst )
      val.get<const T*>({shape[0],shape[1]});
    else
      val.get<T*>({shape[0],shape[1]});  // just check type and discard pointer
  }
  v=val.copy();
}

template <class T>
void DataPassingObjectTyped<T>::setForcePointer( const TypesafePtr & val, const std::vector<std::size_t>& shape ) {
  if( shape.size()==0 ) {
    val.get<T*>();  // just check type and discard pointer
  } else if( shape.size()==1 )
    val.get<T*>({shape[0]});   // just check type and discard pointer
  else if( shape.size()==2 )
    val.get<T*>({shape[0],shape[1]});   // just check type and discard pointer
  f=val.copy();
}

template <class T>
void DataPassingObjectTyped<T>::setData( Value* value ) {
  if( value->getRank()==0 ) {
    *v.template get<T*>() = static_cast<T>(value->get()) / unit;
    return;
  }
  T* pp;
  getPointer( v, value->getShape(), start, stride, pp );
  unsigned nvals=value->getNumberOfValues();
  for(unsigned i=0; i<nvals; ++i) {
    pp[i] = T( value->get(i) );
  }
}

template <class T>
void DataPassingObjectTyped<T>::share_data( const unsigned& j, const unsigned& k, Value* value ) {
  if( value->getRank()==0 ) {
    if( hasbackup ) {
      value->set( unit*bvalue );
    } else {
      value->set( unit*double(v.template get<T>()) );
    }
    return;
  }
  std::vector<std::size_t> s(value->getShape());
  if( s.size()==1 ) {
    s[0]=k-j;
  }
  const T* pp;
  getPointer( v, s, start, stride, pp );
  std::vector<double> & d=value->data;
  #pragma omp parallel for num_threads(value->getGoodNumThreads(j,k))
  for(unsigned i=j; i<k; ++i) {
    d[i]=unit*pp[i*stride];
  }
}

template <class T>
void DataPassingObjectTyped<T>::share_data( std::vector<double>& values ) const {
  std::vector<std::size_t> maxel(1,values.size());
  const T* pp;
  getPointer( v, maxel, start, stride, pp );
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(values))
  for(unsigned i=0; i<values.size(); ++i) {
    values[i]=unit*pp[i*stride];
  }
}

template <class T>
void DataPassingObjectTyped<T>::share_data( const std::vector<AtomNumber>&index, const std::vector<unsigned>& i, Value* value ) {
  plumed_dbg_assert( value->getRank()==1 );
  std::vector<std::size_t> maxel(1,index.size());
#ifndef NDEBUG
// bounds are only checked in debug mode since they require this extra step that is potentially expensive
  maxel[0]=(i.size()>0?*std::max_element(i.begin(),i.end())+1:0);
#else
  maxel[0]=0;
#endif
  const T* pp;
  getPointer( v, maxel, start, stride, pp );
  // cannot be parallelized with omp because access to data is not ordered
  unsigned k=0;
  for(const auto & p : index) {
    value->data[p.index()]=unit*pp[i[k]*stride];
    k++;
  }
}

template <class T>
void DataPassingObjectTyped<T>::add_force( Value* value ) {
  if( value->getRank()==0 ) {
    *f.template get<T*>() += funit*static_cast<T>(value->getForce(0));
    return;
  }
  T* pp;
  getPointer( f, value->getShape(), start, stride, pp );
  unsigned nvals=value->getNumberOfValues();
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(pp,nvals))
  for(unsigned i=0; i<nvals; ++i) {
    pp[i*stride] += funit*T(value->getForce(i));
  }
}

template <class T>
void DataPassingObjectTyped<T>::add_force( const std::vector<int>& index, Value* value ) {
  plumed_assert( value->getRank()==1 );
  std::vector<std::size_t> s(1,index.size());
  T* pp;
  getPointer( f, s, start, stride, pp );
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(pp,index.size()))
  for(unsigned i=0; i<index.size(); ++i) {
    pp[i*stride] += funit*T(value->getForce(index[i]));
  }
}

template <class T>
void DataPassingObjectTyped<T>::add_force( const std::vector<AtomNumber>& index, const std::vector<unsigned>& i, Value* value ) {
  plumed_dbg_assert( value->getRank()==1 );
  std::vector<std::size_t> maxel(1,index.size());
#ifndef NDEBUG
// bounds are only checked in debug mode since they require this extra step that is potentially expensive
  maxel[0]=(i.size()>0?*std::max_element(i.begin(),i.end())+1:0);
#else
  maxel[0]=0;
#endif
  T* pp;
  getPointer( f, maxel, start, stride, pp );
  unsigned k=0;
  for(const auto & p : index) {
    pp[stride*i[k]] += funit*T(value->getForce(p.index()));
    k++;
  }
}

template <class T>
void DataPassingObjectTyped<T>::rescale_force( const unsigned& n, const double& factor, Value* value ) {
  plumed_assert( value->getRank()>0 );
  std::vector<std::size_t> s( value->getShape() );
  if( s.size()==1 ) {
    s[0] = n;
  }
  T* pp;
  getPointer( f, s, start, stride, pp );
  #pragma omp parallel for num_threads(OpenMP::getGoodNumThreads(pp,n))
  for(unsigned i=0; i<n; ++i) {
    pp[i*stride] *= T(factor);
  }
}

}
