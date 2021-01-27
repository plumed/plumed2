/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018 The plumed team
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
#ifndef __PLUMED_tools_TypesafePtr_h
#define __PLUMED_tools_TypesafePtr_h

#include "Exception.h"

#include <memory>
#include <iosfwd>
#include <map>
#include <utility>
#include <mutex>
#include <cstdio>
#include <array>

namespace PLMD {

static inline bool typesafePtrSkipCheck() {
  static const bool ret=std::getenv("PLUMED_TYPESAFE_IGNORE");
  return ret;
}

template<class T>
std::size_t typesafePtrSizeof() {
  return sizeof(T);
}

template<>
inline
std::size_t typesafePtrSizeof<void>() {
  return 0;
}

template<>
inline
std::size_t typesafePtrSizeof<const void>() {
  return 0;
}

/**
\ingroup TOOLBOX
Class to deal with propoagation of typesafe pointers.

*/
class TypesafePtr {
  TypesafePtr(void* ptr, std::size_t nelem, const std::size_t* shape, unsigned long int flags):
    ptr(ptr),
    nelem(nelem),
    flags(flags)
  {
    this->shape[0]=0;
    if(shape) {
      std::size_t nelem_=1;
      unsigned i=0;
      for(i=0; i<this->shape.size(); i++) {
        this->shape[i]=*shape;
        if(*shape==0) break;
        nelem_*=*shape;
        shape++;
      }
      plumed_assert(i<this->shape.size()); // check that last element is actually zero
      if(nelem==0) nelem=nelem_;
      plumed_assert(nelem==nelem_) << "Inconsistent shape/nelem";
    }
  }

public:

  static const unsigned maxrank=4;
  static TypesafePtr fromSafePtr(void* safe);
  static constexpr unsigned short is_integral=3;
  static constexpr unsigned short is_floating_point=4;
  static constexpr unsigned short is_file=5;

  TypesafePtr() {
    shape[0]=0;
  }

  TypesafePtr(const void*ptr) :
    ptr(const_cast<void*>(ptr))
  {}

  ~TypesafePtr() {
  }


  TypesafePtr(const TypesafePtr&other) = delete;

  TypesafePtr & operator=(const TypesafePtr & other) = delete;

  TypesafePtr(TypesafePtr&&other):
    ptr(other.ptr),
    nelem(other.nelem),
    shape(other.shape),
    flags(other.flags)
  {
    other.ptr=nullptr;
  }

  TypesafePtr copy() const;

  TypesafePtr & operator=(TypesafePtr && other) {
    ptr=other.ptr;
    flags=other.flags;
    nelem=other.nelem;
    shape=other.shape;
    other.ptr=nullptr;
    return *this;
  }

  std::string type_str() const {
    auto type=(flags>>16)&0xff;
    if(type==0) return "wildcard";
    if(type==1) return "void";
    if(type==2) return "integral";
    if(type==3) return "integral";
    if(type==4) return "floating point";
    if(type==5) return "FILE";
    return "unknown";
  }

  template<typename T>
  T* get_priv(std::size_t nelem, const std::size_t* shape, bool byvalue) const {

    if(typesafePtrSkipCheck()) return (T*) ptr;
    typedef typename std::remove_const<T>::type T_noconst;
    typedef typename std::remove_pointer<T>::type T_noptr;
    if(flags==0) return (T*) ptr; // no check
    auto size=flags&0xffff;
    auto type=(flags>>16)&0xff;
    // auto unsi=(flags>>24)&0x1; // ignored
    auto cons=(flags>>25)&0x7;

    // type=0: ignore check
    // type>5: undefined yet
    if(type!=0 && type<=5) {
      if(std::is_integral<T_noptr>::value && (type!=is_integral)) {
        throw ExceptionTypeError() <<"This command expects an integer type. Received a " << type_str() << " instead";
      }
      if(std::is_floating_point<T_noptr>::value && (type!=is_floating_point)) {
        throw ExceptionTypeError() <<"This command expects a floating point type. Received a " << type_str() << " instead";
      }
      if(std::is_same<std::remove_const<FILE>,T_noconst>::value && (type!=is_file)) {
        throw ExceptionTypeError() <<"This command expects a FILE. Received a " << type_str() << " instead";
      }
    }

    if(size>0 && typesafePtrSizeof<T_noptr>() >0 && typesafePtrSizeof<T_noptr>()!=size) {
      throw ExceptionTypeError() << "This command expects a type with size " << typesafePtrSizeof<T_noptr>() << ". Received type has size " << size << " instead";
    }

    if(!byvalue) if(cons==1) {
        throw ExceptionTypeError() << "This command is trying to take the address of an argument that was passed by value";
      }

    // cons==1 (by value) is here treated as cons==3 (const type*)
    if(!std::is_pointer<T>::value) {
      if(std::is_void<T>::value) {
        if(cons==1) {
          throw ExceptionTypeError() << "This command expects a void pointer. It received a value instead";
        }
      } else {
        if(cons!=1 && cons!=2 && cons!=3) {
          throw ExceptionTypeError() << "This command expects a pointer or a value. It received a pointer-to-pointer instead";
        }
      }
      if(!std::is_const<T>::value) {
        if(cons==3) {
          throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a non modifiable pointer instead (const T*)";
        } else if(cons==1) {
          throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a value instead (T)";
        }
      }
    } else {
      if(!std::is_const<T>::value) {
        if(cons==1) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a value intead";
        if(cons==2 || cons==3) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a pointer intead";
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4) throw ExceptionTypeError() << "This command expects a modifiable pointer-to-pointer (T**)";
        } else {
          if(cons!=6) throw ExceptionTypeError() << "This command expects a modifiable pointer to unmodifiable pointer (const T**)";
        }
      } else {
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4 && cons!=5) throw ExceptionTypeError() << "This command expects T*const* pointer, and can only receive T**  or T*const* pointers";
        }
      }
    }
    // check full shape, if possible
    if(shape && shape[0] && this->shape[0]) {
      for(unsigned i=0; i<this->shape.size(); i++) {
        if(shape[i]==0 && this->shape[i]!=0) {
          throw ExceptionTypeError() << "Incorrect number of axis (passed greater than requested)";
        }
        if(shape[i]!=0 && this->shape[i]==0) {
          throw ExceptionTypeError() << "Incorrect number of axis (requested greater than passed)";
        }
        if(shape[i]==0) break;
        if(!(shape[i]<=this->shape[i])) {
          throw ExceptionTypeError() << "This command wants to access " << shape[i] << " on axis " << i <<" of this pointer, but only " << this->shape[i] << " have been passed";
        }
      }
    }
    if(nelem==0 && shape && shape[0]>0) {
      nelem=1;
      for(unsigned i=0; i<this->shape.size(); i++) {
        if(shape[i]==0) break;
        nelem*=shape[i];
      }
    }
    // check number of elements
    if(nelem>0 && this->nelem>0) if(!(nelem<=this->nelem)) {
        throw ExceptionTypeError() << "This command wants to access " << nelem << " from this pointer, but only " << this->nelem << " have been passed";
      }
    return (T*) ptr;
  }

  template<typename T>
  T* get(std::size_t nelem=0) const {
    return get_priv<T>(nelem,nullptr,false);
  }

  template<typename T>
  T* get(std::initializer_list<std::size_t> shape) const {
    plumed_assert(shape.size()<=maxrank);
    std::array<std::size_t,maxrank+1> shape_;
    unsigned j=0;
    for(auto i : shape) {
      shape_[j]=i;
      j++;
    }
    shape_[j]=0;
    return get_priv<T>(0,&shape_[0],false);
  }

  template<typename T>
  T getVal() const {
    return *get_priv<const T>(1,nullptr,true);
  }

  operator bool() const noexcept {
    return ptr;
  }

private:
  void* ptr=nullptr;
  std::size_t nelem=0;
  std::array<std::size_t,maxrank+1> shape;
  unsigned long int flags=0;
};

}


#endif
