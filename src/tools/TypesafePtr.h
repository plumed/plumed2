/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2023 The plumed team
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
#include <cstring>
#include <type_traits>
#include <climits>
#include <initializer_list>

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
  inline void init_shape(const std::size_t* shape) {
    this->shape[0]=0;
    if(shape && *shape>0) {
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

  static std::string extra_msg();

public:

  TypesafePtr(void* ptr, std::size_t nelem, const std::size_t* shape, std::size_t flags):
    ptr(ptr),
    nelem(nelem),
    flags(flags)
  {
    buffer[0]='\0';
    init_shape(shape);
  }

  static const unsigned maxrank=4;
  static TypesafePtr fromSafePtr(void* safe);
  static TypesafePtr setNelemAndShape(const TypesafePtr &other, std::size_t nelem, const std::size_t* shape) {
    return TypesafePtr(other.ptr,nelem,shape,other.flags);
  }
  static TypesafePtr unchecked(const void* ptr) {
    return TypesafePtr(const_cast<void*>(ptr),0,nullptr,0);
  }
  static constexpr unsigned short is_integral=3;
  static constexpr unsigned short is_floating_point=4;
  static constexpr unsigned short is_file=5;

  TypesafePtr() {
    shape[0]=0;
    buffer[0]='\0';
  }

  TypesafePtr(std::nullptr_t)
  {
    shape[0]=0;
    buffer[0]='\0';
  }

/// Macro that generate a constructor with given type and flags
#define __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type,type_,flags_) \
  TypesafePtr(type_*ptr, std::size_t nelem=0, const std::size_t* shape=nullptr) : \
    ptr((void*)const_cast<type*>(ptr)), \
    nelem(nelem), \
    flags(flags_) \
  { \
    init_shape(shape); \
    buffer[0]='\0'; \
  }

/// Macro that uses __PLUMED_WRAPPER_TYPESAFEPTR_INNER to generate constructors with
/// all possible pointer-const combinations
#define __PLUMED_WRAPPER_TYPESAFEPTR(type, code,size) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type, type,             size | (0x10000*(code)) | (0x2000000*2)) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type, type const,       size | (0x10000*(code)) | (0x2000000*3)) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type*,type*,            size | (0x10000*(code)) | (0x2000000*4)) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type*,type*const,       size | (0x10000*(code)) | (0x2000000*5)) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type*,type const*,      size | (0x10000*(code)) | (0x2000000*6)) \
  __PLUMED_WRAPPER_TYPESAFEPTR_INNER(type*,type const*const, size | (0x10000*(code)) | (0x2000000*7))

/// Macro that generates the constructors from empy types (those of which sizeof cannot be computed)
#define __PLUMED_WRAPPER_TYPESAFEPTR_EMPTY(type,code) __PLUMED_WRAPPER_TYPESAFEPTR(type,code,0)

/// Macro that generates the constructors from sized types (those of which sizeof can be computed).
/// In addition to generating constructors with all pointer types, it generates a constructor to
/// allow pass-by-value
#define __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(type,code) \
  __PLUMED_WRAPPER_TYPESAFEPTR(type,code,sizeof(type)) \
  TypesafePtr(type val, std::size_t nelem=0, const std::size_t* shape=nullptr): \
      nelem(1), \
      flags(sizeof(type) | (0x10000*(code)) | (0x2000000*1)) \
    { \
    plumed_assert(sizeof(type)<=32); \
    ptr=&buffer[0]; \
    flags=sizeof(type) | (0x10000*(code)) | (0x2000000*1); \
    std::memcpy(&buffer[0],&val,sizeof(type)); \
    init_shape(shape); \
  }

/// Here we create all the required instances
/// 1: void
/// 3: integral
/// 4: floating
/// 5: FILE
/// 0x100: unsigned
  __PLUMED_WRAPPER_TYPESAFEPTR_EMPTY(void,1)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(char,(CHAR_MIN==0)*0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(signed char,3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(unsigned char,0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(short,3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(unsigned short,0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(int,3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(unsigned int,0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(long,3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(unsigned long,0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(long long,3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(unsigned long long,0x100+3)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(float,4)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(double,4)
  __PLUMED_WRAPPER_TYPESAFEPTR_SIZED(long double,4)
  __PLUMED_WRAPPER_TYPESAFEPTR_EMPTY(FILE,5)

  ~TypesafePtr() {
  }


  TypesafePtr(const TypesafePtr&other) = delete;

  TypesafePtr & operator=(const TypesafePtr & other) = delete;

  TypesafePtr(TypesafePtr&&other):
    buffer(other.buffer),
    ptr(other.ptr==&other.buffer[0] ? &buffer[0] : other.ptr),
    nelem(other.nelem),
    shape(other.shape),
    flags(other.flags)
  {
    other.ptr=nullptr;
  }

  TypesafePtr copy() const;

  TypesafePtr & operator=(TypesafePtr && other) {
    ptr=(other.ptr==&other.buffer[0] ? &buffer[0] : other.ptr);
    flags=other.flags;
    nelem=other.nelem;
    shape=other.shape;
    buffer=other.buffer;
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

private:

  template<typename T>
  T* get_priv(std::size_t nelem, const std::size_t* shape, bool byvalue) const {

    if(typesafePtrSkipCheck()) return (T*) ptr;
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
        throw ExceptionTypeError() <<"This command expects an integer type. Received a " << type_str() << " instead"<<extra_msg();
      }
      if(std::is_floating_point<T_noptr>::value && (type!=is_floating_point)) {
        throw ExceptionTypeError() <<"This command expects a floating point type. Received a " << type_str() << " instead"<<extra_msg();
      }
      if(std::is_same<FILE,typename std::remove_const<T_noptr>::type>::value && (type!=is_file)) {
        throw ExceptionTypeError() <<"This command expects a FILE. Received a " << type_str() << " instead"<<extra_msg();
      }
    }

    if(size>0 && typesafePtrSizeof<T_noptr>() >0 && typesafePtrSizeof<T_noptr>()!=size) {
      throw ExceptionTypeError() << "This command expects a type with size " << typesafePtrSizeof<T_noptr>() << ". Received type has size " << size << " instead"<<extra_msg();
    }

    if(!byvalue) if(cons==1) {
        throw ExceptionTypeError() << "This command is trying to take the address of an argument that was passed by value"<<extra_msg();
      }

    // cons==1 (by value) is here treated as cons==3 (const type*)
    if(cons>0) {
      if(!std::is_pointer<T>::value) {
        if(std::is_void<T>::value) {
          if(cons==1) {
            throw ExceptionTypeError() << "This command expects a void pointer. It received a value instead"<<extra_msg();
          }
        } else {
          if(cons!=1 && cons!=2 && cons!=3) {
            throw ExceptionTypeError() << "This command expects a pointer or a value. It received a pointer-to-pointer instead"<<extra_msg();
          }
        }
        if(!std::is_const<T>::value) {
          if(cons==3) {
            throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a non modifiable pointer instead (const T*)"<<extra_msg();
          } else if(cons==1) {
            throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a value instead (T)"<<extra_msg();
          }
        }
      } else {
        if(!std::is_const<T>::value) {
          if(cons==1) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a value intead"<<extra_msg();
          if(cons==2 || cons==3) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a pointer intead"<<extra_msg();
          if(!std::is_const<T_noptr>::value) {
            if(cons!=4) throw ExceptionTypeError() << "This command expects a modifiable pointer-to-pointer (T**)"<<extra_msg();
          } else {
            if(cons!=6) throw ExceptionTypeError() << "This command expects a modifiable pointer to unmodifiable pointer (const T**)"<<extra_msg();
          }
        } else {
          if(!std::is_const<T_noptr>::value) {
            if(cons!=4 && cons!=5) throw ExceptionTypeError() << "This command expects T*const* pointer, and can only receive T**  or T*const* pointers"<<extra_msg();
          }
        }
      }
    }
    // check full shape, if possible
    if(shape && shape[0] && this->shape[0]) {
      for(unsigned i=0; i<this->shape.size(); i++) {
        if(shape[i]==0 && this->shape[i]!=0) {
          throw ExceptionTypeError() << "Incorrect number of axis (passed greater than requested)"<<extra_msg();
        }
        if(shape[i]!=0 && this->shape[i]==0) {
          throw ExceptionTypeError() << "Incorrect number of axis (requested greater than passed)"<<extra_msg();
        }
        if(shape[i]==0) break;
        if(!(shape[i]<=this->shape[i])) {
          throw ExceptionTypeError() << "This command wants to access " << shape[i] << " on axis " << i <<" of this pointer, but only " << this->shape[i] << " have been passed"<<extra_msg();
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
        throw ExceptionTypeError() << "This command wants to access " << nelem << " from this pointer, but only " << this->nelem << " have been passed"<<extra_msg();
      }
    return (T*) ptr;
  }

public:

  template<typename T>
  void set(T val) const {
    *get_priv<T>(0,nullptr,false)=val;
  }

  template<typename T>
  typename std::enable_if<std::is_pointer<T>::value,T>::type get() const {
    typedef typename std::remove_pointer<T>::type T_noptr;
    return get_priv<T_noptr>(0,nullptr,false);
  }

  template<typename T>
  typename std::enable_if<!std::is_pointer<T>::value,T>::type get() const {
    return *get_priv<const T>(1,nullptr,true);
  }

  template<typename T>
  T get(std::size_t nelem) const {
    static_assert(std::is_pointer<T>::value,"only pointer types allowed here");
    typedef typename std::remove_pointer<T>::type T_noptr;
    return get_priv<T_noptr>(nelem,nullptr,false);
  }

  template<typename T>
  T get(std::initializer_list<std::size_t> shape) const {
    static_assert(std::is_pointer<T>::value,"only pointer types allowed here");
    plumed_assert(shape.size()<=maxrank);
    std::array<std::size_t,maxrank+1> shape_;
    typedef typename std::remove_pointer<T>::type T_noptr;
    unsigned j=0;
    for(auto i : shape) {
      shape_[j]=i;
      j++;
    }
    shape_[j]=0;
    return get_priv<T_noptr>(0,&shape_[0],false);
  }

  operator bool() const noexcept {
    return ptr;
  }

  void* getRaw() const noexcept {
    return ptr;
  }

  std::size_t getNelem() const noexcept {
    return nelem;
  }

  const std::size_t* getShape() const noexcept {
    return shape.data();
  }

  std::size_t getFlags() const noexcept {
    return flags;
  }

private:
  std::array<char,32> buffer;
  void* ptr=nullptr;
  std::size_t nelem=0;
  std::array<std::size_t,maxrank+1> shape; // make sure to initialize this!
  std::size_t flags=0;
};

}


#endif
