/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_tools_OpenACC_h
#define __PLUMED_tools_OpenACC_h

#include <vector>
#include "View.h"

#ifdef __PLUMED_HAS_OPENACC
#include <openacc.h>
#include "Exception.h"
#endif //__PLUMED_HAS_OPENACC

namespace PLMD {

namespace OpenACC {

/** @brief this  little tool is a RAII helper to put and remove data on the device

 the captured data will need to have two functions: toACCDevice() and removeFromACCDevice():
 in toACCDevice  you should declare a `#pragma acc enter data` statement with the object
  within your structure to put on the device and eventual calls to toACCDevice of contained objects.

 In the removeFromACCDevice function you should declare a `#pragma acc exit data` statement
 to remove the object from the device but the onbject names should be delcared in the opposite order.fromToDataHelper

 Remember to stat/finish with `this[0:1]`.
 For example:

 @code{c++}
 struct dataContainer {
  int * ptr;
  size_t size;
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1],ptr[0:size], size)
  }
  void removeFromACCDevice() const  {
#pragma acc exit data delete(size,ptr[0:size],this[0:1])
  }
};
 @endcode

 Or, for a slightly more complex example:
 @code{c++}
 struct ParallelActionsInput {
  bool noderiv{false};
  const Pbc* pbc;
  unsigned ncomponents{0};
  unsigned dataSize{0};
  double *inputdata{nullptr};
  ParallelActionsInput( const Pbc& box )
    : pbc(&box) {}
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], noderiv, pbc[0:1],ncomponents, dataSize, inputdata[0:dataSize])
    pbc->toACCDevice();
  }
  void removeFromACCDevice() const  {
    pbc->removeFromACCDevice();
    // assuming dataSize is not changed
#pragma acc exit data delete(inputdata[0:dataSize],dataSize,ncomponents, pbc[0:1],noderiv,this[0:1])
  }
};
 @endcode
*/

template<typename accData>
class fromToDataHelper {
  accData &m;
public:
  fromToDataHelper(accData &d) : m(d) {
    m.toACCDevice();
  }
  ~fromToDataHelper() {
    m.removeFromACCDevice();
  }
};
#ifdef __PLUMED_HAS_OPENACC
///C++ wrapper around acc_malloc, allocates sizeof(T) bytes and return a typed pointer
template <typename T> T *myAccMalloc(size_t size) {
  return reinterpret_cast<T *>(acc_malloc(size * sizeof(T)));
}

///C++ wrapper around acc_free
template <typename T> void myAccFree(T *ptr) {
  acc_free(ptr);
}

///Memory manager for openacc
template <typename T> class memoryManager {
  /// the number of T element stored by the pointer
  size_t size_ {0};
  /// the size of the array stored on the GPU
  size_t stored_ {0};
  /// the device address, not the host address!!!
  void *ptr_ {nullptr};
public:
  /// @brief Initilialize an empty placeholder
  memoryManager()=default;
  /// @brief allocate the memory on the device, and store the pointer
  /// @param sz the number of T that you want to allocate
  memoryManager(size_t sz)
    : size_{sz},
      stored_ {sz} {
    if (size_ > 0) {
      ptr_ = acc_malloc(size_ * sizeof(T));
    }
  }
  /// @brief allocates and copies to the memory of the device the given vector
  memoryManager(const std::vector<T>& data)
    : size_{data.size()},
      stored_{data.size()},
      ptr_{acc_malloc(size_ * sizeof(T))} {
    copyToDevice(data.data());
  }
  /// @brief allocates and copies to the memory of the device the given vector
  template <size_t N>
  memoryManager(const View<T,N>& data)
    : size_{data.size()},
      stored_{data.size()},
      ptr_{acc_malloc(size_ * sizeof(T))} {
    copyToDevice(data.data());
  }
  /// owns device memory: copying would free it twice
  memoryManager(const memoryManager&)=delete;
  memoryManager& operator=(const memoryManager&)=delete;
  memoryManager(memoryManager&& other) noexcept
    : size_{other.size_},
  stored_{other.stored_},
  ptr_{other.ptr_} {
    other.size_=0;
    other.stored_=0;
    other.ptr_=nullptr;
  }
  memoryManager& operator=(memoryManager&& other) noexcept {
    if (this!=&other) {
      acc_free(ptr_);
      size_=other.size_;
      stored_=other.stored_;
      ptr_=other.ptr_;
      other.size_=0;
      other.stored_=0;
      other.ptr_=nullptr;
    }
    return *this;
  }
  /// frees the stored memory
  ~memoryManager() {
    acc_free(ptr_);
  }
  /// @brief gets the device addres
  /// @return gets the device addres, typed
  constexpr T *devicePtr() const {
    return reinterpret_cast<T *>(ptr_);
  }
  /// Resizes the data stored on the device. If the new size is not larger than
  /// the allocated one the memory is not reallocated. Otherwise the old content
  /// is lost and the allocation grows by at least 50%, so that a size that
  /// increases a little at every step does not cause a reallocation at every
  /// step (acc_free usually synchronizes the whole device)
  void resize(size_t sz) {
    size_=sz;
    if (size_ > stored_) {
      const size_t grown=stored_+stored_/2;
      stored_=(size_>grown)?size_:grown;
      acc_free(ptr_);
      ptr_=acc_malloc(stored_ * sizeof(T));
      plumed_massert(ptr_, "acc_malloc failed: not enough memory on the device?");
    }
  }
  /// return the current number of elements stored on the device
  constexpr size_t size() const {
    return size_;
  }
  /// @brief Copies the data from the address to the device
  void copyToDevice(T *data) {
    acc_memcpy_to_device(ptr_, data, sizeof(T) * size_);
  }
  /// @brief Copies the data from the address to the device
  void copyToDevice(T const *data) {
    acc_memcpy_to_device(ptr_, const_cast<T *>(data), sizeof(T) * size_);
  }
  /// @brief copies the data from the device to the host
  void copyFromDevice(T *data) {
    acc_memcpy_from_device(data, ptr_, sizeof(T) * size_);
  }
  /// @brief Queues the copy of the data from the address to the device on the
  /// given async queue: the host data must stay valid until the queue is waited
  void copyToDeviceAsync(T const *data, int queue) {
    if(size_>0) {
      acc_memcpy_to_device_async(ptr_, const_cast<T *>(data), sizeof(T) * size_, queue);
    }
  }
  /// @brief Queues the copy of the data from the device to the host on the
  /// given async queue: the host data is valid only after the queue is waited
  void copyFromDeviceAsync(T *data, int queue) {
    if(size_>0) {
      acc_memcpy_from_device_async(data, ptr_, sizeof(T) * size_, queue);
    }
  }
};
#endif //__PLUMED_HAS_OPENACC
}//namespace OpenACC
}//namespace PLMD

#endif
