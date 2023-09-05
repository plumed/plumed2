#ifndef __PLUMED_cuda_helpers_cuh
#define __PLUMED_cuda_helpers_cuh

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

namespace PLMD {
namespace CUDAHELPERS {

  template<class T>
  __device__ constexpr const T& mymin(const T& a, const T& b) {
    return (b < a) ? b : a;
  }

  template <typename T>
  __device__ T* shared_memory_proxy() {
    // do we need an __align__() here?
    extern __shared__ unsigned char memory[];
    return reinterpret_cast<T*>(memory);
  }

  //onlymoveable helper
  template <typename T>
  class memoryHolder{
    /*
    struct info {
      T* pointer_{nullptr};
      //the dimension in memory
      unsigned dim_{0};
      //the used part of the array
      unsigned usedim_{0};
      //mutex
      //number of implementations
    };
    */
  T* pointer_{nullptr};
  //the dimension in memory
  unsigned dim_{0};
  //the used part of the array
  unsigned usedim_{0};
  public:
    memoryHolder(){}
    memoryHolder(memoryHolder& other) = delete;
    memoryHolder(memoryHolder&& other)
    : pointer_(other.pointer_),dim_(other.dim_),
  usedim_(other.usedim_){other.pointer_=nullptr;}
    memoryHolder(const unsigned newDim)
    :dim_(newDim),usedim_(newDim) {
      cudaMalloc(&pointer_,dim_*sizeof(T));
    }
    memoryHolder& operator=(memoryHolder&& other) {
      pointer_ = other.pointer_;
      dim_ = other.dim_;
      usedim_ = other.usedim_;
      other.pointer_=nullptr;
      other.dim_=0;
      other.usedim_=0;
      return *this;
      }
    
    ~memoryHolder(){cudaFree(pointer_);}
    T* pointer(){return pointer_;}
    unsigned size() const {return usedim_;}
    unsigned reserved() const {return dim_;}
    void/*CUresult*/ resize(const unsigned newDim, bool keepData=false){
      
      if (newDim > dim_){
        T *new_ptr = nullptr;
        cudaMalloc(&new_ptr,newDim*sizeof(T));
        if (keepData && pointer_){
          //copy the old data into the new pointer
          cudaMemcpy(new_ptr, pointer_, usedim_*sizeof(T),cudaMemcpyDeviceToDevice);
        }
        cudaFree(pointer_);
        pointer_= new_ptr;
        dim_=newDim;
      }
      usedim_ = newDim;
    }
    void/*CUresult*/ reserve(const unsigned newDim, bool keepData=false){
      if (newDim > dim_){
        T *new_ptr = nullptr;
        cudaMalloc(&new_ptr,newDim*sizeof(T));
        if (keepData && pointer_){
          //copy the old data into the new pointer
          cudaMemcpy(new_ptr, pointer_, usedim_*sizeof(T));
        }
        cudaFree(pointer_);
        pointer_= new_ptr;
        dim_=newDim;
      }
    }

    void/*CUresult*/ copyToCuda(T* cpupointer){
      cudaMemcpy(
        pointer_,
        cpupointer,
        usedim_ * sizeof(T),
        cudaMemcpyHostToDevice);
    }
    
    void/*CUresult*/ copyToCuda(T* cpupointer,cudaStream_t stream){
      cudaMemcpyAsync(
        pointer_,
        cpupointer,
        usedim_ * sizeof(T),
        cudaMemcpyHostToDevice,
        stream);
    }

    void/*CUresult*/ copyFromCuda(T* cpupointer){
      cudaMemcpy(
        cpupointer,
        pointer_,
        usedim_ * sizeof(T),
        cudaMemcpyDeviceToHost);
    }

    void/*CUresult*/ copyFromCuda(T* cpupointer,cudaStream_t stream){
      cudaMemcpyAsync(
        cpupointer,
        pointer_,
        usedim_ * sizeof(T),
        cudaMemcpyDeviceToHost,
        stream);
    }

    void swap(memoryHolder &other) {
        std::swap(this->pointer_,other.pointer_);
        std::swap(this->dim_,other.dim_);
        std::swap(this->usedim_,other.usedim_);
    }
  };

  template <typename T>
  void swap(memoryHolder<T>& a, memoryHolder<T>& b) {
    //a swap with no no copies :)
    memoryHolder<T> t = std::move(a);
    a = std::move(b);
    b = std::move(t);
  }

} //namespace CUDAHELPERS
} //namespace PLMD
#endif //__PLUMED_cuda_helpers_cuh
