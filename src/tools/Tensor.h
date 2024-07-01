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
#ifndef __PLUMED_tools_Tensor_h
#define __PLUMED_tools_Tensor_h

#include "MatrixSquareBracketsAccess.h"
#include "Vector.h"
#include "LoopUnroller.h"
#include "Exception.h"

#include <array>

namespace PLMD {

//should I add something for calling ssyevr?
/// Small class to contain local utilities.
/// Should not be used outside of the TensorGeneric class.
struct TensorGenericAux {
// local redefinition, just to avoid including lapack.h here
  static void local_dsyevr(const char *jobz, const char *range, const char *uplo, int *n,
                           double *a, int *lda, double *vl, double *vu, int *
                           il, int *iu, double *abstol, int *m, double *w,
                           double *z__, int *ldz, int *isuppz, double *work,
                           int *lwork, int *iwork, int *liwork, int *info);
};

/**
\ingroup TOOLBOX
Class implementing fixed size matrices of Ts

\tparam n The number rows
\tparam m The number columns

This class implements a matrix of Ts with size fixed at
compile time. It is useful for small fixed size objects (e.g.
3x3 tensors) as it does not waste space to store the vector size.
Moreover, as the compiler knows the size, it can be completely
opimized inline.
Most of the loops are explicitly unrolled using PLMD::LoopUnroller class
Matrix elements are initialized to zero by default. Notice that
this means that constructor is a bit slow. This point might change
in future if we find performance issues.
It takes advantage of MatrixSquareBracketsAccess to provide both
() and [] syntax for access.
Several functions are declared as friends even if not necessary so as to
properly appear in Doxygen documentation.

Aliases are defined to simplify common declarations (Tensor, Tensor2d, Tensor3d, Tensor4d).
Also notice that some operations are only available for 3x3 tensors.

Example of usage
\verbatim

#include "Tensor.h"

using namespace PLMD;

int main(){
  Tensor a;
  TensorGeneric<3,2> b;
  TensorGeneric<3,2> c=matmul(a,b);
  return 0;
}

\endverbatim
*/

//there are some "hiccups" with the friend operators, so I declared them inline in the body of tensor
//in that way they are automatically declared templaterd without declaration+implementation
//that despite was the suggested way of doing the implementation on my compiler seems to not work as intended

template<typename T, unsigned n, unsigned m> class TensorGeneric;
//no need to friend-declare these functions
/// \ingroup TOOLBOX
/// returns the determinant of a tensor
template<typename T> T determinant(const TensorGeneric<T,3,3>&);
/// \ingroup TOOLBOX
/// returns the inverse of a 3x3 tensor (same as TensorGeneric<T,3,3>::inverse())
template<typename T> TensorGeneric<T,3,3> inverse(const TensorGeneric<T,3,3>&);
/// \ingroup TOOLBOX
/// returns the determinant of a 3x3 tensor (same as TensorGeneric<T,3,3>::determinant())
template<typename T> T determinant(const TensorGeneric<T,3,3>&);
template<typename T> TensorGeneric<T,3,3> dcrossDv1(const VectorGeneric<T,3>&,const VectorGeneric<T,3>&);
template<typename T> TensorGeneric<T,3,3> dcrossDv2(const VectorGeneric<T,3>&,const VectorGeneric<T,3>&);
template<typename T> TensorGeneric<T,3,3> VcrossTensor(const VectorGeneric<T,3>&,const TensorGeneric<T,3,3>&);
template<typename T> TensorGeneric<T,3,3> VcrossTensor(const TensorGeneric<T,3,3>&,const VectorGeneric<T,3>&);
/// \ingroup TOOLBOX
/// Derivative of a normalized vector
template<typename T> TensorGeneric<T,3,3> deriNorm(const VectorGeneric<T,3>&,const TensorGeneric<T,3,3>&);

template<typename T, unsigned n, unsigned m>
std::ostream & operator<< (std::ostream &os, const TensorGeneric<T,n,m>&);

template<typename T, unsigned n, unsigned m>
class TensorGeneric:
  public MatrixSquareBracketsAccess<TensorGeneric<T,n,m>,T>
{
  std::array<T,n*m> d;
/// Auxiliary private function for constructor
  void auxiliaryConstructor();
/// Auxiliary private function for constructor
  template<typename... Args>
  void auxiliaryConstructor(T first,Args... arg);
public:
/// Constructor accepting n*m T parameters.
/// Can be used as Tensor<2,2>(1.0,2.0,3.0,4.0)
/// In case a wrong number of parameters is given, a static assertion will fail.
  template<typename... Args>
  TensorGeneric(T first,Args... arg);
/// initialize the tensor to zero
  TensorGeneric();
/// initialize a tensor as an external product of two Vector
  TensorGeneric(const VectorGeneric<T,n>&v1,const VectorGeneric<T,m>&v2);
/// set it to zero
  void zero();
/// get the underline pointer to data
  T* data();
/// get the underline pointer to data
  const T* data() const;
/// access element
  T & operator() (unsigned i,unsigned j);
/// access element
  const T & operator() (unsigned i,unsigned j)const;
/// increment
  TensorGeneric& operator +=(const TensorGeneric& b);
/// decrement
  TensorGeneric& operator -=(const TensorGeneric& b);
/// multiply
  TensorGeneric& operator *=(T);
/// divide
  TensorGeneric& operator /=(T);
/// return +t
  TensorGeneric operator +()const;
/// return -t
  TensorGeneric operator -()const;
/// set j-th column
  TensorGeneric& setCol(unsigned j,const VectorGeneric<T,n> & c);
/// set i-th row
  TensorGeneric& setRow(unsigned i,const VectorGeneric<T,m> & r);
/// get j-th column
  VectorGeneric<T,n> getCol(unsigned j)const;
/// get i-th row
  VectorGeneric<T,m> getRow(unsigned i)const;
/// return t1+t2
  friend TensorGeneric operator+ (TensorGeneric t1,const TensorGeneric& t2) {return t1+=t2;}
/// return t1-t2
  friend TensorGeneric operator- (TensorGeneric t1,const TensorGeneric& t2) {return t1-=t2;}
/// scale the tensor by a factor s
  template<typename U>
  friend TensorGeneric operator*(TensorGeneric t1,U s) {return t1*=s;}
/// scale the tensor by a factor s
  template<typename U>
  friend TensorGeneric operator*(U s,TensorGeneric t1) {return t1*=s;}
/// scale the tensor by a factor 1/s
  template<typename U>
  friend TensorGeneric operator/(TensorGeneric t1,U s) {return t1*=(T(1.0)/s);}
/// returns the determinant
  T determinant()const;
/// return an identity tensor
  static TensorGeneric<T,n,n> identity();
/// return the matrix inverse
  TensorGeneric inverse()const;
/// return the transpose matrix
  TensorGeneric<T,m,n> transpose()const;
/// << operator.
/// Allows printing tensor `t` with `std::cout<<t;`
  friend std::ostream & operator<< <>(std::ostream &os, const TensorGeneric&);
/// Diagonalize tensor.
/// Syntax is the same as Matrix::diagMat.
/// In addition, it is possible to call if with m_ smaller than n_. In this case,
/// only the first (smaller) m_ eigenvalues and eigenvectors are retrieved.
/// If case lapack fails (info!=0) it throws an exception.
/// Notice that tensor is assumed to be symmetric!!!
  template<typename U, unsigned n_, unsigned m_>
  friend void diagMatSym(const TensorGeneric<U,n_,n_>&,VectorGeneric<U,m_>&evals,TensorGeneric<U,m_,n_>&evec);
};

template<typename T, unsigned n, unsigned m>
void TensorGeneric<T,n,m>::auxiliaryConstructor()
{}

template<typename T, unsigned n, unsigned m>
template<typename... Args>
void TensorGeneric<T,n,m>::auxiliaryConstructor(T first,Args... arg)
{
  d[n*m-(sizeof...(Args))-1]=first;
  auxiliaryConstructor(arg...);
}

template<typename T, unsigned n, unsigned m>
template<typename... Args>
TensorGeneric<T,n,m>::TensorGeneric(T first,Args... arg)
{
  static_assert((sizeof...(Args))+1==n*m,"you are trying to initialize a Tensor with the wrong number of arguments");
  auxiliaryConstructor(first,arg...);
}

template<typename T, unsigned n, unsigned m>
T* TensorGeneric<T,n,m>::data() {return d.data();}
template<typename T, unsigned n, unsigned m>
const T* TensorGeneric<T,n,m>::data() const {return d.data();}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>::TensorGeneric() {
  LoopUnroller<T,n*m>::_zero(d.data());
}

/* between RVO and compile time this should be faster, but slows down openACC, a lot
template <typename T, unsigned i, unsigned j, unsigned m>
void external_rec(T*const out,const T*const v1, const T*const v2){
      if constexpr (j>0) {
        external_rec<i,j-1,m>(out,v1,v2);
      } else if constexpr (i>0) {
          external_rec<i-1,m-1,m>(out,v1,v2);
      }
      out[i*m+j]=v1[i]*v2[j];
}

template<typename T, unsigned n, unsigned m>
std::array<T,n*m> externaProd(const VectorGeneric<T,n>&v1,const VectorGeneric<T,m>&v2){
std::array<T,n*m> toRet;
external_rec<n-1,m-1,m>(toRet.data(),v1.data(),v2.data());
return toRet;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>::TensorGeneric(const VectorGeneric<T,n>&v1,const VectorGeneric<T,m>&v2)
:d(externaProd(v1,v2)) {}
*/

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>::TensorGeneric(const VectorGeneric<T,n>&v1,const VectorGeneric<T,m>&v2) {
  for(unsigned i=0; i<n; i++) {
    for(unsigned j=0; j<m; j++) {
      d[i*m+j]=v1[i]*v2[j];
    }
  }
}

template<typename T, unsigned n, unsigned m>
void TensorGeneric<T,n,m>::zero() {
  LoopUnroller<T,n*m>::_zero(d.data());
}

template<typename T, unsigned n, unsigned m>
T & TensorGeneric<T,n,m>::operator() (unsigned i,unsigned j) {
#ifdef _GLIBCXX_DEBUG
// index i is implicitly checked by the std::array class
  plumed_assert(j<m);
#endif
  return d[m*i+j];
}

template<typename T, unsigned n, unsigned m>
const T & TensorGeneric<T,n,m>::operator() (unsigned i,unsigned j)const {
#ifdef _GLIBCXX_DEBUG
// index i is implicitly checked by the std::array class
  plumed_assert(j<m);
#endif
  return d[m*i+j];
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::operator +=(const TensorGeneric<T,n,m>& b) {
  LoopUnroller<T,n*m>::_add(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::operator -=(const TensorGeneric<T,n,m>& b) {
  LoopUnroller<T,n*m>::_sub(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::operator *=(T s) {
  LoopUnroller<T,n*m>::_mul(d.data(),s);
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::operator /=(T s) {
  LoopUnroller<T,n*m>::_mul(d.data(),1.0/s);
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m> TensorGeneric<T,n,m>::operator+()const {
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m> TensorGeneric<T,n,m>::operator-()const {
  TensorGeneric<T,n,m> r;
  LoopUnroller<T,n*m>::_neg(r.d.data(),d.data());
  return r;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::setCol(unsigned j,const VectorGeneric<T,n> & c) {
  for(unsigned i=0; i<n; ++i) (*this)(i,j)=c(i);
  return *this;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m>& TensorGeneric<T,n,m>::setRow(unsigned i,const VectorGeneric<T,m> & r) {
  for(unsigned j=0; j<m; ++j) (*this)(i,j)=r(j);
  return *this;
}

template<typename T, unsigned n, unsigned m>
VectorGeneric<T,n> TensorGeneric<T,n,m>::getCol(unsigned j)const {
  VectorGeneric<T,n> v;
  for(unsigned i=0; i<n; ++i) v(i)=(*this)(i,j);
  return v;
}

template<typename T, unsigned n, unsigned m>
VectorGeneric<T,m> TensorGeneric<T,n,m>::getRow(unsigned i)const {
  VectorGeneric<T,m> v;
  for(unsigned j=0; j<m; ++j) v(j)=(*this)(i,j);
  return v;
}

template<typename T, unsigned n, unsigned m>
inline
T TensorGeneric<T,n,m>::determinant()const {
  static_assert(n==3&&m==3,"determinant can be called only for 3x3 Tensors");
  return
    d[0]*d[4]*d[8]
    + d[1]*d[5]*d[6]
    + d[2]*d[3]*d[7]
    - d[0]*d[5]*d[7]
    - d[1]*d[3]*d[8]
    - d[2]*d[4]*d[6];
}

//consider to make this a constexpr function
template<typename T, unsigned n, unsigned m>
inline
TensorGeneric<T,n,n> TensorGeneric<T,n,m>::identity() {
  static_assert(n==m, "Identity matrix can only be called from a square tensor");
  TensorGeneric<T,n,n> t;
  for(unsigned i=0; i<n; i++) t(i,i)=1.0;
  return t;
}

template<typename T, unsigned n, unsigned m>
TensorGeneric<T,m,n> TensorGeneric<T,n,m>::transpose()const {
  TensorGeneric<T,m,n> t;
  for(unsigned i=0; i<m; i++)for(unsigned j=0; j<n; j++) t(i,j)=(*this)(j,i);
  return t;
}

template<typename T, unsigned n, unsigned m>
inline
TensorGeneric<T,n,m> TensorGeneric<T,n,m>::inverse()const {
  static_assert(n==3&&m==3,"inverse can be called only for 3x3 Tensors");
  TensorGeneric t;
  T invdet=1.0/determinant();
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++)
      t(j,i)=invdet*(   (*this)((i+1)%3,(j+1)%3)*(*this)((i+2)%3,(j+2)%3)
                        -(*this)((i+1)%3,(j+2)%3)*(*this)((i+2)%3,(j+1)%3));
  return t;
}

/// matrix-matrix multiplication
template<typename T, unsigned n, unsigned m, unsigned l>
TensorGeneric<T,n,l> matmul(const TensorGeneric<T,n,m>&a,const TensorGeneric<T,m,l>&b) {
  TensorGeneric<T,n,l> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<l; j++) for(unsigned k=0; k<m; k++) {
        t(i,j)+=a(i,k)*b(k,j);
      }
  return t;
}

/// matrix-vector multiplication
template<typename T, unsigned n, unsigned m>
VectorGeneric<T,n> matmul(const TensorGeneric<T,n,m>&a,const VectorGeneric<T,m>&b) {
  VectorGeneric<T,n> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<m; j++) t(i)+=a(i,j)*b(j);
  return t;
}
/// vector-matrix multiplication
template<typename T, unsigned n, unsigned m>
VectorGeneric<T,n> matmul(const VectorGeneric<T,m>&a,const TensorGeneric<T,m,n>&b) {
  VectorGeneric<T,n> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<m; j++) t(i)+=a(j)*b(j,i);
  return t;
}

/// vector-vector multiplication (maps to dotProduct)
template<typename T, unsigned n_>
T matmul(const VectorGeneric<T,n_>&a,const VectorGeneric<T,n_>&b) {
  return dotProduct(a,b);
}

/// matrix-matrix-matrix multiplication
template<typename T, unsigned n, unsigned m, unsigned l, unsigned i>
TensorGeneric<T,n,i> matmul(const TensorGeneric<T,n,m>&a,const TensorGeneric<T,m,l>&b,const TensorGeneric<T,l,i>&c) {
  return matmul(matmul(a,b),c);
}

/// matrix-matrix-vector multiplication
template<typename T, unsigned n, unsigned m, unsigned l>
VectorGeneric<T,n> matmul(const TensorGeneric<T,n,m>&a,const TensorGeneric<T,m,l>&b,const VectorGeneric<T,l>&c) {
  return matmul(matmul(a,b),c);
}

/// vector-matrix-matrix multiplication
template<typename T, unsigned n, unsigned m, unsigned l>
VectorGeneric<T,l> matmul(const VectorGeneric<T,n>&a,const TensorGeneric<T,n,m>&b,const TensorGeneric<T,m,l>&c) {
  return matmul(matmul(a,b),c);
}

/// vector-matrix-vector multiplication
template<typename T, unsigned n, unsigned m>
T matmul(const VectorGeneric<T,n>&a,const TensorGeneric<T,n,m>&b,const VectorGeneric<T,m>&c) {
  return matmul(matmul(a,b),c);
}

template <typename T>
inline
T determinant(const TensorGeneric<T,3,3>&t) {
  return t.determinant();
}

template <typename T>
inline
TensorGeneric<T,3,3> inverse(const TensorGeneric<T,3,3>&t) {
  return t.inverse();
}
/// returns the transpose of a tensor (same as TensorGeneric::transpose())
template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m> transpose(const TensorGeneric<T,m,n>&t) {
  return t.transpose();
}

/// returns the transpose of a tensor (same as TensorGeneric(const VectorGeneric&,const VectorGeneric&))
template<typename T, unsigned n, unsigned m>
TensorGeneric<T,n,m> extProduct(const VectorGeneric<T,n>&v1,const VectorGeneric<T,m>&v2) {
  return TensorGeneric<T,n,m>(v1,v2);
}

template <typename T>
inline
TensorGeneric<T,3,3> dcrossDv1(const VectorGeneric<T,3>&v1,const VectorGeneric<T,3>&v2) {
  (void) v1; // this is to avoid warnings. still the syntax of this function is a bit dummy...
  return TensorGeneric<T,3,3>(
           0.0, v2[2],-v2[1],
           -v2[2],   0.0, v2[0],
           v2[1],-v2[0],   0.0);
}

template <typename T>
inline
TensorGeneric<T,3,3> dcrossDv2(const VectorGeneric<T,3>&v1,const VectorGeneric<T,3>&v2) {
  (void) v2; // this is to avoid warnings. still the syntax of this function is a bit dummy...
  return TensorGeneric<T,3,3>(
           0.0,-v1[2],v1[1],
           v1[2],0.0,-v1[0],
           -v1[1],v1[0],0.0);
}

template<typename T, unsigned n, unsigned m>
std::ostream & operator<<(std::ostream &os, const TensorGeneric<T,n,m>& t) {
  for(unsigned i=0; i<n; i++)for(unsigned j=0; j<m; j++) {
      if(i!=(n-1) || j!=(m-1)) os<<t(i,j)<<" ";
      else os<<t(i,j);
    }
  return os;
}

/// \ingroup TOOLBOX
typedef TensorGeneric<double,1,1> Tensor1d;
/// \ingroup TOOLBOX
typedef TensorGeneric<double,2,2> Tensor2d;
/// \ingroup TOOLBOX
typedef TensorGeneric<double,3,3> Tensor3d;
/// \ingroup TOOLBOX
typedef TensorGeneric<double,4,4> Tensor4d;
/// \ingroup TOOLBOX
typedef TensorGeneric<double,5,5> Tensor5d;
/// \ingroup TOOLBOX
typedef Tensor3d Tensor;

template <typename T>
inline
TensorGeneric<T,3,3> VcrossTensor(const VectorGeneric<T,3>&v1,const TensorGeneric<T,3,3>&v2) {
  TensorGeneric<T,3,3> t;
  for(unsigned i=0; i<3; i++) {
    t.setRow(i,matmul(dcrossDv2(v1,v1),v2.getRow(i)));
  }
  return t;
}

template <typename T>
inline
TensorGeneric<T,3,3> VcrossTensor(const TensorGeneric<T,3,3>&v2,const VectorGeneric<T,3>&v1) {
  TensorGeneric<T,3,3> t;
  for(unsigned i=0; i<3; i++) {
    t.setRow(i,-matmul(dcrossDv2(v1,v1),v2.getRow(i)));
  }
  return t;
}

template <typename T>
inline
TensorGeneric<T,3,3> deriNorm(const VectorGeneric<T,3>&v1,const TensorGeneric<T,3,3>&v2) {
  // delta(v) = delta(v1/v1.norm) = 1/v1.norm*(delta(v1) - (v.delta(v1))cross v;
  T over_norm = 1./v1.modulo();
  return over_norm*(v2 - over_norm*over_norm*(extProduct(matmul(v2,v1),v1)));
}

template<typename T, unsigned n, unsigned m>
void diagMatSym(const TensorGeneric<T,n,n>&mat,VectorGeneric<T,m>&evals,TensorGeneric<T,m,n>&evec) {
  // some guess number to make sure work is large enough.
  // for correctness it should be >=20. However, it is recommended to be the block size.
  // I put some likely exaggerated number
  constexpr int bs=100;
  // temporary data, on stack so as to avoid allocations
  std::array<int,10*n> iwork;
  std::array<T,(6+bs)*n> work;
  std::array<int,2*m> isup;
  // documentation says that "matrix is destroyed" !!!
  auto mat_copy=mat;
  // documentation says this is size n (even if m<n)
  std::array<T,n> evals_tmp;
  int nn=n;              // dimension of matrix
  T vl=0.0, vu=1.0; // ranges - not used
  int one=1,mm=m;        // minimum and maximum index
  T abstol=0.0;     // tolerance
  int mout=0;            // number of eigenvalues found (same as mm)
  int info=0;            // result
  int liwork=iwork.size();
  int lwork=work.size();
  TensorGenericAux::local_dsyevr("V", (n==m?"A":"I"), "U", &nn, mat_copy.data(), &nn, &vl, &vu, &one, &mm,
                                 &abstol, &mout, evals_tmp.data(), evec.data(), &nn,
                                 isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if(info!=0) plumed_error()<<"Error diagonalizing matrix\n"
                              <<"Matrix:\n"<<mat<<"\n"
                              <<"Info: "<<info<<"\n";
  plumed_assert(mout==m);
  for(unsigned i=0; i<m; i++) evals[i]=evals_tmp[i];
  // This changes eigenvectors so that the first non-null element
  // of each of them is positive
  // We can do it because the phase is arbitrary, and helps making
  // the result reproducible
  for(unsigned i=0; i<m; ++i) {
    unsigned j=0;
    for(j=0; j<n; j++) if(evec(i,j)*evec(i,j)>1e-14) break;
    if(j<n) if(evec(i,j)<0.0) for(j=0; j<n; j++) evec(i,j)*=-1;
  }
}

static_assert(sizeof(Tensor)==9*sizeof(double), "code cannot work if this is not satisfied");

} //PLMD

#endif