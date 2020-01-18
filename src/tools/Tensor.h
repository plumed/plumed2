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
#ifndef __PLUMED_tools_Tensor_h
#define __PLUMED_tools_Tensor_h

#include "MatrixSquareBracketsAccess.h"
#include "Vector.h"
#include "LoopUnroller.h"
#include "Exception.h"

#include <array>

namespace PLMD {

/// Small class to contain local utilities.
/// Should not be used outside of the TensorGeneric class.
class TensorGenericAux {
public:
// local redefinition, just to avoid including lapack.h here
  static void local_dsyevr(const char *jobz, const char *range, const char *uplo, int *n,
                           double *a, int *lda, double *vl, double *vu, int *
                           il, int *iu, double *abstol, int *m, double *w,
                           double *z__, int *ldz, int *isuppz, double *work,
                           int *lwork, int *iwork, int *liwork, int *info);
};

/**
\ingroup TOOLBOX
Class implementing fixed size matrices of doubles

\tparam n The number rows
\tparam m The number columns

This class implements a matrix of doubles with size fixed at
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
template <unsigned n,unsigned m>
class TensorGeneric:
  public MatrixSquareBracketsAccess<TensorGeneric<n,m>,double>
{
  std::array<double,n*m> d;
/// Auxiliary private function for constructor
  void auxiliaryConstructor();
/// Auxiliary private function for constructor
  template<typename... Args>
  void auxiliaryConstructor(double first,Args... arg);
public:
/// Constructor accepting n*m double parameters.
/// Can be used as Tensor<2,2>(1.0,2.0,3.0,4.0)
/// In case a wrong number of parameters is given, a static assertion will fail.
  template<typename... Args>
  TensorGeneric(double first,Args... arg);
/// initialize the tensor to zero
  TensorGeneric();
/// initialize a tensor as an external product of two Vector
  TensorGeneric(const VectorGeneric<n>&v1,const VectorGeneric<m>&v2);
/// set it to zero
  void zero();
/// access element
  double & operator() (unsigned i,unsigned j);
/// access element
  const double & operator() (unsigned i,unsigned j)const;
/// increment
  TensorGeneric& operator +=(const TensorGeneric<n,m>& b);
/// decrement
  TensorGeneric& operator -=(const TensorGeneric<n,m>& b);
/// multiply
  TensorGeneric& operator *=(double);
/// divide
  TensorGeneric& operator /=(double);
/// return +t
  TensorGeneric operator +()const;
/// return -t
  TensorGeneric operator -()const;
/// set j-th column
  TensorGeneric& setCol(unsigned j,const VectorGeneric<n> & c);
/// set i-th row
  TensorGeneric& setRow(unsigned i,const VectorGeneric<m> & r);
/// get j-th column
  VectorGeneric<n> getCol(unsigned j)const;
/// get i-th row
  VectorGeneric<m> getRow(unsigned i)const;
/// return t1+t2
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> operator+(const TensorGeneric<n_,m_>&,const TensorGeneric<n_,m_>&);
/// return t1+t2
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> operator-(const TensorGeneric<n_,m_>&,const TensorGeneric<n_,m_>&);
/// scale the tensor by a factor s
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> operator*(double,const TensorGeneric<n_,m_>&);
/// scale the tensor by a factor s
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> operator*(const TensorGeneric<n_,m_>&,double s);
/// scale the tensor by a factor 1/s
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> operator/(const TensorGeneric<n_,m_>&,double s);
/// returns the determinant
  double determinant()const;
/// return an identity tensor
  static TensorGeneric<n,n> identity();
/// return the matrix inverse
  TensorGeneric inverse()const;
/// return the transpose matrix
  TensorGeneric<m,n> transpose()const;
/// matrix-matrix multiplication
  template<unsigned n_,unsigned m_,unsigned l_>
  friend TensorGeneric<n_,l_> matmul(const TensorGeneric<n_,m_>&,const TensorGeneric<m_,l_>&);
/// matrix-vector multiplication
  template<unsigned n_,unsigned m_>
  friend VectorGeneric<n_> matmul(const TensorGeneric<n_,m_>&,const VectorGeneric<m_>&);
/// vector-matrix multiplication
  template<unsigned n_,unsigned m_>
  friend VectorGeneric<n_> matmul(const VectorGeneric<m_>&,const TensorGeneric<m_,n_>&);
/// vector-vector multiplication (maps to dotProduct)
  template<unsigned n_>
  friend double matmul(const VectorGeneric<n_>&,const VectorGeneric<n_>&);
/// matrix-matrix-matrix multiplication
  template<unsigned n_,unsigned m_,unsigned l_,unsigned i_>
  friend TensorGeneric<n_,i_> matmul(const TensorGeneric<n_,m_>&,const TensorGeneric<m_,l_>&,const TensorGeneric<l_,i_>&);
/// matrix-matrix-vector multiplication
  template<unsigned n_,unsigned m_,unsigned l_>
  friend VectorGeneric<n_> matmul(const TensorGeneric<n_,m_>&,const TensorGeneric<m_,l_>&,const VectorGeneric<l_>&);
/// vector-matrix-matrix multiplication
  template<unsigned n_,unsigned m_,unsigned l_>
  friend VectorGeneric<l_> matmul(const VectorGeneric<n_>&,const TensorGeneric<n_,m_>&,const TensorGeneric<m_,l_>&);
/// vector-matrix-vector multiplication
  template<unsigned n_,unsigned m_>
  friend double matmul(const VectorGeneric<n_>&,const TensorGeneric<n_,m_>&,const VectorGeneric<m_>&);
/// returns the determinant of a tensor
  friend double determinant(const TensorGeneric<3,3>&);
/// returns the inverse of a tensor (same as inverse())
  friend TensorGeneric<3,3> inverse(const TensorGeneric<3,3>&);
/// returns the transpose of a tensor (same as transpose())
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> transpose(const TensorGeneric<m_,n_>&);
/// returns the transpose of a tensor (same as TensorGeneric(const VectorGeneric&,const VectorGeneric&))
  template<unsigned n_,unsigned m_>
  friend TensorGeneric<n_,m_> extProduct(const VectorGeneric<n>&,const VectorGeneric<m>&);
  friend TensorGeneric<3,3> dcrossDv1(const VectorGeneric<3>&,const VectorGeneric<3>&);
  friend TensorGeneric<3,3> dcrossDv2(const VectorGeneric<3>&,const VectorGeneric<3>&);
  friend TensorGeneric<3,3> VcrossTensor(const VectorGeneric<3>&,const TensorGeneric<3,3>&);
  friend TensorGeneric<3,3> VcrossTensor(const TensorGeneric<3,3>&,const VectorGeneric<3>&);
/// Derivative of a normalized vector
  friend TensorGeneric<3,3> deriNorm(const VectorGeneric<3>&,const TensorGeneric<3,3>&);
/// << operator.
/// Allows printing tensor `t` with `std::cout<<t;`
  template<unsigned n_,unsigned m_>
  friend std::ostream & operator<<(std::ostream &os, const TensorGeneric<n_,m_>&);
/// Diagonalize tensor.
/// Syntax is the same as Matrix::diagMat.
/// In addition, it is possible to call if with m_ smaller than n_. In this case,
/// only the first (smaller) m_ eigenvalues and eigenvectors are retrieved.
/// If case lapack fails (info!=0) it throws an exception.
/// Notice that tensor is assumed to be symmetric!!!
  template<unsigned n_,unsigned m_>
  friend void diagMatSym(const TensorGeneric<n_,n_>&,VectorGeneric<m_>&evals,TensorGeneric<m_,n_>&evec);
};

template <unsigned n,unsigned m>
void TensorGeneric<n,m>::auxiliaryConstructor()
{}

template <unsigned n,unsigned m>
template<typename... Args>
void TensorGeneric<n,m>::auxiliaryConstructor(double first,Args... arg)
{
  d[n*m-(sizeof...(Args))-1]=first;
  auxiliaryConstructor(arg...);
}

template <unsigned n,unsigned m>
template<typename... Args>
TensorGeneric<n,m>::TensorGeneric(double first,Args... arg)
{
  static_assert((sizeof...(Args))+1==n*m,"you are trying to initialize a Tensor with the wrong number of arguments");
  auxiliaryConstructor(first,arg...);
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>::TensorGeneric() {
  LoopUnroller<n*m>::_zero(d.data());
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>::TensorGeneric(const VectorGeneric<n>&v1,const VectorGeneric<m>&v2) {
  for(unsigned i=0; i<n; i++)for(unsigned j=0; j<m; j++)d[i*m+j]=v1[i]*v2[j];
}

template<unsigned n,unsigned m>
void TensorGeneric<n,m>::zero() {
  LoopUnroller<n*m>::_zero(d.data());
}

template<unsigned n,unsigned m>
double & TensorGeneric<n,m>::operator() (unsigned i,unsigned j) {
#ifdef _GLIBCXX_DEBUG
// index i is implicitly checked by the std::array class
  plumed_assert(j<m);
#endif
  return d[m*i+j];
}

template<unsigned n,unsigned m>
const double & TensorGeneric<n,m>::operator() (unsigned i,unsigned j)const {
#ifdef _GLIBCXX_DEBUG
// index i is implicitly checked by the std::array class
  plumed_assert(j<m);
#endif
  return d[m*i+j];
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::operator +=(const TensorGeneric<n,m>& b) {
  LoopUnroller<n*m>::_add(d.data(),b.d.data());
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::operator -=(const TensorGeneric<n,m>& b) {
  LoopUnroller<n*m>::_sub(d.data(),b.d.data());
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::operator *=(double s) {
  LoopUnroller<n*m>::_mul(d.data(),s);
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::operator /=(double s) {
  LoopUnroller<n*m>::_mul(d.data(),1.0/s);
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> TensorGeneric<n,m>::operator+()const {
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> TensorGeneric<n,m>::operator-()const {
  TensorGeneric<n,m> r;
  LoopUnroller<n*m>::_neg(r.d.data(),d.data());
  return r;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::setCol(unsigned j,const VectorGeneric<n> & c) {
  for(unsigned i=0; i<n; ++i) (*this)(i,j)=c(i);
  return *this;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m>& TensorGeneric<n,m>::setRow(unsigned i,const VectorGeneric<m> & r) {
  for(unsigned j=0; j<m; ++j) (*this)(i,j)=r(j);
  return *this;
}

template<unsigned n,unsigned m>
VectorGeneric<n> TensorGeneric<n,m>::getCol(unsigned j)const {
  VectorGeneric<n> v;
  for(unsigned i=0; i<n; ++i) v(i)=(*this)(i,j);
  return v;
}

template<unsigned n,unsigned m>
VectorGeneric<m> TensorGeneric<n,m>::getRow(unsigned i)const {
  VectorGeneric<m> v;
  for(unsigned j=0; j<m; ++j) v(j)=(*this)(i,j);
  return v;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> operator+(const TensorGeneric<n,m>&t1,const TensorGeneric<n,m>&t2) {
  TensorGeneric<n,m> t(t1);
  t+=t2;
  return t;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> operator-(const TensorGeneric<n,m>&t1,const TensorGeneric<n,m>&t2) {
  TensorGeneric<n,m> t(t1);
  t-=t2;
  return t;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> operator*(const TensorGeneric<n,m>&t1,double s) {
  TensorGeneric<n,m> t(t1);
  t*=s;
  return t;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> operator*(double s,const TensorGeneric<n,m>&t1) {
  return t1*s;
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> operator/(const TensorGeneric<n,m>&t1,double s) {
  return t1*(1.0/s);
}

template<>
inline
double TensorGeneric<3,3>::determinant()const {
  return
    d[0]*d[4]*d[8]
    + d[1]*d[5]*d[6]
    + d[2]*d[3]*d[7]
    - d[0]*d[5]*d[7]
    - d[1]*d[3]*d[8]
    - d[2]*d[4]*d[6];
}

template<unsigned n,unsigned m>
inline
TensorGeneric<n,n> TensorGeneric<n,m>::identity() {
  TensorGeneric<n,n> t;
  for(unsigned i=0; i<n; i++) t(i,i)=1.0;
  return t;
}

template<unsigned n,unsigned m>
TensorGeneric<m,n> TensorGeneric<n,m>::transpose()const {
  TensorGeneric<m,n> t;
  for(unsigned i=0; i<m; i++)for(unsigned j=0; j<n; j++) t(i,j)=(*this)(j,i);
  return t;
}

template<>
inline
TensorGeneric<3,3> TensorGeneric<3,3>::inverse()const {
  TensorGeneric t;
  double invdet=1.0/determinant();
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++)
      t(j,i)=invdet*(   (*this)((i+1)%3,(j+1)%3)*(*this)((i+2)%3,(j+2)%3)
                        -(*this)((i+1)%3,(j+2)%3)*(*this)((i+2)%3,(j+1)%3));
  return t;
}

template<unsigned n,unsigned m,unsigned l>
TensorGeneric<n,l> matmul(const TensorGeneric<n,m>&a,const TensorGeneric<m,l>&b) {
  TensorGeneric<n,l> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<l; j++) for(unsigned k=0; k<m; k++) {
        t(i,j)+=a(i,k)*b(k,j);
      }
  return t;
}

template<unsigned n,unsigned m>
VectorGeneric<n> matmul(const TensorGeneric<n,m>&a,const VectorGeneric<m>&b) {
  VectorGeneric<n> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<m; j++) t(i)+=a(i,j)*b(j);
  return t;
}

template<unsigned n,unsigned m>
VectorGeneric<n> matmul(const VectorGeneric<m>&a,const TensorGeneric<m,n>&b) {
  VectorGeneric<n> t;
  for(unsigned i=0; i<n; i++) for(unsigned j=0; j<m; j++) t(i)+=a(j)*b(j,i);
  return t;
}

template<unsigned n_>
double matmul(const VectorGeneric<n_>&a,const VectorGeneric<n_>&b) {
  return dotProduct(a,b);
}

template<unsigned n,unsigned m,unsigned l,unsigned i>
TensorGeneric<n,i> matmul(const TensorGeneric<n,m>&a,const TensorGeneric<m,l>&b,const TensorGeneric<l,i>&c) {
  return matmul(matmul(a,b),c);
}

template<unsigned n,unsigned m,unsigned l>
VectorGeneric<n> matmul(const TensorGeneric<n,m>&a,const TensorGeneric<m,l>&b,const VectorGeneric<l>&c) {
  return matmul(matmul(a,b),c);
}

template<unsigned n,unsigned m,unsigned l>
VectorGeneric<l> matmul(const VectorGeneric<n>&a,const TensorGeneric<n,m>&b,const TensorGeneric<m,l>&c) {
  return matmul(matmul(a,b),c);
}

template<unsigned n,unsigned m>
double matmul(const VectorGeneric<n>&a,const TensorGeneric<n,m>&b,const VectorGeneric<m>&c) {
  return matmul(matmul(a,b),c);
}

inline
double determinant(const TensorGeneric<3,3>&t) {
  return t.determinant();
}

inline
TensorGeneric<3,3> inverse(const TensorGeneric<3,3>&t) {
  return t.inverse();
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> transpose(const TensorGeneric<m,n>&t) {
  return t.transpose();
}

template<unsigned n,unsigned m>
TensorGeneric<n,m> extProduct(const VectorGeneric<n>&v1,const VectorGeneric<m>&v2) {
  return TensorGeneric<n,m>(v1,v2);
}

inline
TensorGeneric<3,3> dcrossDv1(const VectorGeneric<3>&v1,const VectorGeneric<3>&v2) {
  (void) v1; // this is to avoid warnings. still the syntax of this function is a bit dummy...
  return TensorGeneric<3,3>(
           0.0, v2[2],-v2[1],
           -v2[2],   0.0, v2[0],
           v2[1],-v2[0],   0.0);
}

inline
TensorGeneric<3,3> dcrossDv2(const VectorGeneric<3>&v1,const VectorGeneric<3>&v2) {
  (void) v2; // this is to avoid warnings. still the syntax of this function is a bit dummy...
  return TensorGeneric<3,3>(
           0.0,-v1[2],v1[1],
           v1[2],0.0,-v1[0],
           -v1[1],v1[0],0.0);
}

template<unsigned n,unsigned m>
std::ostream & operator<<(std::ostream &os, const TensorGeneric<n,m>& t) {
  for(unsigned i=0; i<n; i++)for(unsigned j=0; j<m; j++) {
      if(i!=(n-1) || j!=(m-1)) os<<t(i,j)<<" ";
      else os<<t(i,j);
    }
  return os;
}

/// \ingroup TOOLBOX
typedef TensorGeneric<1,1> Tensor1d;
/// \ingroup TOOLBOX
typedef TensorGeneric<2,2> Tensor2d;
/// \ingroup TOOLBOX
typedef TensorGeneric<3,3> Tensor3d;
/// \ingroup TOOLBOX
typedef TensorGeneric<4,4> Tensor4d;
/// \ingroup TOOLBOX
typedef TensorGeneric<5,5> Tensor5d;
/// \ingroup TOOLBOX
typedef Tensor3d Tensor;

inline
TensorGeneric<3,3> VcrossTensor(const VectorGeneric<3>&v1,const TensorGeneric<3,3>&v2) {

  TensorGeneric<3,3> t;
  for(unsigned i=0; i<3; i++) {
    t.setRow(i,matmul(dcrossDv2(v1,v1),v2.getRow(i)));
  }
  return t;
}

inline
TensorGeneric<3,3> VcrossTensor(const TensorGeneric<3,3>&v2,const VectorGeneric<3>&v1) {
  TensorGeneric<3,3> t;
  for(unsigned i=0; i<3; i++) {
    t.setRow(i,-matmul(dcrossDv2(v1,v1),v2.getRow(i)));
  }
  return t;
}


inline
TensorGeneric<3,3> deriNorm(const VectorGeneric<3>&v1,const TensorGeneric<3,3>&v2) {
  // delta(v) = delta(v1/v1.norm) = 1/v1.norm*(delta(v1) - (v.delta(v1))cross v;
  double over_norm = 1./v1.modulo();
  return over_norm*(v2 - over_norm*over_norm*(extProduct(matmul(v2,v1),v1)));
}

template<unsigned n,unsigned m>
void diagMatSym(const TensorGeneric<n,n>&mat,VectorGeneric<m>&evals,TensorGeneric<m,n>&evec) {
  // some guess number to make sure work is large enough.
  // for correctness it should be >=20. However, it is recommended to be the block size.
  // I put some likely exaggerated number
  constexpr int bs=100;
  // temporary data, on stack so as to avoid allocations
  std::array<int,10*n> iwork;
  std::array<double,(6+bs)*n> work;
  std::array<int,2*m> isup;
  int nn=n;              // dimension of matrix
  double vl=0.0, vu=1.0; // ranges - not used
  int one=1,mm=m;        // minimun and maximum index
  double abstol=0.0;     // tolerance
  int mout=0;            // number of eigenvalues found (same as mm)
  int info=0;            // result
  int liwork=iwork.size();
  int lwork=work.size();
  TensorGenericAux::local_dsyevr("V", (n==m?"A":"I"), "U", &nn, const_cast<double*>(&mat[0][0]), &nn, &vl, &vu, &one, &mm,
                                 &abstol, &mout, &evals[0], &evec[0][0], &nn,
                                 isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if(info!=0) plumed_error()<<"Error diagonalizing matrix\n"
                              <<"Matrix:\n"<<mat<<"\n"
                              <<"Info: "<<info<<"\n";
  plumed_assert(mout==m);
  // This changes eigenvectors so that the first non-null element
  // of each of them is positive
  // We can do it because the phase is arbitrary, and helps making
  // the result reproducible
  for(int i=0; i<m; ++i) {
    int j=0;
    for(j=0; j<n; j++) if(evec(i,j)*evec(i,j)>1e-14) break;
    if(j<n) if(evec(i,j)<0.0) for(j=0; j<n; j++) evec(i,j)*=-1;
  }
}


}

#endif

