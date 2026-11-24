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
#ifndef __PLUMED_tools_Matrix_h
#define __PLUMED_tools_Matrix_h
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include "Exception.h"
#include "MatrixSquareBracketsAccess.h"
#include "Tools.h"
#include "Log.h"
#include "lapack/lapack.h"

namespace PLMD {

/// Calculate the dot product between two vectors
template <typename T> T dotProduct( const std::vector<T>& A, const std::vector<T>& B ) {
  plumed_assert( A.size()==B.size() );
  T val{};
  for(unsigned i=0; i<A.size(); ++i) {
    val+=A[i]*B[i];
  }
  return val;
}

/// Calculate the dot product between a vector and itself
template <typename T> T norm( const std::vector<T>& A ) {
  T val{};
  for(unsigned i=0; i<A.size(); ++i) {
    val+=A[i]*A[i];
  }
  return val;
}

/// This class stores a full matrix and allows one to do some simple matrix operations
template <typename T>
class Matrix:
  public MatrixSquareBracketsAccess<Matrix<T>,T> {
  /// Multiply matrix by scalar
  template <typename U> friend Matrix<U> operator*(U&, const Matrix<U>& );
  /// Matrix matrix multiply
  template <typename U> friend void mult( const Matrix<U>&, const Matrix<U>&, Matrix<U>& );
  /// Matrix times a std::vector
  template <typename U> friend void mult( const Matrix<U>&, const std::vector<U>&, std::vector<U>& );
  /// std::vector times a Matrix
  template <typename U> friend void mult( const std::vector<U>&, const Matrix<U>&, std::vector<U>& );
  /// Matrix transpose
  template <typename U> friend void transpose( const Matrix<U>&, Matrix<U>& );
  /// Output the entire matrix on a single line
  template <typename U> friend Log& operator<<(Log&, const Matrix<U>& );
  /// Output the Matrix in matrix form
  template <typename U> friend void matrixOut( Log&, const Matrix<U>& );
  /// Diagonalize a symmetric matrix - returns zero if diagonalization worked
  template <typename U> friend int diagMat( const Matrix<U>&, std::vector<U>&, Matrix<U>& );
  /// Diagonalize a symmetric matrix computing only the eigenpairs from first to last (0-based, ascending eigenvalues) - returns zero if diagonalization worked
  template <typename U> friend int diagMatRange( const Matrix<U>&, unsigned, unsigned, std::vector<U>&, Matrix<U>& );
  /// Calculate the Moore-Penrose Pseudoinverse of a matrix
  template <typename U> friend int pseudoInvert( const Matrix<U>&, Matrix<U>& );
  /// Calculate the logarithm of the determinant of a symmetric matrix - returns zero if successful
  template <typename U> friend int logdet( const Matrix<U>&, U& );
  /// Invert a matrix (works for both symmetric and asymmetric matrices) - returns zero if sucesfull
  template <typename U> friend int Invert( const Matrix<U>&, Matrix<U>& );
  /// Do a cholesky decomposition of a matrix
  template <typename U> friend void cholesky( const Matrix<U>&, Matrix<U>& );
  /// Solve a system of equations using the cholesky decomposition
  template <typename U> friend void chol_elsolve( const Matrix<U>&, const std::vector<U>&, std::vector<U>& );
private:
  /// Number of elements in matrix (nrows*ncols)
  unsigned sz;
  /// Number of rows in matrix
  unsigned rw;
  /// Number of columns in matrix
  unsigned cl;
  /// The data in the matrix
  std::vector<T> data;
  template <typename Y>
  struct  openAccVectorHelper {
    using pointer = Y*;
    size_t sz=0;
    pointer data=nullptr;
    void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],sz,data[0:sz])
    }
    void removeFromACCDevice() const {
#pragma acc exit data delete(data[0:sz],sz,this[0:1])
    }
    void update(std::vector<Y>&x) {
      data = x.data();
      sz = x.size();
    }
  };
  openAccVectorHelper<T> helper{};

public:
  explicit Matrix(const unsigned nr=0, const unsigned nc=0 )  : sz(nr*nc), rw(nr), cl(nc), data(nr*nc) {
    helper.update(data);
  }
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    helper.toACCDevice();
  }
  void removeFromACCDevice() const {
    helper.removeFromACCDevice();
    // and delete
#pragma acc exit data delete( this[0:1])
  }
  Matrix(const Matrix<T>& t) : sz(t.sz), rw(t.rw), cl(t.cl), data(t.data) {
    helper.update(data);
  }
  /// Resize the matrix
  void resize( const unsigned nr, const unsigned nc ) {
    rw=nr;
    cl=nc;
    sz=nr*nc;
    data.resize(sz);
    helper.update(data);
  }
  /// Return the number of rows
  inline unsigned nrows() const {
    return rw;
  }
  /// Return the number of columns
  inline unsigned ncols() const {
    return cl;
  }
  /// Return the contents of the matrix as a vector of length rw*cl
  inline std::vector<T>& getVector() {
    return data;
  }
  /// Set the matrix from a vector input
  inline void setFromVector( const std::vector<T>& vecin ) {
    plumed_assert( vecin.size()==sz );
    for(unsigned i=0; i<sz; ++i) {
      data[i]=vecin[i];
    }
  }
  /// Return element i,j of the matrix
  inline const T& operator () (const unsigned& i, const unsigned& j) const {
    return helper.data[j+i*cl];
  }
  /// Return a referenre to element i,j of the matrix
  inline T& operator () (const unsigned& i, const unsigned& j)      {
    return helper.data[j+i*cl];
  }
  /// Set all elements of the matrix equal to the value of v
  Matrix<T>& operator=(const T& v) {
    for(unsigned i=0; i<sz; ++i) {
      data[i]=v;
    }
    return *this;
  }
  /// Set the Matrix equal to another Matrix
  Matrix<T>& operator=(const Matrix<T>& m) {
    sz=m.sz;
    rw=m.rw;
    cl=m.cl;
    data=m.data;
    helper.update(data);
    return *this;
  }
  /// Set the Matrix equal to the value of a standard vector - used for readin
  Matrix<T>& operator=(const std::vector<T>& v) {
    plumed_dbg_assert( v.size()==sz );
    for(unsigned i=0; i<sz; ++i) {
      data[i]=v[i];
    }
    return *this;
  }
  /// Add v to all elements of the Matrix
  Matrix<T> operator+=(const T& v) {
    for(unsigned i=0; i<sz; ++i) {
      data[i]+=v;
    }
    return *this;
  }
  /// Multiply all elements by v
  Matrix<T> operator*=(const T& v) {
    for(unsigned i=0; i<sz; ++i) {
      data[i]*=v;
    }
    return *this;
  }
  /// Matrix addition
  Matrix<T>& operator+=(const Matrix<T>& m) {
    plumed_dbg_assert( m.rw==rw && m.cl==cl );
    data+=m.data;
    return *this;
  }
  /// Subtract v from all elements of the Matrix
  Matrix<T> operator-=(const T& v) {
    for(unsigned i=0; i<sz; ++i) {
      data-=v;
    }
    return *this;
  }
  /// Matrix subtraction
  Matrix<T>& operator-=(const Matrix<T>& m) {
    plumed_dbg_assert( m.rw==rw && m.cl==cl );
    data-=m.data;
    return *this;
  }
  /// Test if the matrix is symmetric or not
  unsigned isSymmetric() const {
    if (rw!=cl) {
      return 0;
    }
    unsigned sym=1;
    for(unsigned i=1; i<rw; ++i)
      for(unsigned j=0; j<i; ++j)
        if( std::fabs(data[i+j*cl]-data[j+i*cl])>1.e-10 ) {
          sym=0;
          break;
        }
    return sym;
  }
};

/// Multiply matrix by scalar
template <typename T> Matrix<T> operator*(T& v, const Matrix<T>& m ) {
  Matrix<T> new_m(m);
  new_m*=v;
  return new_m;
}

template <typename T> void mult( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ) {
  plumed_assert(A.cl==B.rw);
  if( A.rw !=C.rw  || B.cl !=C.cl ) {
    C.resize( A.rw, B.cl );
  }
  C=static_cast<T>( 0 );
  // i-k-j loop order: the innermost loop runs over contiguous elements of B and C,
  // so it is cache friendly and can be vectorized. For each element of C the
  // products are summed in the same order as with the i-j-k loop order.
  const T* a=A.data.data();
  const T* b=B.data.data();
  T* c=C.data.data();
  for(unsigned i=0; i<A.rw; ++i) {
    T* ci=c+i*B.cl;
    for(unsigned k=0; k<A.cl; ++k) {
      const T aik=a[i*A.cl+k];
      const T* bk=b+k*B.cl;
      for(unsigned j=0; j<B.cl; ++j) {
        ci[j]+=aik*bk[j];
      }
    }
  }
}

template <typename T> void mult( const Matrix<T>& A, const std::vector<T>& B, std::vector<T>& C) {
  plumed_assert( A.cl==B.size() );
  if( C.size()!=A.rw  ) {
    C.resize(A.rw);
  }
  for(unsigned i=0; i<A.rw; ++i) {
    C[i]= static_cast<T>( 0 );
  }
  for(unsigned i=0; i<A.rw; ++i)
    for(unsigned k=0; k<A.cl; ++k) {
      C[i]+=A(i,k)*B[k] ;
    }
}

template <typename T> void mult( const std::vector<T>& A, const Matrix<T>& B, std::vector<T>& C) {
  plumed_assert( B.rw==A.size() );
  if( C.size()!=B.cl ) {
    C.resize( B.cl );
  }
  for(unsigned i=0; i<B.cl; ++i) {
    C[i]=static_cast<T>( 0 );
  }
  // loop on the rows of B in the outer loop, so that the inner loop is on contiguous elements
  for(unsigned k=0; k<B.rw; ++k)
    for(unsigned i=0; i<B.cl; ++i) {
      C[i]+=A[k]*B(k,i);
    }
}

template <typename T> void transpose( const Matrix<T>& A, Matrix<T>& AT ) {
  if( A.rw!=AT.cl || A.cl!=AT.rw ) {
    AT.resize( A.cl, A.rw );
  }
  for(unsigned i=0; i<A.cl; ++i)
    for(unsigned j=0; j<A.rw; ++j) {
      AT(i,j)=A(j,i);
    }
}

template <typename T> Log& operator<<(Log& ostr, const Matrix<T>& mat) {
  for(unsigned i=0; i<mat.sz; ++i) {
    ostr<<mat.data[i]<<" ";
  }
  return ostr;
}

template <typename T> void matrixOut( Log& ostr, const Matrix<T>& mat) {
  for(unsigned i=0; i<mat.rw; ++i) {
    for(unsigned j=0; j<mat.cl; ++j) {
      ostr<<mat(i,j)<<" ";
    }
    ostr<<"\n";
  }
  return;
}

/// Cholesky factorization A=L*L^T of a symmetric positive definite n x n matrix.
/// Only the upper triangle of A is read (as done by diagMat). L is returned as a lower
/// triangular matrix stored row by row in a vector of size n*n.
/// Returns false if A is not positive definite (or contains NaNs).
template <typename T> bool choleskyFactorization( const Matrix<T>& A, std::vector<T>& L ) {
  const unsigned n=A.nrows();
  plumed_assert( A.ncols()==n );
  // Right-looking factorization A=U^T*U working on the rows of the upper triangle,
  // so that the innermost loop is a vectorizable update of contiguous elements
  std::vector<T> U(static_cast<std::size_t>(n)*n,static_cast<T>(0));
  for(unsigned j=0; j<n; ++j)
    for(unsigned i=j; i<n; ++i) {
      U[static_cast<std::size_t>(j)*n+i]=A(j,i);
    }
  for(unsigned j=0; j<n; ++j) {
    T* Uj=U.data()+static_cast<std::size_t>(j)*n;
    if( !(Uj[j]>0) ) {
      return false;
    }
    Uj[j]=std::sqrt(Uj[j]);
    const T invd=static_cast<T>(1)/Uj[j];
    for(unsigned i=j+1; i<n; ++i) {
      Uj[i]*=invd;
    }
    for(unsigned i=j+1; i<n; ++i) {
      const T uji=Uj[i];
      T* Ui=U.data()+static_cast<std::size_t>(i)*n;
      for(unsigned k=i; k<n; ++k) {
        Ui[k]-=uji*Uj[k];
      }
    }
  }
  L.assign(static_cast<std::size_t>(n)*n,static_cast<T>(0));
  for(unsigned i=0; i<n; ++i)
    for(unsigned j=0; j<=i; ++j) {
      L[static_cast<std::size_t>(i)*n+j]=U[static_cast<std::size_t>(j)*n+i];
    }
  return true;
}

/// Inverse of a symmetric positive definite matrix from its Cholesky factorization.
/// Returns false (and leaves inverse untouched) if A is not positive definite.
template <typename T> bool choleskyInvert( const Matrix<T>& A, Matrix<T>& inverse ) {
  const unsigned n=A.nrows();
  std::vector<T> L;
  if( !choleskyFactorization( A, L ) ) {
    return false;
  }
  // Compute M=L^{-1} (lower triangular) row by row: M_i = (e_i - sum_{k<i} L_ik M_k) / L_ii
  std::vector<T> M(static_cast<std::size_t>(n)*n,static_cast<T>(0));
  for(unsigned i=0; i<n; ++i) {
    const T* Li=L.data()+static_cast<std::size_t>(i)*n;
    T* Mi=M.data()+static_cast<std::size_t>(i)*n;
    Mi[i]=static_cast<T>(1);
    for(unsigned k=0; k<i; ++k) {
      const T lik=Li[k];
      const T* Mk=M.data()+static_cast<std::size_t>(k)*n;
      for(unsigned j=0; j<=k; ++j) {
        Mi[j]-=lik*Mk[j];
      }
    }
    const T invd=static_cast<T>(1)/Li[i];
    for(unsigned j=0; j<=i; ++j) {
      Mi[j]*=invd;
    }
  }
  // A^{-1} = M^T M = sum_k (row k of M)^T (row k of M); compute the lower triangle
  if( inverse.nrows()!=n || inverse.ncols()!=n ) {
    inverse.resize(n,n);
  }
  inverse=static_cast<T>(0);
  for(unsigned k=0; k<n; ++k) {
    const T* Mk=M.data()+static_cast<std::size_t>(k)*n;
    for(unsigned p=0; p<=k; ++p) {
      const T mkp=Mk[p];
      for(unsigned q=0; q<=p; ++q) {
        inverse(p,q)+=mkp*Mk[q];
      }
    }
  }
  for(unsigned p=0; p<n; ++p)
    for(unsigned q=0; q<p; ++q) {
      inverse(q,p)=inverse(p,q);
    }
  return true;
}

template <typename T> int diagMatRange( const Matrix<T>& A, unsigned first, unsigned last, std::vector<T>& eigenvals, Matrix<T>& eigenvecs ) {

  // Check matrix is square and symmetric
  plumed_assert( A.rw==A.cl );
  plumed_assert( A.isSymmetric()==1 );
  if( eigenvals.size()!=A.cl ) {
    eigenvals.resize( A.cl );
  }
  if( eigenvecs.rw!=A.rw || eigenvecs.cl!=A.cl ) {
    eigenvecs.resize( A.rw, A.cl );
  }
  if( A.cl==0 ) {
    return 0;
  }
  plumed_assert( first<=last && last<A.cl );
  std::vector<T> da(A.sz);
  unsigned k=0;
  std::vector<T> evals(A.cl);
  // Transfer the matrix to the local array
  for (unsigned i=0; i<A.cl; ++i)
    for (unsigned j=0; j<A.rw; ++j) {
      da[k++]=static_cast<T>( A(j,i) );
    }

  int n=A.cl;
  int lwork=-1;
  int liwork=-1;
  int m;
  int info=0;
#ifdef __PLUMED_HAS_EXTERNAL_LAPACK
  const bool partial=true;
#else
  // The internal dsyevr reads out of bounds when only a subset of the eigenpairs is
  // requested, so all of them are computed and then the requested ones are copied
  const bool partial=false;
#endif
  // dsyevr numbers the eigenvalues from 1, in ascending order
  const unsigned offset=( partial ? first : 0 );
  int il=offset+1;
  int iu=( partial ? last+1 : A.cl );
  const unsigned nev=iu-il+1;
  std::vector<T> work(A.cl);
  std::vector<int> iwork(A.cl);
  T vl, vu, abstol=0.0;
  std::vector<int> isup(2*nev);
  std::vector<T> evecs(static_cast<std::size_t>(A.cl)*nev);
  if constexpr (std::is_same_v<T,double>) {
    plumed_lapack_dsyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_ssyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "diagMat is only compatible with float and double");
  }
  if (info!=0) {
    return info;
  }

  // Retrieve correct sizes for work and iwork then reallocate
  liwork=iwork[0];
  iwork.resize(liwork);
  lwork=static_cast<int>( work[0] );
  work.resize(lwork);

  if constexpr (std::is_same_v<T,double>) {
    plumed_lapack_dsyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_ssyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &il, &iu,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "diagMat is only compatible with float and double");
  }
  if (info!=0) {
    return info;
  }

  k=(first-offset)*A.rw;
  for(unsigned i=first; i<=last; ++i) {
    eigenvals[i]=evals[i-offset];
    // N.B. For ease of producing projectors we store the eigenvectors
    // ROW-WISE in the eigenvectors matrix.  The first index is the
    // eigenvector number and the second the component
    for(unsigned j=0; j<A.rw; ++j) {
      eigenvecs(i,j)=evecs[k++];
    }
  }

  // This changes eigenvectors so that the first non-null element
  // of each of them is positive
  // We can do it because the phase is arbitrary, and helps making
  // the result reproducible
  for(unsigned i=first; i<=last; ++i) {
    int j;
    for(j=0; j<n; j++)
      if(eigenvecs(i,j)*eigenvecs(i,j)>1e-14) {
        break;
      }
    if(j<n)
      if(eigenvecs(i,j)<0.0)
        for(j=0; j<n; j++) {
          eigenvecs(i,j)*=-1;
        }
  }
  return 0;
}

template <typename T> int diagMat( const Matrix<T>& A, std::vector<T>& eigenvals, Matrix<T>& eigenvecs ) {
  if( A.cl==0 ) {
    return diagMatRange( A, 0, 0, eigenvals, eigenvecs );
  }
  return diagMatRange( A, 0, A.cl-1, eigenvals, eigenvecs );
}

template <typename T> int pseudoInvert( const Matrix<T>& A, Matrix<T>& pseudoinverse ) {
  const unsigned nrows=A.rw, ncols=A.cl;
  const unsigned nsv=std::min(nrows,ncols);
  if( pseudoinverse.rw!=ncols || pseudoinverse.cl!=nrows ) {
    pseudoinverse.resize( ncols, nrows );
  }
  pseudoinverse=static_cast<T>( 0 );
  if( nsv==0 ) {
    return 0;
  }

  std::vector<T> da(A.sz);
  unsigned k=0;
  // Transfer the matrix to the local array
  for (unsigned i=0; i<A.cl; ++i)
    for (unsigned j=0; j<A.rw; ++j) {
      da[k++]=static_cast<T>( A(j,i) );
    }

  int info=0;
  int nr=nrows, nc=ncols;

  // Create some containers for stuff from the singular value decomposition A=U*S*VT
  // (U is nrows x nrows and VT is ncols x ncols, both stored column by column).
  // The thin decomposition (JOBZ='S') would be enough, but the internal dgesdd
  // writes out of bounds with JOBZ='S', so the full one is used.
  std::vector<T> S(nsv);
  std::vector<T> U(static_cast<std::size_t>(nrows)*nrows);
  std::vector<T> VT(static_cast<std::size_t>(ncols)*ncols);
  std::vector<int> iwork(8*nsv);

  // This optimizes the size of the work array used in lapack singular value decomposition
  int lwork=-1;
  std::vector<T> work(1);
  if constexpr (std::is_same_v<T,double> ) {
    plumed_lapack_dgesdd( "A", &nr, &nc, da.data(), &nr, S.data(),
                          U.data(), &nr, VT.data(), &nc, work.data(),
                          &lwork, iwork.data(), &info );
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_sgesdd( "A", &nr, &nc, da.data(), &nr, S.data(),
                          U.data(), &nr, VT.data(), &nc, work.data(),
                          &lwork, iwork.data(), &info );
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "pseudoInvert works only for float and double");
  }
  if(info!=0) {
    return info;
  }

  // Retrieve correct sizes for work and rellocate
  lwork=(int) work[0];
  work.resize(lwork);

  // This does the singular value decomposition
  if constexpr (std::is_same_v<T,double> ) {
    plumed_lapack_dgesdd( "A", &nr, &nc, da.data(), &nr, S.data(),
                          U.data(), &nr, VT.data(), &nc, work.data(),
                          &lwork, iwork.data(), &info );
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_sgesdd( "A", &nr, &nc, da.data(), &nr, S.data(),
                          U.data(), &nr, VT.data(), &nc, work.data(),
                          &lwork, iwork.data(), &info );
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "pseudoInvert works only for float and double");
  }
  if(info!=0) {
    return info;
  }

  // Compute the tolerance on the singular values ( machine epsilon * number of singular values * maximum singular value )
  T tol;
  tol=S[0];
  for(unsigned i=1; i<nsv; ++i) {
    if( S[i]>tol ) {
      tol=S[i];
    }
  }
  tol*=nsv*epsilon;

  // pseudoinverse = V * S^{-1} * U^T = sum_k (column k of V) (column k of U)^T / S_k,
  // summing only over the singular values larger than the tolerance
  for(unsigned isv=0; isv<nsv; ++isv) {
    if( !(S[isv]>tol) ) {
      continue;
    }
    const T inv=static_cast<T>(1)/S[isv];
    const T* Uk=U.data()+static_cast<std::size_t>(isv)*nrows;
    for(unsigned i=0; i<ncols; ++i) {
      const T a=inv*VT[isv+static_cast<std::size_t>(i)*ncols];
      for(unsigned j=0; j<nrows; ++j) {
        pseudoinverse(i,j)+=a*Uk[j];
      }
    }
  }

  return 0;
}

template <typename T> int Invert( const Matrix<T>& A, Matrix<T>& inverse ) {

  if( A.isSymmetric()==1 ) {
    // Symmetric positive definite matrices (e.g. covariance and metric matrices) are inverted
    // with a Cholesky factorization, which is several times faster than diagonalizing them.
    if( choleskyInvert( A, inverse ) ) {
      return 0;
    }
    // GAT -- I only ever use symmetric matrices so I can invert them like this.
    // I choose to do this as I have had problems with the more general way of doing this that
    // is implemented below.
    std::vector<T> eval(A.rw);
    Matrix<T> evec(A.rw,A.cl), tevec(A.rw,A.cl);
    int err;
    err=diagMat( A, eval, evec );
    if(err!=0) {
      return err;
    }
    for (unsigned i=0; i<A.rw; ++i)
      for (unsigned j=0; j<A.cl; ++j) {
        tevec(i,j)=evec(j,i)/eval[j];
      }
    mult(tevec,evec,inverse);
  } else {
    std::vector<T> da(A.sz);
    std::vector<int> ipiv(A.cl);
    unsigned k=0;
    int n=A.rw;
    int info=0;
    for(unsigned i=0; i<A.cl; ++i)
      for(unsigned j=0; j<A.rw; ++j) {
        da[k++]=static_cast<T>( A(j,i) );
      }
    if constexpr (std::is_same_v<T,double> ) {
      plumed_lapack_dgetrf(&n,&n,da.data(),&n,ipiv.data(),&info);
    } else if constexpr (std::is_same_v<T,float>) {
      plumed_lapack_sgetrf(&n,&n,da.data(),&n,ipiv.data(),&info);
    } else {
      static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "Invert works only for float and double");
    }
    if(info!=0) {
      return info;
    }

    int lwork=-1;
    std::vector<T> work(A.cl);
    if constexpr (std::is_same_v<T,double> ) {
      plumed_lapack_dgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    } else if constexpr (std::is_same_v<T,float>) {
      plumed_lapack_sgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    } else {
      static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "Invert works only for float and double");
    }
    if(info!=0) {
      return info;
    }

    lwork=static_cast<int>( work[0] );
    work.resize(lwork);
    if constexpr (std::is_same_v<T,double> ) {
      plumed_lapack_dgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    } else if constexpr (std::is_same_v<T,float>) {
      plumed_lapack_sgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    } else {
      static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "Invert works only for float and double");
    }
    if(info!=0) {
      return info;
    }

    if( inverse.cl!=A.cl || inverse.rw!=A.rw ) {
      inverse.resize(A.rw,A.cl);
    }
    k=0;
    for(unsigned i=0; i<A.rw; ++i)
      for(unsigned j=0; j<A.cl; ++j) {
        inverse(j,i)=da[k++];
      }
  }

  return 0;
}

template <typename T> void cholesky( const Matrix<T>& A, Matrix<T>& B ) {

  plumed_assert( A.rw==A.cl && A.isSymmetric() );
  Matrix<T> L(A.rw,A.cl);
  L=0.;
  std::vector<T> D(A.rw,0.);
  for(unsigned i=0; i<A.rw; ++i) {
    L(i,i)=static_cast<T>( 1 );
    for (unsigned j=0; j<i; ++j) {
      L(i,j)=A(i,j);
      for (unsigned k=0; k<j; ++k) {
        L(i,j)-=L(i,k)*L(j,k)*D[k];
      }
      if (D[j]!=0.) {
        L(i,j)/=D[j];
      } else {
        L(i,j)=static_cast<T>( 0 );
      }
    }
    D[i]=A(i,i);
    for (unsigned k=0; k<i; ++k) {
      D[i]-=L(i,k)*L(i,k)*D[k];
    }
  }

  for(unsigned i=0; i<A.rw; ++i) {
    D[i]=(D[i]>0.?std::sqrt(D[i]):0.);
  }
  if( B.rw!=A.rw || B.cl!=A.cl ) {
    B.resize( A.rw, A.cl);
  }
  B=0.;
  for(unsigned i=0; i<A.rw; ++i)
    for(unsigned j=0; j<=i; ++j) {
      B(i,j)+=L(i,j)*D[j];
    }
}

template <typename T> void chol_elsolve( const Matrix<T>& M, const std::vector<T>& b, std::vector<T>& y ) {

  plumed_assert( M.rw==M.cl && M(0,1)==0.0 && b.size()==M.rw );
  if( y.size()!=M.rw ) {
    y.resize( M.rw );
  }
  for(unsigned i=0; i<M.rw; ++i) {
    y[i]=b[i];
    for(unsigned j=0; j<i; ++j) {
      y[i]-=M(i,j)*y[j];
    }
    y[i]*=1.0/M(i,i);
  }
}

template <typename T> int logdet( const Matrix<T>& M, T& ldet ) {
  // Check matrix is square and symmetric
  plumed_assert( M.rw==M.cl || M.isSymmetric() );

  // For a positive definite matrix, log(det(M)) = 2 sum_i log(L_ii) with M=L*L^T
  {
    std::vector<T> L;
    if( choleskyFactorization( M, L ) ) {
      ldet=0;
      for(unsigned i=0; i<M.rw; ++i) {
        ldet+=2*std::log(L[static_cast<std::size_t>(i)*M.rw+i]);
      }
      return 0;
    }
  }

  std::vector<T> da(M.sz);
  unsigned k=0;
  std::vector<T> evals(M.cl);
  // Transfer the matrix to the local array
  for (unsigned i=0; i<M.rw; ++i)
    for (unsigned j=0; j<M.cl; ++j) {
      da[k++]=static_cast<T>( M(j,i) );
    }

  int n=M.cl;
  int lwork=-1;
  int liwork=-1;
  int info=0;
  int m;
  int one=1;
  std::vector<T> work(M.rw);
  std::vector<int> iwork(M.rw);
  T vl, vu, abstol=0.0;
  std::vector<int> isup(2*M.rw);
  std::vector<T> evecs(M.sz);
  if constexpr (std::is_same_v<T,double> ) {
    plumed_lapack_dsyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_ssyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "logdet works only for float and double");
  }
  if (info!=0) {
    return info;
  }

  // Retrieve correct sizes for work and iwork then reallocate
  lwork=static_cast<int>( work[0] );
  work.resize(lwork);
  liwork=iwork[0];
  iwork.resize(liwork);
  if constexpr (std::is_same_v<T,double> ) {
    plumed_lapack_dsyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else if constexpr (std::is_same_v<T,float>) {
    plumed_lapack_ssyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                         &abstol, &m, evals.data(), evecs.data(), &n,
                         isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  } else {
    static_assert(std::is_same_v<T,double> || std::is_same_v<T,float>, "logdet works only for float and double");
  }
  if (info!=0) {
    return info;
  }

  // Transfer the eigenvalues and eigenvectors to the output
  ldet=0;
  for(unsigned i=0; i<M.cl; i++) {
    ldet+=log(evals[i]);
  }

  return 0;
}



}
#endif
