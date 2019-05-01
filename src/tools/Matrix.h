/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
  T val; for(unsigned i=0; i<A.size(); ++i) { val+=A[i]*B[i]; }
  return val;
}

/// Calculate the dot product between a vector and itself
template <typename T> T norm( const std::vector<T>& A ) {
  T val; for(unsigned i=0; i<A.size(); ++i) { val+=A[i]*A[i]; }
  return val;
}

/// This class stores a full matrix and allows one to do some simple matrix operations
template <typename T>
class Matrix:
  public MatrixSquareBracketsAccess<Matrix<T>,T>
{
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
  template <typename U> friend int diagMat( const Matrix<U>&, std::vector<double>&, Matrix<double>& );
  /// Calculate the Moore-Penrose Pseudoinverse of a matrix
  template <typename U> friend int pseudoInvert( const Matrix<U>&, Matrix<double>& );
  /// Calculate the logarithm of the determinant of a symmetric matrix - returns zero if succesfull
  template <typename U> friend int logdet( const Matrix<U>&, double& );
  /// Invert a matrix (works for both symmetric and assymetric matrices) - returns zero if sucesfull
  template <typename U> friend int Invert( const Matrix<U>&, Matrix<double>& );
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
public:
  Matrix(const unsigned nr=0, const unsigned nc=0 )  : sz(nr*nc), rw(nr), cl(nc), data(nr*nc) {}
  Matrix(const Matrix<T>& t) : sz(t.sz), rw(t.rw), cl(t.cl), data(t.data) {}
  /// Resize the matrix
  void resize( const unsigned nr, const unsigned nc ) { rw=nr; cl=nc; sz=nr*nc; data.resize(sz); }
  /// Return the number of rows
  inline unsigned nrows() const { return rw; }
  /// Return the number of columns
  inline unsigned ncols() const { return cl; }
  /// Return the contents of the matrix as a vector of length rw*cl
  inline std::vector<T>& getVector() { return data; }
  /// Set the matrix from a vector input
  inline void setFromVector( const std::vector<T>& vecin ) { plumed_assert( vecin.size()==sz ); for(unsigned i=0; i<sz; ++i) data[i]=vecin[i]; }
  /// Return element i,j of the matrix
  inline T operator () (const unsigned& i, const unsigned& j) const { return data[j+i*cl]; }
  /// Return a referenre to element i,j of the matrix
  inline T& operator () (const unsigned& i, const unsigned& j)      { return data[j+i*cl]; }
  /// Set all elements of the matrix equal to the value of v
  Matrix<T>& operator=(const T& v) {
    for(unsigned i=0; i<sz; ++i) { data[i]=v; }
    return *this;
  }
  /// Set the Matrix equal to another Matrix
  Matrix<T>& operator=(const Matrix<T>& m) {
    sz=m.sz;
    rw=m.rw;
    cl=m.cl;
    data=m.data;
    return *this;
  }
  /// Set the Matrix equal to the value of a standard vector - used for readin
  Matrix<T>& operator=(const std::vector<T>& v) {
    plumed_dbg_assert( v.size()==sz );
    for(unsigned i=0; i<sz; ++i) { data[i]=v[i]; }
    return *this;
  }
  /// Add v to all elements of the Matrix
  Matrix<T> operator+=(const T& v) {
    for(unsigned i=0; i<sz; ++i) { data[i]+=v; }
    return *this;
  }
  /// Multiply all elements by v
  Matrix<T> operator*=(const T& v) {
    for(unsigned i=0; i<sz; ++i) { data[i]*=v; }
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
    for(unsigned i=0; i<sz; ++i) { data-=v; }
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
    if (rw!=cl) { return 0; }
    unsigned sym=1;
    for(unsigned i=1; i<rw; ++i) for(unsigned j=0; j<i; ++j) if( std::fabs(data[i+j*cl]-data[j+i*cl])>1.e-10 ) { sym=0; break; }
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
  if( A.rw !=C.rw  || B.cl !=C.cl ) { C.resize( A.rw, B.cl ); } C=static_cast<T>( 0 );
  for(unsigned i=0; i<A.rw; ++i) for(unsigned j=0; j<B.cl; ++j) for (unsigned k=0; k<A.cl; ++k) C(i,j)+=A(i,k)*B(k,j);
}

template <typename T> void mult( const Matrix<T>& A, const std::vector<T>& B, std::vector<T>& C) {
  plumed_assert( A.cl==B.size() );
  if( C.size()!=A.rw  ) { C.resize(A.rw); }
  for(unsigned i=0; i<A.rw; ++i) { C[i]= static_cast<T>( 0 ); }
  for(unsigned i=0; i<A.rw; ++i) for(unsigned k=0; k<A.cl; ++k) C[i]+=A(i,k)*B[k] ;
}

template <typename T> void mult( const std::vector<T>& A, const Matrix<T>& B, std::vector<T>& C) {
  plumed_assert( B.rw==A.size() );
  if( C.size()!=B.cl ) {C.resize( B.cl );}
  for(unsigned i=0; i<B.cl; ++i) { C[i]=static_cast<T>( 0 ); }
  for(unsigned i=0; i<B.cl; ++i) for(unsigned k=0; k<B.rw; ++k) C[i]+=A[k]*B(k,i);
}

template <typename T> void transpose( const Matrix<T>& A, Matrix<T>& AT ) {
  if( A.rw!=AT.cl || A.cl!=AT.rw ) AT.resize( A.cl, A.rw );
  for(unsigned i=0; i<A.cl; ++i) for(unsigned j=0; j<A.rw; ++j) AT(i,j)=A(j,i);
}

template <typename T> Log& operator<<(Log& ostr, const Matrix<T>& mat) {
  for(unsigned i=0; i<mat.sz; ++i) ostr<<mat.data[i]<<" ";
  return ostr;
}

template <typename T> void matrixOut( Log& ostr, const Matrix<T>& mat) {
  for(unsigned i=0; i<mat.rw; ++i) {
    for(unsigned j=0; j<mat.cl; ++j) { ostr<<mat(i,j)<<" "; }
    ostr<<"\n";
  }
  return;
}

template <typename T> int diagMat( const Matrix<T>& A, std::vector<double>& eigenvals, Matrix<double>& eigenvecs ) {

  // Check matrix is square and symmetric
  plumed_assert( A.rw==A.cl ); plumed_assert( A.isSymmetric()==1 );
  std::vector<double> da(A.sz);
  unsigned k=0;
  std::vector<double> evals(A.cl);
  // Transfer the matrix to the local array
  for (unsigned i=0; i<A.cl; ++i) for (unsigned j=0; j<A.rw; ++j) da[k++]=static_cast<double>( A(j,i) );

  int n=A.cl; int lwork=-1, liwork=-1, m, info, one=1;
  std::vector<double> work(A.cl);
  std::vector<int> iwork(A.cl);
  double vl, vu, abstol=0.0;
  std::vector<int> isup(2*A.cl);
  std::vector<double> evecs(A.sz);

  plumed_lapack_dsyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                       &abstol, &m, evals.data(), evecs.data(), &n,
                       isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info!=0) return info;

  // Retrieve correct sizes for work and iwork then reallocate
  liwork=iwork[0]; iwork.resize(liwork);
  lwork=static_cast<int>( work[0] ); work.resize(lwork);

  plumed_lapack_dsyevr("V", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                       &abstol, &m, evals.data(), evecs.data(), &n,
                       isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info!=0) return info;

  if( eigenvals.size()!=A.cl ) { eigenvals.resize( A.cl ); }
  if( eigenvecs.rw!=A.rw || eigenvecs.cl!=A.cl ) { eigenvecs.resize( A.rw, A.cl ); }
  k=0;
  for(unsigned i=0; i<A.cl; ++i) {
    eigenvals[i]=evals[i];
    // N.B. For ease of producing projectors we store the eigenvectors
    // ROW-WISE in the eigenvectors matrix.  The first index is the
    // eigenvector number and the second the component
    for(unsigned j=0; j<A.rw; ++j) { eigenvecs(i,j)=evecs[k++]; }
  }

  // This changes eigenvectors so that the first non-null element
  // of each of them is positive
  // We can do it because the phase is arbitrary, and helps making
  // the result reproducible
  for(int i=0; i<n; ++i) {
    int j;
    for(j=0; j<n; j++) if(eigenvecs(i,j)*eigenvecs(i,j)>1e-14) break;
    if(j<n) if(eigenvecs(i,j)<0.0) for(j=0; j<n; j++) eigenvecs(i,j)*=-1;
  }
  return 0;
}

template <typename T> int pseudoInvert( const Matrix<T>& A, Matrix<double>& pseudoinverse ) {
  std::vector<double> da(A.sz);
  unsigned k=0;
  // Transfer the matrix to the local array
  for (unsigned i=0; i<A.cl; ++i) for (unsigned j=0; j<A.rw; ++j) da[k++]=static_cast<double>( A(j,i) );

  int nsv, info, nrows=A.rw, ncols=A.cl;
  if(A.rw>A.cl) {nsv=A.cl;} else {nsv=A.rw;}

  // Create some containers for stuff from single value decomposition
  std::vector<double> S(nsv);
  std::vector<double> U(nrows*nrows);
  std::vector<double> VT(ncols*ncols);
  std::vector<int> iwork(8*nsv);

  // This optimizes the size of the work array used in lapack singular value decomposition
  int lwork=-1;
  std::vector<double> work(1);
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  if(info!=0) return info;

  // Retrieve correct sizes for work and rellocate
  lwork=(int) work[0]; work.resize(lwork);

  // This does the singular value decomposition
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  if(info!=0) return info;

  // Compute the tolerance on the singular values ( machine epsilon * number of singular values * maximum singular value )
  double tol; tol=S[0]; for(int i=1; i<nsv; ++i) { if( S[i]>tol ) { tol=S[i]; } } tol*=nsv*epsilon;

  // Get the inverses of the singlular values
  Matrix<double> Si( ncols, nrows ); Si=0.0;
  for(int i=0; i<nsv; ++i) { if( S[i]>tol ) { Si(i,i)=1./S[i]; } else { Si(i,i)=0.0; } }

  // Now extract matrices for pseudoinverse
  Matrix<double> V( ncols, ncols ), UT( nrows, nrows ), tmp( ncols, nrows );
  k=0; for(int i=0; i<nrows; ++i) { for(int j=0; j<nrows; ++j) { UT(i,j)=U[k++]; } }
  k=0; for(int i=0; i<ncols; ++i) { for(int j=0; j<ncols; ++j) { V(i,j)=VT[k++]; } }

  // And do matrix algebra to construct the pseudoinverse
  if( pseudoinverse.rw!=ncols || pseudoinverse.cl!=nrows ) pseudoinverse.resize( ncols, nrows );
  mult( V, Si, tmp ); mult( tmp, UT, pseudoinverse );

  return 0;
}

template <typename T> int Invert( const Matrix<T>& A, Matrix<double>& inverse ) {

  if( A.isSymmetric()==1 ) {
    // GAT -- I only ever use symmetric matrices so I can invert them like this.
    // I choose to do this as I have had problems with the more general way of doing this that
    // is implemented below.
    std::vector<double> eval(A.rw); Matrix<double> evec(A.rw,A.cl), tevec(A.rw,A.cl);
    int err; err=diagMat( A, eval, evec );
    if(err!=0) return err;
    for (unsigned i=0; i<A.rw; ++i) for (unsigned j=0; j<A.cl; ++j) tevec(i,j)=evec(j,i)/eval[j];
    mult(tevec,evec,inverse);
  } else {
    std::vector<double> da(A.sz);
    std::vector<int> ipiv(A.cl);
    unsigned k=0; int n=A.rw, info;
    for(unsigned i=0; i<A.cl; ++i) for(unsigned j=0; j<A.rw; ++j) da[k++]=static_cast<double>( A(j,i) );

    plumed_lapack_dgetrf(&n,&n,da.data(),&n,ipiv.data(),&info);
    if(info!=0) return info;

    int lwork=-1;
    std::vector<double> work(A.cl);
    plumed_lapack_dgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    if(info!=0) return info;

    lwork=static_cast<int>( work[0] ); work.resize(lwork);
    plumed_lapack_dgetri(&n,da.data(),&n,ipiv.data(),work.data(),&lwork,&info);
    if(info!=0) return info;

    if( inverse.cl!=A.cl || inverse.rw!=A.rw ) { inverse.resize(A.rw,A.cl); }
    k=0; for(unsigned i=0; i<A.rw; ++i) for(unsigned j=0; j<A.cl; ++j) inverse(j,i)=da[k++];
  }

  return 0;
}

template <typename T> void cholesky( const Matrix<T>& A, Matrix<T>& B ) {

  plumed_assert( A.rw==A.cl && A.isSymmetric() );
  Matrix<T> L(A.rw,A.cl); L=0.;
  std::vector<T> D(A.rw,0.);
  for(unsigned i=0; i<A.rw; ++i) {
    L(i,i)=static_cast<T>( 1 );
    for (unsigned j=0; j<i; ++j) {
      L(i,j)=A(i,j);
      for (unsigned k=0; k<j; ++k) L(i,j)-=L(i,k)*L(j,k)*D[k];
      if (D[j]!=0.) L(i,j)/=D[j]; else L(i,j)=static_cast<T>( 0 );
    }
    D[i]=A(i,i);
    for (unsigned k=0; k<i; ++k) D[i]-=L(i,k)*L(i,k)*D[k];
  }

  for(unsigned i=0; i<A.rw; ++i) D[i]=(D[i]>0.?sqrt(D[i]):0.);
  if( B.rw!=A.rw || B.cl!=A.cl ) { B.resize( A.rw, A.cl); }
  B=0.; for(unsigned i=0; i<A.rw; ++i) for(unsigned j=0; j<=i; ++j) B(i,j)+=L(i,j)*D[j];
}

template <typename T> void chol_elsolve( const Matrix<T>& M, const std::vector<T>& b, std::vector<T>& y ) {

  plumed_assert( M.rw==M.cl && M(0,1)==0.0 && b.size()==M.rw );
  if( y.size()!=M.rw ) { y.resize( M.rw ); }
  for(unsigned i=0; i<M.rw; ++i) {
    y[i]=b[i];
    for(unsigned j=0; j<i; ++j) y[i]-=M(i,j)*y[j];
    y[i]*=1.0/M(i,i);
  }
}

template <typename T> int logdet( const Matrix<T>& M, double& ldet ) {
  // Check matrix is square and symmetric
  plumed_assert( M.rw==M.cl || M.isSymmetric() );

  std::vector<double> da(M.sz);
  unsigned k=0;
  std::vector<double> evals(M.cl);
  // Transfer the matrix to the local array
  for (unsigned i=0; i<M.rw; ++i) for (unsigned j=0; j<M.cl; ++j) da[k++]=static_cast<double>( M(j,i) );

  int n=M.cl; int lwork=-1, liwork=-1, info, m, one=1;
  std::vector<double> work(M.rw);
  std::vector<int> iwork(M.rw);
  double vl, vu, abstol=0.0;
  std::vector<int> isup(2*M.rw);
  std::vector<double> evecs(M.sz);
  plumed_lapack_dsyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                       &abstol, &m, evals.data(), evecs.data(), &n,
                       isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info!=0) return info;

  // Retrieve correct sizes for work and iwork then reallocate
  lwork=static_cast<int>( work[0] ); work.resize(lwork);
  liwork=iwork[0]; iwork.resize(liwork);

  plumed_lapack_dsyevr("N", "I", "U", &n, da.data(), &n, &vl, &vu, &one, &n,
                       &abstol, &m, evals.data(), evecs.data(), &n,
                       isup.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
  if (info!=0) return info;

  // Transfer the eigenvalues and eigenvectors to the output
  ldet=0; for(unsigned i=0; i<M.cl; i++) { ldet+=log(evals[i]); }

  return 0;
}



}
#endif
