#ifndef __PLUMED_MatrixSquareBracketsAccess_h
#define __PLUMED_MatrixSquareBracketsAccess_h

namespace PLMD{

/**
Utility class to add [][] access

It implements the trick described in C++ FAQ 13.12 to allow [][] access
to matrix-like classes, based on the (,) syntax, thus translating
[i][j] into (i,j).  In practice, one only needs to implement the (,) syntax and to inherit from
MatrixSquareBracketsAccess.
The first template parameter (T) should be the
class itself, the second (C) is the type of the returned value,
and the third (I) is the type of the index (unsigned by default).
As everything is inlined, no overhead is expected.

\verbatim

class MyMatrixClass:
  public MatrixSquareBracketsAccess<MyMatrixClass,double>
{
  double data[16];
public:
  double & operator ()(unsigned i,unsigned j){
    return data[4*i+j];
  }
  const double & operator ()(unsigned i,unsigned j)const{
    return data[4*i+j];
  }
};

int main(){
  MyMatrixClass m;
  m[0][1]=3.0;
  return 0;
}

*/

template<class T,class C,class I=unsigned>
class MatrixSquareBracketsAccess{
/// Small utility class which just contains a pointer to the T and the row number
  class Const_row{
    friend class MatrixSquareBracketsAccess; // this so as to allow only T to instantiate Const_row
                    // the user should not manipulate it directly
    const MatrixSquareBracketsAccess& t;
    const I i;
    Const_row(const MatrixSquareBracketsAccess&t,I i); // constructor is private and cannot be manipulated by the user
  public:
  /// access element
    const C & operator[] (I j)const;
  };
/// Small utility class which just contains a pointer to the T and the row number
  class Row{
    friend class MatrixSquareBracketsAccess; // this so as to allow only T to instantiate Const_row
                         // the user should not manipulate it directly
    MatrixSquareBracketsAccess& t;
    const I i;
    Row(MatrixSquareBracketsAccess&t,I i); // constructor is private and cannot be manipulated by the user
  public:
  /// access element
    C & operator[] (I j);
  };
public:
/// access element (with [][] syntax)
  Row operator[] (I i);
/// access element (with [][] syntax)
  Const_row operator[] (I i)const;
};

template<class T,class C,class I>
MatrixSquareBracketsAccess<T,C,I>::Const_row::Const_row(const MatrixSquareBracketsAccess&t,I i):
  t(t),i(i){}

template<class T,class C,class I>
MatrixSquareBracketsAccess<T,C,I>::Row::Row(MatrixSquareBracketsAccess&t,I i):
  t(t),i(i){}

template<class T,class C,class I>
const C & MatrixSquareBracketsAccess<T,C,I>::Const_row::operator[] (I j)const{
  return (*static_cast<T*>(&t))(i,j);
}

template<class T,class C,class I>
C & MatrixSquareBracketsAccess<T,C,I>::Row::operator[] (I j){
  return (*static_cast<T*>(&t))(i,j);
}

template<class T,class C,class I>
typename MatrixSquareBracketsAccess<T,C,I>::Row MatrixSquareBracketsAccess<T,C,I>::operator[] (I i){
  return Row(*this,i);
}

template<class T,class C,class I>
typename MatrixSquareBracketsAccess<T,C,I>::Const_row MatrixSquareBracketsAccess<T,C,I>::operator[] (I i)const{
  return Const_row(*this,i);
}

}


#endif


