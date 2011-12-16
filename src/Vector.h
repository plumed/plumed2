#ifndef __PLUMED_Vector_h
#define __PLUMED_Vector_h

namespace PLMD{

/// 3d vector of double.
/// Useful to simplify syntax. All the methods are inlined for better optimization.
/// Accepts both [] and () syntax for access.
class Vector{
  double d[3];
public:
/// return v2-v1
  friend Vector delta(const Vector&v1,const Vector&v2);
/// return s*v
  friend Vector operator*(double,const Vector&);
/// return v*s
  friend Vector operator*(const Vector&,double);
/// return v1+v2
  friend Vector operator+(const Vector&,const Vector&);
/// return -v;
  friend Vector operator-(const Vector&);
/// return v1 .scalar. v2
  friend double dotProduct(const Vector&,const Vector&);
/// return v1 .vector. v2
  friend Vector crossProduct(const Vector&,const Vector&);
/// array-like access [i]
  double & operator[](int i);
/// array-like access [i]
  const double & operator[](int i)const;
/// parenthesis access (i)
  double & operator()(int i);
/// parenthesis access (i)
  const double & operator()(int i)const;
/// set it to zero
  void clear();
/// create it with components x0, x1 and x2
  Vector(double x0,double x1,double x2);
/// create it null
  Vector();
/// compute the squared modulo
  double modulo2()const;
/// compute the modulo
  double modulo()const;
/// scale the vector by a factor s
  void scale(double s);
/// increment
  Vector& operator +=(const Vector& b);
};

inline
double & Vector::operator[](int i){
return d[i];
}

inline
const double & Vector::operator[](int i)const{
  return d[i];
}

inline
double & Vector::operator()(int i){
  return d[i];
}

inline
const double & Vector::operator()(int i)const{
  return d[i];
}

inline
Vector delta(const Vector&v1,const Vector&v2){
  return Vector(v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]);
}

inline
Vector operator*(double s,const Vector&v){
  return Vector(s*v[0],s*v[1],s*v[2]);
}

inline
Vector operator*(const Vector&v,double s){
  return Vector(s*v[0],s*v[1],s*v[2]);
}

inline
void Vector::clear(){
  d[0]=0.0;
  d[1]=0.0;
  d[2]=0.0;
}

inline
Vector::Vector(double x0,double x1,double x2){
  d[0]=x0;
  d[1]=x1;
  d[2]=x2;
}

inline
Vector::Vector(){
  d[0]=0.0;
  d[1]=0.0;
  d[2]=0.0;
}

inline
double Vector::modulo2()const{
  return d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
}

inline
void Vector::scale(double s){
  d[0]*=s;d[1]*=s;d[2]*=s;
}

inline
Vector operator+(const Vector&v1,const Vector&v2){
  return Vector(v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]);
}

inline
Vector& Vector::operator +=(const Vector& b){
  d[0]+=b(0);
  d[1]+=b(1);
  d[2]+=b(2);
  return *this;
}

inline
Vector operator-(const Vector&v){
  return Vector(-v[0],-v[1],-v[2]);
}

inline
double dotProduct(const Vector& v1,const Vector& v2){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

inline
Vector crossProduct(const Vector& v1,const Vector& v2){
  return Vector(
    v1[1]*v2[2]-v1[2]*v2[1],
    v1[2]*v2[0]-v1[0]*v2[2],
    v1[0]*v2[1]-v1[1]*v2[0]);
}

}


#endif

