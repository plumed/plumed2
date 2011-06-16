#ifndef __PLUMED_MDAtoms_h
#define __PLUMED_MDAtoms_h

#include "Tensor.h"
#include "Vector.h"
#include <vector>

namespace PLMD {

class MDAtomsBase
{
protected:
  double scalep,scalef;
  double scaleb,scalev;
  int stride;
public:
  MDAtomsBase();
  virtual ~MDAtomsBase(){};
  virtual void setm(void*m)=0;
  virtual void setc(void*m)=0;
  virtual void setBox(void*)=0;
  virtual void setp(void*p)=0;
  virtual void setVirial(void*)=0;
  virtual void setf(void*f)=0;
  virtual void setp(void*p,int i)=0;
  virtual void setf(void*f,int i)=0;
          void setUnits(double,double);
  virtual void MD2double(const void*,double&)const=0;
  virtual void double2MD(const double&,void*)const=0;
  virtual void getBox(Tensor &)const=0;
  virtual void getPositions(const std::vector<int>&index,std::vector<Vector>&)const=0;
  virtual void getMasses(const std::vector<int>&index,std::vector<double>&)const=0;
  virtual void getCharges(const std::vector<int>&index,std::vector<double>&)const=0;
  virtual void updateVirial(const Tensor&)const=0;
  virtual void updateForces(const std::vector<int>&index,const std::vector<Vector>&)=0;
  virtual void rescaleForces(int first,int last,double factor)=0;
  virtual int  getRealPrecision()const=0;
  static MDAtomsBase* create(int);
};

}


#endif
