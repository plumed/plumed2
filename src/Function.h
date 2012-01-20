#ifndef __PLUMED_Function_h
#define __PLUMED_Function_h

#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD{

class Keywords;

/// Action representing a function of other actions
class Function:
  public ActionWithValue,
  public ActionWithArguments
{
public:
  Function(const ActionOptions&);
  virtual ~Function(){};
  void apply();
  void setDerivatives(int,double);
  void setDerivatives(Value*,int,double);
  static void registerKeywords(Keywords&);
};

inline
void Function::setDerivatives(Value*v,int i,double d){
  v->setDerivatives(i,d);
}

inline
void Function::setDerivatives(int i,double d){
  setDerivatives(getValue(0),i,d);
}

}

#endif

