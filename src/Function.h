#ifndef __PLUMED_Function_h
#define __PLUMED_Function_h

#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD{

/**
\ingroup INHERIT
Inherit from here if you are implementing a function of CVs
*/

/// Action representing a function of other actions
class Function:
  public ActionWithValue,
  public ActionWithArguments
{
protected:
  void setDerivative(int,double);
  void setDerivative(Value*,int,double);
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name ); 
public:
  Function(const ActionOptions&);
  virtual ~Function(){};
  void apply();
  static void registerKeywords(Keywords&);
};

inline
void Function::setDerivative(Value*v,int i,double d){
  v->addDerivative(i,d);
}

inline
void Function::setDerivative(int i,double d){
  setDerivative(getPntrToValue(),i,d);
}

}

#endif

