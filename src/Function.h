#ifndef __PLUMED_Function_h
#define __PLUMED_Function_h

#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD{

//+DEVELDOC INHERIT Function
/**
Inherit from here if you are implementing a function of a set of CVs
*/
//+ENDDEVELDOC

class Keywords;

/// Action representing a function of other actions
class Function:
  public ActionWithValue,
  public ActionWithArguments
{
protected:
  void setDerivatives(int,double);
  void setDerivatives(Value*,int,double);
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name ); 
public:
  Function(const ActionOptions&);
  virtual ~Function(){};
  void apply();
  static void registerKeywords(Keywords&);
};

inline
void Function::setDerivatives(Value*v,int i,double d){
  v->addDerivative(i,d);
}

inline
void Function::setDerivatives(int i,double d){
  setDerivatives(getPntrToValue(),i,d);
}

}

#endif

