#ifndef __PLUMED_Value_h
#define __PLUMED_Value_h

#include <vector>
#include <string>
#include <cassert>


namespace PLMD{

class ActionWithValue;

/// Class containing a value which can be addressed by PLMD::ActionWithArguments.
/// It also contains the derivative of this value with respect to
/// an arbitrary number of parameters.
/// Typically, an object of type PLMD::ActionWithValue will contain one or
/// more objects of type PLUMD::Value, one per component.
class Value{
  ActionWithValue&action;
  double value;
  double inputForce;
  std::vector<double> derivatives;
  std::string name;
  bool deriv;
public:
  void set(double v){value=v;}
  void setValue(double v){value=v;}
  double getValue()const{return value;}
  const std::string& getName()const{return name;}
  const std::string getFullName()const;
  Value(ActionWithValue&action,const std::string& name):action(action),value(0.0),name(name),deriv(false){};
  void enableDerivatives();
  bool hasDerivatives(){return deriv;};
  void setNumberOfParameters(int n){if(deriv)derivatives.resize(n);}
  void setDerivatives(int i,double d){derivatives[i]=d;}
  void clearInputForce(){inputForce=0.0;}
  void clearDerivatives(){derivatives.assign(derivatives.size(),0.0);}
  double  getForce()const{return inputForce;}
  void  addForce(double f){assert(hasDerivatives());inputForce+=f;}
  const std::vector<double> &  getDerivatives()const{return derivatives;}
  ActionWithValue& getAction(){return action;};
};

}

#endif

