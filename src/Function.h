#ifndef __PLUMED_Function_h
#define __PLUMED_Function_h

//#include "ActionWithValue.h"
#include "ActionWithArguments.h"

namespace PLMD{

/// Action representing a function of other actions
class Function : public ActionWithArguments {
private:
  std::vector<double> arguments, derivatives;
protected:
  void readFunction();
public:
  Function(const ActionOptions&);
  virtual ~Function(){};
  void calculate();
  void apply();
  virtual double compute( const std::vector<double>& arguments, std::vector<double> derivatives )=0;
};

}

#endif

