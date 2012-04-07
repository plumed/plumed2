#ifndef __PLUMED_ActionWithField_h
#define __PLUMED_ActionWithField_h

#include "Action.h"
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "PlumedException.h"
#include "ActionWithDistribution.h"
#include "Grid.h"
#include "Field.h"

namespace PLMD{

class ActionWithField : 
  public ActionWithValue,
  public ActionPilot
  {
private:
  bool serial;
  Action* apply_action;
  Field* myfield;
  Grid* bias;
  std::vector<double> buffer;
  std::vector<unsigned> blocks;
  std::vector<double> derivatives;
protected:
  Grid* getPntrToBias();
  void clearBias();
  void addFieldToBias( const double& hh );
public:
  static void registerKeywords(Keywords& keys);
  ActionWithField(const ActionOptions&ao);
  ~ActionWithField();
  void calculate();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void apply();
};

inline
Grid* ActionWithField::getPntrToBias(){
  return bias;
}

}

#endif
