#ifndef __PLUMED_FieldBias_h
#define __PLUMED_FieldBias_h

#include "Action.h"
#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "PlumedException.h"
#include "ActionWithDistribution.h"
#include "Grid.h"
#include "Field.h"

namespace PLMD{

class FieldBias : 
  public ActionWithValue,
  public ActionPilot
  {
private:
  bool serial, debug;
  ActionWithValue* apply_action;
  Field* myfield;
  Grid* bias;
  double norm;
  std::vector<Value*> f_arg;
  std::vector<double> buffer;
  std::vector<unsigned> blocks;
  std::vector<double> derivatives;
protected:
  Grid* getPntrToBias();
  double get_normalizer() const ;
  std::vector<double>& get_buffer();
  void clearBias();
public:
  static void registerKeywords(Keywords& keys);
  FieldBias(const ActionOptions&ao);
  ~FieldBias();
  void calculate();
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void apply();
};

inline
Grid* FieldBias::getPntrToBias(){
  return bias;
}

inline
double FieldBias::get_normalizer() const {
  return pow( buffer[0], 1./static_cast<double>(norm) );
}

inline
std::vector<double>& FieldBias::get_buffer(){
  return buffer;
}

}

#endif
