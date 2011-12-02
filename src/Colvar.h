#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include <string>
#include <cassert>
#include <vector>
#include "ActionAtomistic.h"
#include "ActionWithValue.h"

namespace PLMD {

class Colvar : public ActionAtomistic {
private:
/// The forces on the atoms and on the virial
  std::vector<double> forces;
/// The forces on the atoms
  std::vector<Vector> f;
protected:
  void readActionColvar();
public:
  Colvar(const ActionOptions&);
  ~Colvar(){};
  void apply();
};

} 

#endif

