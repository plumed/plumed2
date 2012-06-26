#ifndef __PLUMED_TargetDist_h
#define __PLUMED_TargetDist_h

#include "Value.h"
#include "ActionWithValue.h"
#include "PDB.h"
#include <vector>
#include <string>

namespace PLMD{

class Log;
class PDB;

class TargetDist {
private:
  std::vector<Value*> args; 
  std::vector<double> target;
  Log &log;
public:
  TargetDist(Log& log) : log(log) {};
  void read( const PDB& pdb, std::vector<Value*> args ); 
  void read( const std::vector<double>& targ, std::vector<Value*> ar );
  double calculate( std::vector<double>& derivs );
};

}

#endif
