#ifndef __PLUMED_ColvarWithModifiers_h
#define __PLUMED_ColvarWithModifiers_h

#include <string>
#include <cassert>
#include <vector>
#include "Colvar.h"

namespace PLMD {

class ColvarModifier;

class ColvarWithModifiers : public Colvar {
friend class ColvarModifier;
private:
  bool readInput;
  int maxatoms;
  std::vector< std::pair<std::string,std::string> > mod_params;
  std::vector<ColvarModifier*> modifiers;
  std::vector<Vector> derivatives; 
protected:
  void readAtomsKeyword( int& natoms );
  void readGroupsKeyword( int& maxgroups );
  void finishColvarSetup( const double& min, const double& max );
public:
  ColvarWithModifiers(const ActionOptions&);
  void calculate();
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial )=0; 
  virtual void interpretGroupsKeyword( const std::vector<unsigned>& boundaries, unsigned& maxatoms );
};
  
}

#endif
