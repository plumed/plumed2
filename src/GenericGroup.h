#ifndef __PLUMED_GenericGroup_h
#define __PLUMED_GenericGroup_h

#include "ActionAtomistic.h"
#include "PlumedMain.h"

namespace PLMD{

class GenericGroup : public ActionAtomistic {

public:
  GenericGroup(const ActionOptions&ao);
  ~GenericGroup();
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){ assert(false); }
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
  virtual void updateAtomSelection( std::vector<bool>& skips ){}
  void updateNeighbourList( const double& cutoff, std::vector<bool>& skips ){ assert(false); }
  void calculate(){}
  void apply(){}
};

}

#endif
