#ifndef __PLUMED_MolInfo_h
#define __PLUMED_MolInfo_h

#include "ActionSetup.h"
#include "PlumedMain.h"
#include "PlumedException.h"
#include "PDB.h"

namespace PLMD {

class MolInfo : public virtual ActionSetup {
private:
  PDB pdb;
  std::vector< std::vector<AtomNumber> > read_backbone;
public:
  static void registerKeywords( Keywords& keys );
  MolInfo(const ActionOptions&ao);
  void getBackbone( std::vector<std::string>& resstrings, const std::vector<std::string>& atnames, std::vector< std::vector<AtomNumber> >& backbone );
};

}

#endif
