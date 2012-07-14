#ifndef __PLUMED_ColvarRMSD_h
#define __PLUMED_ColvarRMSD_h

#include "Colvar.h"
#include "PlumedMain.h"
#include "ActionRegister.h"
#include "PDB.h"
#include "RMSD.h"
#include "Atoms.h"


using namespace std;

namespace PLMD{
   
class ColvarRMSD : public Colvar {
	
  RMSD rmsd;
	
  bool squared; 

  vector<Vector> derivs;

public:
  ColvarRMSD(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
  friend class FunctionPathMSD;
};

}
#endif
