#ifndef __PLUMED_MultiColvarSecondaryStructureRMSD_h
#define __PLUMED_MultiColvarSecondaryStructureRMSD_h

#include "MultiColvar.h"
#include "DRMSD.h"
#include "RMSD.h"
#include <vector>

namespace PLMD {

/// Base action for calculating things like AlphRMSD, AntibetaRMSD, etc

class MultiColvarSecondaryStructureRMSD : public MultiColvar {
private:
  RMSD* secondary_rmsd;
  DRMSD* secondary_drmsd;
protected:
  void setSecondaryStructure( std::vector<Vector>& structure, double bondlength, double units );
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarSecondaryStructureRMSD(const ActionOptions&);
  virtual ~MultiColvarSecondaryStructureRMSD();
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(const unsigned nn){ return false; }
};

}

#endif
