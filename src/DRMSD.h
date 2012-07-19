#ifndef __PLUMED_DRMSD_h
#define __PLUMED_DRMSD_h

#include "Tensor.h"
#include "Vector.h"
#include "Matrix.h"
#include "Pbc.h"
#include <vector>
#include <limits>
#include <map>

namespace PLMD{

class PDB;

/// A class that implements DRMSD calculations
class DRMSD {
  std::map< std::pair <unsigned,unsigned> , double> targets;
  unsigned natoms;
  public:
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void setFromPDB(const PDB&, double lbound=0.0, double ubound=std::numeric_limits<double>::max( ));
/// set reference coordinates
  void setReference(const std::vector<Vector> & reference, double lbound=0.0, double ubound=std::numeric_limits<double>::max( ));
/// Compute drmsd ( no pbc )
  double calculate(const std::vector<Vector> & positions,
                   std::vector<Vector> &derivatives, Tensor& virial) const ;
/// Compute drmsd ( with pbc )
  double calculate(const std::vector<Vector>& positions, const Pbc& pbc,
                   std::vector<Vector> &derivatives, Tensor& virial, bool do_pbc=true) const ;
};

}

#endif

