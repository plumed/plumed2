#ifndef __PLUMED_DRMSD_h
#define __PLUMED_DRMSD_h

#include "Tensor.h"
#include "Vector.h"
#include "Matrix.h"
#include "Pbc.h"
#include <vector>

namespace PLMD{

class PDB;

/// A class that implements DRMSD calculations
class DRMSD {
  unsigned natoms, npairs;
  Matrix<double> targets;
public:
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void setFromPDB( const double& bondc, const PDB& );
/// set reference coordinates
  void setReference( const double& bondc, const std::vector<Vector> & reference );
/// Compute drmsd ( no pbc )
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, Tensor& virial) const ;
/// Compute drmsd ( with pbc )
  double calculate(const std::vector<Vector>& positions, const Pbc& pbc, std::vector<Vector> &derivatives, Tensor& virial) const ;
};

}

#endif

