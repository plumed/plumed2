#ifndef __PLUMED_RMSD_h
#define __PLUMED_RMSD_h

#include "Tensor.h"
#include "Vector.h"
#include <vector>

namespace PLMD{

class PDB;

/// A class that implements RMSD calculations
class RMSD
{
  std::vector<Vector> reference;
  std::vector<double> align;
  std::vector<double> displace;
public:
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void setFromPDB(const PDB&);
/// set reference coordinates
  void setReference(const std::vector<Vector> & reference);
/// set weights
  void setAlign(const std::vector<double> & align);
/// set align
  void setDisplace(const std::vector<double> & displace);
/// Compute rmsd
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, Tensor& virial ) const;
};

}

#endif

