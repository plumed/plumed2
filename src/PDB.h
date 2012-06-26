#ifndef __PLUMED_PDB_h
#define __PLUMED_PDB_h

#include "AtomNumber.h"
#include "Vector.h"
#include <vector>
#include <string>

namespace PLMD{

/// Minimalistic pdb parser.
/// Contain positions, atomic indexes, occupancy and beta.
/// We should also add other info (e.g. residue name etc).
class PDB{
  std::vector<Vector> positions;
  std::vector<double> occupancy;
  std::vector<double> beta;
  std::vector<AtomNumber> numbers;
public:
/// Read the pdb from a file, scaling positions by a factor scale
  void read(const std::string&file,bool naturalUnits,double scale);
/// Access to the position array
  const std::vector<Vector>     & getPositions()const;
/// Access to the occupancy array
  const std::vector<double>     & getOccupancy()const;
/// Access to the beta array
  const std::vector<double>     & getBeta()const;
/// Access to the indexes
  const std::vector<AtomNumber> & getAtomNumbers()const;
/// Returns the number of atoms
  unsigned                        size()const;
};

}
#endif
