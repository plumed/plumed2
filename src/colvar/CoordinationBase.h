/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_colvar_CoordinationBase_h
#define __PLUMED_colvar_CoordinationBase_h
#include "Colvar.h"
#include <memory>
#include <vector>
#include "tools/Tensor.h"

namespace PLMD {

class NeighborList;

namespace colvar {

class CoordinationBase : public Colvar {
  bool pbc;
  bool serial;
/// The neighbor list, only with NLIST, NLISTCELLS or PAIR (otherwise the pairs come from the groups)
  std::unique_ptr<NeighborList> nl;
  bool invalidateList;
  bool firsttime;
/// Without a neighbor list: the pairs are the rows times the columns
/// (all the pairs i<j in GROUPA when GROUPB is empty).
/// The rows are the larger group, split in blocks of rowsPerBlock rows.
  unsigned nrows=0;
  unsigned ncols=0;
/// Index in the list of atoms of the first row and of the first column
  unsigned rowBase=0;
  unsigned colBase=0;
/// The atom of the row is the first of the pair (the one in GROUPA)
  bool rowFirst=true;
/// The pairs are i<j in a single group
  bool singleGroup=false;
  unsigned rowsPerBlock=1;
/// A pair can be made of the same atom twice (the groups share some atoms)
  bool samePairsPossible=false;
/// Contributions of each block to the derivatives of the columns, and the value and the virial of each block
  std::vector<Vector> blockColumns;
  std::vector<double> blockValue;
  std::vector<Tensor> blockVirial;
/// Buffers of the threads for the neighbor list
  std::vector<std::vector<Vector>> threadDeriv;
  std::vector<double> threadValue;
  std::vector<Tensor> threadVirial;
  std::vector<char> threadUsed;
/// Derivatives summed over the MPI processes
  std::vector<Vector> mpiDeriv;

  template<bool usePbc>
  void calculateBlocks();
  template<bool usePbc>
  void calculateList();
/// Add d to the derivatives of atom i: to the value, or to mpiDeriv when the processes have to sum them
  void addAtomDerivatives(bool toBuffer,unsigned i,const Vector& d) {
    if(toBuffer) {
      mpiDeriv[i]+=d;
    } else {
      setAtomsDerivatives(i,d);
    }
  }
/// Sum the derivatives of mpiDeriv over the processes and add them to the value
  void sumAndSetDerivatives();

public:
  explicit CoordinationBase(const ActionOptions&);
  ~CoordinationBase();
// active methods:
  void calculate() override;
  void prepare() override;
  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const=0;
/// Square of the distance beyond which pairing() is exactly zero, with zero derivative:
/// the pairs farther than this are skipped. By default all the pairs are computed.
  virtual double pairCutoff2() const;
  static void registerKeywords( Keywords& keys );
};

}
}
#endif
