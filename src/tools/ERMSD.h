/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#ifndef __PLUMED_tools_ERMSD_h
#define __PLUMED_tools_ERMSD_h

#include "Tensor.h"
#include "Vector.h"
#include <vector>
#include <limits>
#include <string>
#include <map>
#include <utility>

namespace PLMD {

class PDB;
class Pbc;

/// A class that implements ERMSD calculations
class ERMSD {
  //std::map< std::pair <unsigned,unsigned> , double> targets;
  //unsigned natoms;
  std::vector<Vector4d> reference_mat;
  unsigned natoms;
  unsigned nresidues;
  std::vector<std::pair <unsigned,unsigned> > pairs;
  double cutoff;

public:
/// Constructor
  ERMSD(): natoms(0),nresidues(0), cutoff(0.0) {}

/// clear the structure
  void clear();

  bool inPair(unsigned i, unsigned j);

/// set reference coordinates
  void setReference(const std::vector<Vector> & reference, const std::vector<unsigned> & pairs_vec,double mycutoff=0.24);

  void calcMat(const std::vector<Vector> & positions, const Pbc& pbc,std::vector<Vector4d> &mat,std::vector<TensorGeneric<4,3> > & Gderivatives);

/// Compute ermsd ( no pbc )
//  double calculate(const std::vector<Vector> & positions,
//                   std::vector<Vector> &derivatives, Tensor& virial) const ;
/// Compute ermsd ( with pbc )
  double calculate(const std::vector<Vector>& positions, const Pbc& pbc,std::vector<Vector> &derivatives, Tensor& virial);
};

}

#endif

