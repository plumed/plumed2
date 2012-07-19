/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
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

