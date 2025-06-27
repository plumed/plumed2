/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_adjmat_ContactMatrix_h
#define __PLUMED_adjmat_ContactMatrix_h

#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC

#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"
namespace PLMD {
namespace adjmat {

struct ContactMatrix {
#ifdef __PLUMED_USE_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<ContactMatrix>* action );
  static void calculateWeight( const ContactMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

}
}
#endif
