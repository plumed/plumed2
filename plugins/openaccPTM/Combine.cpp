/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "plumed/function/Combine.h"
#include "plumed/function/FunctionOfVector.h"
#include "plumed/function/FunctionOfMatrix.h"
#include "plumed/core/ActionRegister.h"

#include "ACCParallelTaskManager.h"

namespace PLMD {
namespace function {

typedef FunctionOfVector<Combine,PLMD::ACCPTM> VectorCombine;
PLUMED_REGISTER_ACTION(VectorCombine,"COMBINE_VECTORACC")
typedef FunctionOfMatrix<Combine,PLMD::ACCPTM> MatrixCombine;
PLUMED_REGISTER_ACTION(MatrixCombine,"COMBINE_MATRIXACC")

} // namespace function
} // namespace PLMD
