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
#include "plumed/core/ActionRegister.h"
#include "plumed/gridtools/KDE.h"
#include "plumed/gridtools/KDEGridTools.h"
#include "plumed/gridtools/SphericalKDEGridTools.h"

#include "ACCParallelTaskManager.h"

namespace PLMD {
namespace gridtools {
typedef KDE<DiagonalKernelParams,DiscreteKernel,KDEGridTools<DiagonalKernelParams,DiscreteKernel>,PLMD::ACCPTM> discretekde;
PLUMED_REGISTER_ACTION(discretekde,"KDE_DISCRETEACC")
// typedef KDE<DiagonalKernelParams,HistogramBeadKernel,KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>,PLMD::ACCPTM> beadkde;
// PLUMED_REGISTER_ACTION(beadkde,"KDE_BEADS")
typedef KDE<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>,PLMD::ACCPTM> flatkde;
PLUMED_REGISTER_ACTION(flatkde,"KDE_KERNELSACC")
typedef KDE<VonMissesKernelParams,UniversalVonMisses,SphericalKDEGridTools,PLMD::ACCPTM> sphericalkde;
PLUMED_REGISTER_ACTION(sphericalkde,"SPHERICAL_KDEACC")
typedef KDE<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>,PLMD::ACCPTM> flatfkde;
PLUMED_REGISTER_ACTION(flatfkde,"KDE_FULLCOVARACC")
} // namespace colvar
} // namespace PLMD
