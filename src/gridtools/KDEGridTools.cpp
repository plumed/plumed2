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
#include "KDEGridTools.h"

namespace PLMD {
namespace gridtools {

template <>
void KDEGridTools<DiagonalKernelParams,DiscreteKernel>::readBandwidthAndHeight( const DiscreteKernel& params, ActionWithArguments* action ) {
  std::size_t nargs = action->getNumberOfArguments();
  KDEGridTools<DiagonalKernelParams,DiscreteKernel>::readHeightKeyword( false, nargs, std::vector<std::string>(), action );
}

template <>
void KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>::readBandwidthAndHeight( const RegularKernel<DiagonalKernelParams>& params, ActionWithArguments* action ) {
  std::vector<Value*> bwargs;
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidth( nargs, action, bw );
  std::string volstr;
  action->parse("VOLUMES",volstr);
  if( volstr.length()>0 ) {
    if( !params.canusevol ) {
      action->error("cannot use normalized kernels with selected kernel type");
    }
    // Check if we are using Gaussian kernels
    KDEHelper<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>>::readKernelParameters( volstr, action, "_volumes", false );
    convertHeightsToVolumes(nargs, bw, volstr, action);
  } else {
    KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>::readHeightKeyword( params.canusevol, nargs, bw, action );
  }
}

template <>
void KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>::readBandwidthAndHeight( const HistogramBeadKernel& params, ActionWithArguments* action ) {
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidth( nargs, action, bw );
  KDEGridTools<DiagonalKernelParams,HistogramBeadKernel>::readHeightKeyword( false, nargs, bw, action );
}

template <>
void KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>::readBandwidthAndHeight( const RegularKernel<NonDiagonalKernelParams>& params, ActionWithArguments* action ) {
  std::vector<Value*> bwargs;
  std::vector<std::string> bw;
  std::size_t nargs = action->getNumberOfArguments();
  readBandwidthKeyword( nargs, action, bw, bwargs );
  if( nargs>1 && bw.size()==1 ) {
    if( bwargs[0]->getRank()!=2 || bwargs[0]->getShape()[0]!=nargs || bwargs[0]->getShape()[1]!=nargs  ) {
      action->error("invalid input for bandwidth parameter");
    }
    std::string str_i, str_j;
    bw.resize( nargs*nargs );
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      for(unsigned j=0; j<nargs; ++j) {
        Tools::convert( j+1, str_j );
        action->plumed.readInputWords( Tools::getWords(action->getLabel() + "_scalar_bw" + str_i + "_" + str_j + ": SELECT_COMPONENTS ARG=" + bwargs[0]->getName() + " COMPONENTS=" + str_i + "." + str_j ), false );
        bw[i*nargs+j] = action->getLabel() + "_bw" + str_i + "_" + str_j;
        action->plumed.readInputWords( Tools::getWords( bw[i*nargs+j] + ": CUSTOM ARG=" + action->getLabel() + "_bwones," + action->getLabel() + "_scalar_bw" + str_i + "_" + str_j + " FUNC=x*y PERIODIC=NO"), false );
      }
    }
  } else if( bw.size()==nargs*nargs ) {
    std::string str_i, str_j;
    for(unsigned i=0; i<nargs; ++i) {
      Tools::convert( i+1, str_i );
      for(unsigned j=0; j<nargs; ++j) {
        Tools::convert( j+1, str_j );
        KDEHelper<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>>::readKernelParameters( bw[i*nargs+j], action, "_bw" + str_i + "_" + str_j, true );
      }
    }
  } else {
    action->error("wrong number of arguments specified in input to bandwidth parameter");
  }
  std::string volstr;
  action->parse("VOLUMES",volstr);
  if( volstr.length()>0 ) {
    if( !params.canusevol ) {
      action->error("cannot use normalized kernels with selected kernel type");
    }
    // Check if we are using Gaussian kernels
    KDEHelper<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>,KDEGridTools<DiagonalKernelParams,RegularKernel<DiagonalKernelParams>>>::readKernelParameters( volstr, action, "_volumes", false );
    convertHeightsToVolumes(nargs, bw, volstr, action);
  } else {
    KDEGridTools<NonDiagonalKernelParams,RegularKernel<NonDiagonalKernelParams>>::readHeightKeyword( params.canusevol, nargs, bw, action );
  }
}

}
}
