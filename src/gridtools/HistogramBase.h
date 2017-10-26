/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#ifndef __PLUMED_gridtools_HistogramBase_h
#define __PLUMED_gridtools_HistogramBase_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "GridCoordinatesObject.h"

namespace PLMD {
namespace gridtools {

class HistogramBase : 
public ActionWithValue,
public ActionWithArguments
{
protected:
  double norm;
  unsigned heights_index;
  bool one_kernel_at_a_time, unorm;
  GridCoordinatesObject gridobject;
  void addValueWithDerivatives( const std::vector<unsigned>& shape );
public:
  static void registerKeywords( Keywords& keys );
  explicit HistogramBase(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  virtual void buildSingleKernel( std::vector<unsigned>& tflags, const double& height, std::vector<double>& args ) = 0;
  virtual double calculateValueOfSingleKernel( const std::vector<double>& args, std::vector<double>& der ) const = 0;
  virtual void addKernelToGrid( const double& height, const std::vector<double>& args, const unsigned& bufstart, std::vector<double>& buffer ) const = 0;
  void calculate();
  virtual void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void gatherGridAccumulators( const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void apply();
};

} 
}
#endif
