/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_vesselbase_FunctionVessel_h
#define __PLUMED_vesselbase_FunctionVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "ValueVessel.h"
#include "core/Value.h"

namespace PLMD {
namespace vesselbase {

/**
\ingroup TOOLBOX
Objects that inherit from FunctionVessel can be used (in tandem with PLMD::vesselbase::ActionWithVessel) to calculate
functions of the form \f$\prod_k H_k[ \sum_j \prod_i g_i(x) ]\f$.  They should take in a series of values
and return one single value.
*/

class FunctionVessel : public ValueVessel {
protected:
/// Are the derivatives differentiable
  bool diffweight;
/// Are we normalising by the weight
  bool norm;
/// Are we using the tolerance
  bool usetol;
public:
  static void registerKeywords( Keywords& keys );
  explicit FunctionVessel( const VesselOptions& );
/// This does the resizing of the buffer
  void resize() override;
/// Do the calcualtion
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const override;
/// Do any transformations of the value that are required
  virtual double calcTransform( const double& val, double& df ) const ;
/// Finish the calculation of the quantity
  void finish( const std::vector<double>& buffer ) override;
/// Finish with any transforms required
  virtual double finalTransform( const double& val, double& dv );
};

}
}
#endif
