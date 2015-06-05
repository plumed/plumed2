/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Vessel.h"
#include "core/Value.h"

namespace PLMD {
namespace vesselbase{

/**
\ingroup TOOLBOX
Objects that inherit from FunctionVessel can be used (in tandem with PLMD::vesselbase::ActionWithVessel) to calculate
functions of the form \f$\prod_k H_k[ \sum_j \prod_i g_i(x) ]\f$.  They should take in a series of values
and return one single value.   
*/

class FunctionVessel : public Vessel {
private:
/// This is the pointer to the value we are creating
  Value* final_value;
protected:
/// The number of derivatives
  unsigned nderivatives;
/// Are the derivatives differentiable
  bool diffweight;
/// Are we normalising by the weight
  bool norm;
/// Are we using the tolerance
  bool usetol;
/// Set the final value
  void setOutputValue( const double& val );
/// Resize the vector containing the derivatives
  void setNumberOfDerivatives( const unsigned& nder );
/// Return a pointer to the final value
  void addDerivativeToFinalValue( const unsigned& j, const double& der  );
public:
  static void registerKeywords( Keywords& keys );
  FunctionVessel( const VesselOptions& );
/// This does the resizing of the buffer
  virtual void resize();
/// This applies all the forces
  bool applyForce( std::vector<double>& forces );
/// The description for the log
  std::string description();
/// The rest of the description of what we are calculating
  virtual std::string function_description()=0;
/// Do the calcualtion
  virtual bool calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const ;
/// Do any transformations of the value that are required
  virtual double calcTransform( const double& val, double& df ) const ;
/// Finish the calculation of the quantity
  virtual void finish( const std::vector<double>& buffer );
/// Finish with any transforms required
  virtual double finalTransform( const double& val, double& dv );
};

inline
void FunctionVessel::setOutputValue( const double& val ){
  final_value->set( val );
}

inline
void FunctionVessel::addDerivativeToFinalValue( const unsigned& j, const double& der ){
  plumed_dbg_assert( j<nderivatives );
  final_value->addDerivative( j, der );
}

}
}
#endif
