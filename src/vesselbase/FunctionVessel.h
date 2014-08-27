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
/// The number of derivatives
  unsigned nderivatives;
/// This is the pointer to the value we are creating
  Value* final_value;
protected:
/// Are the derivatives differentiable
  bool diffweight;
/// Add some value to the accumulator if it is greater than tolerance
  bool addValueUsingTolerance( const unsigned& jval, const double& val );
/// Add some value to the accumulator and ignore the tolerance
  void addValueIgnoringTolerance( const unsigned& jval, const double& val );
/// Set the final value
  void setOutputValue( const double& val );
/// Get the nth value in the distribution
  double getFinalValue( const unsigned& j );
/// This does a combination of the product and chain rules
  void mergeFinalDerivatives( const std::vector<double>& df );
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
};

inline
bool FunctionVessel::addValueUsingTolerance( const unsigned& jval, const double& val ){
  if( fabs(val)<getTolerance() ) return false; 
  addToBufferElement( (nderivatives+1)*jval, val );
  return true;
}

inline
void FunctionVessel::addValueIgnoringTolerance( const unsigned& jval, const double& val ){
  addToBufferElement( (nderivatives+1)*jval, val );
}

inline
double FunctionVessel::getFinalValue(const unsigned& j){
  return getBufferElement( (nderivatives+1)*j );
}

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
