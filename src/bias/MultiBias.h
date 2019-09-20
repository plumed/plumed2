/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_bias_MultiBias_h
#define __PLUMED_bias_MultiBias_h

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace bias {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new simulation biases, within it there is
information as to how to go about implementing a new bias.
*/

class MultiBias :
  public ActionPilot,
  public ActionWithValue,
  public ActionWithArguments
{
private:
  unsigned nderivatives;
  std::vector<double> forcesToApply;
protected:
  void setBias( const double& bb, MultiValue& myvals ) const ;
  void setNonBiasComponent( const unsigned& ival, const double val, MultiValue& myvals ) const ;
  void addBiasDerivative( const unsigned& jder, const double& der, MultiValue& myvals ) const ;
public:
  static void registerKeywords(Keywords&);
  explicit MultiBias(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void calculate();
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  virtual void calculateBias( const std::vector<double>& args, MultiValue& myvals ) const = 0;
  void apply();
};

inline
unsigned MultiBias::getNumberOfDerivatives() const {
  return nderivatives;
}

inline
void MultiBias::setBias( const double& bb, MultiValue& myvals ) const {
  myvals.addValue( getPntrToOutput(0)->getPositionInStream(), bb );
}

inline
void MultiBias::setNonBiasComponent( const unsigned& ival, const double val, MultiValue& myvals ) const {
  plumed_dbg_assert( ival+1<getNumberOfComponents() );
  myvals.addValue( getPntrToOutput(ival+1)->getPositionInStream(), val );
}

inline
void MultiBias::addBiasDerivative( const unsigned& jder, const double& der, MultiValue& myvals ) const {
  plumed_dbg_assert( jder<getNumberOfArguments() );
  if( actionInChain() ) {
    unsigned istrn = getArgumentPositionInStream(jder, myvals);
    unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
    for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
      unsigned kind=myvals.getActiveIndex(istrn,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[jder] + kind, der*myvals.getDerivative( istrn, kind ) );
    }
    return;
  }
  if( arg_ends.size()>0 ) {
    unsigned base=0;
    for(unsigned i=0; i<jder; ++i) {
      for(unsigned j=arg_ends[i]; j<arg_ends[i+1]; ++j) base += getPntrToArgument(j)->getNumberOfValues( getLabel() );
    }

    if( arg_ends[jder+1]==(arg_ends[jder]+1) && getPntrToArgument(arg_ends[jder])->getRank()==0 ) { 
        myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), base, der );
    } else {
       myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), base + myvals.getTaskIndex(), der );
    }
    return;
  } 
  myvals.addDerivative( getPntrToOutput(0)->getPositionInStream(), jder, der );
}

}
}

#endif

