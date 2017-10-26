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
#ifndef __PLUMED_function_Function_h
#define __PLUMED_function_Function_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"

namespace PLMD {
namespace function {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new CV function, within it there is
\ref AddingAFunction "information" as to how to go about implementing a new function.
*/

class Function:
  public ActionWithValue,
  public ActionWithArguments
{
private:
  bool firststep;
  unsigned nderivatives;
  std::vector<double> forcesToApply;
  bool hasGridOutput() const ;
  std::vector<unsigned> getShape();
  void evaluateAllFunctions();
protected:
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name );
  void addValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const ;
public:
  static void registerKeywords(Keywords&);
  explicit Function(const ActionOptions&);
  virtual ~Function() {}
  virtual void calculate();
  void buildCurrentTaskList( std::vector<unsigned>& tflags );
  void getInfoForGridHeader( std::vector<std::string>& argn, std::vector<std::string>& min, 
                             std::vector<std::string>& max, std::vector<unsigned>& nbin, 
                             std::vector<double>& spacing, std::vector<bool>& pbc ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void gatherGridAccumulators( const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
  virtual void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const = 0;
  void apply();
  void update();
  void runFinalJobs();
  unsigned getNumberOfDerivatives() const ;
};

inline
unsigned Function::getNumberOfDerivatives() const {
  return nderivatives;
}

inline
void Function::addValue( const unsigned& ival, const double& val, MultiValue& myvals ) const {
  myvals.addValue( getPntrToOutput(ival)->getPositionInStream(), val ); 
}

inline
void Function::addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() && !(getPntrToOutput(ival)->getRank()>0 && getPntrToOutput(ival)->hasDerivatives()) ) return ;

  if( actionInChain() ){ 
      unsigned istrn = getPntrToArgument(jder)->getPositionInStream();
      unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
      for(unsigned k=0;k<myvals.getNumberActive(istrn);++k){ 
          unsigned kind=myvals.getActiveIndex(istrn,k); 
          myvals.addDerivative( ostrn, arg_deriv_starts[jder] + kind, der*myvals.getDerivative( istrn, kind ) );
      }
      return; 
  }
  if( getPntrToOutput(ival)->getRank()>0 && getPntrToOutput(ival)->hasDerivatives() ) {
      unsigned np = myvals.getTaskIndex(), ostrn = getPntrToOutput(ival)->getPositionInStream();
      for(unsigned i=0;i<getPntrToArgument(jder)->getRank();++i) {
          myvals.addDerivative( ostrn, i, der*getPntrToArgument(jder)->getGridDerivative( np, i ) ); 
      }
      return;
  }
  if( getPntrToArgument(0)->getRank()>0 ) { 
      unsigned base=0; for(unsigned i=0;i<jder;++i) base += getPntrToArgument(i)->getSize();
      myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), base + myvals.getTaskIndex(), der ); 
      return; 
  }
  myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), arg_ends[jder] + myvals.getTaskIndex(), der ); 
}

}
}

#endif

