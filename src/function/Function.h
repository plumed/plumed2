/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
  bool matinp, matout;
  unsigned nderivatives;
  std::vector<double> forcesToApply;
  bool hasGridOutput() const ;
  std::vector<unsigned> getShape();
  void evaluateAllFunctions();
  void fixTimeSeries();
protected:
  bool getPeriodFromArg;
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name );
  void addValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addDerivative( const unsigned& ival, const unsigned& jder, const double& der, MultiValue& myvals ) const ;
public:
  static void registerKeywords(Keywords&);
  explicit Function(const ActionOptions&);
  virtual ~Function() {}
  virtual void calculate() override; 
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
  virtual void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
  virtual void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const = 0;
  void apply() override;
  void update() override;
  void runFinalJobs();
  unsigned getNumberOfDerivatives() const override;
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

  if( actionInChain() ) {
    unsigned istrn = getArgumentPositionInStream(jder, myvals);
    unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
    for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
      unsigned kind=myvals.getActiveIndex(istrn,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[jder] + kind, der*myvals.getDerivative( istrn, kind ) );
    }
    return;
  }
  if( getPntrToOutput(ival)->getRank()>0 && getPntrToOutput(ival)->hasDerivatives() ) {
    if( getPntrToArgument(jder)->getRank()==0 ) {
      myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), getPntrToOutput(ival)->getRank()+jder, der );
    } else {
      unsigned np = myvals.getTaskIndex(), ostrn = getPntrToOutput(ival)->getPositionInStream();
      for(unsigned i=0; i<getPntrToArgument(jder)->getRank(); ++i) {
        myvals.addDerivative( ostrn, i, der*getPntrToArgument(jder)->getGridDerivative( np, i ) );
      } 
      if( nderivatives>getPntrToArgument(jder)->getRank() ) myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), getPntrToOutput(ival)->getRank()+jder, der );
    }
    return;
  }
  if( arg_ends.size()>0 ) {
    unsigned base=0;
    for(unsigned i=0; i<jder; ++i) {
      for(unsigned j=arg_ends[i]; j<arg_ends[i+1]; ++j) base += getPntrToArgument(j)->getNumberOfValues( getLabel() );
    }
    if( arg_ends[jder+1]==(arg_ends[jder]+1) && getPntrToArgument(arg_ends[jder])->getRank()==0 ) {
      myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), base, der );
    } else {
      myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), base + myvals.getTaskIndex(), der );
    }
    return;
  }
  myvals.addDerivative( getPntrToOutput(ival)->getPositionInStream(), jder, der );
}

}
}

#endif

