/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/Angle.h"

namespace PLMD {
namespace adjmat {

class Neighbors :
  public ActionWithValue,
  public ActionWithArguments
{
private:
  bool lowest;
  unsigned number;
public:
  static void registerKeywords( Keywords& keys );
  explicit Neighbors(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void calculate() { if( !actionInChain() ) plumed_error(); }
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void apply() {}
};

PLUMED_REGISTER_ACTION(Neighbors,"NEIGHBORS")

void Neighbors::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","NLOWEST","0","in each row of the output matrix set the elements that correspond to the n lowest elements in each row of the input matrix equal to one");
  keys.add("compulsory","NHIGHEST","0","in each row of the output matrix set the elements that correspond to the n highest elements in each row of the input matrix equal to one");
}

Neighbors::Neighbors(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only input one argument to Neighbors");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input to neighbors should be a matrix");
  // Add this to the chain
  std::vector<std::string> alabels(1); alabels[0]=(getPntrToArgument(0)->getPntrToAction())->getLabel();
  (getPntrToArgument(0)->getPntrToAction())->addActionToChain( alabels, this );
  // Now create a value
  std::vector<unsigned> shape(2); getPntrToArgument(0)->buildDataStore( getLabel() );
  shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(0)->getShape()[1];
  for(unsigned i=0; i<shape[0]; ++i) addTaskToList( i );
  addValue( shape );

  unsigned nlow; parse("NLOWEST",nlow);
  unsigned nhigh; parse("NHIGHEST",nhigh);
  if( nlow==0 && nhigh==0 ) error("missing NLOWEST or NHIGHEST keyword one of these two keywords must be set in input");
  if( nlow>0 && nhigh>0 ) error("should only be one of NLOWEST or NHIGHEST set in input");
  if( nlow>0 ) {
    number=nlow; lowest=true;
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d lowest elements in each row of the input matrix\n",number);
  }
  if( nhigh>0 ) {
    number=nhigh; lowest=false;
    log.printf("  output matrix will have non-zero values for elements that correpsond to the %d highest elements in each row of the input matrix\n",number);
  }
  checkRead();
}

unsigned Neighbors::getNumberOfDerivatives() const {
  return 0;
}

void Neighbors::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( !myvals.inVectorCall() ) return ;
  // Work out how many things we have to work with
  unsigned nind=0, matind = getPntrToArgument(0)->getPositionInMatrixStash();
  for(unsigned i=0; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i);
    double weighti=myvals.getStashedMatrixElement( matind, iind );
    if( !lowest && weighti<epsilon ) continue ;
    nind++;
  }
  if( number>nind ) error("not enough matrix elements were stored");
  // Now build vectors for doing sorting
  std::vector<std::pair<double,unsigned> > rows( nind ); unsigned n=0;
  for(unsigned i=0; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i);
    double weighti=myvals.getStashedMatrixElement( matind, iind );
    if( !lowest && weighti<epsilon ) continue ;
    rows[n].first=weighti; rows[n].second=iind; n++;
  }
  // Now do the sort
  std::sort( rows.begin(), rows.end() );
  // And do everything that follows me for only the relevant elements of the matrix
  unsigned matout = getPntrToOutput(0)->getPositionInStream();
  ActionWithValue* av = (getPntrToArgument(0)->getPntrToAction())->getActionThatCalculates();
  for(unsigned i=0; i<number; ++i) {
    myvals.setValue( matout, 1.0 );
    unsigned colno = rows[nind-1-i].second; if( lowest ) colno = rows[i].second;
    av->runTask( av->getLabel(), myvals.getTaskIndex(), current, myvals.getNumberOfIndicesInFirstBlock()+colno, myvals );
    av->clearMatrixElements( myvals );
  }
}

}
}
