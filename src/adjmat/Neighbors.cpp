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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR NEIGHBORS
/*
Build a matrix with ones in for the N nearest neighbours of an atom

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class Neighbors : public ActionWithMatrix {
  bool lowest;
  unsigned number;
public:
  static void registerKeywords( Keywords& keys );
  explicit Neighbors(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  unsigned getNumberOfColumns() const override { return number; }
  bool canBeAfterInChain( ActionWithVector* av ) { return av->getLabel()!=(getPntrToArgument(0)->getPntrToAction())->getLabel(); }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override {}
  void turnOnDerivatives() override ;
};

PLUMED_REGISTER_ACTION(Neighbors,"NEIGHBORS")

void Neighbors::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","NLOWEST","0","in each row of the output matrix set the elements that correspond to the n lowest elements in each row of the input matrix equal to one");
  keys.add("compulsory","NHIGHEST","0","in each row of the output matrix set the elements that correspond to the n highest elements in each row of the input matrix equal to one");
  keys.setValueDescription("a matrix in which the ij element is one if the ij-element of the input matrix is one of the NLOWEST/NHIGHEST elements on that row of the input matrix and zero otherwise");
}

Neighbors::Neighbors(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao)
{
  if( getNumberOfArguments()!=1 ) error("found wrong number of arguments in input");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument should be a matrix");
  getPntrToArgument(0)->buildDataStore();

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

  // And get the shape
  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() );
  addValue( shape ); setNotPeriodic();
  checkRead();
}

void Neighbors::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  warning("think about whether your symmetry functions are continuous. If the symmetry function can be calculated from distances only then you can use NEIGHBORS. If you calculate angles between vectors or use the vectors directly then the symmetry function computed using NEIGHBORS is not continuous.  It does not make sense to use such CVs when biasing");
}

unsigned Neighbors::getNumberOfDerivatives() {
  return 0;
}

void Neighbors::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  const Value* wval = getPntrToArgument(0); unsigned nbonds = wval->getRowLength( task_index ), ncols = wval->getShape()[1];
  if( indices.size()!=1+number ) indices.resize( 1 + number );
  myvals.setSplitIndex(1+number);

  unsigned nind=0;
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned ipos = ncols*task_index + wval->getRowIndex( task_index, i );
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) continue ;
    nind++;
  }
  if( number>nind ) plumed_merror("not enough matrix elements were stored");

  // Now build vectors for doing sorting
  std::vector<std::pair<double,unsigned> > rows( nind ); unsigned n=0;
  for(unsigned i=0; i<nbonds; ++i) {
    unsigned iind = wval->getRowIndex( task_index, i );
    unsigned ipos = ncols*task_index + iind;
    double weighti = wval->get( ipos );
    if( weighti<epsilon ) continue ;
    rows[n].first=weighti; rows[n].second=iind; n++;
  }

  // Now do the sort and clear all the stored values ready for recompute
  std::sort( rows.begin(), rows.end() );
  // This is to make this action consistent with what in other matrix actions
  unsigned start_n = getPntrToArgument(0)->getShape()[0];
  // And setup the lowest indices, which are the ones we need to calculate
  for(unsigned i=0; i<number; ++i) {
    indices[i+1] = start_n + rows[nind-1-i].second;
    if( lowest ) indices[i+1] = start_n + rows[i].second;
  }
}

void Neighbors::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  myvals.addValue( getConstPntrToComponent(0)->getPositionInStream(), 1.0 );
}

}
}
