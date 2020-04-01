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

class Voronoi :
  public ActionWithValue,
  public ActionWithArguments
{
private:
  bool highest;
public:
  static void registerKeywords( Keywords& keys );
  explicit Voronoi(const ActionOptions&);
  unsigned getNumberOfDerivatives() const ;
  void calculate();
  void performTask( const unsigned& current, MultiValue& myvals ) const {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void apply() {}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(Voronoi,"VORONOI")

void Voronoi::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.addFlag("HIGHEST",false,"input matrix is such that atoms are close if the matrix element is large");
}

Voronoi::Voronoi(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only input one argument");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input should be a matrix");
  // Now create a value
  std::vector<unsigned> shape(2); getPntrToArgument(0)->buildDataStore( getLabel() );
  shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(0)->getShape()[1];
  for(unsigned i=0; i<shape[1]; ++i) addTaskToList( i );
  addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues();

  parseFlag("HIGHEST",highest);
  if( highest ) log.printf("  connecting each element to largest matrix element \n");
  else log.printf("  connecting each element to smallest matrix element\n");
  checkRead();
}

unsigned Voronoi::getNumberOfDerivatives() const {
  return 0;
}

void Voronoi::calculate() {
  if( skipCalculate() ) return;
  plumed_dbg_assert( !actionInChain() );
  runAllTasks();
}

void Voronoi::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  runAllTasks();
}

void Voronoi::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() && getFullNumberOfTasks()==0 );
  std::vector<unsigned> shape( getPntrToArgument(0)->getShape() ); 
  for(unsigned i=0; i<shape[1]; ++i) addTaskToList( i );
  getPntrToOutput(0)->setShape( shape ); runAllTasks();
} 

void Voronoi::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                 const unsigned& bufstart, std::vector<double>& buffer ) const {
  Value* arg0 = getPntrToArgument(0); unsigned nv = 0; double minmax = arg0->get( code ); 
  for(unsigned i=0;i<arg0->getShape()[0];++i) {
      double value = arg0->get( i*arg0->getShape()[1] + code );
      if( highest && value>minmax ) { minmax = value; nv = i; }
      else if( !highest && value<minmax ) { minmax = value; nv = i; } 
  }
  buffer[bufstart + nv*arg0->getShape()[1] + code] = 1;
}

}
}
