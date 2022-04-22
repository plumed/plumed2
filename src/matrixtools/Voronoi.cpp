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
#include "ActionWithInputMatrices.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class Voronoi : public ActionWithInputMatrices {
private:
  bool highest;
  unsigned npairs;
public:
  static void registerKeywords( Keywords& keys );
  explicit Voronoi(const ActionOptions&);
  unsigned getNumberOfColumns() const override;
  void performTask( const unsigned& current, MultiValue& myvals ) const override {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
  void completeMatrixOperations() override;
  void apply() override;
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& krow ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(Voronoi,"VORONOI")

void Voronoi::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys ); 
  keys.addFlag("HIGHEST",false,"input matrix is such that atoms are close if the matrix element is large");
}

Voronoi::Voronoi(const ActionOptions&ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only input one argument");
  // Now create a value
  std::vector<unsigned> shape(2); 
  shape[0]=getPntrToArgument(0)->getShape()[0]; 
  shape[1]=getPntrToArgument(0)->getShape()[1];
  addValue( shape ); getPntrToOutput(0)->setNumberOfTasks( shape[1] );

  parseFlag("HIGHEST",highest);
  if( highest ) log.printf("  connecting each element to largest matrix element \n");
  else log.printf("  connecting each element to smallest matrix element\n");
  checkRead();
}

unsigned Voronoi::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[1];
}

void Voronoi::completeMatrixOperations() {
  runAllTasks();
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

void Voronoi::apply() {
  if( getPntrToOutput(0)->forcesWereAdded() ) error("forces on voronoi matrix cannot work as voroni matrix is not differentiable");
}

}
}
