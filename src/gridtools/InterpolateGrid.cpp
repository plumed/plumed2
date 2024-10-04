/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "EvaluateGridFunction.h"
#include "ActionWithGrid.h"

//+PLUMEDOC GRIDANALYSIS INTERPOLATE_GRID
/*
Interpolate a smooth function stored on a grid onto a grid with a smaller grid spacing.

This action takes a function evaluated on a grid as input and can be used to interpolate the values of that
function on to a finer grained grid.  The interpolation within this algorithm is done using splines.

\par Examples

The input below can be used to post process a trajectory.  It calculates a \ref HISTOGRAM as a function the
distance between atoms 1 and 2 using kernel density estimation.  During the calculation the values of the kernels
are evaluated at 100 points on a uniform grid between 0.0 and 3.0.  Prior to outputting this function at the end of the
simulation this function is interpolated onto a finer grid of 200 points between 0.0 and 3.0.

\plumedfile
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ii: INTERPOLATE_GRID GRID=hA1 GRID_BIN=200
DUMPGRID GRID=ii FILE=histo.dat
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class InterpolateGrid : public ActionWithGrid {
private:
  bool midpoints;
  std::vector<unsigned> nbin;
  std::vector<double> gspacing;
  EvaluateGridFunction input_grid;
  GridCoordinatesObject output_grid;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  void setupOnFirstStep( const bool incalc ) override ;
  unsigned getNumberOfDerivatives() override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const override ;
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.use("ARG");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.addFlag("MIDPOINTS",false,"interpolate the values of the function at the midpoints of the grid coordinates of the input grid");
  EvaluateGridFunction ii;
  ii.registerKeywords( keys );
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) {
    error("input to this action should be a grid");
  }

  parseFlag("MIDPOINTS",midpoints);
  parseVector("GRID_BIN",nbin);
  parseVector("GRID_SPACING",gspacing);
  unsigned dimension = getPntrToArgument(0)->getRank();
  if( !midpoints && nbin.size()!=dimension && gspacing.size()!=dimension ) {
    error("MIDPOINTS, GRID_BIN or GRID_SPACING must be set");
  }
  if( midpoints ) {
    log.printf("  evaluating function at midpoints of cells in input grid\n");
  } else if( nbin.size()==dimension ) {
    log.printf("  number of bins in grid %d", nbin[0]);
    for(unsigned i=1; i<nbin.size(); ++i) {
      log.printf(", %d", nbin[i]);
    }
    log.printf("\n");
  } else if( gspacing.size()==dimension ) {
    log.printf("  spacing for bins in grid %f", gspacing[0]);
    for(unsigned i=1; i<gspacing.size(); ++i) {
      log.printf(", %d", gspacing[i]);
    }
    log.printf("\n");
  }
  // Create the input grid
  input_grid.read( this );
  // Need this for creation of tasks
  output_grid.setup( "flat", input_grid.getPbc(), 0, 0.0 );

  // Now add a value
  std::vector<unsigned> shape( dimension, 0 );
  addValueWithDerivatives( shape );

  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  setupOnFirstStep( false );
}

void InterpolateGrid::setupOnFirstStep( const bool incalc ) {
  input_grid.setup( this );
  ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  const GridCoordinatesObject& mygrid = ag->getGridCoordinatesObject();
  if( midpoints ) {
    double min, max;
    nbin.resize( getPntrToComponent(0)->getRank() );
    std::vector<std::string> str_min( input_grid.getMin() ), str_max(input_grid.getMax() );
    for(unsigned i=0; i<nbin.size(); ++i) {
      if( incalc ) {
        Tools::convert( str_min[i], min );
        Tools::convert( str_max[i], max );
        min += 0.5*input_grid.getGridSpacing()[i];
      }
      if( input_grid.getPbc()[i] ) {
        nbin[i] = input_grid.getNbin()[i];
        if( incalc ) {
          max += 0.5*input_grid.getGridSpacing()[i];
        }
      } else {
        nbin[i] = input_grid.getNbin()[i] - 1;
        if( incalc ) {
          max -= 0.5*input_grid.getGridSpacing()[i];
        }
      }
      if( incalc ) {
        Tools::convert( min, str_min[i] );
        Tools::convert( max, str_max[i] );
      }
    }
    output_grid.setBounds( str_min, str_max, nbin,  gspacing );
  } else {
    output_grid.setBounds( mygrid.getMin(), mygrid.getMax(), nbin, gspacing );
  }
  getPntrToComponent(0)->setShape( output_grid.getNbin(true) );
  if( !incalc ) {
    gspacing.resize(0);
  }
}

unsigned InterpolateGrid::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getRank();
}

const GridCoordinatesObject& InterpolateGrid::getGridCoordinatesObject() const {
  return output_grid;
}

std::vector<std::string> InterpolateGrid::getGridCoordinateNames() const {
  ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  return ag->getGridCoordinateNames();
}

void InterpolateGrid::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( output_grid.getDimension() );
  output_grid.getGridPointCoordinates( current, pos );
  std::vector<double> val(1);
  Matrix<double> der( 1, output_grid.getDimension() );
  input_grid.calc( this, pos, val, der );
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  myvals.setValue( ostrn, val[0] );
  for(unsigned i=0; i<output_grid.getDimension(); ++i) {
    myvals.addDerivative( ostrn, i, der(0,i) );
    myvals.updateIndex( ostrn, i );
  }
}

void InterpolateGrid::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
    const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 );
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  unsigned istart = bufstart + (1+output_grid.getDimension())*code;
  buffer[istart] += myvals.get( ostrn );
  for(unsigned i=0; i<output_grid.getDimension(); ++i) {
    buffer[istart+1+i] += myvals.getDerivative( ostrn, i );
  }
}

void InterpolateGrid::gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  std::vector<double> pos(output_grid.getDimension());
  double ff = myval->getForce(itask);
  output_grid.getGridPointCoordinates( itask, pos );
  input_grid.applyForce( this, pos, ff, forces );
}


}
}
