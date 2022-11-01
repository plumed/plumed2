/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "EvaluateGridFunction.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

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

class InterpolateGrid : 
public ActionWithValue, 
public ActionWithArguments {
private:
  bool firststep, midpoints;
  std::vector<unsigned> nbin;
  std::vector<double> gspacing;
  std::vector<double> forcesToApply;
  EvaluateGridFunction input_grid;
  GridCoordinatesObject output_grid;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  unsigned getNumberOfDerivatives() const override;
  std::vector<unsigned> getValueShapeFromArguments() override;
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
  void apply() override;
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const ;
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.addFlag("MIDPOINTS",false,"interpolate the values of the function at the midpoints of the grid coordinates of the input grid");
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  EvaluateGridFunction ii; ii.registerKeywords( keys );
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  firststep(true)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument to this action");
  if( getPntrToArgument(0)->getRank()==0 || (!getPntrToArgument(0)->hasDerivatives() && !getPntrToArgument(0)->isTimeSeries()) ) error("input to this action should be a grid");

  parseFlag("MIDPOINTS",midpoints); parseVector("GRID_BIN",nbin); parseVector("GRID_SPACING",gspacing); unsigned dimension = getPntrToArgument(0)->getRank();
  if( !midpoints && nbin.size()!=dimension && gspacing.size()!=dimension ) error("MIDPOINTS, GRID_BIN or GRID_SPACING must be set");
  if( midpoints ) {
    log.printf("  evaluating function at midpoints of cells in input grid\n");
  } else if( nbin.size()==dimension ) {
    log.printf("  number of bins in grid %d", nbin[0]);
    for(unsigned i=1; i<nbin.size(); ++i) log.printf(", %d", nbin[i]);
    log.printf("\n");
  } else if( gspacing.size()==dimension ) {
    log.printf("  spacing for bins in grid %f", gspacing[0]);
    for(unsigned i=1; i<gspacing.size(); ++i) log.printf(", %d", gspacing[i]);
    log.printf("\n");
  }

  // Create the input grid
  input_grid.read( this );
  // Need this for creation of tasks
  output_grid.setup( "flat", input_grid.getPbc(), 0, 0.0 );

  // Now add a value
  std::vector<unsigned> shape( dimension, 1 ); 
  if( getPntrToArgument(0)->isTimeSeries() ) addValue( shape );  
  else addValueWithDerivatives( shape );

  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max; getPntrToArgument(0)->getDomain( min, max ); setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  getPntrToOutput(0)->alwaysStoreValues();
}

unsigned InterpolateGrid::getNumberOfDerivatives() const {
  return getPntrToArgument(0)->getRank();
}

void InterpolateGrid::calculate() {
  if( firststep ) {
      input_grid.setup(this); 
      if( midpoints ) {
          nbin.resize( getPntrToOutput(0)->getRank() );
          std::vector<std::string> str_min( input_grid.getMin() ), str_max(input_grid.getMax() ); 
          for(unsigned i=0; i<nbin.size(); ++i) {
              double max, min; Tools::convert( str_min[i], min ); Tools::convert( str_max[i], max ); 
              min += 0.5*input_grid.getGridSpacing()[i];
              if( input_grid.getPbc()[i] ) {
                  nbin[i] = input_grid.getNbin()[i]; max += 0.5*input_grid.getGridSpacing()[i];
              } else {
                  nbin[i] = input_grid.getNbin()[i] - 1; max -= 0.5*input_grid.getGridSpacing()[i];
              }
              Tools::convert( min, str_min[i] ); Tools::convert( max, str_max[i] );
          }
          output_grid.setBounds( str_min, str_max, nbin,  gspacing );
      } else output_grid.setBounds( input_grid.getMin(), input_grid.getMax(), nbin, gspacing );
      getPntrToOutput(0)->setShape( output_grid.getNbin(true) );
      forcesToApply.resize( getPntrToArgument(0)->getNumberOfValues() ); firststep=false;
  }
  plumed_assert( !actionInChain() ); runAllTasks();
}

std::vector<unsigned> InterpolateGrid::getValueShapeFromArguments() {
  return getPntrToOutput(0)->getShape();
}  

void InterpolateGrid::update() {
  if( skipUpdate() ) return;
  calculate();
}
  
void InterpolateGrid::runFinalJobs() {
  if( skipUpdate() ) return;
  calculate();
}

void InterpolateGrid::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
    std::vector<std::string>& max, std::vector<unsigned>& nbin,
    std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  getPntrToArgument(0)->getPntrToAction()->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, dumpcube );
  bool isdists=dumpcube; double units=1.0; gtype="flat";
  for(unsigned i=0; i<argn.size(); ++i) {
    if( argn[i].find(".")==std::string::npos ) { isdists=false; break; }
    std::size_t dot = argn[i].find("."); std::string name = argn[i].substr(dot+1);
    if( name!="x" && name!="y" && name!="z" ) { isdists=false; break; }
  }
  if( isdists ) {
    if( plumed.usingNaturalUnits() ) units = 1.0/0.5292;
    else units = plumed.getAtoms().getUnits().getLength()/.05929;
  }
  if( !firststep ) { 
    std::vector<unsigned> nn( output_grid.getNbin( false ) );
    for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
      if( output_grid.getGridSpacing().size()>0 ) spacing[i]=units*output_grid.getGridSpacing()[i];
      nbin[i]=nn[i];
    }
  }
}

void InterpolateGrid::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  output_grid.getGridPointCoordinates( ind, indices, coords );
}

void InterpolateGrid::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( output_grid.getDimension() ); output_grid.getGridPointCoordinates( current, pos );
  std::vector<double> val(1); Matrix<double> der( 1, output_grid.getDimension() ); input_grid.calc( this, pos, val, der );
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream(); myvals.setValue( ostrn, val[0] );
  if( getPntrToOutput(0)->isTimeSeries() ) return ;
  for(unsigned i=0; i<output_grid.getDimension(); ++i) { myvals.addDerivative( ostrn, i, der(0,i) ); myvals.updateIndex( ostrn, i ); }
}

void InterpolateGrid::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                         const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 ); unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  unsigned istart = bufstart + (1+getNumberOfDerivatives())*code; 
  if( getPntrToOutput(0)->isTimeSeries() ) istart = bufstart + code;
  buffer[istart] += myvals.get( ostrn ); if( getPntrToOutput(0)->isTimeSeries() ) return ;
  for(unsigned i=0; i<output_grid.getDimension(); ++i) buffer[istart+1+i] += myvals.getDerivative( ostrn, i );
}

void InterpolateGrid::apply() {
  if( doNotCalculateDerivatives() ) return; 

  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss; 
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

void InterpolateGrid::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  std::vector<double> pos(output_grid.getDimension()); double ff = getPntrToOutput(0)->getForce(itask); 
  output_grid.getGridPointCoordinates( itask, pos ); input_grid.applyForce( 0, pos, ff, forces );
}

}
}
