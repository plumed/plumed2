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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"

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

class InterpolateGrid : public ActionWithInputGrid {
private:
  std::vector<unsigned> nbin;
  std::vector<double> gspacing;
  GridCoordinatesObject gridcoords;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void finishOutputSetup();
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ;
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ) {
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithInputGrid(ao)
{
  parseVector("GRID_BIN",nbin); parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=gridobject.getDimension() && gspacing.size()!=gridobject.getDimension() ) error("GRID_BIN or GRID_SPACING must be set");
  if( nbin.size()==gridobject.getDimension() ) {
    log.printf("  number of bins in grid %d", nbin[0]);
    for(unsigned i=1; i<nbin.size(); ++i) log.printf(", %d", nbin[i]);
    log.printf("\n");
  } else if( gspacing.size()==gridobject.getDimension() ) {
    log.printf("  spacing for bins in grid %f", gspacing[0]);
    for(unsigned i=1; i<gspacing.size(); ++i) log.printf(", %d", gspacing[i]);
    log.printf("\n");
  }

  // Need this for creation of tasks
  std::vector<bool> ipbc( gridobject.getDimension() ); for(unsigned i=0; i<ipbc.size(); ++i) ipbc[i] = gridobject.isPeriodic(i);
  gridcoords.setup( "flat", ipbc, 0, 0.0 );

  // Now add a value
  std::vector<unsigned> shape( gridobject.getDimension() ); addValueWithDerivatives( shape );
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max; getPntrToArgument(0)->getDomain( min, max ); setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  getPntrToOutput(0)->alwaysStoreValues();
}

void InterpolateGrid::finishOutputSetup() {
  gridcoords.setBounds( gridobject.getMin(), gridobject.getMax(), nbin, gspacing );
  getPntrToOutput(0)->setShape( gridcoords.getNbin(true) );
  for(unsigned i=0; i<gridcoords.getNumberOfPoints(); ++i) addTaskToList(i);
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
    if( plumed.getAtoms().usingNaturalUnits() ) units = 1.0/0.5292;
    else units = plumed.getAtoms().getUnits().getLength()/.05929;
  }
  std::vector<unsigned> nn( gridcoords.getNbin( false ) );
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    if( gridcoords.getGridSpacing().size()>0 ) spacing[i]=units*gridcoords.getGridSpacing()[i];
    nbin[i]=nn[i];
  }
}

void InterpolateGrid::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridcoords.getGridPointCoordinates( ind, indices, coords );
}

void InterpolateGrid::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( gridcoords.getDimension() ); gridcoords.getGridPointCoordinates( current, pos );
  std::vector<double> der( gridcoords.getDimension() ); double val = getFunctionValueAndDerivatives( pos, der );
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream(); myvals.setValue( ostrn, val );
  for(unsigned i=0; i<gridcoords.getDimension(); ++i) { myvals.addDerivative( ostrn, i, der[i] ); myvals.updateIndex( ostrn, i ); }
}

void InterpolateGrid::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                         const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 ); unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  unsigned istart = bufstart + (1+getNumberOfDerivatives())*code; buffer[istart] += myvals.get( ostrn );
  for(unsigned i=0; i<gridcoords.getDimension(); ++i) buffer[istart+1+i] += myvals.getDerivative( ostrn, i );
}

}
}
