/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "ActionWithInputGrid.h"

//+PLUMEDOC GRIDANALYSIS CUMULATIVE_INTEGRAL
/*
Calculate a cumulative integral of a function using the trapesium rule.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class CumulativeIntegral : public ActionWithInputGrid {
private:
  double xspace;
public:
  static void registerKeywords( Keywords& keys );
  explicit CumulativeIntegral(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() const ;
  void finishOutputSetup();
  void performTask( const unsigned& current, MultiValue& myvals ) const {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ; 
};

PLUMED_REGISTER_ACTION(CumulativeIntegral,"CUMULATIVE_INTEGRAL")

void CumulativeIntegral::registerKeywords( Keywords& keys ) {
  ActionWithInputGrid::registerKeywords( keys );
}

CumulativeIntegral::CumulativeIntegral(const ActionOptions&ao):
  Action(ao),
  ActionWithInputGrid(ao)
{
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) error("input should be a function on a grid");
  // Now create the shape
  std::vector<unsigned> shape( getNumberOfArguments() );
  for(unsigned i=0;i<shape.size();++i) shape[i] = getPntrToArgument(0)->getShape()[i];
  // Retrieve information about the grid
  std::vector<std::string> argn( shape.size() ), min( shape.size() ), max( shape.size() ); std::string gtype;
  std::vector<unsigned> nbin( shape.size() ); std::vector<double> spacing( shape.size() ); std::vector<bool> ipbc( shape.size() );
  (getPntrToArgument(0)->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, ipbc, false ); 
  if( gtype!="flat" ) error("grid should be flat - does not work with spherical grids");
  for(unsigned i=1;i<getNumberOfArguments();++i) {
      std::vector<std::string> iargn( shape.size() ), imin( shape.size() ), imax( shape.size() ); 
      std::vector<unsigned> inbin( shape.size() ); std::vector<double> ispacing( shape.size() ); std::vector<bool> iipbc( shape.size() );
      (getPntrToArgument(i)->getPntrToAction())->getInfoForGridHeader( gtype, iargn, imin, imax, inbin, ispacing, iipbc, false );
      if( gtype!="flat" ) error("grid should be flat - does not work with spherical grids");
      for(unsigned j=0;j<shape.size();++j) {
          if( iargn[j]!=argn[j] || imin[j]!=min[j] || imax[j]!=max[j] || inbin[j]!=nbin[j] || ispacing[j]!=spacing[j] || iipbc[j]!=ipbc[j] ) {
              error("mismatch between grid containing components of vectors for integration"); 
          }
      }
  }
  plumed_assert( shape.size()==1 );   // Higher dimensions not implemented yet
  // And create the value 
  addValueWithDerivatives( shape ); xspace=spacing[0];
  // Turn off all parallelization
  runInSerial(); runWithoutOpenMP();
}

void CumulativeIntegral::finishOutputSetup() {
  for(unsigned i=0; i<getPntrToArgument(0)->getNumberOfValues(); ++i) addTaskToList(i);
}

unsigned CumulativeIntegral::getNumberOfDerivatives() const {
  return (getPntrToArgument(0)->getPntrToAction())->getNumberOfDerivatives();
}

void CumulativeIntegral::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                            const unsigned& bufstart, std::vector<double>& buffer ) const { 
  if( code==0 ) return;
  Value* argv=getPntrToArgument(0); plumed_dbg_assert( bufstart + (1+getNumberOfDerivatives())*code < buffer.size() ); 
  buffer[ bufstart + (1+getNumberOfDerivatives())*code ] = buffer[ bufstart + (1+getNumberOfDerivatives())*(code-1) ] + (argv->get(code-1) + argv->get(code))*xspace / 2.; 
}

}
}
