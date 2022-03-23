/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "CreateGridInSetup.h"

using namespace std;

namespace PLMD {
namespace gridtools {

void CreateGridInSetup::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.remove("NUMERICAL_DERIVATIVES"); keys.remove("SERIAL"); keys.remove("TIMINGS");
  keys.add( "hidden", "LABEL", "a label for the action so that its output can be referenced in the input to other actions.  Actions with scalar output are referenced using their label only.  Actions with vector output must have a separate label for every component.  Individual componets are then refered to using label.component" );
}

CreateGridInSetup::CreateGridInSetup(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao),
  ActionWithValue(ao)
{
}

void CreateGridInSetup::createGridAndValue( const std::string& gtype, const std::vector<bool>& ipbc, const unsigned& nfermi, 
                                            const std::vector<std::string>& gmin, const std::vector<std::string>& gmax,
                                            const std::vector<unsigned>& gbin ) {
  gridobject.setup( gtype, ipbc, nfermi, 0.0 ); std::vector<double> gspacing;
  if( gtype=="flat" ) {
      gridobject.setBounds( gmin, gmax, gbin, gspacing );
      // Now create the value
      std::vector<unsigned> shape( gridobject.getNbin(true) );
      ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  } else {
      std::vector<unsigned> shape( 3 ); shape[0]=gbin[0]; shape[1]=shape[2]=1;
      ActionWithValue::addValueWithDerivatives( shape ); setNotPeriodic();
  }
  for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
}

unsigned CreateGridInSetup::getNumberOfDerivatives() const {
  return labels.size();
}

void CreateGridInSetup::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                            std::vector<std::string>& max, std::vector<unsigned>& out_nbin,
                                            std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  plumed_assert( !dumpcube ); gtype=gridobject.getGridType(); 
  if( gtype=="fibonacci" ) { 
    out_nbin[0]=getPntrToOutput(0)->getShape()[0];  
    for(unsigned i=0;i<labels.size();++i) argn[i] = labels[i];
  } else {
    std::vector<unsigned> nbin( gridobject.getNbin(false) );
    for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
      argn[i]=labels[i]; double gmin, gmax;
      if( gridobject.getMin().size()>0 ) {
        min[i] = gridobject.getMin()[i]; max[i] = gridobject.getMax()[i];
      }
      if( nbin.size()>0 ) out_nbin[i]=nbin[i];
      if( spacing.size()>0 ) spacing[i]=gridobject.getGridSpacing()[i];
      pbc[i]=gridobject.isPeriodic(i);
    }
  }
}

void CreateGridInSetup::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridobject.getGridPointCoordinates( ind, indices, coords );
}

void CreateGridInSetup::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  if( setlength ) gridobject.putCoordinateAtValue( ind, getPntrToOutput(0)->get(ind), coords );
  else  gridobject.putCoordinateAtValue( ind, 1.0, coords );
}

}
}

