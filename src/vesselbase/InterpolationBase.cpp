/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
  
   See http://www.plumed-code.org for more information.
  
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
#include "InterpolationBase.h"
#include "FunctionOnGrid.h"
#include "tools/Exception.h"
#include "CubicInterpolation.h"
#include "BicubicInterpolation.h"

namespace PLMD {
namespace vesselbase {

InterpolationBase::InterpolationBase( GridVesselBase* gg, const unsigned dstart ):
tablesAreSet(false),
dimension(gg->dimension),
mygrid(gg)
{
  mygrid->interpolating=true;
  data = dstart*(dimension + 1); ddx.resize( dimension );
  FunctionOnGrid* noint = dynamic_cast<FunctionOnGrid*>( gg );
  plumed_massert( !noint, "Cannot interpolate grid without derivatives");
  plumed_massert( data<gg->nper && gg->nper%(gg->dimension+1)==0 , "Inappropriate data stored on grid"); 
}

void InterpolationBase::getGridPointSpacing( std::vector<double>& spacing ) const {
  plumed_dbg_assert( spacing.size()==dimension );
  for(unsigned i=0;i<dimension;++i) spacing[i]=mygrid->dx[i];
}

void InterpolationBase::getNumberOfGridPoints( std::vector<unsigned>& bin ) const {
  plumed_dbg_assert( bin.size()==dimension );
  for(unsigned i=0;i<dimension;++i) bin[i]=mygrid->nbin[i];
}

double InterpolationBase::getFunctionValue( const std::vector<double>& pos ){
  plumed_dbg_assert( interpolationTablesWereSet() && pos.size()==dimension );
  unsigned box = mygrid->getLocationOnGrid( pos, ddx );
  return interpolateFunction( box, ddx );
}

unsigned InterpolationBase::getGridStride( const unsigned& i ) const {
  return mygrid->stride[i];
}

unsigned InterpolationBase::getBoxIndex( const std::vector<unsigned>& pn ) const {
  return mygrid->getIndex( pn );
}

double InterpolationBase::getValue( const std::vector<unsigned>& pn ){
  unsigned mybox = mygrid->getIndex( pn );
  return mygrid->getGridElement( mybox, data );
}

double InterpolationBase::getValue( const unsigned& mybox ){
  return mygrid->getGridElement( mybox, data );
}

double InterpolationBase::getValueAndDerivatives( const std::vector<unsigned>& pn, std::vector<double>& der ){
  plumed_dbg_assert( pn.size()==dimension && der.size()==dimension );
  unsigned mybox = mygrid->getIndex( pn );
  for(unsigned i=0;i<dimension;++i) der[i]=mygrid->getGridElement( mybox, data + 1 + i ); 
  return mygrid->getGridElement( mybox, data );
}

void InterpolationBase::set_table(){
  tablesAreSet=true; mygrid->dataHasChangedSinceInterpol=false;
  setInterpolationTables();
}

InterpolationBase* InterpolationBase::createCubicInterpolator( GridVesselBase* gg, const unsigned dstart ){
  InterpolationBase* my_interp;
  if( gg->dimension==1 ) my_interp = new CubicInterpolation( gg, dstart );
  else if( gg->dimension==2 ) my_interp = new BicubicInterpolation( gg, dstart );
  else plumed_merror("Cannot use cubic interpolation with 3 dimensional functions");
  return my_interp;
}

}
}
