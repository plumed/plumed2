/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "GridFunction.h"
#include "ActionWithInputGrid.h"

namespace PLMD {
namespace gridtools {

void GridFunction::registerKeywords( Keywords& keys ){
  GridVessel::registerKeywords( keys );
}

GridFunction::GridFunction( const vesselbase::VesselOptions& da ):
GridVessel(da)
{
}

void GridFunction::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_dbg_assert( myvals.getNumberOfValues()==(nper+1) );
  for(unsigned i=0;i<myvals.getNumberOfValues();++i) buffer[bufstart + nper*current + i] += myvals.get(i+1);
  return;
}


void GridFunction::finish( const std::vector<double>& buffer ){
  for(unsigned i=0;i<data.size();++i) data[i]+=buffer[bufstart + i];
}

void GridFunction::incorporateRestartDataIntoGrid( const double& old_norm, std::vector<double>& indata ){
  ActionWithInputGrid* myfunc = dynamic_cast<ActionWithInputGrid*>( getAction() );
  std::vector<double> pin( nper ), pout( nper ); 
  for(unsigned i=0;i<getNumberOfPoints();++i){
      for(unsigned j=0;j<nper;++j) pin[j]=indata[i*nper+j];
      myfunc->invertTask( pin, pout ); 
      for(unsigned j=0;j<nper;++j) indata[i*nper+j]=pout[j]; 
  } 
  (myfunc->mygrid)->incorporateRestartDataIntoGrid( old_norm, indata );
}

}
}
