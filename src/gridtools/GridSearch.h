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
#ifndef __PLUMED_gridtools_GridSearch_h
#define __PLUMED_gridtools_GridSearch_h

#include "tools/MinimiseBase.h"
#include "GridVessel.h"
#include <iostream>
#include <memory>

namespace PLMD {
namespace gridtools {

template <class FCLASS>
class GridSearch {
private:
/// This is the pointer to the member function in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
  FCLASS* myclass_func;
  std::unique_ptr<GridVessel> mygrid;
  std::unique_ptr<GridVessel> myfgrid;
public:
  GridSearch( const std::vector<double>& mmin, const std::vector<double>& mmax, const std::vector<unsigned>& ng, const std::vector<unsigned>& nfg, FCLASS* funcc ) :
    myclass_func( funcc )
  {
    // Create the grid objects
    std::string nstr, vstring="COMPONENTS=func COORDINATES=x1";
    for(unsigned i=1; i<mmin.size(); ++i) { Tools::convert(i+1,nstr); vstring += ",x" + nstr; }
    vstring += " PBC=F"; for(unsigned i=1; i<mmin.size(); ++i) vstring += ",F";
    vesselbase::VesselOptions da("mygrid","",-1,vstring,NULL);
    Keywords keys; gridtools::GridVessel::registerKeywords( keys );
    vesselbase::VesselOptions dar( da, keys );
    mygrid=Tools::make_unique<GridVessel>(dar);
    if( nfg[0]>0 ) myfgrid=Tools::make_unique<GridVessel>(dar);

    // Now setup the min and max values for the grid
    std::vector<std::string> gmin( nfg.size() ), gmax( nfg.size() ); std::vector<double> dummy_spacing;
    for(unsigned i=0; i<nfg.size(); ++i) { Tools::convert(mmin[i],gmin[i]); Tools::convert(mmax[i],gmax[i]); }
    mygrid->setBounds( gmin, gmax, ng, dummy_spacing ); mygrid->resize();
    if( myfgrid ) myfgrid->setBounds( gmin, gmax, nfg, dummy_spacing );
  }
  bool minimise( std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
bool GridSearch<FCLASS>::minimise( std::vector<double>& p, engf_pointer myfunc ) {
  std::vector<double> der( p.size() ); std::vector<double> coords( p.size() );
  double initial_eng = (myclass_func->*myfunc)( p, der );
  mygrid->getGridPointCoordinates( 0, coords );
  double emin=(myclass_func->*myfunc)( coords, der );
  mygrid->setValueAndDerivatives( 0, 0, emin, der ); unsigned pmin=0;
  for(unsigned i=1; i<mygrid->getNumberOfPoints(); ++i) {
    mygrid->getGridPointCoordinates( i, coords );
    double eng = (myclass_func->*myfunc)( coords, der );
    mygrid->setValueAndDerivatives( i, 0, eng, der );
    if( eng<emin ) { emin=eng; pmin=i; }
  }
  // This prevents division by zero
  mygrid->setNorm( 1.0 );

  if( myfgrid ) {
    myfgrid->getGridPointCoordinates( 0, coords ); pmin=0;
    double emin=mygrid->getValueAndDerivatives( coords, 0, der );
    for(unsigned i=1; i<myfgrid->getNumberOfPoints(); ++i) {
      myfgrid->getGridPointCoordinates( i, coords );
      double eng = mygrid->getValueAndDerivatives( coords, 0, der );
      if( eng<emin ) { emin=eng; pmin=i; }
    }
    myfgrid->getGridPointCoordinates( pmin, coords );
    double checkEng = (myclass_func->*myfunc)( coords, der );
    if( checkEng<initial_eng ) {
      myfgrid->getGridPointCoordinates( pmin, p );
      return true;
    } else {
      return false;
    }
  }

  if( emin<initial_eng ) {
    mygrid->getGridPointCoordinates( pmin, p );
    return true;
  } else {
    return false;
  }
}

}
}
#endif

