/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2024 The plumed team
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

#include "core/Value.h"
#include "tools/MinimiseBase.h"
#include "GridCoordinatesObject.h"
#include "Interpolator.h"
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
  bool using_fgrid;
  std::unique_ptr<Value> value;
  std::unique_ptr<Interpolator> myinterp;
  GridCoordinatesObject mygrid;
  GridCoordinatesObject myfgrid;
public:
  GridSearch( const std::vector<double>& mmin, const std::vector<double>& mmax, const std::vector<std::size_t>& ng, const std::vector<std::size_t>& nfg, FCLASS* funcc ) :
    myclass_func( funcc ),
    using_fgrid(false) {
    // Setup the min and max values for the grid
    std::vector<std::string> gmin( nfg.size() ), gmax( nfg.size() );
    std::vector<bool> pbc( nfg.size(), false );
    std::vector<double> dummy_spacing;
    for(unsigned i=0; i<nfg.size(); ++i) {
      Tools::convert(mmin[i],gmin[i]);
      Tools::convert(mmax[i],gmax[i]);
    }
    // Create the coarse grid grid objects
    mygrid.setup( "flat", pbc, 0, 0.0 );
    mygrid.setBounds( gmin, gmax, ng, dummy_spacing );
    // Setup the fine grid object
    if( nfg[0]>0 ) {
      using_fgrid=true;
      myfgrid.setup("flat", pbc, 0, 0.0 );
      dummy_spacing.resize(0);
      myfgrid.setBounds( gmin, gmax, nfg, dummy_spacing );
    }
    value.reset( new Value( NULL, "gval", true, mygrid.getNbin(true)) );
    myinterp.reset( new Interpolator( value.get(), mygrid ) );
  }
  void setGridElement( const unsigned& ind, const double& emin, const std::vector<double>& der );
  bool minimise( std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
void GridSearch<FCLASS>::setGridElement( const unsigned& ind, const double& emin, const std::vector<double>& der ) {
  value->set( ind, emin );
  for(unsigned j=0; j<der.size(); ++j) {
    value->setGridDerivatives( ind, j, der[j] );
  }
}

template <class FCLASS>
bool GridSearch<FCLASS>::minimise( std::vector<double>& p, engf_pointer myfunc ) {
  std::vector<double> der( p.size() );
  std::vector<double> coords( p.size() );
  double initial_eng = (myclass_func->*myfunc)( p, der );
  mygrid.getGridPointCoordinates( 0, coords );
  unsigned pmin=0;
  double emin=(myclass_func->*myfunc)( coords, der );
  setGridElement( 0, emin, der );
  for(unsigned i=1; i<mygrid.getNumberOfPoints(); ++i) {
    mygrid.getGridPointCoordinates( i, coords );
    double eng = (myclass_func->*myfunc)( coords, der );
    setGridElement( i, eng, der );
    if( eng<emin ) {
      emin=eng;
      pmin=i;
    }
  }

  if( using_fgrid ) {
    myfgrid.getGridPointCoordinates( 0, coords );
    pmin=0;
    double eminVal = myinterp->splineInterpolation( coords, der );
    for(unsigned i=1; i<myfgrid.getNumberOfPoints(); ++i) {
      myfgrid.getGridPointCoordinates( i, coords );
      const double eng = myinterp->splineInterpolation( coords, der );
      if( eng<eminVal ) {
        eminVal=eng;
        pmin=i;
      }
    }
    myfgrid.getGridPointCoordinates( pmin, coords );
    double checkEng = (myclass_func->*myfunc)( coords, der );
    if( checkEng<initial_eng ) {
      myfgrid.getGridPointCoordinates( pmin, p );
      return true;
    } else {
      return false;
    }
  }

  if( emin<initial_eng ) {
    mygrid.getGridPointCoordinates( pmin, p );
    return true;
  }
  return false;
}

}
}
#endif

