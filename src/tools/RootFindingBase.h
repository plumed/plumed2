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
#ifndef __PLUMED_tools_RootFindingBase_h
#define __PLUMED_tools_RootFindingBase_h

#include "MinimiseBase.h"
#include "Brent1DRootSearch.h"

namespace PLMD {

template <class FCLASS>
class RootFindingBase {
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
  typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der ) const ;
  typedef double(FCLASS::*engfnc_pointer)( const std::vector<double>& p, std::vector<double>& der ) ;
/// The class that calculates the energy given a position
  FCLASS* myclass_func;
/// This actually does the search for a root
  void doSearch( const std::vector<double>& dir, std::vector<double>& p, F1dim<FCLASS>& f1dim  ) const ;
public:
  explicit RootFindingBase( FCLASS* funcc ) : myclass_func(funcc) {}
/// This is the line minimiser
  void linesearch( const std::vector<double>& dir, std::vector<double>& p, engf_pointer myfunc ) const ;
  void lsearch( const std::vector<double>& dir, std::vector<double>& p, engfnc_pointer myfunc ) const ;
};

template <class FCLASS>
void RootFindingBase<FCLASS>::doSearch( const std::vector<double>& dir, std::vector<double>& p, F1dim<FCLASS>& f1dim ) const {
  // Construct an object that will do the line search for the minimum
  Brent1DRootSearch<F1dim<FCLASS> > bb(f1dim);

  // This does the actual search for the root
  double ax=0.0, xx=1.0;
  bb.bracket( ax, xx, &F1dim<FCLASS>::getEng );
  double xmin=bb.search( &F1dim<FCLASS>::getEng );
  for(unsigned i=0; i<p.size(); ++i) p[i] += xmin*dir[i];
}

template <class FCLASS>
void RootFindingBase<FCLASS>::linesearch( const std::vector<double>& dir, std::vector<double>& p, engf_pointer myfunc ) const {
  // Construct the object that turns points on a line into vectors
  F1dim<FCLASS> f1dim( p, dir, myclass_func, myfunc, NULL );
  // Actually do the search
  doSearch( dir, p, f1dim );
}

template <class FCLASS>
void RootFindingBase<FCLASS>::lsearch( const std::vector<double>& dir, std::vector<double>& p, engfnc_pointer myfunc ) const {
  // Construct the object that turns points on a line into vectors
  F1dim<FCLASS> f1dim( p, dir, myclass_func, NULL, myfunc );
  // Actually do the search
  doSearch( dir, p, f1dim );
}

}
#endif
