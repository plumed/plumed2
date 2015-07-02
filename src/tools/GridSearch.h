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
#ifndef __PLUMED_tools_GridSearch_h
#define __PLUMED_tools_GridSearch_h

#include "MinimiseBase.h"
#include <iostream>
#include <math.h>
namespace PLMD{

template <class FCLASS>
class GridSearch { 
private:
/// This is the pointer to the member funciton in the energy
/// calculating class that calculates the energy
   typedef double(FCLASS::*engf_pointer)( const std::vector<double>& p, std::vector<double>& der );
   FCLASS* myclass_func;
   std::vector<unsigned> ngrid;
   std::vector<double> min, delr;
   unsigned np;
public:
   GridSearch( const std::vector<double>& mmin, const std::vector<double>& mmax, const std::vector<unsigned>& ng, FCLASS* funcc ) : 
     myclass_func( funcc ), 
     ngrid(ng),
     min(mmin),
     delr(mmin.size())
     {
        // Work out the stride
        for(unsigned i=0;i<delr.size();++i) delr[i] = ( mmax[i] - mmin[i] ) / static_cast<double>( ng[i] ); 
        // And the number of poitns in the grid
        np=1; for(unsigned i=0;i<ngrid.size();++i) np*=ngrid[i];
     }
     void minimise( std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
void GridSearch<FCLASS>::minimise( std::vector<double>& p, engf_pointer myfunc ){

   std::vector<double> fake_der( p.size() ); std::vector<unsigned> ind( p.size() ); 
   unsigned pmin=0; double emin=(myclass_func->*myfunc)( min, fake_der );
   for(unsigned i=1;i<np;++i){
      unsigned myind=i; ind[0]=myind%ngrid[0];
      for(unsigned j=1;j<p.size()-1;++j){
          myind=(myind-ind[j-1])/ngrid[j-1];
          ind[j]=myind%ngrid[j];
      }
      if( p.size()>=2 ) ind[p.size()-1] = (myind - ind[p.size()-2])/ngrid[p.size()-2];
      for(unsigned j=0;j<p.size();++j) p[j] = min[j] + ind[j]*delr[j];

      double eng = (myclass_func->*myfunc)( p, fake_der );
      if( eng<emin ){ emin=eng; pmin=i; }
   }

   // Now recover lowest energy point
   ind[0]=pmin%ngrid[0];
   for(unsigned j=1;j<p.size()-1;++j){
       pmin=(pmin-ind[j-1])/ngrid[j-1];
       ind[j]=pmin%ngrid[j];
   }
   if( p.size()>=2 ) ind[p.size()-1] = (pmin - ind[p.size()-2])/ngrid[p.size()-2];
   for(unsigned j=0;j<p.size();++j) p[j] = min[j] + ind[j]*delr[j];
}
  
}
#endif

