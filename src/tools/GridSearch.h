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
#include "Grid.h"
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
   Grid* mygrid; 
   Grid* myfgrid;
public:
   GridSearch( const std::vector<double>& mmin, const std::vector<double>& mmax, const std::vector<unsigned>& ng, const std::vector<unsigned>& nfg, FCLASS* funcc ) : 
     myclass_func( funcc ), 
     myfgrid(NULL)
     {
        std::vector<std::string> fake_args( nfg.size() ), gmin( nfg.size() ), gmax( nfg.size() );
        for(unsigned i=0;i<nfg.size();++i){
            Tools::convert(i+1,fake_args[i]); 
            Tools::convert(mmin[i],gmin[i]); 
            Tools::convert(mmax[i],gmax[i]);
        }  
        std::vector<bool> isperiodic( nfg.size(), false );
        mygrid = new Grid("searcher",fake_args,gmin,gmax,ng,true,true,true,isperiodic,gmin,gmax);  
        if( nfg[0]>0 ) myfgrid = new Grid("searcher",fake_args,gmin,gmax,nfg,false,false,true,isperiodic,gmin,gmax);
     }
     ~GridSearch(){ delete mygrid; if(myfgrid) delete myfgrid; }
     bool minimise( std::vector<double>& p, engf_pointer myfunc );
};

template <class FCLASS>
bool GridSearch<FCLASS>::minimise( std::vector<double>& p, engf_pointer myfunc ){
   std::vector<double> der( p.size() ); 
   double initial_eng = (myclass_func->*myfunc)( p, der );

   double emin=(myclass_func->*myfunc)( mygrid->getPoint(0), der );
   mygrid->setValueAndDerivatives( 0, emin, der ); unsigned pmin=0;
   for(unsigned i=1;i<mygrid->getSize();++i){
      double eng = (myclass_func->*myfunc)( mygrid->getPoint(i), der );
      mygrid->setValueAndDerivatives( i, eng, der );
      if( eng<emin ){ emin=eng; pmin=i; }
   }

   if( myfgrid ){
       pmin=0; double emin=mygrid->getValueAndDerivatives( myfgrid->getPoint(0), der );
       for(unsigned i=1;i<myfgrid->getSize();++i){
           double eng = mygrid->getValueAndDerivatives( myfgrid->getPoint(i), der );
           if( eng<emin ){ emin=eng; pmin=i; }
       } 
       double checkEng = (myclass_func->*myfunc)( myfgrid->getPoint(pmin), der );
       if( checkEng<initial_eng ){
           p=myfgrid->getPoint(pmin);
           return true; 
       } else {
           return false;
       }
   }
  
   if( emin<initial_eng ){
      p=mygrid->getPoint(pmin);
      return true;
   } else {
      return false;
   }
}
  
}
#endif

