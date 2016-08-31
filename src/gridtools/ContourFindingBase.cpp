/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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
#include "ContourFindingBase.h"
#include "vesselbase/StoreDataVessel.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace gridtools {

void ContourFindingBase::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","CONTOUR","the value we would like to draw the contour at in the space");
  keys.add("compulsory","FILE","file on which to output coordinates");
  keys.add("compulsory","UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
  keys.remove("KERNEL"); keys.remove("BANDWIDTH");
}

ContourFindingBase::ContourFindingBase(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
mymin(this),
mydata(NULL)
{
  if( ingrid->noDerivatives() ) error("cannot find contours if input grid has no derivatives");
  parse("CONTOUR",contour); 
  log.printf("  calculating dividing surface along which function equals %f \n", contour);

  if( keywords.exists("FILE") ){
      std::string file; parse("FILE",file);
      if(file.length()==0 && keywords.style("FILE","compulsory") ) error("name out output file was not specified");
      else if( file.length()>0 ){
         std::string type=Tools::extension(file);
         log<<"  file name "<<file<<"\n";
         if(type!="xyz") error("can only print xyz file type with contour finding");

         fmt_xyz="%f";
         std::string precision; parse("PRECISION",precision);
         if(precision.length()>0){
           int p; Tools::convert(precision,p);
           log<<"  with precision "<<p<<"\n";
           std::string a,b;
           Tools::convert(p+5,a);
           Tools::convert(p,b);
           fmt_xyz="%"+a+"."+b+"f";
         }

         std::string unitname; parse("UNITS",unitname);
         if(unitname!="PLUMED"){
           Units myunit; myunit.setLength(unitname);
           lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
         }
         else lenunit=1.0; 

         of.link(*this); of.open(file); 

         // Now create a store data vessel to hold the points
         mydata=buildDataStashes( NULL );
      }
  }
}

void ContourFindingBase::finishAveraging(){
  // And this looks after the output
  if( mydata ){
      std::vector<double> point( 1 + ingrid->getDimension() );
      of.printf("%u\n",mydata->getNumberOfStoredValues());  
      of.printf("Points found on isocontour\n");
      for(unsigned i=0;i<mydata->getNumberOfStoredValues();++i){
          mydata->retrieveSequentialValue( i, false, point ); of.printf("X");
          for(unsigned j=0;j<ingrid->getDimension();++j) of.printf( (" " + fmt_xyz).c_str(), lenunit*point[1+j] );
          of.printf("\n");
      }   
  }
}

}
}
