/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "GridPrintingBase.h"
#include "core/ActionRegister.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDANALYSIS DUMPGRID
/*
Output the function on the grid to a file with the PLUMED grid format.

\par Examples

*/
//+ENDPLUMEDOC

class DumpGrid : public GridPrintingBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpGrid(const ActionOptions&ao); 
  void printGrid( OFile& ofile ) const ;
};

PLUMED_REGISTER_ACTION(DumpGrid,"DUMPGRID")

void DumpGrid::registerKeywords( Keywords& keys ){
  GridPrintingBase::registerKeywords( keys );
}

DumpGrid::DumpGrid(const ActionOptions&ao):
Action(ao),
GridPrintingBase(ao)
{
  fmt = " " + fmt; checkRead();
}

void DumpGrid::printGrid( OFile& ofile ) const {
  ofile.addConstantField("normalisation");
  for(unsigned i=0;i<ingrid->getDimension();++i){
     ofile.addConstantField("min_" + ingrid->getComponentName(i) );
     ofile.addConstantField("max_" + ingrid->getComponentName(i) );
     ofile.addConstantField("nbins_" + ingrid->getComponentName(i) );
     ofile.addConstantField("periodic_" + ingrid->getComponentName(i) );
  }

  std::vector<double> xx( ingrid->getDimension() );
  std::vector<unsigned> ind( ingrid->getDimension() );
  for(unsigned i=0;i<ingrid->getNumberOfPoints();++i){
     ingrid->getIndices( i, ind );
     if(i>0 && ingrid->getDimension()==2 && ind[ingrid->getDimension()-2]==0) ofile.printf("\n");
     ofile.fmtField(fmt); ofile.printField("normalisation", ingrid->getNorm() );
     for(unsigned j=0;j<ingrid->getDimension();++j){
         ofile.printField("min_" + ingrid->getComponentName(j), ingrid->getMin()[j] );
         ofile.printField("max_" + ingrid->getComponentName(j), ingrid->getMax()[j] );
         ofile.printField("nbins_" + ingrid->getComponentName(j), static_cast<int>(ingrid->getNbin()[j]) );
         if( ingrid->isPeriodic(j) ) ofile.printField("periodic_" + ingrid->getComponentName(j), "true" );
         else          ofile.printField("periodic_" + ingrid->getComponentName(j), "false" );
     }
     // Retrieve and print the grid coordinates
     ingrid->getGridPointCoordinates(i, xx ); 
     for(unsigned j=0;j<ingrid->getDimension();++j){ ofile.fmtField(fmt); ofile.printField(ingrid->getComponentName(j),xx[j]); }
     for(unsigned j=0;j<ingrid->getNumberOfQuantities();++j){
        ofile.fmtField(fmt); ofile.printField(ingrid->arg_names[ingrid->dimension+j], ingrid->getGridElement( i, j ) );
     }
     ofile.printField();
  }
}

}
}
