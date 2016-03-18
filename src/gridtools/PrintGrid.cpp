/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "ActionWithInputGrid.h"
#include "core/ActionRegister.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

class PrintGrid : public ActionWithInputGrid {
private:
  std::string fmt, filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit PrintGrid(const ActionOptions&ao); 
  void performOperationsWithGrid( const bool& from_update );
  unsigned getNumberOfDerivatives(){ return 0; }
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const {}
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(PrintGrid,"PRINT_GRID")

void PrintGrid::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","FILE","density","the file on which to write the grid."); 
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

PrintGrid::PrintGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
fmt("%f")
{
  parse("FMT",fmt); fmt=" "+fmt; parse("FILE",filename); 
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  outputting grid to file named %s with format %s \n",filename.c_str(), fmt.c_str() );

  checkRead();
}

void PrintGrid::performOperationsWithGrid( const bool& from_update ){
  if( !from_update && !single_run ) return ;
  OFile ofile; ofile.link(*this);
  if( from_update ) ofile.setBackupString("analysis");
  ofile.open( filename );

  ofile.addConstantField("normalisation");
  for(unsigned i=0;i<mygrid->getDimension();++i){
     ofile.addConstantField("min_" + mygrid->getComponentName(i) );
     ofile.addConstantField("max_" + mygrid->getComponentName(i) );
     ofile.addConstantField("nbins_" + mygrid->getComponentName(i) );
     ofile.addConstantField("periodic_" + mygrid->getComponentName(i) );
  }

  double norm = mygrid->getNorm();
  std::vector<double> xx( mygrid->getDimension() );
  std::vector<unsigned> ind( mygrid->getDimension() );
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
     mygrid->getIndices( i, ind );
     if(i>0 && mygrid->getDimension()==2 && ind[mygrid->getDimension()-2]==0) ofile.printf("\n");
     ofile.fmtField(fmt); ofile.printField("normalisation", norm );
     for(unsigned j=0;j<mygrid->getDimension();++j){
         ofile.printField("min_" + mygrid->getComponentName(j), mygrid->getMin()[j] );
         ofile.printField("max_" + mygrid->getComponentName(j), mygrid->getMax()[j] );
         ofile.printField("nbins_" + mygrid->getComponentName(j), static_cast<int>(mygrid->getNbin()[j]) );
         if( mygrid->isPeriodic(j) ) ofile.printField("periodic_" + mygrid->getComponentName(j), "true" );
         else          ofile.printField("periodic_" + mygrid->getComponentName(j), "false" );
     }
     // Retrieve and print the grid coordinates
     mygrid->getGridPointCoordinates(i, xx ); 
     for(unsigned j=0;j<mygrid->getDimension();++j){ ofile.fmtField(fmt); ofile.printField(mygrid->getComponentName(j),xx[j]); }
     for(unsigned j=0;j<mygrid->getNumberOfQuantities();++j){ 
        ofile.fmtField(fmt); ofile.printField(mygrid->getComponentName(mygrid->getDimension()+j), getGridElement( i, j ) ); 
     }
     ofile.printField();
  }

  ofile.close();
  // Clear the grid ready for next time
  if( from_update ) mygrid->reset();
}

}
}
