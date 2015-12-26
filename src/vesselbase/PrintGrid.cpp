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
#include "core/ActionPilot.h"
#include "ActionWithInputVessel.h"
#include "GridVessel.h"
#include "core/ActionRegister.h"
#include "tools/OFile.h"

namespace PLMD {
namespace vesselbase {

class PrintGrid : 
public ActionWithInputVessel,
public ActionPilot
{
private:
  bool single_run;
  GridVessel* mygrid;
  std::string fmt, filename;
  void printGrid( OFile& ofile );
public:
  static void registerKeywords( Keywords& keys );
  explicit PrintGrid(const ActionOptions&ao); 
  void calculate(){}
  void apply(){}
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(PrintGrid,"PRINT_GRID")

void PrintGrid::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.use("GRID"); keys.remove("DATA");
  keys.add("compulsory","STRIDE","the frequency with which to output the grid");
  keys.add("compulsory","FILE","density","the file on which to write the grid."); 
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.addFlag("USE_ALL_DATA",false,"use the data from the entire trajectory to perform the analysis");
}

PrintGrid::PrintGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputVessel(ao),
ActionPilot(ao),
mygrid(NULL),
fmt("%f")
{
  readArgument( "grid" );
  mygrid=dynamic_cast<GridVessel*>( getPntrToArgument() );

  parse("FMT",fmt); fmt=" "+fmt;
  parse("FILE",filename); parseFlag("USE_ALL_DATA",single_run);
  if(filename.length()==0) error("name out output file was not specified");

  if( !single_run ) log.printf("  outputting grid every %u steps to file named %s \n", getStride(), filename.c_str() );
  else log.printf("  outputting grid to file named %s with format %s \n",filename.c_str(), fmt.c_str() );

  checkRead();
}

void PrintGrid::update(){
  OFile gridfile; gridfile.link(*this); 
  gridfile.setBackupString("analysis");
  gridfile.open( filename );
  printGrid( gridfile );
  gridfile.close();

  // Clear the grid ready for next time
  mygrid->clear();
}

void PrintGrid::runFinalJobs(){
  if( !single_run ) return ;
  OFile gridfile; gridfile.link(*this); 
  gridfile.open( filename );
  printGrid( gridfile );
  gridfile.close();
}

void PrintGrid::printGrid( OFile& ofile ){
  for(unsigned i=0;i<mygrid->getDimension();++i){
     ofile.addConstantField("min_" + mygrid->getComponentName(i) );
     ofile.addConstantField("max_" + mygrid->getComponentName(i) );
     ofile.addConstantField("nbins_" + mygrid->getComponentName(i) );
     ofile.addConstantField("periodic_" + mygrid->getComponentName(i) );
  }

  std::vector<double> xx( mygrid->getDimension() );
  std::vector<unsigned> ind( mygrid->getDimension() );
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
     mygrid->getIndices( i, ind );
     if(i>0 && mygrid->getDimension()==2 && ind[mygrid->getDimension()-2]==0) ofile.printf("\n");
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
        ofile.fmtField(fmt); ofile.printField(mygrid->getComponentName(mygrid->getDimension()+j), mygrid->getGridElement( i, j ) ); 
     }
     ofile.printField();
  }
}

}
}
