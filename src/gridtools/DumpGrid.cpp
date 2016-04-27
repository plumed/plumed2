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
#include "GridFunction.h"
#include "AverageOnGrid.h"
#include "tools/IFile.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDANALYSIS DUMPGRID
/*
Output the function on the grid to a file with the PLUMED grid format.

\par Examples

*/
//+ENDPLUMEDOC

class DumpGrid : public ActionWithInputGrid {
private:
  bool printav;
  std::string fmt, filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpGrid(const ActionOptions&ao); 
  void performOperationsWithGrid( const bool& from_update );
  unsigned getNumberOfDerivatives(){ return 0; }
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const {}
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(DumpGrid,"DUMPGRID")

void DumpGrid::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","FILE","density","the file on which to write the grid."); 
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.addFlag("PRINT_AVERAGE",false,"if your input is a \\ref MULTICOLVARDENS output the average in the grid file.");
}

DumpGrid::DumpGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
printav(false),
fmt("%f")
{
  // This ensures that restarting of grids works with memory
  mygrid->foundprint=true;
  GridFunction* mygf = dynamic_cast<GridFunction*> ( mygrid );
  if( mygf ){
      ActionWithInputGrid* ming = dynamic_cast<ActionWithInputGrid*>( mygf->getAction() );
      plumed_assert( ming ); (ming->mygrid)->foundprint=true; 
  }
  AverageOnGrid* myav = dynamic_cast<AverageOnGrid*>( mygrid );
  if( myav ){
     parseFlag("PRINT_AVERAGE",printav);
     if( printav && !mygrid->nomemory && getStride()>0 ) error("cannot use PRINT_AVERAGE flag if you are outputting with a stride");
   }

  parse("FMT",fmt); fmt=" "+fmt; parse("FILE",filename); 
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  outputting grid to file named %s with format %s \n",filename.c_str(), fmt.c_str() );
  checkRead();
}

void DumpGrid::performOperationsWithGrid( const bool& from_update ){
  if( !from_update && getStride()>0 ) return ;

  OFile ofile; ofile.link(*this);
  if ( from_update ) ofile.setBackupString("analysis");
  ofile.open( filename );

  ofile.addConstantField("normalisation");
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
     ofile.fmtField(fmt); ofile.printField("normalisation", mygrid->norm );
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
     if( printav ){
        unsigned nnorm=mygrid->dimension+1; if( mygrid->noderiv ) nnorm=1;
        for(unsigned j=0;j<mygrid->getNumberOfQuantities()-nnorm;++j){
           ofile.fmtField(fmt); ofile.printField(mygrid->arg_names[mygrid->dimension+j], mygrid->getGridElement( i, j ) );
        }
     } else {
        for(unsigned j=0;j<mygrid->getNumberOfQuantities();++j){ 
           ofile.fmtField(fmt); ofile.printField(mygrid->arg_names[mygrid->dimension+j], mygrid->getGridElementForPrint( i, j ) );  
        }
     }
     ofile.printField();
  }
  ofile.close();
}

}
}
