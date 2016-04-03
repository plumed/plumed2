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
#include "tools/IFile.h"
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
  bool isGridPrint() const { return true; }
  void readGridFile( IFile& ifile );
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
  // This ensures that restarting of grids works with NOMEMORY
  mygrid->foundprint=true;
  parse("FMT",fmt); fmt=" "+fmt; parse("FILE",filename); 
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  outputting grid to file named %s with format %s \n",filename.c_str(), fmt.c_str() );

  checkRead();
}

void PrintGrid::performOperationsWithGrid( const bool& from_update ){
  if( !from_update && !single_run ) return ;

  // Read in the old grid and ensure that it is considered
  if( from_update && !mygrid->nomemory ){
     IFile oldf; oldf.link(*this);
     if( oldf.FileExist(filename) ){
         oldf.open(filename);
         readGridFile( oldf );
         oldf.close();
     }
  }

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
     for(unsigned j=0;j<mygrid->getNumberOfQuantities();++j){ 
        ofile.fmtField(fmt); ofile.printField(mygrid->arg_names[mygrid->dimension+j], mygrid->getGridElement( i, j ) );  
     }
     ofile.printField();
  }

  ofile.close();
  // Clear the grid ready for next time
  if( from_update ) mygrid->reset();
}

void PrintGrid::readGridFile( IFile& ifile ){
  // First check that everything in the file is as expected
  std::string min; int nb; std::vector<std::string> fieldnames; ifile.scanFieldList( fieldnames );
  for(unsigned i=0;i<mygrid->dimension;++i){
      ifile.scanField("min_" + mygrid->arg_names[i], min );
      if( min!=mygrid->str_min[i] ) error("minimum of grid in restart file does not match stored minimum");
      ifile.scanField("max_" + mygrid->arg_names[i], min );
      if( min!=mygrid->str_max[i] ) error("maximum of grid in restart file does not match stored minimum");
      ifile.scanField("periodic_" + mygrid->arg_names[i], min );
      if( mygrid->pbc[i] && min!="true" ) error("periodicity of grid in restart file does not match stored periodicity");
      else if( mygrid->pbc[i] && min!="false" ) error("periodicity of grid in restart file does not match stored periodicity");
      ifile.scanField("nbins_" + mygrid->arg_names[i], nb );
      if( mygrid->pbc[i] && nb!=mygrid->nbin[i] ) error("number of bins in restart file does not match stored number of bins");
      else if( !mygrid->pbc[i] && (nb+1)!=mygrid->nbin[i] ) error("number of bins in restart file does not match stored number of bins");
      bool found_field=false;
      for(unsigned j=0;j<fieldnames.size();++j){
         if( fieldnames[j]==mygrid->arg_names[i] ){ found_field=true; break; }
      }
      if( !found_field ) error("missing field " + mygrid->arg_names[i] + " in restart file");
  }
  // Check if all the fields we need are there
  for(unsigned i=0;i<mygrid->arg_names.size();++i){
      if( !ifile.FieldExist( mygrid->arg_names[i] ) ) error("an expected field is missing from the input grid");
  }
  // Now we read in the grid
  unsigned np=0; std::vector<double> indata( mygrid->npoints*mygrid->nper ); double f, old_norm;
  while( ifile.scanField(mygrid->arg_names[0],f) ){
     // Read in coordinates
     for(unsigned i=1;i<mygrid->dimension;++i) ifile.scanField(mygrid->arg_names[i],f);
     // Read header 
     ifile.scanField("normalisation", old_norm );
     for(unsigned i=0;i<mygrid->dimension;++i){
        ifile.scanField("min_" + mygrid->arg_names[i], min );
        ifile.scanField("max_" + mygrid->arg_names[i], min );
        ifile.scanField("periodic_" + mygrid->arg_names[i], min );
        ifile.scanField("nbins_" + mygrid->arg_names[i], nb );
     }
     // Now read data of interest
     for(unsigned i=mygrid->dimension;i<mygrid->arg_names.size();++i){ ifile.scanField(mygrid->arg_names[i],indata[np]); np++; }
     ifile.scanField();
  }
  mygrid->incorporateRestartDataIntoGrid( old_norm, indata );
}

}
}
