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
#include "core/ActionRegister.h"
#include "ActionWithInputGrid.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

class PrintCube : public ActionWithInputGrid {
private:
  std::string fmt, filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit PrintCube(const ActionOptions&ao); 
  void performOperationsWithGrid( const bool& from_update );
};

PLUMED_REGISTER_ACTION(PrintCube,"PRINT_CUBE")

void PrintCube::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("compulsory","FILE","density","the file on which to write the grid."); 
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

PrintCube::PrintCube(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao),
fmt("%f")
{
  if( mygrid->getDimension()!=3 ) error("cannot print cube file if grid does not contain three dimensional data");
  if( mygrid->getNumberOfComponents()!=1 ) error("cannot print cube file if data on grid is a vector field");

  parse("FMT",fmt); fmt=" "+fmt; parse("FILE",filename); 
  if(filename.length()==0) error("name out output file was not specified");
  log.printf("  outputting grid to file named %s with format %s \n",filename.c_str(), fmt.c_str() );

  checkRead();
}

void PrintCube::performOperationsWithGrid( const bool& from_update ){
  if( !from_update && !single_run ) return ;
  OFile ofile; ofile.link(*this);
  if( from_update ) ofile.setBackupString("analysis");
  ofile.open( filename );

  double lunit = mygrid->getCubeUnits();

  ofile.printf("PLUMED CUBE FILE\n");
  ofile.printf("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  // Number of atoms followed by position of origin (origin set so that center of grid is in center of cell)
  ofile.printf("%d %f %f %f\n",1,-0.5*lunit*mygrid->getGridExtent(0),-0.5*lunit*mygrid->getGridExtent(1),-0.5*lunit*mygrid->getGridExtent(2));
  ofile.printf("%u %f %f %f\n",mygrid->getNbin()[0],lunit*mygrid->getGridSpacing()[0],0.0,0.0);  // Number of bins in each direction followed by 
  ofile.printf("%u %f %f %f\n",mygrid->getNbin()[1],0.0,lunit*mygrid->getGridSpacing()[1],0.0);  // shape of voxel
  ofile.printf("%u %f %f %f\n",mygrid->getNbin()[2],0.0,0.0,lunit*mygrid->getGridSpacing()[2]);
  ofile.printf("%u %f %f %f\n",1,0.0,0.0,0.0); // Fake atom otherwise VMD doesn't work
  std::vector<unsigned> pp(3); std::vector<unsigned> nbin( mygrid->getNbin() );
  for(pp[0]=0;pp[0]<nbin[0];++pp[0]){
      for(pp[1]=0;pp[1]<nbin[1];++pp[1]){
          for(pp[2]=0;pp[2]<nbin[2];++pp[2]){
              ofile.printf("%f ",getGridElement( pp, 0 ) );
              if(pp[2]%6==5) ofile.printf("\n");
          }
          ofile.printf("\n");
     }
  }

  ofile.close();
  // Clear the grid ready for next time
  if( from_update ) mygrid->clear();
}

}
}
