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
#include "core/ActionRegister.h"
#include "GridPrintingBase.h"
#include "tools/OFile.h"

//+PLUMEDOC GRIDANALYSIS DUMPCUBE
/*
Output a three dimensional grid using the Gaussian cube file format.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class DumpCube : public GridPrintingBase {
private:
  unsigned mycomp;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpCube(const ActionOptions&ao); 
  void printGrid( OFile& ofile ) const ;
};

PLUMED_REGISTER_ACTION(DumpCube,"DUMPCUBE")

void DumpCube::registerKeywords( Keywords& keys ){
  GridPrintingBase::registerKeywords( keys );
  keys.add("optional","COMPONENT","if your input is a vector field use this to specifiy the component of the input vector field for which you wish to output");
}

DumpCube::DumpCube(const ActionOptions&ao):
Action(ao),
GridPrintingBase(ao)
{
  fmt = fmt + " ";
  if( ingrid->getDimension()!=3 ) error("cannot print cube file if grid does not contain three dimensional data");

  if( ingrid->getNumberOfComponents()==1 ){
     mycomp=0;
  } else {
     int tcomp=-1; parse("COMPONENT",tcomp);
     if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
     mycomp=tcomp*(1+ingrid->getDimension()); if( ingrid->noDerivatives() ) mycomp=tcomp;
     log.printf("  using %dth component of grid \n",tcomp );
  }

  checkRead();
}

void DumpCube::printGrid( OFile& ofile ) const {
  double lunit = ingrid->getCubeUnits();

  ofile.printf("PLUMED CUBE FILE\n");
  ofile.printf("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  // Number of atoms followed by position of origin (origin set so that center of grid is in center of cell)
  std::string ostr = "%d " + fmt + fmt + fmt + "\n";
  ofile.printf(ostr.c_str(),1,-0.5*lunit*ingrid->getGridExtent(0),-0.5*lunit*ingrid->getGridExtent(1),-0.5*lunit*ingrid->getGridExtent(2));
  ofile.printf(ostr.c_str(),ingrid->getNbin()[0],lunit*ingrid->getGridSpacing()[0],0.0,0.0);  // Number of bins in each direction followed by 
  ofile.printf(ostr.c_str(),ingrid->getNbin()[1],0.0,lunit*ingrid->getGridSpacing()[1],0.0);  // shape of voxel
  ofile.printf(ostr.c_str(),ingrid->getNbin()[2],0.0,0.0,lunit*ingrid->getGridSpacing()[2]);
  ofile.printf(ostr.c_str(),1,0.0,0.0,0.0); // Fake atom otherwise VMD doesn't work
  std::vector<unsigned> pp(3); std::vector<unsigned> nbin( ingrid->getNbin() );
  for(pp[0]=0;pp[0]<nbin[0];++pp[0]){
      for(pp[1]=0;pp[1]<nbin[1];++pp[1]){
          for(pp[2]=0;pp[2]<nbin[2];++pp[2]){
              ofile.printf(fmt.c_str(), ingrid->getGridElement( ingrid->getIndex(pp), mycomp ) );
              if(pp[2]%6==5) ofile.printf("\n");
          }
          ofile.printf("\n");
     }
  }
}

}
}
