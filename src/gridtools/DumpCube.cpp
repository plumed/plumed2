/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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

Suppose you have calculated the value of a function on a three dimensional grid.
This function might be a \ref HISTOGRAM or it might be a free energy energy surface
that was calculated from this histogram by using \ref CONVERT_TO_FES.  Alternatively,
your function might be a phase-field that was calculated using \ref MULTICOLVARDENS.
Whatever the function is, however, you obviously cannot show it using a typical contour
plotting program such as gnuplot as you have three input variables.

Tools like VMD have nice features for plotting these types of three dimensional functions
but typically you are required to use a Gaussian cube file format to input the data.  This
action thus allows you to output a function evaluated on a grid to a Gaussian cube file format.

\par Examples

The input below can be used to post process a trajectory.  A histogram as a function of the distance
between atoms 1 and 2, the distance between atom 1 and 3 and the angle between the vector connecting atoms 1 and
2 and 1 and 3 is computed using kernel density estimation.  Once all the data contained in the trajectory has been read in and
all the kernels have been added the resulting histogram is output to a file called histoA1.cube.  This file has the
Gaussian cube file format.  The histogram can thus be visualized using tools such as VMD.

\plumedfile
x1: DISTANCE ATOMS=1,2
x2: DISTANCE ATOMS=1,3
x3: ANGLE ATOMS=1,2,3

hA1: HISTOGRAM ARG=x1,x2,x3 GRID_MIN=0.0,0.0,0.0 GRID_MAX=3.0,3.0,3.0 GRID_BIN=10,10,10 BANDWIDTH=1.0,1.0,1.0
DUMPCUBE GRID=hA1 FILE=histoA1.cube
\endplumedfile

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
  void printGrid( OFile& ofile ) const override;
};

PLUMED_REGISTER_ACTION(DumpCube,"DUMPCUBE")

void DumpCube::registerKeywords( Keywords& keys ) {
  GridPrintingBase::registerKeywords( keys );
  keys.add("optional","COMPONENT","if your input is a vector field use this to specify the component of the input vector field for which you wish to output");
}

DumpCube::DumpCube(const ActionOptions&ao):
  Action(ao),
  GridPrintingBase(ao)
{
  fmt = fmt + " ";
  if( ingrid->getType()!="flat" ) error("cannot dump grid of type " + ingrid->getType() + " using DUMPCUBE");
  if( ingrid->getDimension()!=3 ) error("cannot print cube file if grid does not contain three dimensional data");

  if( ingrid->getNumberOfComponents()==1 ) {
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
  for(pp[0]=0; pp[0]<nbin[0]; ++pp[0]) {
    for(pp[1]=0; pp[1]<nbin[1]; ++pp[1]) {
      for(pp[2]=0; pp[2]<nbin[2]; ++pp[2]) {
        ofile.printf(fmt.c_str(), ingrid->getGridElement( ingrid->getIndex(pp), mycomp ) );
        if(pp[2]%6==5) ofile.printf("\n");
      }
      ofile.printf("\n");
    }
  }
}

}
}
