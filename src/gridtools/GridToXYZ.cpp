/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "tools/OFile.h"
#include "core/Atoms.h"

namespace PLMD {
namespace gridtools {

//+PLUMEDOC GRIDANALYSIS GRID_TO_XYZ
/*
Output the function on the grid to an xyz file

\par Examples

*/
//+ENDPLUMEDOC

class GridToXYZ : public GridPrintingBase {
private:
  double lenunit;
  unsigned mycomp;
public:
  static void registerKeywords( Keywords& keys );
  explicit GridToXYZ(const ActionOptions&ao);
  void printGrid( OFile& ofile ) const override;
};

PLUMED_REGISTER_ACTION(GridToXYZ,"GRID_TO_XYZ")

void GridToXYZ::registerKeywords( Keywords& keys ) {
  GridPrintingBase::registerKeywords( keys );
  keys.add("optional","COMPONENT","if your input is a vector field use this to specify the component of the input vector field for which you wish to output");
  keys.add("compulsory","UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
  keys.remove("FMT");
}

GridToXYZ::GridToXYZ(const ActionOptions&ao):
  Action(ao),
  GridPrintingBase(ao)
{
  if( ingrid->getDimension()!=3 ) error("cannot print grid xyz file if grid does not contain three dimensional data");
  fmt = " " + fmt;

  if( ingrid->getNumberOfComponents()==1 ) {
    mycomp=0;
  } else {
    int tcomp=-1; parse("COMPONENT",tcomp);
    if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
    mycomp=tcomp*(1+ingrid->getDimension()); if( ingrid->noDerivatives() ) mycomp=tcomp;
    log.printf("  using %dth component of grid \n",tcomp );
  }
  fmt="%f";
  std::string precision; parse("PRECISION",precision);
  if(precision.length()>0) {
    int p; Tools::convert(precision,p);
    log<<"  with precision "<<p<<"\n";
    std::string a,b;
    Tools::convert(p+5,a);
    Tools::convert(p,b);
    fmt="%"+a+"."+b+"f";
  }
  std::string unitname; parse("UNITS",unitname);
  if(unitname!="PLUMED") {
    Units myunit; myunit.setLength(unitname);
    lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
  }
  else lenunit=1.0;
  checkRead();
}

void GridToXYZ::printGrid( OFile& ofile ) const {
  std::vector<double> point( 3 );
  ofile.printf("%u\n",ingrid->getNumberOfPoints());
  ofile.printf("Grid converted to xyz file \n");
  for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) {
    ingrid->getGridPointCoordinates( i, point );
    ofile.printf("X");
    double val;
    if( ingrid->getType()=="flat" ) val=1.0;
    else val=ingrid->getGridElement( i, 0 );
    for(unsigned j=0; j<3; ++j) { ofile.printf( (" " + fmt).c_str(), val*lenunit*point[j] ); }
    ofile.printf("\n");
  }
}

}
}
