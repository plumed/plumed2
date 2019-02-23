/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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

PLUMED provides a number of actions that calculate the values of functions on grids.
For instance, whenever you calculate a free energy as a function of a collective variable using
\ref HISTOGRAM and \ref CONVERT_TO_FES you will generally want to output the value of the free energy at a number of points on a
discrete grid that covers the CV space uniformly.  Alternatively you may want to calculate
what value some symmetry function takes at different points inside your simulation cell using \ref MULTICOLVARDENS.

This action allows you to output these functions calculated on a grid using a format that can be read in using gnuplot
and other such plotting programs.  The file output using this action will have a header that contains some essential
information about the function plotted and that looks something like this:

\verbatim
#! FIELDS x y hA1 dhA1_x dhA1_x
#! SET normalisation    2.0000
#! SET min_x 0.0
#! SET max_x 3.0
#! SET nbins_x  100
#! SET periodic_x false
#! SET min_y 0.0
#! SET max_y 3.0
#! SET nbins_y  100
#! SET periodic_y false
\endverbatim

The header shown here tells us that we have grid showing the values that a function with two arguments x and y
takes at various points in our cell.  The lines beneath the first line then tell us a little bit about these two
input arguments.

The remaining lines of the file give us information on the positions of our grid points and the value the function and
its partial derivatives with respect to x and y.  If the header is as above a list of values of the function that have
x=0 and 100 values of y between 0.0 and 3.0 will be provided.  This block of data will be followed with a blank line.
There will then be a second block of values which will all have been evaluated the same value of x and all possible values
for y.  This block is then followed by a blank line again and this pattern continues until all points of the grid have been covered.

\par Examples

The following input monitors two torsional angles during a simulation
and outputs a continuous histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a discrete histogram as a function of them at the end of the simulation.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  KERNEL=DISCRETE
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs the histogram accumulated thus far every 100000 steps.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

The following input monitors two torsional angles during a simulation
and outputs a separate histogram for each 100000 steps worth of trajectory.
Notice how the CLEAR keyword is used here and how it is not used in the
previous example.

\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
HISTOGRAM ...
  ARG=r1,r2 CLEAR=100000
  GRID_MIN=-3.14,-3.14
  GRID_MAX=3.14,3.14
  GRID_BIN=200,200
  BANDWIDTH=0.05,0.05
  GRID_WFILE=histo
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

*/
//+ENDPLUMEDOC

class DumpGrid : public GridPrintingBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpGrid(const ActionOptions&ao);
  void printGrid( OFile& ofile ) const ;
};

PLUMED_REGISTER_ACTION(DumpGrid,"DUMPGRID")

void DumpGrid::registerKeywords( Keywords& keys ) {
  GridPrintingBase::registerKeywords( keys );
}

DumpGrid::DumpGrid(const ActionOptions&ao):
  Action(ao),
  GridPrintingBase(ao)
{
  if( ingrid->getType()!="flat" ) error("cannot dump grid of type " + ingrid->getType() + " using DUMPGRID");
  fmt = " " + fmt; checkRead();
}

void DumpGrid::printGrid( OFile& ofile ) const {
  ofile.addConstantField("normalisation");
  for(unsigned i=0; i<ingrid->getDimension(); ++i) {
    ofile.addConstantField("min_" + ingrid->getComponentName(i) );
    ofile.addConstantField("max_" + ingrid->getComponentName(i) );
    ofile.addConstantField("nbins_" + ingrid->getComponentName(i) );
    ofile.addConstantField("periodic_" + ingrid->getComponentName(i) );
  }

  std::vector<double> xx( ingrid->getDimension() );
  std::vector<unsigned> ind( ingrid->getDimension() );
  for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) {
    ingrid->getIndices( i, ind );
    if(i>0 && ingrid->getDimension()==2 && ind[ingrid->getDimension()-2]==0) ofile.printf("\n");
    ofile.fmtField(fmt); ofile.printField("normalisation", ingrid->getNorm() );
    for(unsigned j=0; j<ingrid->getDimension(); ++j) {
      ofile.printField("min_" + ingrid->getComponentName(j), ingrid->getMin()[j] );
      ofile.printField("max_" + ingrid->getComponentName(j), ingrid->getMax()[j] );
      ofile.printField("nbins_" + ingrid->getComponentName(j), static_cast<int>(ingrid->getNbin()[j]) );
      if( ingrid->isPeriodic(j) ) ofile.printField("periodic_" + ingrid->getComponentName(j), "true" );
      else          ofile.printField("periodic_" + ingrid->getComponentName(j), "false" );
    }
    // Retrieve and print the grid coordinates
    ingrid->getGridPointCoordinates(i, xx );
    for(unsigned j=0; j<ingrid->getDimension(); ++j) { ofile.fmtField(fmt); ofile.printField(ingrid->getComponentName(j),xx[j]); }
    for(unsigned j=0; j<ingrid->getNumberOfQuantities(); ++j) {
      ofile.fmtField(fmt); ofile.printField(ingrid->arg_names[ingrid->dimension+j], ingrid->getGridElement( i, j ) );
    }
    ofile.printField();
  }
}

}
}
