/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "ActionWithGrid.h"
#include "tools/OFile.h"

namespace PLMD {
namespace gridtools {

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
  LABEL=hh
... HISTOGRAM

DUMPGRID GRID=hh FILE=histo STRIDE=100000
\endplumedfile

*/
//+ENDPLUMEDOC

class DumpGrid :
  public ActionWithArguments,
  public ActionPilot {
private:
  std::string fmt, filename;
  bool onefile, xyzfile;
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpGrid(const ActionOptions&ao);
  ~DumpGrid() {}
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpGrid,"DUMPCUBE")
PLUMED_REGISTER_ACTION(DumpGrid,"DUMPGRID")

void DumpGrid::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.use("ARG");
  keys.add("optional","GRID","the grid you would like to print (can also use ARG for specifying what is being printed)");
  keys.add("compulsory","STRIDE","0","the frequency with which the grid should be output to the file.  Default of zero means dump at end of calculation");
  keys.add("compulsory","FILE","density","the file on which to write the grid.");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.addFlag("PRINT_XYZ",false,"output coordinates on fibonacci grid to xyz file");
  keys.addFlag("PRINT_ONE_FILE",false,"output grids one after the other in a single file");
}

DumpGrid::DumpGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  fmt("%f") {
  if( getNumberOfArguments()==0 ) {
    std::vector<Value*> grids;
    parseArgumentList("GRID",grids);
    requestArguments(grids);
  }
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument");
  }
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) {
    error("input should be a grid");
  }
  if( getName()=="DUMPCUBE" && getPntrToArgument(0)->getRank()!=3 ) {
    error("input should be a three dimensional grid");
  }
  parse("FILE",filename);
  if(filename.length()==0) {
    error("name out output file was not specified");
  }
  ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  if( !ag ) {
    error( getPntrToArgument(0)->getName() + " is not grid");
  }

  log.printf("  outputting grid with label %s to file named %s",getPntrToArgument(0)->getName().c_str(), filename.c_str() );
  parse("FMT",fmt);
  log.printf(" with format %s \n", fmt.c_str() );
  if( getName()=="DUMPGRID" ) {
    fmt = " " + fmt;
  } else if( getName()=="DUMPCUBE" ) {
    fmt = fmt + " ";
  }
  parseFlag("PRINT_ONE_FILE", onefile);
  parseFlag("PRINT_XYZ",xyzfile);
  if( xyzfile ) {
    if( getName()=="DUMPCUBE" ) {
      error("PRINT_XYZ flag not compatible with DUMPCUBE");
    }
    if( ag->getGridCoordinatesObject().getGridType()=="flat" ) {
      error("can only use PRINT_XYZ option for fibonacci grids");
    }
    log.printf("  outputting grid to xyzfile\n");
  }
  if( onefile ) {
    log.printf("  printing all grids on a single file \n");
  } else {
    log.printf("  printing all grids on separate files \n");
  }
}

void DumpGrid::update() {
  OFile ofile;
  ofile.link(*this);
  if( onefile ) {
    ofile.enforceRestart();
    ofile.open( filename );
  } else {
    ofile.setBackupString("analysis");
    ofile.open( filename );
  }

  ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  const GridCoordinatesObject & mygrid = ag->getGridCoordinatesObject();

  if( getName()=="DUMPCUBE" ) {
    double lunit=1.0;
    std::vector<std::string> argnames( ag->getGridCoordinateNames() );

    ofile.printf("PLUMED CUBE FILE\n");
    ofile.printf("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
    // Number of atoms followed by position of origin (origin set so that center of grid is in center of cell)
    bool isdists=true;
    std::vector<double> extent(3);
    std::vector<std::string> min( mygrid.getMin() ), max( mygrid.getMax() );
    for(unsigned j=0; j<3; ++j) {
      if( argnames[j].find(".")!=std::string::npos ) {
        std::size_t dot = argnames[j].find(".");
        std::string name = argnames[j].substr(dot+1);
        if( name!="x" && name!="y" && name!="z" ) {
          isdists=false;
        }
      } else {
        isdists=false;
      }

      double mind, maxd;
      Tools::convert( min[j], mind );
      Tools::convert( max[j], maxd );
      if( mygrid.isPeriodic(j) ) {
        extent[j]=maxd-mind;
      } else {
        extent[j]=maxd-mind+mygrid.getGridSpacing()[j];
      }
    }
    if( isdists ) {
      if( plumed.usingNaturalUnits() ) {
        lunit = 1.0/0.5292;
      } else {
        lunit = plumed.getUnits().getLength()/.05929;
      }
    }
    std::string ostr = "%d " + fmt + fmt + fmt + "\n";
    Value* gval=getPntrToArgument(0);
    ofile.printf(ostr.c_str(),1,-0.5*lunit*extent[0],-0.5*lunit*extent[1],-0.5*lunit*extent[2]);
    ofile.printf(ostr.c_str(),mygrid.getNbin(true)[0],lunit*mygrid.getGridSpacing()[0],0.0,0.0);  // Number of bins in each direction followed by
    ofile.printf(ostr.c_str(),mygrid.getNbin(true)[1],0.0,lunit*mygrid.getGridSpacing()[1],0.0);  // shape of voxel
    ofile.printf(ostr.c_str(),mygrid.getNbin(true)[2],0.0,0.0,lunit*mygrid.getGridSpacing()[2]);
    ofile.printf(ostr.c_str(),1,0.0,0.0,0.0); // Fake atom otherwise VMD doesn't work
    std::vector<unsigned> pp(3);
    std::vector<unsigned> nbin( mygrid.getNbin(true) );
    for(pp[0]=0; pp[0]<nbin[0]; ++pp[0]) {
      for(pp[1]=0; pp[1]<nbin[1]; ++pp[1]) {
        for(pp[2]=0; pp[2]<nbin[2]; ++pp[2]) {
          unsigned ival=pp[pp.size()-1];
          for(unsigned i=pp.size()-1; i>0; --i) {
            ival=ival*nbin[i-1]+pp[i-1];
          }
          ofile.printf(fmt.c_str(), gval->get(ival) );
          if(pp[2]%6==5) {
            ofile.printf("\n");
          }
        }
        ofile.printf("\n");
      }
    }
  } else if( xyzfile ) {
    std::vector<double> coords(3);
    Value* myarg = getPntrToArgument(0);
    unsigned nvals = myarg->getNumberOfValues();
    ofile.printf("%d\n\n", nvals);
    for(unsigned i=0; i<nvals; ++i) {
      double val = myarg->get(i);
      mygrid.getGridPointCoordinates( i, coords );
      ofile.printf("X");
      for(unsigned j=0; j<3; ++j) {
        ofile.printf((" " + fmt).c_str(), val*coords[j] );
      }
      ofile.printf("\n");
    }
  } else {
    Value* gval=getPntrToArgument(0);
    std::vector<double> xx( gval->getRank() );
    std::vector<unsigned> ind( gval->getRank() );
    std::vector<std::string> argn( ag->getGridCoordinateNames() );
    if( mygrid.getGridType()=="fibonacci" ) {
      ofile.addConstantField("nbins");
    } else {
      plumed_assert( mygrid.getGridType()=="flat" );
      for(unsigned i=0; i<gval->getRank(); ++i) {
        ofile.addConstantField("min_" + argn[i] );
        ofile.addConstantField("max_" + argn[i] );
        ofile.addConstantField("nbins_" + argn[i] );
        ofile.addConstantField("periodic_" + argn[i] );
      }
    }

    for(unsigned i=0; i<gval->getNumberOfValues(); ++i) {
      // Retrieve and print the grid coordinates
      mygrid.getGridPointCoordinates( i, ind, xx );
      if(i>0 && gval->getRank()==2 && ind[gval->getRank()-2]==0) {
        ofile.printf("\n");
      }
      ofile.fmtField(fmt);
      if( mygrid.getGridType()=="fibonacci" ) {
        ofile.printField("nbins", static_cast<int>(gval->getNumberOfValues()) );
      } else {
        for(unsigned j=0; j<gval->getRank(); ++j) {
          ofile.printField("min_" + argn[j], mygrid.getMin()[j] );
          ofile.printField("max_" + argn[j], mygrid.getMax()[j] );
          ofile.printField("nbins_" + argn[j], static_cast<int>(mygrid.getNbin(false)[j]) );
          if( mygrid.isPeriodic(j) ) {
            ofile.printField("periodic_" + argn[j], "true" );
          } else {
            ofile.printField("periodic_" + argn[j], "false" );
          }
        }
      }
      // Print the grid coordinates
      for(unsigned j=0; j<gval->getRank(); ++j) {
        ofile.fmtField(fmt);
        ofile.printField(argn[j],xx[j]);
      }
      // Print value
      ofile.fmtField(fmt);
      ofile.printField( gval->getName(), gval->get(i) );
      // Print the derivatives
      for(unsigned j=0; j<gval->getRank(); ++j) {
        ofile.fmtField(fmt);
        ofile.printField( "d" + gval->getName() + "_" + argn[j], gval->getGridDerivative(i,j) );
      }
      ofile.printField();
    }
  }
}


}
}
