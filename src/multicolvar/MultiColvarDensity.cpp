/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "tools/Pbc.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include "core/ActionSet.h"
#include "MultiColvarBase.h"
#include "gridtools/ActionWithGrid.h"

using namespace std;

namespace PLMD
{
namespace multicolvar {

//+PLUMEDOC GRIDCALC MULTICOLVARDENS
/*
Evaluate the average value of a multicolvar on a grid.

This keyword allows one to construct a phase field representation for a symmetry function from
an atomistic description.  If each atom has an associated order parameter, \f$\phi_i\f$ then a
smooth phase field function \f$\phi(r)\f$ can be computed using:

\f[
\phi(\mathbf{r}) = \frac{\sum_i K(\mathbf{r}-\mathbf{r}_i) \phi_i }{ \sum_i K(\mathbf{r} - \mathbf{r}_i )}
\f]

where \f$\mathbf{r}_i\f$ is the position of atom \f$i\f$, the sums run over all the atoms input
and \f$K(\mathbf{r} - \mathbf{r}_i)\f$ is one of the \ref kernelfunctions implemented in plumed.
This action calculates the above function on a grid, which can then be used in the input to further
actions.

\par Examples

The following example shows perhaps the simplest way in which this action can be used.  The following
input computes the density of atoms at each point on the grid and outputs this quantity to a file.  In
other words this input instructs plumed to calculate \f$\rho(\mathbf{r}) = \sum_i K(\mathbf{r} - \mathbf{r}_i )\f$

\plumedfile
dens: DENSITY SPECIES=1-100
grid: MULTICOLVARDENS DATA=dens ORIGIN=1 DIR=xyz NBINS=100,100,100 BANDWIDTH=0.05,0.05,0.05 STRIDE=1
DUMPGRID GRID=grid STRIDE=500 FILE=density
\endplumedfile

In the above example density is added to the grid on every step.  The PRINT_GRID instruction thus tells PLUMED to
output the average density at each point on the grid every 500 steps of simulation.  Notice that the that grid output
on step 1000 is an average over all 1000 frames of the trajectory.  If you would like to analyze these two blocks
of data separately you must use the CLEAR flag.

This second example computes an order parameter (in this case \ref FCCUBIC) and constructs a phase field model
for this order parameter using the equation above.

\plumedfile
fcc: FCCUBIC SPECIES=1-5184 SWITCH={CUBIC D_0=1.2 D_MAX=1.5} ALPHA=27
dens: MULTICOLVARDENS DATA=fcc ORIGIN=1 DIR=xyz NBINS=14,14,28 BANDWIDTH=1.0,1.0,1.0 STRIDE=1 CLEAR=1
DUMPCUBE GRID=dens STRIDE=1 FILE=dens.cube
\endplumedfile

In this example the phase field model is computed and output to a file on every step of the simulation.  Furthermore,
because the CLEAR=1 keyword is set on the MULTICOLVARDENS line each Gaussian cube file output is a phase field
model for a particular trajectory frame. The average value accumulated thus far is cleared at the start of every single
timestep and there is no averaging over trajectory frames in this case.

*/
//+ENDPLUMEDOC

class MultiColvarDensity : public gridtools::ActionWithGrid {
  bool fractional;
  MultiColvarBase* mycolv;
  std::vector<unsigned> nbins;
  std::vector<double> gspacing;
  std::vector<bool> confined;
  std::vector<double> cmin, cmax;
  vesselbase::StoreDataVessel* stash;
  Vector origin;
  std::vector<unsigned> directions;
public:
  explicit MultiColvarDensity(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfQuantities() const ;
  bool isPeriodic() { return false; }
  void clearAverage();
  void prepareForAveraging();
  void compute( const unsigned&, MultiValue& ) const ;
  void apply() {}
};

PLUMED_REGISTER_ACTION(MultiColvarDensity,"MULTICOLVARDENS")

void MultiColvarDensity::registerKeywords( Keywords& keys ) {
  gridtools::ActionWithGrid::registerKeywords( keys );
  keys.add("atoms","ORIGIN","we will use the position of this atom as the origin");
  keys.add("compulsory","DATA","the multicolvar which you would like to calculate the density profile for");
  keys.add("compulsory","DIR","the direction in which to calculate the density profile");
  keys.add("optional","NBINS","the number of bins to use to represent the density profile");
  keys.add("optional","SPACING","the approximate grid spacing (to be used as an alternative or together with NBINS)");
  keys.addFlag("FRACTIONAL",false,"use fractional coordinates for the various axes");
  keys.addFlag("XREDUCED",false,"limit the calculation of the density/average to a portion of the z-axis only");
  keys.add("optional","XLOWER","this is required if you are using XREDUCED. It specifies the lower bound for the region of the x-axis that for "
           "which you are calculating the density/average");
  keys.add("optional","XUPPER","this is required if you are using XREDUCED. It specifies the upper bound for the region of the x-axis that for "
           "which you are calculating the density/average");
  keys.addFlag("YREDUCED",false,"limit the calculation of the density/average to a portion of the y-axis only");
  keys.add("optional","YLOWER","this is required if you are using YREDUCED. It specifies the lower bound for the region of the y-axis that for "
           "which you are calculating the density/average");
  keys.add("optional","YUPPER","this is required if you are using YREDUCED. It specifies the upper bound for the region of the y-axis that for "
           "which you are calculating the density/average");
  keys.addFlag("ZREDUCED",false,"limit the calculation of the density/average to a portion of the z-axis only");
  keys.add("optional","ZLOWER","this is required if you are using ZREDUCED. It specifies the lower bound for the region of the z-axis that for "
           "which you are calculating the density/average");
  keys.add("optional","ZUPPER","this is required if you are using ZREDUCED. It specifies the upper bound for the region of the z-axis that for "
           "which you are calculating the density/average");
}

MultiColvarDensity::MultiColvarDensity(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao)
{
  std::vector<AtomNumber> atom;
  parseAtomList("ORIGIN",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  origin is at position of atom : %d\n",atom[0].serial() );

  std::string mlab; parse("DATA",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlab);
  if(!mycolv) error("action labelled " +  mlab + " does not exist or is not a MultiColvar");
  stash = mycolv->buildDataStashes( NULL );

  parseFlag("FRACTIONAL",fractional);
  std::string direction; parse("DIR",direction);
  log.printf("  calculating for %s density profile along ", mycolv->getLabel().c_str() );
  if( direction=="x" ) {
    log.printf("x axis");
    directions.resize(1); directions[0]=0;
  } else if( direction=="y" ) {
    log.printf("y axis");
    directions.resize(1); directions[0]=1;
  } else if( direction=="z" ) {
    log.printf("z axis");
    directions.resize(1); directions[0]=2;
  } else if( direction=="xy" ) {
    log.printf("x and y axes");
    directions.resize(2); directions[0]=0; directions[1]=1;
  } else if( direction=="xz" ) {
    log.printf("x and z axes");
    directions.resize(2); directions[0]=0; directions[1]=2;
  } else if( direction=="yz" ) {
    log.printf("y and z axis");
    directions.resize(2); directions[0]=1; directions[1]=2;
  } else if( direction=="xyz" ) {
    log.printf("x, y and z axes");
    directions.resize(3); directions[0]=0; directions[1]=1; directions[2]=2;
  } else {
    error( direction + " is not valid gradient direction");
  }
  log.printf(" for colvars calculated by action %s \n",mycolv->getLabel().c_str() );
  parseVector("NBINS",nbins); parseVector("SPACING",gspacing);
  if( nbins.size()!=directions.size() && gspacing.size()!=directions.size() ) error("NBINS or SPACING must be set");

  confined.resize( directions.size() ); cmin.resize( directions.size(), 0 ); cmax.resize( directions.size(), 0 );
  for(unsigned i=0; i<directions.size(); ++i) {
    if( directions[i]==0 ) {
      bool tflag; parseFlag("XREDUCED",tflag); confined[i]=tflag;
      if( confined[i] ) {
        cmin[i]=cmax[i]=0.0; parse("XLOWER",cmin[i]); parse("XUPPER",cmax[i]);
        if( fractional ) error("XREDUCED is incompatible with FRACTIONAL");
        if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for x axis makes no sense");
        log.printf("  confining calculation in x direction to between %f and %f \n",cmin[i],cmax[i]);
      }
    } else if( directions[i]==1 ) {
      bool tflag; parseFlag("YREDUCED",tflag); confined[i]=tflag;
      if( confined[i] ) {
        cmin[i]=cmax[i]=0.0; parse("YLOWER",cmin[i]); parse("YUPPER",cmax[i]);
        if( fractional ) error("YREDUCED is incompatible with FRACTIONAL");
        if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for y axis makes no sense");
        log.printf("  confining calculation in y direction to between %f and %f \n",cmin[i],cmax[i]);
      }
    } else if( directions[i]==2 ) {
      bool tflag; parseFlag("ZREDUCED",tflag); confined[i]=tflag;
      if( confined[i] ) {
        cmin[i]=cmax[i]=0.0; parse("ZLOWER",cmin[i]); parse("ZUPPER",cmax[i]);
        if( fractional ) error("ZREDUCED is incompatible with FRACTIONAL");
        if( fabs(cmin[i]-cmax[i])<epsilon ) error("range set for z axis search makes no sense");
        log.printf("  confining calculation in z direction to between %f and %f \n",cmin[i],cmax[i]);
      }
    }
  }

  std::string vstring;
  if( confined[0] ) vstring +="PBC=F";
  else vstring += " PBC=T";
  for(unsigned i=1; i<directions.size(); ++i) {
    if( confined[i] ) vstring += ",F";
    else vstring += ",T";
  }
  vstring +=" COMPONENTS=" + mycolv->getLabel() + ".dens";
  vstring +=" COORDINATES=";
  if( directions[0]==0 ) vstring+="x";
  else if( directions[0]==1 ) vstring+="y";
  else if( directions[0]==2 ) vstring+="z";
  for(unsigned i=1; i<directions.size(); ++i) {
    if( directions[i]==0 ) vstring+=",x";
    else if( directions[i]==1 ) vstring+=",y";
    else if( directions[i]==2 ) vstring+=",z";
  }
  // Create a task list
  for(unsigned i=0; i<mycolv->getFullNumberOfTasks(); ++i) addTaskToList(i);
  // And create the grid
  std::unique_ptr<gridtools::GridVessel> grid;
  if( mycolv->isDensity() ) grid=createGrid( "histogram", vstring );
  else grid=createGrid( "average", vstring );
  mygrid=grid.get();
  // And finish the grid setup
  setAveragingAction( std::move(grid), true );

  // Enusre units for cube files are set correctly
  if( !fractional ) {
    if( plumed.getAtoms().usingNaturalUnits() ) mygrid->setCubeUnits( 1.0/0.5292 );
    else mygrid->setCubeUnits( plumed.getAtoms().getUnits().getLength()/.05929 );
  }

  checkRead(); requestAtoms(atom);
  // Stupid dependencies cleared by requestAtoms - why GBussi why? That's got me so many times
  addDependency( mycolv );
}

unsigned MultiColvarDensity::getNumberOfQuantities() const {
  return directions.size() + 2;
}

void MultiColvarDensity::clearAverage() {
  std::vector<double> min(directions.size()), max(directions.size());
  std::vector<std::string> gmin(directions.size()), gmax(directions.size());;
  for(unsigned i=0; i<directions.size(); ++i) { min[i]=-0.5; max[i]=0.5; }
  if( !fractional ) {
    if( !mycolv->getPbc().isOrthorombic() ) {
      error("I think that density profiles with non-orthorhombic cells don't work.  If you want it have a look and see if you can work it out");
    }

    for(unsigned i=0; i<directions.size(); ++i) {
      if( !confined[i] ) {
        min[i]*=mycolv->getBox()(directions[i],directions[i]);
        max[i]*=mycolv->getBox()(directions[i],directions[i]);
      } else {
        min[i]=cmin[i]; max[i]=cmax[i];
      }
    }
  }
  for(unsigned i=0; i<directions.size(); ++i) { Tools::convert(min[i],gmin[i]); Tools::convert(max[i],gmax[i]); }
  ActionWithAveraging::clearAverage();
  mygrid->setBounds( gmin, gmax, nbins, gspacing ); resizeFunctions();
}

void MultiColvarDensity::prepareForAveraging() {
  for(unsigned i=0; i<directions.size(); ++i) {
    if( confined[i] ) continue;
    std::string max; Tools::convert( 0.5*mycolv->getBox()(directions[i],directions[i]), max );
    if( max!=mygrid->getMax()[i] ) error("box size should be fixed.  Use FRACTIONAL");
  }
  // Ensure we only work with active multicolvars
  deactivateAllTasks();
  for(unsigned i=0; i<stash->getNumberOfStoredValues(); ++i) taskFlags[i]=1;
  lockContributors();
  // Retrieve the origin
  origin = getPosition(0);
}

void MultiColvarDensity::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> cvals( mycolv->getNumberOfQuantities() ); stash->retrieveSequentialValue( current, false, cvals );
  Vector fpos, apos = pbcDistance( origin, mycolv->getCentralAtomPos( mycolv->getPositionInFullTaskList(current) ) );
  if( fractional ) { fpos = getPbc().realToScaled( apos ); } else { fpos=apos; }

  myvals.setValue( 0, cweight*cvals[0] ); for(unsigned j=0; j<directions.size(); ++j) myvals.setValue( 1+j, fpos[ directions[j] ] );
  myvals.setValue( 1+directions.size(), cvals[1] );
}

}
}
