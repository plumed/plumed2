/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "RDF.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS RDF
/*
Calculate the radial distribution function

This shortcut can be used to evaluate the [radial distribution function](https://en.wikipedia.org/wiki/Radial_distribution_function) for
a collection of atoms.  The following example illustrates how it can be used to compute the average radial distribution from a
trajectory:

```plumed
rdf: RDF GROUP=1-108 MAXR=2.5 GRID_BIN=25 KERNEL=DISCRETE STRIDE=1 CLEAR=0
DUMPGRID ARG=rdf FILE=rdf.dat
```

If you expand the shortcut in this input you can see how the radial distribution function is computed using [KDE](KDE.md) actions and that the
average is obtained by using [ACCUMULATE](ACCUMULATE.md) actions.  By using expanded versions of this shortcut you can thus calculate and print
the instantaneous value of the radial distribution function.  Alternatively, if you use the CLEAR option in the [ACCUMULATE](ACCUMULATE.md) commands
you can calculate radial distribution functions from various parts of the trajectory.

The following example shows how you can use this function to calculate the radial distribution for the atoms in GROUPB around the atoms in GROUPA.
Notice that the DISCRETE keyword has been removed here so we are using [kernel density estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation)
to compute the radial distribution function.

```plumed
rdf: RDF GROUPA=1-108 GROUPB=109-300 MAXR=2.5 GRID_BIN=25 BANDWIDTH=0.01 DENSITY=1 STRIDE=1 CLEAR=0
DUMPGRID ARG=rdf FILE=rdf.dat
```

We have also used the `DENSITY` keyword to set the background density that is used when normalizing the radial distribution function explicity to 1 atom$/nm^{3}$.
When this keyword is not used, this density is calculated by dividing the number of atoms by the volume of the box as you can see if you expand the shortcut in the
first input above.

Notice that you can also use the `NO_AVERAGE` keyword as shown below to calculate an independent RDF around each of the atoms:

```plumed
rdf: RDF GROUPA=1-3 GROUPB=4-108 MAXR=2.5 GRID_BIN=25 BANDWIDTH=0.01 DENSITY=1 STRIDE=1 NO_AVERAGE CLEAR=0
DUMPGRID ARG=rdf FILE=rdf.dat
```

This command will output a function on a $25 \times 3$ grid.  The first coordinate is the distance between atoms.  The second coordinate is then then an index that tells
you, which central atom.  In other words, the input above will output three radial distribution functions.  The first of these is the radial distribution function around atom 1,
the second is the radial distribution function around atom 2 and the last is the radial distribution function around atom 3.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

PLUMED_REGISTER_ACTION(RDF,"RDF")

void RDF::getDistanceMatrixShape( const std::string& lab, ActionShortcut* action, std::vector<std::string>& shape_str ) {
  std::vector<Value*> args;
  std::vector<std::string> arglabs(1);
  arglabs[0] = lab;
  ActionWithArguments::interpretArgumentList( arglabs, action->plumed.getActionSet(), action, args );
  Tools::convert( args[0]->getShape()[0], shape_str[0] );
  Tools::convert( args[0]->getShape()[1], shape_str[1] );
}

void RDF::createX2ReferenceObject( const std::string& lab, const std::string& grid_setup, const bool& calc_dens, const bool& no_average, ActionShortcut* action ) {
  // Create grid with normalizing function
  if( no_average ) {
    action->readInputLine( lab  + "_x2: REFERENCE_GRID PERIODIC=NO,NO FUNC=x*x+0*y " + grid_setup );
  } else {
    action->readInputLine( lab  + "_x2: REFERENCE_GRID PERIODIC=NO FUNC=x*x " + grid_setup );
  }
  // Compute density if required
  if( calc_dens ) {
    action->readInputLine( lab + "_vol: VOLUME" );
  }
}

void RDF::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","GROUP","the atoms that are being used to calculate the RDF");
  keys.addDeprecatedKeyword("ATOMS","GROUP");
  keys.add("atoms-2","GROUPA","the atoms that you would like to compute the RDF about.  Must be used with GROUPB.");
  keys.add("atoms-2","GROUPB","the atoms that you would like to to use when computing the RDF around the atoms that were specified with GROUPA");
  keys.add("compulsory","GRID_BIN","the number of bins to use when computing the RDF");
  keys.add("compulsory","KERNEL","GAUSSIAN","the type of kernel to use for computing the histograms for the RDF");
  keys.add("compulsory","CUTOFF","6.25","the cutoff at which to stop evaluating the kernel functions is set equal to sqrt(2*x)*bandwidth in each direction where x is this number");
  keys.add("compulsory","MAXR","the maximum distance to use for the rdf");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  if( keys.getDisplayName()=="RDF") {
    keys.add("compulsory","CLEAR","1","the frequency with which to clear the estimate of the rdf.  Set equal to 0 if you want to compute an rdf over the whole trajectory");
    keys.add("compulsory","STRIDE","1","the frequency with which to compute the rdf and accumulate averages");
  }
  keys.add("optional","DENSITY","the reference density to use when normalizing the RDF");
  keys.addFlag("NO_AVERAGE",false,"output a two-dimensional grid with separate RDF functions for each of the central atoms");
  keys.add("hidden","REFERENCE","this is the label of the reference objects");
  keys.setValueDescription("grid","the radial distribution function");
  keys.needsAction("REFERENCE_GRID");
  keys.needsAction("VOLUME");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("CUSTOM");
  keys.needsAction("KDE");
  keys.needsAction("ACCUMULATE");
  keys.needsAction("CONSTANT");
  keys.needsAction("ONES");
  keys.needsAction("OUTER_PRODUCT");
}

RDF::RDF(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  bool noaverage;
  parseFlag("NO_AVERAGE",noaverage);
  // Read in grid extent and number of bins
  std::string maxr, nbins, dens;
  parse("MAXR",maxr);
  parse("GRID_BIN",nbins);
  parse("DENSITY",dens);
  // Read input to histogram
  std::string cutoff, kernel, bandwidth, kernel_data;
  parse("KERNEL",kernel);
  if( kernel=="DISCRETE" ) {
    cutoff = maxr;
    kernel_data="KERNEL=DISCRETE";
    warning("rdf is normalised by dividing by the surface area at the grid value and not by the volume of the bin as it should be with discrete kernels");
  } else {
    parse("BANDWIDTH",bandwidth);
    double bw;
    Tools::convert( bandwidth, bw );
    if( noaverage ) {
      bandwidth += ",0";
    }
    double rcut;
    parse("CUTOFF",rcut);
    kernel_data="KERNEL=" + kernel + " BANDWIDTH=" + bandwidth;
    double fcut;
    Tools::convert( maxr, fcut );
    Tools::convert( fcut + sqrt(2.0*rcut)*bw, cutoff );
  }

  // Create contact matrix
  std::string atom_str, group_str, groupa_str, groupb_str;
  parse("GROUP",group_str);
  if( group_str.size()==0 ) {
    parse("ATOMS",group_str);
  }
  if( group_str.length()>0 ) {
    atom_str="GROUP=" + group_str;
  } else {
    parse("GROUPA",groupa_str);
    parse("GROUPB",groupb_str);
    atom_str="GROUPA=" + groupa_str + " GROUPB=" + groupb_str;
  }

  // Calculate all the distances
  readInputLine( getShortcutLabel() + "_mat: DISTANCE_MATRIX CUTOFF=" + cutoff + " " + atom_str);
  // Get the shape of the matrix
  std::vector<std::string> shape_str(2);
  getDistanceMatrixShape( getShortcutLabel() + "_mat", this, shape_str );

  // Setup the grid
  std::string grid_setup = "GRID_MIN=0 GRID_MAX=" + maxr + " GRID_BIN=" + nbins;
  if( noaverage ) {
    int iatoms;
    std::string num, str_values;
    Tools::convert( 0, str_values );
    Tools::convert( shape_str[0], iatoms );
    for(int i=1; i<iatoms; ++i) {
      Tools::convert( i, num );
      str_values += "," + num;
    }
    Tools::convert( iatoms-1, num );
    readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + shape_str[1] );
    readInputLine( getShortcutLabel() + "_index: CONSTANT VALUES=" + str_values );
    readInputLine( getShortcutLabel() + "_atno: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_index," + getShortcutLabel() + "_ones" );
    grid_setup = "GRID_MIN=0,0 GRID_MAX=" + maxr + "," + num + " GRID_BIN=" + nbins + "," + num;
    shape_str[0] = "1";
  }
  // Create grid with normalizing function on it
  std::string refstr;
  parse("REFERENCE",refstr);
  if( refstr.length()==0 ) {
    if( noaverage ) {
      readInputLine( getShortcutLabel() + "_x2: REFERENCE_GRID PERIODIC=NO,NO FUNC=x*x+0*y " + grid_setup );
    } else {
      readInputLine( getShortcutLabel() + "_x2: REFERENCE_GRID PERIODIC=NO FUNC=x*x " + grid_setup );
    }
    refstr = getShortcutLabel();
  }
  if( dens.length()==0 ) {
    readInputLine( getShortcutLabel() + "_vol: VOLUME" );
  }
  // Calculate weights of distances
  readInputLine( getShortcutLabel() + "_wmat: CUSTOM ARG=" + getShortcutLabel() + "_mat FUNC=step(" + cutoff + "-x) PERIODIC=NO");
  // Now create a histogram from the contact matrix
  unsigned clear, stride;
  parse("CLEAR",clear);
  parse("STRIDE",stride);
  if( clear==1 ) {
    if( noaverage ) {
      readInputLine( getShortcutLabel() + "_kde: KDE ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_atno VOLUMES=" + getShortcutLabel() + "_wmat " + grid_setup + " " + kernel_data);
    } else {
      readInputLine( getShortcutLabel() + "_kde: KDE ARG=" + getShortcutLabel() + "_mat VOLUMES=" + getShortcutLabel() + "_wmat " + grid_setup + " " + kernel_data);
    }
  } else {
    std::string stridestr, clearstr;
    Tools::convert( stride, stridestr );
    Tools::convert( clear, clearstr );
    if( noaverage ) {
      readInputLine( getShortcutLabel() + "_okde: KDE ARG=" + getShortcutLabel() + "_mat," + getShortcutLabel() + "_atno HEIGHTS=" + getShortcutLabel() + "_wmat " + grid_setup + " " + kernel_data);
    } else {
      readInputLine( getShortcutLabel() + "_okde: KDE ARG=" + getShortcutLabel() + "_mat HEIGHTS=" + getShortcutLabel() + "_wmat " + grid_setup + " " + kernel_data);
    }
    readInputLine( getShortcutLabel() + "_kde: ACCUMULATE ARG=" + getShortcutLabel() + "_okde STRIDE=" + stridestr + " CLEAR=" + clearstr );
    readInputLine( getShortcutLabel() + "_one: CONSTANT VALUE=1");
    readInputLine( getShortcutLabel() + "_norm: ACCUMULATE ARG=" + getShortcutLabel() + "_one STRIDE=" + stridestr + " CLEAR=" + clearstr );
  }
  // Transform the histogram by normalizing factor for rdf
  readInputLine( getShortcutLabel() + "_vrdf: CUSTOM ARG=" + getShortcutLabel() + "_kde," + refstr + "_x2 FUNC=x/(4*pi*y) PERIODIC=NO");
  // And normalize by density and number of atoms (separated from above to avoid nans)
  std::string func_str = "PERIODIC=NO ARG=" + getShortcutLabel() + "_vrdf";
  if( dens.length()>0 ) {
    if( clear==1 ) {
      func_str += " FUNC=x/(" + dens + "*" + shape_str[0] + ")";
    } else {
      func_str += "," + getShortcutLabel() + "_norm FUNC=x/(y*" + dens + "*" + shape_str[0] + ")";
    }
  } else {
    if( clear==1 ) {
      func_str += "," + refstr + "_vol FUNC=x*y/(" + shape_str[1] + "*" + shape_str[0] + ")";
    } else {
      func_str += "," + refstr + "_vol," + getShortcutLabel() + "_norm FUNC=x*y/(z*" + shape_str[1] + "*" + shape_str[0] + ")";
    }
  }
  readInputLine( getShortcutLabel() + ": CUSTOM " + func_str);
}

}
}
