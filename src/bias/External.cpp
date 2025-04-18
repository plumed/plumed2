/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "Bias.h"
#include "core/ActionRegister.h"
#include "tools/Grid.h"
#include "tools/Exception.h"
#include "tools/File.h"

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS EXTERNAL
/*
Calculate a restraint that is defined on a grid that is read during start up

In this action you define an external potential that you would like to use as a bias in a file.
This file contains the value of your bias potential on a grid of points.  The instantaneous value
of the bias potential is then computed via interpolation.


The following input illustrates how this works in practice.  Here the external potential is
defined in the file `extras/1d_bias.dat`. This potential acts on the distance between atoms 3 and 5.

```plumed
#SETTINGS INPUTFILES=extras/1d_bias.grid
d1: DISTANCE ATOMS=3,5
m: EXTERNAL ARG=d1 FILE=extras/1d_bias.grid
```

As you can see from the example above, the file that is input to the EXTERNAL command
contains the value of the function and its derivative on a grid of points.  This file
has the grid format that is discussed in the documentation for [DUMPGRID](DUMPGRID.md).

You can also include grids that are a function of more than one collective
variable.  For instance the following would be the input for an external
potential acting on two torsional angles:

```plumed
#SETTINGS INPUTFILES=extras/2d_bias.grid
t1: TORSION ATOMS=4,5,6,7
t2: TORSION ATOMS=6,7,8,9
ext: EXTERNAL ARG=t1,t2 FILE=extras/2d_bias.grid
```

Please note the order that the order of arguments in the plumed.dat file must be the same as
the order of arguments in the header of the grid file.
*/
//+ENDPLUMEDOC

class External : public Bias {

private:
  std::unique_ptr<GridBase> BiasGrid_;
  double scale_;

public:
  explicit External(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(External,"EXTERNAL")

void External::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.add("compulsory","FILE","the name of the file containing the external potential.");
  keys.addFlag("NOSPLINE",false,"specifies that no spline interpolation is to be used when calculating the energy and forces due to the external potential");
  keys.addFlag("SPARSE",false,"specifies that the external potential uses a sparse grid");
  keys.add("compulsory","SCALE","1.0","a factor that multiplies the external potential, useful to invert free energies");
}

External::External(const ActionOptions& ao):
  PLUMED_BIAS_INIT(ao) {
  std::string filename;
  parse("FILE",filename);
  if( filename.length()==0 ) {
    error("No external potential file was specified");
  }
  bool sparsegrid=false;
  parseFlag("SPARSE",sparsegrid);
  bool nospline=false;
  parseFlag("NOSPLINE",nospline);
  bool spline=!nospline;
  parse("SCALE",scale_);

  checkRead();

  log.printf("  External potential from file %s\n",filename.c_str());
  log.printf("  Multiplied by %lf\n",scale_);
  if(spline) {
    log.printf("  External potential uses spline interpolation\n");
  }
  if(sparsegrid) {
    log.printf("  External potential uses sparse grid\n");
  }

// read grid
  IFile gridfile;
  gridfile.open(filename);
  std::string funcl=getLabel() + ".bias";
  BiasGrid_=GridBase::create(funcl,getArguments(),gridfile,sparsegrid,spline,true);
  if(BiasGrid_->getDimension()!=getNumberOfArguments()) {
    error("mismatch between dimensionality of input grid and number of arguments");
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isPeriodic()!=BiasGrid_->getIsPeriodic()[i] ) {
      error("periodicity mismatch between arguments and input bias");
    }
  }
}

void External::calculate() {
  unsigned ncv=getNumberOfArguments();
  std::vector<double> cv(ncv), der(ncv);

  for(unsigned i=0; i<ncv; ++i) {
    cv[i]=getArgument(i);
  }

  double ene=scale_*BiasGrid_->getValueAndDerivatives(cv,der);

  setBias(ene);

  for(unsigned i=0; i<ncv; ++i) {
    const double f=-scale_*der[i];
    setOutputForce(i,f);
  }
}

}
}
