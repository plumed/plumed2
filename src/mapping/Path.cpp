/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "Path.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/PDB.h"

//+PLUMEDOC COLVAR PATH
/*
Path collective variables with a more flexible framework for the distance metric being used.

The Path Collective Variables developed by Branduardi and that are described in the first paper that is cited below alow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional
path.  The progress along the path (s) is computed using:

$$
s = \frac{ \sum_{i=1}^N i \exp( -\lambda R[X - X_i] ) }{ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) }
$$

while the distance from the path (z) is measured using:

$$
z = -\frac{1}{\lambda} \ln\left[ \sum_{i=1}^N \exp( -\lambda R[X - X_i] ) \right]
$$

In these expressions $N$ high-dimensional frames ($X_i$) are used to describe the path in the high-dimensional
space. The two expressions above are then functions of the distances from each of the high-dimensional frames $R[X - X_i]$.
Within PLUMED there are multiple ways to define the distance from a high-dimensional configuration.  You could calculate
the RMSD distance or you could calculate the amount by which a set of collective variables change.  As such this shortcut
of the path CV allows one to compute the distances from the paths in a variety of different ways. This is as opposed to
the alternative implementation of path ([PATHMSD](PATHMSD.md)) which is a bit faster but which only allows one to use the [RMSD](RMSD.md) distance.

## Examples

In the example below the path is defined using RMSD distance from frames.

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1: PATH REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL LAMBDA=69087
PRINT ARG=p1.spath,p1.zpath STRIDE=1 FILE=colvar
```

In the example below the path is defined using the values of two torsional angles (t1 and t2).  You can specify these arguments in the PLUMED input
as illustrated below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: PATH ARG=t1,t2 REFERENCE=regtest/mapping/rt-tpath/epath.pdb LAMBDA=1.0
PRINT ARG=pp.* FILE=colvar
```

However, this is not strictly necessary as PLUMED can also get the names of the arguments from the input pdb file directly as illustrated in the following
input.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: PATH REFERENCE=regtest/mapping/rt-tpath/epath.pdb LAMBDA=1.0
PRINT ARG=pp.* FILE=colvar
```

If you look at the pdb file for the input above you can see that the remarks tell PLUMED the labels that are being used to define the position in the
high dimensional space and the values that these arguments have at each point on the path.

## Controlling what is calculated

If you only want to calculate the $s$ coordinate you can use the NOZPATH flag as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1: PATH REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL NOZPATH LAMBDA=69087
PRINT ARG=p1.spath STRIDE=1 FILE=colvar
```

Similarly, if you only want to calculate the $z$ coordinate you use the NOSPATH flag as shown below

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1: PATH REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL NOSPATH LAMBDA=69087
PRINT ARG=p1.zpath STRIDE=1 FILE=colvar
```

You can also tell PLUMED to calculate the distance along and distance from the path using the formulas that are given in the documentation for the
[GPATH](GPATH.md) shortcut as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
p1: PATH REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL GPATH LAMBDA=69087
PRINT ARG=p1.spath,p1.zpath,p1.gspath,p1.gzpath STRIDE=1 FILE=colvar
```

The PATH command in this input calculates four quantities:

- p1.spath - the s coordinate evaluated using the expression above
- p1.zpath - the z coordinate evluated using the expression above
- p1.gspath - the s coordinate evaluated using the expression in the documentation for [GPATH](GPATH.md)
- p1.gzpath - the z coordinate evaluated using the expression in the documentation for [GPATH](GPATH.md)

## Paths using normalized Euclidean distances

If you are defining the positions of the frames in the path using arguments you may wish to calculate the distances between the instantaneous values of the
arguments and the values of these arguments on the path using the [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) rather than the
[EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md). In other words, instead of calculating the distance between the vectors coordinates of coordinates $x$ and $y$ using:

$$
D^2 = \sum_{i=1}^N (x_i - y_i)^2
$$

where the sum runs over the number of arguments, you can calculate this distance as:

$$
D^2 = \sum_{i=1}^N C_i (x_i - y_i)^2
$$

where the $C_i$ are a set of coefficients that are specified using the COEFFICIENTS keyword as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: PATH ...
  ARG=t1,t2 COEFFICIENTS=0.2,0.3
  REFERENCE=regtest/mapping/rt-tpath/epath.pdb
  LAMBDA=1.0
...
PRINT ARG=pp.* FILE=colvar
```

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR GPROPERTYMAP
/*
Property maps but with a more flexible framework for the distance metric being used.

This colvar calculates a property map using the formalism developed by Spiwok that is referenced in the second paper cited below.
In essence if you have the value of some property, $X_i$, that it takes at a set of high-dimensional
positions then you calculate the value of the property at some arbitrary point in the high-dimensional space
using:

$$
X=\frac{\sum_i X_i*\exp(-\lambda D_i(x))}{\sum_i  \exp(-\lambda D_i(x))}
$$

In these expressions the value of the property $X$ is given at the position of $N$ high-dimensional frames ($X_i$). The distances, $D_i$, between
these frames and our instaneous position are then used to compute the value of the property at the instaneous position in the high dimensional space.
As illustrated in the examples below, this implementation of the property map allows you to calculate these distances in various different ways.
This is as opposed to the alternative implementation [PROPERTYMAP](PROPERTYMAP.md) which is a bit faster but which only allows one to use the RMSD distance.

## Examples

The input shown below can be used to calculate the interpolated values of two properties called X and Y based on the values
that these properties take at a set of reference configurations and using the formula above.  For this input the distances
between the reference configurations and the instantaneous configurations are calculated using the OPTIMAL metric that is
discussed at length in the manual pages on [RMSD](RMSD.md).

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/allv.pdb
p2: GPROPERTYMAP REFERENCE=regtest/trajectories/path_msd/allv.pdb PROPERTY=X,Y LAMBDA=69087
PRINT ARG=p2_X,p2_Y,p2_zpath STRIDE=1 FILE=colvar
```

Notice that the REFERENCE file here gives the values of the properties at each of the points of interest.
In this second input the value of the property at each point of interest is also given in the REFERENCE file.
Here as well the REFERENCE file also tells us that each of the reference points is defined based on the values of the
two torsions `t1` and `t2`.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: GPROPERTYMAP ...
   ARG=t1,t2 REFERENCE=regtest/mapping/rt-tpath/epath.pdb
   PROPERTY=X LAMBDA=1.0
...
PRINT ARG=pp.* FILE=colvar
```

Notice that you can use the ARG keyword to specify the arguments are used in the definition of the PATH.  However, you do not have to
use this keyword.  The following input also works as PLUMED gets the input arguments from the file that contains the definition of the
path.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: GPROPERTYMAP ...
   REFERENCE=regtest/mapping/rt-tpath/epath.pdb
   PROPERTY=X LAMBDA=1.0
...
PRINT ARG=pp.* FILE=colvar

## Controlling what is calculated

If you dont want to bother calculating the z coordinate you can use the NOZPATH flag as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/trajectories/path_msd/allv.pdb
p2: GPROPERTYMAP ...
   REFERENCE=regtest/trajectories/path_msd/allv.pdb
   PROPERTY=X,Y LAMBDA=69087 NOZPATH
...
PRINT ARG=p2_X,p2_Y STRIDE=1 FILE=colvar
```

## Property maps using normalized Euclidean distances

If you are defining the positions of the frames in the path using arguments you may wish to calculate the distances between the instantaneous values of the
arguments and the values of these arguments on the path using the [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) rather than the
[EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md). In other words, instead of calculating the distance between the vectors coordinates of coordinates $x$ and $y$ using:

$$
D^2 = \sum_{i=1}^N (x_i - y_i)^2
$$

where the sum runs over the number of arguments, you can calculate this distance as:

$$
D^2 = \sum_{i=1}^N C_i (x_i - y_i)^2
$$

where the $C_i$ are a set of coefficients that are specified using the COEFFICIENTS keyword as shown below.

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-tpath/epath.pdb
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17
pp: GPROPERTYMAP ...
   REFERENCE=regtest/mapping/rt-tpath/epath.pdb
   PROPERTY=X COEFFICIENTS=0.2,0.3 LAMBDA=1.0
...
PRINT ARG=pp.* FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

PLUMED_REGISTER_ACTION(Path,"PATH")
PLUMED_REGISTER_ACTION(Path,"GPROPERTYMAP")

void Path::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  Path::registerInputFileKeywords( keys );
  if( keys.getDisplayName()=="GPROPERTYMAP" ) {
    keys.add("optional","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
  }
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.addFlag("NOZPATH",false,"do not calculate the zpath CV");
  keys.add("optional","COEFFICIENTS","the coefficients of the displacements along each argument that should be used when calculating the euclidean distance");
  if( keys.getDisplayName()=="PATH" ) {
    keys.addFlag("NOSPATH",false,"do not calculate the spath CV");
    keys.addFlag("GPATH",false,"calculate the geometric path");
    keys.addOutputComponent("gspath","GPATH","scalar","the position along the path calculated using the geometric formula");
    keys.addOutputComponent("gzpath","GPATH","scalar","the distance from the path calculated using the geometric formula");
  }
  keys.addOutputComponent("spath","default","scalar","the position along the path calculated");
  keys.addOutputComponent("zpath","default","scalar","the distance from the path calculated");
  keys.addDOI("10.1063/1.2432340");
  keys.addDOI("10.1063/1.3660208");
}

void Path::registerInputFileKeywords( Keywords& keys ) {
  keys.add("compulsory","REFERENCE","a pdb file containing the set of reference configurations");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
           "metrics that are available in PLUMED can be found in the section of the manual on "
           "\\ref dists");
  keys.addInputKeyword("optional","ARG","scalar","the list of arguments you would like to use in your definition of the path");
  keys.needsAction("DRMSD");
  keys.needsAction("RMSD");
  keys.needsAction("LOWEST");
  keys.needsAction("GPATH");
  keys.needsAction("EUCLIDEAN_DISTANCE");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
  keys.needsAction("COMBINE");
  keys.needsAction("NORMALIZED_EUCLIDEAN_DISTANCE");
  keys.needsAction("PDB2CONSTANT");
  keys.needsAction("CONSTANT");
}

Path::Path( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {
  // Now create all other parts of the calculation
  std::string lambda;
  parse("LAMBDA",lambda);
  std::string reference;
  parse("REFERENCE",reference);
  bool nospath=false, nozpath=false, gpath=false;
  if( getName()=="PATH" ) {
    parseFlag("NOSPATH",nospath);
    parseFlag("GPATH",gpath);
  }
  parseFlag("NOZPATH",nozpath);
  if( gpath ) {
    readInputLine( getShortcutLabel() + "_gpath: GPATH REFERENCE=" + reference );
    readInputLine( getShortcutLabel() + "_gspath: COMBINE ARG=" + getShortcutLabel() + "_gpath.s PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_gzpath: COMBINE ARG=" + getShortcutLabel() + "_gpath.z PERIODIC=NO");
  }
  if( nospath && nozpath ) {
    return;
  }
  // Setup the properties
  std::vector<std::string> properties, pnames;
  if( getName()=="PATH") {
    properties.resize(1);
  } else {
    parseVector("PROPERTY",pnames);
    properties.resize( pnames.size() );
  }
  std::string type, reference_data;
  std::vector<std::string> argnames;
  parseVector("ARG",argnames);
  parse("TYPE",type);
  if( type.find("DRMSD")!=std::string::npos ) {
    readInputLine( getShortcutLabel() + "_data: DRMSD SQUARED TYPE=" + type + " REFERENCE=" + reference );
  } else {
    readInputFrames( reference, type, argnames, false, this, reference_data );
  }
  // Find the shortest distance to the frames
  readInputLine( getShortcutLabel() + "_mindist: LOWEST ARG=" + getShortcutLabel() + "_data");
  // Now create MATHEVAL object to compute exponential functions
  readInputLine( getShortcutLabel() + "_weights: CUSTOM ARG=" + getShortcutLabel() + "_data," + getShortcutLabel() + "_mindist FUNC=exp(-(x-y)*" + lambda + ") PERIODIC=NO" );
  // Create denominator
  readInputLine( getShortcutLabel() + "_denom: SUM ARG=" + getShortcutLabel() + "_weights PERIODIC=NO");
  // Now compte zpath variable
  if( !nozpath ) {
    readInputLine( getShortcutLabel() + "_z: CUSTOM ARG=" + getShortcutLabel() + "_denom," + getShortcutLabel() + "_mindist FUNC=y-log(x)/" + lambda + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zpath: COMBINE ARG=" + getShortcutLabel() + "_z PERIODIC=NO");
  }
  // Now get coefficients for properies for spath
  readPropertyInformation( pnames, getShortcutLabel(), reference, this );
  // Now create COMBINE objects to compute numerator of path
  for(unsigned i=0; i<properties.size(); ++i) {
    if( pnames.size()>0 ) {
      readInputLine( pnames[i] + "_numer_prod: CUSTOM ARG=" + getShortcutLabel() + "_weights," + pnames[i] + "_ref FUNC=x*y PERIODIC=NO");
      readInputLine( pnames[i] + "_numer: SUM ARG=" + pnames[i]  + "_numer_prod PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_" + pnames[i] + ": CUSTOM ARG=" + pnames[i]  + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    } else if( !nospath ) {
      readInputLine( getShortcutLabel() + "_s_prod: CUSTOM ARG=" + getShortcutLabel() + "_weights," + getShortcutLabel() + "_ind FUNC=x*y PERIODIC=NO");
      readInputLine( getShortcutLabel()  + "_numer: SUM ARG=" + getShortcutLabel() + "_s_prod PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_s: CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_spath: COMBINE ARG=" + getShortcutLabel() + "_s PERIODIC=NO");
    }
  }
}

std::string Path::fixArgumentName( const std::string& argin ) {
  std::string argout = argin;
  std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) {
    argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  }
  return argout;
}

void Path::readInputFrames( const std::string& reference, const std::string& type, std::vector<std::string>& argnames, const bool& displacements, ActionShortcut* action, std::string& reference_data ) {
  FILE* fp=std::fopen(reference.c_str(),"r");
  PDB pdb;
  if(!fp) {
    action->error("could not open reference file " + reference );
  }
  bool do_read=pdb.readFromFilepointer(fp,false,0.1);
  if( !do_read ) {
    action->error("missing file " + reference );
  }
  if( pdb.getPositions().size()!=0 && argnames.size()==0 ) {
    reference_data = action->getShortcutLabel() + "_data_ref";
    if( displacements ) {
      action->readInputLine( action->getShortcutLabel() + "_data: RMSD DISPLACEMENT SQUARED REFERENCE=" + reference + " TYPE=" + type );
    } else {
      action->readInputLine( action->getShortcutLabel() + "_data: RMSD SQUARED REFERENCE=" + reference + " TYPE=" + type );
    }
  } else if( pdb.getPositions().size()!=0 ) {
    reference_data = action->getShortcutLabel() + "_atomdata_ref";
    if( displacements ) {
      action->readInputLine( action->getShortcutLabel() + "_atomdata: RMSD DISPLACEMENT SQUARED REFERENCE=" + reference + " TYPE=" + type );
    } else {
      action->readInputLine( action->getShortcutLabel() + "_atomdata: RMSD SQUARED REFERENCE=" + reference + " TYPE=" + type );
    }
  } else if( argnames.size()==0 ) {
    argnames.resize( pdb.getArgumentNames().size() );
    for(unsigned i=0; i<argnames.size(); ++i) {
      argnames[i] = pdb.getArgumentNames()[i];
    }
  }
  std::vector<Value*> theargs;
  if( argnames.size()>0 ) {
    ActionWithArguments::interpretArgumentList( argnames, action->plumed.getActionSet(), action, theargs );
  }

  if( theargs.size()>0 ) {
    std::string instargs, refargs;
    for(unsigned i=0; i<theargs.size(); ++i) {
      std::string iargn = fixArgumentName( theargs[i]->getName() );
      action->readInputLine( action->getShortcutLabel() + "_ref_" + iargn + ": PDB2CONSTANT REFERENCE=" + reference + " ARG=" + theargs[i]->getName() );
      if( i==0 ) {
        instargs=" ARG2=" + theargs[i]->getName();
        refargs=" ARG1=" + action->getShortcutLabel() + "_ref_" + iargn;
      } else {
        instargs +="," + theargs[i]->getName();
        refargs +="," + action->getShortcutLabel() + "_ref_" + iargn;
      }
      if( pdb.getPositions().size()==0 && i==0 ) {
        reference_data = action->getShortcutLabel() + "_ref_" + iargn;
      } else {
        reference_data += "," + action->getShortcutLabel() + "_ref_" + iargn;
      }
    }
    std::string comname="EUCLIDEAN_DISTANCE SQUARED";
    std::string coeffstr;
    if( action->keywords.exists("COEFFICIENTS") ) {
      action->parse("COEFFICIENTS",coeffstr);
    }
    if( coeffstr.length()>0 ) {
      if( displacements ) {
        action->error("cannot use COEFFICIENTS arguments with GEOMETRIC PATH");
      }
      action->readInputLine( action->getShortcutLabel() + "_coeff: CONSTANT VALUES=" + coeffstr );
      action->readInputLine( action->getShortcutLabel() + "_coeff2: CUSTOM ARG=" + action->getShortcutLabel() + "_coeff FUNC=x*x PERIODIC=NO");
      comname = "NORMALIZED_EUCLIDEAN_DISTANCE SQUARED METRIC=" + action->getShortcutLabel() + "_coeff2";
    } else if( displacements ) {
      comname = "DISPLACEMENT";
    }

    if( pdb.getPositions().size()==0 ) {
      if( displacements ) {
        action->readInputLine( action->getShortcutLabel() + "_dataP: " + comname + instargs + refargs );
        action->readInputLine( action->getShortcutLabel() + "_data: CUSTOM ARG=" + action->getShortcutLabel() + "_dataP FUNC=-x PERIODIC=NO");
      } else {
        action->readInputLine( action->getShortcutLabel() + "_data: " + comname + instargs + refargs );
      }
    } else {
      action->readInputLine( action->getShortcutLabel() + "_argdata: " + comname + instargs + refargs );
    }
  }
}

void Path::readPropertyInformation( const std::vector<std::string>& pnames, const std::string& lab, const std::string& refname, ActionShortcut* action ) {
  if( pnames.size()>0 ) {
    for(unsigned i=0; i<pnames.size(); ++i) {
      action->readInputLine( pnames[i] + ": CONSTANT VALUE=1");
      action->readInputLine( pnames[i] + "_ref: PDB2CONSTANT REFERENCE=" + refname + " ARG=" + pnames[i] );
    }
  } else {
    ActionWithValue* av=action->plumed.getActionSet().selectWithLabel<ActionWithValue*>( lab + "_data" );
    unsigned nfram = av->copyOutput(0)->getShape()[0];
    std::string indices = "VALUES=1";
    for(unsigned i=1; i<nfram; ++i) {
      std::string num;
      Tools::convert( i+1, num );
      indices += "," + num;
    }
    action->readInputLine( lab + "_ind: CONSTANT " + indices );
  }
}

}
}
