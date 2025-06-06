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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "tools/PDB.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR PDB2CONSTANT
/*
Create a constant value from a PDB input file

This shortcut converts the contents of a PDB file to one or more [CONSTANT](CONSTANT.md) actions.
Converting PDB files to Values in this way is useful because it means that when we implement methods, like those
in the [refdist](module_refdist.md) or [mapping](module_mapping.md) modules or the [RMSD](RMSD.md) action, that calculate the distance between two
configurations those two configurations are both stored in PLUMED values. The same code can thus be used to
calculate the difference between the instantaneous configuration
and a constant reference configuration that was read from a file or between two reference configuration.

The following example illustrates how this action can be used to read a set of reference atomic positions

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt19/test0.pdb
ref: PDB2CONSTANT REFERENCE=regtest/basic/rt19/test0.pdb
```

You can see how the reference positions are converted to [CONSTANT](CONSTANT.md) action that outputs a vector by expanding the shortcut.

You can also use this command to read in multiple reference positions as illustrated below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pathtools-3/all.pdb
ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-3/all.pdb
```

The [CONSTANT](CONSTANT.md) that is created by this action is a matrix. Each row of the output matrix contains one set of reference positions.
Notice also that if you have a PDB input which contains multiple reference configurations you can create a vector constant by using the `NUMBER`
keyword to specify the particular configuration that you would like to read in as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pathtools-3/all.pdb
ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-3/all.pdb NUMBER=4
```

The input above will reads in the fourth configuration in the input PDB file only.

## The PDB file format

PLUMED uses the PDB file format here and in several other places

- To read the molecular structure ([MOLINFO](MOLINFO.md)).
- To read reference conformations ([RMSD](RMSD.md), but also in methods such as [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md), etc).

The implemented PDB reader expects a file formatted correctly according to the
[PDB standard](http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html).
In particular, the following columns are read from ATOM records

````
columns | content
1-6     | record name (ATOM or HETATM)
7-11    | serial number of the atom (starting from 1)
13-16   | atom name
18-20   | residue name
22      | chain id
23-26   | residue number
31-38   | x coordinate
39-46   | y coordinate
47-54   | z coordinate
55-60   | occupancy
61-66   | beta factor
````

The PLUMED parser is slightly more permissive than the official PDB format
in the fact that the format of real numbers is not fixed. In other words,
any real number that can be parsed is OK and the dot can be placed anywhere. However,
__columns are interpret strictly__. A sample PDB should look like the following

````
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
````

Notice that serial numbers need not to be consecutive. In the three-line example above,
only the coordinates of three atoms are provided. This is perfectly legal and indicates to PLUMED
that information about these atoms only is available. This could be both for structural
information in [MOLINFO](MOLINFO.md), where the other atoms would have no name assigned, and for
reference structures used in [RMSD](RMSD.md), where only the provided atoms would be used to compute RMSD.

## Including arguments in PDB files

If you wish to specify reference values for PLUMED Values in the REMARKS of a PLUMED input file like this:

````
REMARK t1=-4.3345
REMARK t2=3.4725
END
````

You can read in these reference values by using the PDB2CONSTANT command as follows:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pathtools-4/epath.pdb
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=5,6,7,8
t1_ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-4/epath.pdb ARG=t1
t2_ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-4/epath.pdb ARG=t2
```

In this case the input must define values with the labels that are being read in from the reference file
and separate PDB2CONSTANT commands are required for reading in `t1` and `t2`.  Furthermore, because the
input PDB file contains multiple frames vectors containing all the values for `t1` and `t2` are output from
the constant commands that are created by the shortcuts in the above input.  If you want to read only one of the
configurations in the input PDB file you can use a pdb with a single frame or the `NUMBER` keyword described above.

If, for any reason, you want to read data from a PDB file that is not a reference value for one of the values defined in
your PLUMED input file you use the NOARGS flag as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-pathtools-4/epath.pdb
t1_ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-4/epath.pdb NOARGS ARG=t1
t2_ref: PDB2CONSTANT REFERENCE=regtest/mapping/rt-pathtools-4/epath.pdb NOARGS ARG=t2
```

## Occupancy and beta factors

PLUMED also reads the occupancy and beta factors from the input PDB files. However, these columns of data are
given a very special meaning.
In cases where the PDB structure is used as a reference for an alignment (that's the case
for instance in [RMSD](RMSD.md) and in [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md)), the occupancy column is used
to provide the weight of each atom in the alignment. In cases where, perhaps after alignment,
the displacement between running coordinates and the provided PDB is computed, the beta factors
are used as weight for the displacement.
Since setting the weights to zero is the same as __not__ including an atom in the alignment or
displacement calculation, the two following reference files would be equivalent when used in an [RMSD](RMSD.md)
calculation. First file:

````
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  0.00  0.00
````

Second file:

````
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
````

However notice that many extra atoms with zero weight might slow down the calculation, so
removing lines is better than setting their weights to zero.
In addition, weights for alignment need not to be equivalent to weights for displacement.
Starting with PLUMED 2.7, if all the weights are set to zero they will be normalized to be equal to the
inverse of the number of involved atoms. This means that it will be possible to use files with
the weight columns set to zero obtaining a meaningful result. In previous PLUMED versions,
setting all weights to zero was resulting in an error instead.

## Systems with more than 100k atoms

Notice that it very likely does not make any sense to compute the [RMSD](RMSD.md) or any other structural
deviation __using__ many atoms. However, if the protein for which you want to compute [RMSD](RMSD.md)
has atoms with large serial numbers (e.g. because it is located __after__ solvent in the sorted list of atoms)
you might end up with troubles with the limitations of the PDB format. Indeed, since there are 5
columns available for atom serial number, this number cannot be larger than 99999.
In addition, providing [MOLINFO](MOLINFO.md) with names associated to atoms with a serial larger than 99999 would be impossible.

Since PLUMED 2.4 we allow the [hybrid 36](http://cci.lbl.gov/hybrid_36/) format to be used to specify atom numbers.
This format is not particularly widespread, but has the nice feature that it provides a one-to-one mapping
between numbers up to approximately 80 millions and strings with 5 characters, plus it is backward compatible
for numbers smaller than 100000. This is not true for notations like the hex notation exported by VMD.
Using the hybrid 36 format, the ATOM records for atom ranging from 99997 to 100002 would read like these:

````
ATOM  99997  Ar      X   1      45.349  38.631  15.116  1.00  1.00
ATOM  99998  Ar      X   1      46.189  38.631  15.956  1.00  1.00
ATOM  99999  Ar      X   1      46.189  39.471  15.116  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  39.471  15.956  1.00  1.00
ATOM  A0000  Ar      X   1      45.349  38.631  16.796  1.00  1.00
ATOM  A0001  Ar      X   1      46.189  38.631  17.636  1.00  1.00
````

There are tools that can be found to translate from integers to strings and back using hybrid 36 format
(a simple python script can be found [here](https://sourceforge.net/p/cctbx/code/HEAD/tree/trunk/iotbx/pdb/hybrid_36.py)).
In addition, as of PLUMED 2.5, we provide a command line tool that can be used to renumber atoms in a PDB file.

*/
//+ENDPLUMEDOC

class PDB2Constant : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit PDB2Constant(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PDB2Constant,"PDB2CONSTANT")

void PDB2Constant::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure");
  keys.add("compulsory","NUMBER","0","if there are multiple structures in the pdb file you can specify that you want the RMSD from a specific structure by specifying its place in the file here. If NUMBER=0 then the RMSD from all structures are computed");
  keys.addFlag("NOARGS",false,"the arguments that are being read from the PDB file are not in the plumed input");
  keys.addInputKeyword("optional","ARG","scalar/vector","read this single argument from the input rather than the atomic structure");
  keys.setValueDescription("scalar/vector/matrix","a value that is constructed from the information in the PDB file");
  keys.needsAction("CONSTANT");
}

PDB2Constant::PDB2Constant(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string input;
  parse("REFERENCE",input);
  unsigned frame;
  parse("NUMBER",frame);
  bool noargs=false;
  std::vector<std::string> argn;
  parseVector("ARG",argn);
  std::vector<Value*> theargs;
  if( argn.size()>0 ) {
    parseFlag("NOARGS",noargs);
    if( !noargs ) {
      ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, theargs );
    } else if( argn.size()>1 ) {
      error("can only read one argument at a time from input pdb file");
    } else {
      log.printf("  reading argument %s from file \n", argn[0].c_str() );
    }
  }
  if( theargs.size()>1 ) {
    error("can only read one argument at a time from input pdb file");
  }

  FILE* fp=std::fopen(input.c_str(),"r");
  bool do_read=true;
  std::vector<double> vals;
  if(!fp) {
    plumed_merror("could not open reference file " + input);
  }
  unsigned natoms=0, nframes=0;

  while ( do_read ) {
    PDB mypdb;
    do_read=mypdb.readFromFilepointer(fp,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
    if( !do_read && nframes>0 ) {
      break ;
    }

    if( natoms==0 ) {
      natoms = mypdb.getPositions().size();
    } else if( mypdb.getPositions().size()!=natoms ) {
      plumed_merror("mismatch between sizes of reference configurations");
    }

    if( nframes+1==frame || frame==0 ) {
      std::vector<double> align( mypdb.getOccupancy() );
      double asum=0;
      for(unsigned i=0; i<align.size(); ++i) {
        asum += align[i];
      }
      if( asum>epsilon ) {
        double iasum = 1 / asum;
        for(unsigned i=0; i<align.size(); ++i) {
          align[i] *= iasum;
        }
      } else if( mypdb.size()>0 ) {
        double iasum = 1 / mypdb.size();
        for(unsigned i=0; i<align.size(); ++i) {
          align[i] = iasum;
        }
      }
      Vector center;
      center.zero();
      for(unsigned i=0; i<mypdb.getPositions().size(); ++i) {
        center += align[i]*mypdb.getPositions()[i];
      }

      if( theargs.size()==0 && argn.size()==0 ) {
        for(unsigned j=0; j<3; ++j) {
          for(unsigned i=0; i<mypdb.getPositions().size(); ++i) {
            vals.push_back( mypdb.getPositions()[i][j] - center[j] );
          }
        }
      } else if( noargs ) {
        std::vector<double> argvals( 1 );
        if( !mypdb.getArgumentValue(argn[0], argvals ) ) {
          error("argument " + argn[0] + " was not set in pdb input");
        }
        vals.push_back( argvals[0] );
      } else {
        std::vector<double> argvals( theargs[0]->getNumberOfValues() );
        if( !mypdb.getArgumentValue(theargs[0]->getName(), argvals ) ) {
          error("argument " + theargs[0]->getName() + " was not set in pdb input");
        }
        for(unsigned i=0; i<argvals.size(); ++i) {
          vals.push_back( argvals[i] );
        }
      }
    }
    nframes++;
  }
  if( frame>0 ) {
    nframes=1;
  }
  std::fclose(fp);
  std::string rnum;
  plumed_assert( vals.size()>0 );
  Tools::convert( vals[0], rnum );
  std::string valstr = " VALUES=" + rnum;
  for(unsigned i=1; i<vals.size(); ++i) {
    Tools::convert( vals[i], rnum );
    valstr += "," + rnum;
  }
  if( vals.size()>nframes ) {
    std::string nc, nr;
    Tools::convert( nframes, nr );
    Tools::convert( vals.size()/nframes, nc );
    readInputLine( getShortcutLabel() + ": CONSTANT NROWS=" + nr + " NCOLS=" + nc + valstr );
  } else {
    readInputLine( getShortcutLabel() + ": CONSTANT" + valstr );
  }
}

}
}
