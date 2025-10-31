/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "tools/Units.h"
#include "tools/CheckInRange.h"
#include "core/GenericMolInfo.h"
#include "core/ActionSet.h"
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"


namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPATOMS
/*
Dump selected atoms on a file.

This command can be used to output the positions of a particular set of atoms.
For example, the following input instructs PLUMED to print out the positions
of atoms 1-10 together with the position of the center of mass of atoms 11-20 every
10 steps to a file called file.xyz.

```plumed
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
```

By default, the coordinates in the output xyz file are in nm but you can change these units
by using the `UNITS` keyword as shown below:

```plumed
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1 UNITS=A
```

or by using the [UNITS](UNITS.md) action as shown below:

```plumed
UNITS LENGTH=A
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
```

Notice, however, that if you use the second option here all the quantitities with units of length in your input
file must be provided in Angstrom and not nm.

## DUMPATOMS and WHOLEMOLECULES

The commands [WHOLEMOLECULES](WHOLEMOLECULES.md), [WRAPAROUND](WRAPAROUND.md), [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md)
and [RESET_CELL](RESET_CELL.md) all edit the global positions of the atoms.  If you use an input like this one:

```plumed
DUMPATOMS ATOMS=1-10 FILE=file.xyz
WHOLEMOLECULES ENTITY0=1-10
```

then the positions of the atoms that were passed to PLUMED by the MD code are output.  However, if you use an input
like this one:

```plumed
WHOLEMOLECULES ENTITY0=1-10
DUMPATOMS ATOMS=1-10 FILE=file.xyz
```

the positions outputted by the DUMPATOMS command will have been editted by the [WHOLEMOLECULES](WHOLEMOLECULES.md) command.

## Outputting other file types

The extension that is given to the file specified using the `FILE` keyword determines the output file type. Hence,
the following example will output a gro file rather than an xyz file:

```plumed
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.gro ATOMS=1-10,c1
```

You can also enforce the output file type by using the `TYPE` keyword as shown below:

```plumed
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1 TYPE=gro
FLUSH STRIDE=1
```

Notice that DUMPATOMS command here outputs the atoms in the gro-file format even though the author of this input has used the xyz extension.
Further note that by using the [FLUSH](FLUSH.md) we force PLUMED to output flush all the open files every step and not to store output
data in a buffer before printing it to the output files.

Outputting the atomic positions using the gro file format is particularly advantageous if you also have a [MOLINFO](MOLINFO.md) command in
your input file as shown below:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
# this is required to have proper atom names:
MOLINFO STRUCTURE=regtest/basic/rt32/helix.pdb
# if omitted, atoms will have "X" name...

c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.gro ATOMS=1-10,c1
# notice that last atom is a virtual one and will not have
# a correct name in the resulting gro file
```

The reason that using the gro file format is advantageous in this case is that PLUMED will also output the atom and residue names
for the non-virtual atoms.  PLUMED is able to do this in this case as it is able to use the information that was read in from the
pdb file that was provided to the [MOLINFO](MOLINFO.md) command.

If PLUMED has been compiled with xdrfile support, then PLUMED
can output xtc and trr files as well as xyz and gro files.  If you want to use these output types you should install the xdrfile
library by following the instructions [here](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library).
If the xdrfile library is installed properly the PLUMED configure script should be able to
detect it and enable it.  The following example shows how you can use DUMPATOMS to output an xtc file:

```plumed
c1: COM ATOMS=11-20
DUMPATOMS STRIDE=10 FILE=file.xtc ATOMS=1-10,c1
```

The xtc file that is output by this command will be significantly smaller than a gro and xyz file.

Finally, notice that gro and xtc file store coordinates with limited precision set by the
`PRECISION` keyword. The default value is 3, which means "3 digits after dot" in nm (1/1000 of a nm).
The following will write a larger xtc file with high resolution coordinates:

```plumed
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xtc ATOMS=1-10,c1 PRECISION=7
```

## Outputting atomic positions and vectors

The atoms section of an xyz file normally contains four columns of data - a symbol that tells you the atom type
and then three columns containing the atom's $x$, $y$ and $z$ coordinates.  PLUMED allows you to output more columns
of data in this file.  For example, the following input outputs five columns of data.  The first four columns of data
are the usual ones you would expect in an xyz file, while the fifth contains the coordination numbers for each of the
atom that have been calculated using PLUMED

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
DUMPATOMS ATOMS=1-100 ARG=cc FILE=file.xyz
```

This command is used in the shortcut that recovered the old [DUMPMULTICOLVAR](DUMPMULTICOLVAR.md) command. This new version of
the command is better, however, as you can output more than one vector of symmetry functions at once as is demonstrated by the following
input that outputs the coordination numbers and the values that were obtained when the coordination numbers were all transformed by a
switching function:

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
fcc: LESS_THAN ARG=cc SWITCH={RATIONAL R_0=4}
DUMPATOMS ATOMS=1-100 ARG=cc,fcc FILE=file.xyz
```

Notice that when we use an `ARG` keyword we can also use DUMPATOMS to only print those atoms whose corresponding element in the
the input vector satisfies a certain criteria.  For example the input file below only prints the positions (and coordination numbers) of atoms that
have a coordination number that is greater than or equal to 4.

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
DUMPATOMS ATOMS=1-100 ARG=cc GREATER_THAN_OR_EQUAL=4 FILE=file.xyz
```

Alternatively, the following input allows us to output those atoms that have a coordination number that is between 4 and 6:

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
DUMPATOMS ...
  ATOMS=1-100 ARG=cc
  GREATER_THAN_OR_EQUAL=4
  LESS_THAN_OR_EQUAL=6
  FILE=file.xyz
...
```

Commands like these are useful if you want to print the coordinates of the atom that are in a paricular cluster that has been identified using
the [DFSCLUSTERING](DFSCLUSTERING.md) command.

__You can only use the ARG keyword if you are outputting an xyz file__

## PRINT and RESTART

If you run a calculation with the following input:

```plumed
DUMPATOMS ATOMS=1-100 FILE=file.xyz
```

and a file called `file.xyz` is already present in the directory where the calculation is running, the existing file is backed up
and renamed to `bck.0.file.xyz` so that new data can be output to a new file called `file.xyz`.  If you would like to append to the
existing file you can use the RESTART command as shown below:

```plumed
DUMPATOMS ATOMS=1-100 FILE=file.xyz RESTART=YES
```

You can achieve the same result by using the [RESTART](RESTART.md) action as shown below:

```plumed
RESTART
DUMPATOMS ATOMS=1-100 FILE=file.xyz
```

However, the advantage of using the RESTART keyword is that you can apped to some files and back up others as illustrated below:

```plumed
DUMPATOMS ATOMS=1-100 FILE=file1.xyz
DUMPATOMS ATOMS=1-100 FILE=file2.xyz RESTART=YES
```

If you use the input above the file `file1.xyz` is backed up, while new data will be appended to the file `file2.xyz`.  If you use the
[RESTART](RESTART.md) action instead data will be appended to both colvar files.

## Switching printing on and off

You can use the UPDATE_FROM and UPDATE_UNTIL flags to make the DUMPATOMS command only output data at certain points during the trajectory.
To see how this works consider the following example:

```plumed
DUMPATOMS ATOMS=1-100 FILE=file.xyz UPDATE_FROM=100 UPDATE_UNTIL=500 STRIDE=1
```

During the first 100 ps of a simulation with this input the atomic positions are not output to the file called file.xyz.
The positions are instead only output after the first 100 ps of trajectory have elapsed.  Furthermore, output of the positions stops
once the trajectory is longer than 500 ps. In other words, the positions are only output during the 400 ps time interval after the first
100 ps of the simulation.

*/
//+ENDPLUMEDOC

class DumpAtoms:
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionPilot {
  OFile of;
  double lenunit;
  int iprecision;
  CheckInRange bounds;
  std::vector<std::string> names;
  std::vector<unsigned>    residueNumbers;
  std::vector<std::string> residueNames;
  std::string type;
  std::string fmt_gro_pos;
  std::string fmt_gro_box;
  std::string fmt_xyz;
  xdrfile::XDRFILE* xd;
public:
  explicit DumpAtoms(const ActionOptions&);
  ~DumpAtoms();
  static void registerKeywords( Keywords& keys );
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;
  bool actionHasForces() override {
    return false;
  }
  void lockRequests() override;
  void unlockRequests() override;
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpAtoms,"DUMPATOMS")

void DumpAtoms::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("optional","ARG","vector","the labels of vectors that should be output in the xyz file. The number of elements in the vector should equal the number of atoms output");
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates; extension is automatically detected");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
  keys.add("optional", "TYPE","file type, either xyz, gro, xtc, or trr, can override an automatically detected file extension");
  keys.add("optional","LESS_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value less than or equal to this value");
  keys.add("optional","GREATER_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value greater than or equal to this value");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

DumpAtoms::DumpAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  iprecision(3) {
  std::vector<AtomNumber> atoms;
  std::string file;
  parse("FILE",file);
  if(file.length()==0) {
    error("name out output file was not specified");
  }
  type=Tools::extension(file);
  log<<"  file name "<<file<<"\n";
  if(type=="gro" || type=="xyz" || type=="xtc" || type=="trr") {
    log<<"  file extension indicates a "<<type<<" file\n";
  } else {
    log<<"  file extension not detected, assuming xyz\n";
    type="xyz";
  }
  std::string ntype;
  parse("TYPE",ntype);
  if(ntype.length()>0) {
    if(ntype!="xyz" && ntype!="gro" && ntype!="xtc" && ntype!="trr"
      ) {
      error("TYPE cannot be understood");
    }
    log<<"  file type enforced to be "<<ntype<<"\n";
    type=ntype;
  }

  fmt_gro_pos="%8.3f";
  fmt_gro_box="%12.7f";
  fmt_xyz="%f";

  std::string precision;
  parse("PRECISION",precision);
  if(precision.length()>0) {
    Tools::convert(precision,iprecision);
    log<<"  with precision "<<iprecision<<"\n";
    std::string a,b;
    Tools::convert(iprecision+5,a);
    Tools::convert(iprecision,b);
    fmt_gro_pos="%"+a+"."+b+"f";
    fmt_gro_box=fmt_gro_pos;
    fmt_xyz=fmt_gro_box;
  }

  parseAtomList("ATOMS",atoms);

  std::string unitname;
  parse("UNITS",unitname);
  if(unitname!="PLUMED") {
    Units myunit;
    myunit.setLength(unitname);
    if(myunit.getLength()!=1.0 && type=="gro") {
      error("gro files should be in nm");
    }
    if(myunit.getLength()!=1.0 && type=="xtc") {
      error("xtc files should be in nm");
    }
    if(myunit.getLength()!=1.0 && type=="trr") {
      error("trr files should be in nm");
    }
    lenunit=getUnits().getLength()/myunit.getLength();
  } else if(type=="gro" || type=="xtc" || type=="trr") {
    lenunit=getUnits().getLength();
  } else {
    lenunit=1.0;
  }

  of.link(*this);
  of.open(file);
  std::string path=of.getPath();
  log<<"  Writing on file "<<path<<"\n";
  std::string mode=of.getMode();
  if(type=="xtc") {
    of.close();
    xd=xdrfile::xdrfile_open(path.c_str(),mode.c_str());
  } else if(type=="trr") {
    of.close();
    xd=xdrfile::xdrfile_open(path.c_str(),mode.c_str());
  }
  log.printf("  printing the following atoms in %s :", unitname.c_str() );
  for(unsigned i=0; i<atoms.size(); ++i) {
    log.printf(" %d",atoms[i].serial() );
  }
  log.printf("\n");

  if( getNumberOfArguments()>0 ) {
    if( type!="xyz" ) {
      error("can only print atomic properties when outputting xyz files");
    }

    std::vector<std::string> argnames;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()!=1 ) {
        error("arguments for xyz output should be vectors");
      }
      if( getPntrToArgument(i)->getNumberOfValues()!=atoms.size() ) {
        error("number of elements in vector " + getPntrToArgument(i)->getName() + " is not equal to number of atoms output");
      }
      argnames.push_back( getPntrToArgument(i)->getName() );
    }
    std::vector<std::string> str_upper, str_lower;
    std::string errors;
    parseVector("LESS_THAN_OR_EQUAL",str_upper);
    parseVector("GREATER_THAN_OR_EQUAL",str_lower);
    if( !bounds.setBounds( getNumberOfArguments(), str_lower, str_upper, errors ) ) {
      error( errors );
    }
    if( bounds.wereSet() ) {
      log.printf("  %s \n", bounds.report( argnames ).c_str() );
    }
  }

  requestAtoms(atoms, false);
  auto* infomoldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( infomoldat ) {
    log<<"  MOLINFO DATA found with label " <<infomoldat->getLabel()<<", using proper atom names\n";
    names.resize(atoms.size());
    for(unsigned i=0; i<atoms.size(); i++)
      if(atoms[i].index()<infomoldat->getPDBsize()) {
        names[i]=infomoldat->getAtomName(atoms[i]);
      }
    residueNumbers.resize(atoms.size());
    for(unsigned i=0; i<residueNumbers.size(); ++i)
      if(atoms[i].index()<infomoldat->getPDBsize()) {
        residueNumbers[i]=infomoldat->getResidueNumber(atoms[i]);
      }
    residueNames.resize(atoms.size());
    for(unsigned i=0; i<residueNames.size(); ++i)
      if(atoms[i].index()<infomoldat->getPDBsize()) {
        residueNames[i]=infomoldat->getResidueName(atoms[i]);
      }
  }
}

void DumpAtoms::calculateNumericalDerivatives( ActionWithValue* a ) {
  plumed_merror("this should never be called");
}

void DumpAtoms::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void DumpAtoms::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void DumpAtoms::update() {
  if(type=="xyz") {
    unsigned nat=0;
    std::vector<double> args( getNumberOfArguments() );
    for(unsigned i=0; i<getNumberOfAtoms(); ++i)  {
      for(unsigned j=0; j<getNumberOfArguments(); ++j) {
        args[j] = getPntrToArgument(j)->get(i);
      }
      if( bounds.check( args ) ) {
        nat++;
      }
    }
    of.printf("%d\n",nat);
    const Tensor & t(getPbc().getBox());
    if(getPbc().isOrthorombic()) {
      of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2));
    } else {
      of.printf((" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),
                lenunit*t(0,0),lenunit*t(0,1),lenunit*t(0,2),
                lenunit*t(1,0),lenunit*t(1,1),lenunit*t(1,2),
                lenunit*t(2,0),lenunit*t(2,1),lenunit*t(2,2)
               );
    }
    const std::string posFormatString="%s "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz;
    const std::string extraFormatString=" "+fmt_xyz;
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      for(unsigned j=0; j<getNumberOfArguments(); ++j) {
        args[j] = getPntrToArgument(j)->get(i);
      }
      if( !bounds.check(args) ) {
        continue;
      }
      const char* defname="X";
      const char* atomName=defname;
      if(names.size()>0) {
        if(names[i].length()>0) {
          atomName=names[i].c_str();
        }
      }
      of.printf(posFormatString.c_str(),
                atomName,
                lenunit*getPosition(i)(0),
                lenunit*getPosition(i)(1),
                lenunit*getPosition(i)(2));
      for(unsigned j=0; j<getNumberOfArguments(); ++j) {
        of.printf(extraFormatString.c_str(), getPntrToArgument(j)->get(i) );
      }
      of.printf("\n");
    }
  } else if(type=="gro") {
    const std::string posFormatString="%5u%-5s%5s%5d"+fmt_gro_pos+fmt_gro_pos+fmt_gro_pos+"\n";
    std::string extraFormatString=" "+fmt_xyz;
    const Tensor & t(getPbc().getBox());
    of.printf("Made with PLUMED t=%f\n",getTime()/getUnits().getTime());
    of.printf("%d\n",getNumberOfAtoms());
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      const char* defname="X";
      const char* atomName=defname;
      unsigned residueNumber=0;
      if(names.size()>0)
        if(names[i].length()>0) {
          atomName=names[i].c_str();
        }
      if(residueNumbers.size()>0) {
        residueNumber=residueNumbers[i];
      }
      std::string resname="";
      if(residueNames.size()>0) {
        resname=residueNames[i];
      }
      of.printf(posFormatString.c_str(),
                residueNumber%100000,resname.c_str(),atomName,getAbsoluteIndex(i).serial()%100000,
                lenunit*getPosition(i)(0),lenunit*getPosition(i)(1),lenunit*getPosition(i)(2));
    }
    of.printf((fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+"\n").c_str(),
              lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2),
              lenunit*t(0,1),lenunit*t(0,2),lenunit*t(1,0),
              lenunit*t(1,2),lenunit*t(2,0),lenunit*t(2,1));
  } else if(type=="xtc" || type=="trr") {
    xdrfile::matrix box;
    const Tensor & t(getPbc().getBox());
    int natoms=getNumberOfAtoms();
    int step=getStep();
    float time=getTime()/getUnits().getTime();
    float precision=Tools::fastpow(10.0,iprecision);
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++) {
        box[i][j]=lenunit*t(i,j);
      }
// here we cannot use a std::vector<rvec> since it does not compile.
// we thus use a std::unique_ptr<rvec[]>
    auto pos = Tools::make_unique<xdrfile::rvec[]>(natoms);
    for(int i=0; i<natoms; i++)
      for(int j=0; j<3; j++) {
        pos[i][j]=lenunit*getPosition(i)(j);
      }
    if(type=="xtc") {
      write_xtc(xd,natoms,step,time,box,&pos[0],precision);
    } else if(type=="trr") {
      write_trr(xd,natoms,step,time,0.0,box,&pos[0],NULL,NULL);
    }
  } else {
    plumed_merror("unknown file type "+type);
  }
}

DumpAtoms::~DumpAtoms() {
  if(type=="xtc") {
    xdrfile_close(xd);
  } else if(type=="trr") {
    xdrfile_close(xd);
  }
}


}
}
