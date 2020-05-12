/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include <memory>
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"

#if defined(__PLUMED_HAS_XDRFILE)
#include <xdrfile/xdrfile_xtc.h>
#include <xdrfile/xdrfile_trr.h>
#endif


using namespace std;

namespace PLMD
{
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPATOMS
/*
Dump selected atoms on a file.

This command can be used to output the positions of a particular set of atoms.
The atoms required are output in a xyz or gro formatted file.
If PLUMED has been compiled with xdrfile support, then also xtc and trr files can be written.
To this aim one should install xdrfile library (http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library).
If the xdrfile library is installed properly the PLUMED configure script should be able to
detect it and enable it.
The type of file is automatically detected from the file extension, but can be also
enforced with TYPE.
Importantly, if your
input file contains actions that edit the atoms position (e.g. \ref WHOLEMOLECULES)
and the DUMPATOMS command appears after this instruction, then the edited
atom positions are output.
You can control the buffering of output using the \ref FLUSH keyword on a separate line.

Units of the printed file can be controlled with the UNITS keyword. By default PLUMED units as
controlled in the \ref UNITS command are used, but one can override it e.g. with UNITS=A.
Notice that gro/xtc/trr files can only contain coordinates in nm.

\par Examples

The following input instructs plumed to print out the positions of atoms
1-10 together with the position of the center of mass of atoms 11-20 every
10 steps to a file called file.xyz.
\plumedfile
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1
\endplumedfile
Notice that the coordinates in the xyz file will be expressed in nm, since these
are the defaults units in PLUMED. If you want the xyz file to be expressed in A, you should use the
following input
\plumedfile
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xyz ATOMS=1-10,c1 UNITS=A
\endplumedfile
As an alternative, you might want to set all the length used by PLUMED to Angstrom using the \ref UNITS
action. However, this latter choice will affect all your input and output.

The following input is very similar but dumps a .gro (gromacs) file,
which also contains atom and residue names.
\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
# this is required to have proper atom names:
MOLINFO STRUCTURE=reference.pdb
# if omitted, atoms will have "X" name...

COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.gro ATOMS=1-10,c1
# notice that last atom is a virtual one and will not have
# a correct name in the resulting gro file
\endplumedfile

The `file.gro` will contain coordinates expressed in nm, since this is the convention for gro files.

In case you have compiled PLUMED with `xdrfile` library, you might even write xtc or trr files as follows
\plumedfile
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xtc ATOMS=1-10,c1
\endplumedfile
Notice that xtc files are significantly smaller than gro and xyz files.

Finally, consider that gro and xtc file store coordinates with limited precision set by the
`PRECISION` keyword. Default value is 3, which means "3 digits after dot" in nm (1/1000 of a nm).
The following will write a larger xtc file with high resolution coordinates:
\plumedfile
COM ATOMS=11-20 LABEL=c1
DUMPATOMS STRIDE=10 FILE=file.xtc ATOMS=1-10,c1 PRECISION=7
\endplumedfile



*/
//+ENDPLUMEDOC

class DumpAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  OFile of;
  double lenunit;
  int iprecision;
  std::vector<std::string> names;
  std::vector<unsigned>    residueNumbers;
  std::vector<std::string> residueNames;
  std::string type;
  std::string fmt_gro_pos;
  std::string fmt_gro_box;
  std::string fmt_xyz;
#if defined(__PLUMED_HAS_XDRFILE)
  XDRFILE* xd;
#endif
public:
  explicit DumpAtoms(const ActionOptions&);
  ~DumpAtoms();
  static void registerKeywords( Keywords& keys );
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpAtoms,"DUMPATOMS")

void DumpAtoms::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "FILE", "file on which to output coordinates; extension is automatically detected");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
#if defined(__PLUMED_HAS_XDRFILE)
  keys.add("optional", "TYPE","file type, either xyz, gro, xtc, or trr, can override an automatically detected file extension");
#else
  keys.add("optional", "TYPE","file type, either xyz or gro, can override an automatically detected file extension");
#endif
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

DumpAtoms::DumpAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao),
  iprecision(3)
{
  vector<AtomNumber> atoms;
  string file;
  parse("FILE",file);
  if(file.length()==0) error("name out output file was not specified");
  type=Tools::extension(file);
  log<<"  file name "<<file<<"\n";
  if(type=="gro" || type=="xyz" || type=="xtc" || type=="trr") {
    log<<"  file extension indicates a "<<type<<" file\n";
  } else {
    log<<"  file extension not detected, assuming xyz\n";
    type="xyz";
  }
  string ntype;
  parse("TYPE",ntype);
  if(ntype.length()>0) {
    if(ntype!="xyz" && ntype!="gro" && ntype!="xtc" && ntype!="trr"
      ) error("TYPE cannot be understood");
    log<<"  file type enforced to be "<<ntype<<"\n";
    type=ntype;
  }
#ifndef __PLUMED_HAS_XDRFILE
  if(type=="xtc" || type=="trr") error("types xtc and trr require PLUMED to be linked with the xdrfile library. Please install it and recompile PLUMED.");
#endif

  fmt_gro_pos="%8.3f";
  fmt_gro_box="%12.7f";
  fmt_xyz="%f";

  string precision;
  parse("PRECISION",precision);
  if(precision.length()>0) {
    Tools::convert(precision,iprecision);
    log<<"  with precision "<<iprecision<<"\n";
    string a,b;
    Tools::convert(iprecision+5,a);
    Tools::convert(iprecision,b);
    fmt_gro_pos="%"+a+"."+b+"f";
    fmt_gro_box=fmt_gro_pos;
    fmt_xyz=fmt_gro_box;
  }

  parseAtomList("ATOMS",atoms);

  std::string unitname; parse("UNITS",unitname);
  if(unitname!="PLUMED") {
    Units myunit; myunit.setLength(unitname);
    if(myunit.getLength()!=1.0 && type=="gro") error("gro files should be in nm");
    if(myunit.getLength()!=1.0 && type=="xtc") error("xtc files should be in nm");
    if(myunit.getLength()!=1.0 && type=="trr") error("trr files should be in nm");
    lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
  } else if(type=="gro" || type=="xtc" || type=="trr") lenunit=plumed.getAtoms().getUnits().getLength();
  else lenunit=1.0;

  checkRead();
  of.link(*this);
  of.open(file);
  std::string path=of.getPath();
  log<<"  Writing on file "<<path<<"\n";
#ifdef __PLUMED_HAS_XDRFILE
  std::string mode=of.getMode();
  if(type=="xtc") {
    of.close();
    xd=xdrfile_open(path.c_str(),mode.c_str());
  } else if(type=="trr") {
    of.close();
    xd=xdrfile_open(path.c_str(),mode.c_str());
  }
#endif
  log.printf("  printing the following atoms in %s :", unitname.c_str() );
  for(unsigned i=0; i<atoms.size(); ++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
  std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()==1 ) {
    log<<"  MOLINFO DATA found, using proper atom names \n";
    names.resize(atoms.size(),"");
    for(unsigned i=0; i<atoms.size(); i++) if(atoms[i].index()<moldat[0]->getPDBsize()) names[i]=moldat[0]->getAtomName(atoms[i]);
    residueNumbers.resize(atoms.size());
    for(unsigned i=0; i<atoms.size(); ++i) if(atoms[i].index()<moldat[0]->getPDBsize()) residueNumbers[i]=moldat[0]->getResidueNumber(atoms[i]);
    residueNames.resize(atoms.size());
    for(unsigned i=0; i<atoms.size(); ++i) if(atoms[i].index()<moldat[0]->getPDBsize()) residueNames[i]=moldat[0]->getResidueName(atoms[i]);
  }
}

void DumpAtoms::update() {
  if(type=="xyz") {
    of.printf("%d\n",getNumberOfAtoms());
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
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      const char* defname="X";
      const char* name=defname;
      if(names.size()>0) if(names[i].length()>0) name=names[i].c_str();
      of.printf(("%s "+fmt_xyz+" "+fmt_xyz+" "+fmt_xyz+"\n").c_str(),name,lenunit*getPosition(i)(0),lenunit*getPosition(i)(1),lenunit*getPosition(i)(2));
    }
  } else if(type=="gro") {
    const Tensor & t(getPbc().getBox());
    of.printf("Made with PLUMED t=%f\n",getTime()/plumed.getAtoms().getUnits().getTime());
    of.printf("%d\n",getNumberOfAtoms());
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      const char* defname="X";
      const char* name=defname;
      unsigned residueNumber=0;
      if(names.size()>0) if(names[i].length()>0) name=names[i].c_str();
      if(residueNumbers.size()>0) residueNumber=residueNumbers[i];
      std::string resname="";
      if(residueNames.size()>0) resname=residueNames[i];
      of.printf(("%5u%-5s%5s%5d"+fmt_gro_pos+fmt_gro_pos+fmt_gro_pos+"\n").c_str(),
                residueNumber%100000,resname.c_str(),name,getAbsoluteIndex(i).serial()%100000,
                lenunit*getPosition(i)(0),lenunit*getPosition(i)(1),lenunit*getPosition(i)(2));
    }
    of.printf((fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+" "+fmt_gro_box+"\n").c_str(),
              lenunit*t(0,0),lenunit*t(1,1),lenunit*t(2,2),
              lenunit*t(0,1),lenunit*t(0,2),lenunit*t(1,0),
              lenunit*t(1,2),lenunit*t(2,0),lenunit*t(2,1));
#if defined(__PLUMED_HAS_XDRFILE)
  } else if(type=="xtc" || type=="trr") {
    matrix box;
    const Tensor & t(getPbc().getBox());
    int natoms=getNumberOfAtoms();
    int step=getStep();
    float time=getTime()/plumed.getAtoms().getUnits().getTime();
    float precision=Tools::fastpow(10.0,iprecision);
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) box[i][j]=lenunit*t(i,j);
    std::unique_ptr<rvec[]> pos(new rvec [natoms]);
    for(int i=0; i<natoms; i++) for(int j=0; j<3; j++) pos[i][j]=lenunit*getPosition(i)(j);
    if(type=="xtc") {
      write_xtc(xd,natoms,step,time,box,&pos[0],precision);
    } else if(type=="trr") {
      write_trr(xd,natoms,step,time,0.0,box,&pos[0],NULL,NULL);
    }
#endif
  } else plumed_merror("unknown file type "+type);
}

DumpAtoms::~DumpAtoms() {
#ifdef __PLUMED_HAS_XDRFILE
  if(type=="xtc") {
    xdrfile_close(xd);
  } else if(type=="trr") {
    xdrfile_close(xd);
  }
#endif
}


}
}
