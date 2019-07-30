/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2017 of Pipolo Silvio and Fabio Pietrucci.

The piv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The piv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

#include <string>
#include <cmath>
#include <iostream>

using namespace std;

namespace PLMD
{
namespace piv
{

//+PLUMEDOC PIVMOD_COLVAR PIV
/*
Calculates the PIV-distance.

PIV distance is the squared Cartesian distance between the PIV \cite gallet2013structural \cite pipolo2017navigating
associated to the configuration of the system during the dynamics and a reference configuration provided
as input (PDB file format).
PIV can be used together with \ref FUNCPATHMSD to define a path in the PIV space.

\par Examples

The following example calculates PIV-distances from three reference configurations in Ref1.pdb, Ref2.pdb and Ref3.pdb
and prints the results in a file named colvar.
Three atoms (PIVATOMS=3) with names (pdb file) A B and C are used to construct the PIV and all PIV blocks (AA, BB, CC, AB, AC, BC) are considered.
SFACTOR is a scaling factor that multiplies the contribution to the PIV-distance given by the single PIV block.
NLIST sets the use of neighbor lists for calculating atom-atom distances.
The SWITCH keyword specifies the parameters of the switching function that transforms atom-atom distances.
SORT=1 means that the PIV block elements are sorted (SORT=0 no sorting.)
Values for SORT, SFACTOR and the neighbor list parameters have to be specified for each block.
The order is the following: AA,BB,CC,AB,AC,BC. If ONLYDIRECT (ONLYCROSS) is used the order is AA,BB,CC (AB,AC,BC).
The sorting operation within each PIV block is performed using the counting sort algorithm, PRECISION specifies the size of the counting array.

\plumedfile
PIV ...
LABEL=Pivd1
PRECISION=1000
NLIST
REF_FILE=Ref1.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV
PIV ...
LABEL=Pivd2
PRECISION=1000
NLIST
REF_FILE=Ref2.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV
PIV ...
LABEL=Pivd3
PRECISION=1000
NLIST
REF_FILE=Ref3.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV

PRINT ARG=Pivd1,Pivd2,Pivd3 FILE=colvar
\endplumedfile

WARNING:
Both the "CRYST" and "ATOM" lines of the PDB files must conform precisely to the official pdb format, including the width of each alphanumerical field:

\verbatim
CRYST1   31.028   36.957   23.143  89.93  92.31  89.99 P 1           1
ATOM      1  OW1 wate    1      15.630  19.750   1.520  1.00  0.00
\endverbatim

In each pdb frame, atoms must be numbered in the same order and with the same element symbol as in the input of the MD program.

The following example calculates the PIV-distances from two reference configurations Ref1.pdb and Ref2.pdb
and uses PIV-distances to define a Path Collective Variable (\ref FUNCPATHMSD) with only two references (Ref1.pdb and Ref2.pdb).
With the VOLUME keyword one scales the atom-atom distances by the cubic root of the ratio between the specified value and the box volume of the initial step of the trajectory file.

\plumedfile
PIV ...
LABEL=c1
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref1.pdb
PIVATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.5 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV
PIV ...
LABEL=c2
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref2.pdb
PIVATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV

p1: FUNCPATHMSD ARG=c1,c2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=c1,c2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIV please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

class PIV      : public Colvar
{
private:
  bool pbc, serial, timer;
  ForwardDecl<Stopwatch> stopwatch_fwd;
  Stopwatch& stopwatch=*stopwatch_fwd;
  int updatePIV;
  unsigned Nprec,Natm,Nlist,NLsize;
  double Fvol,Vol0,m_PIVdistance;
  std::string ref_file;
  NeighborList *nlall;
  std::vector<SwitchingFunction> sfs;
  std::vector<std:: vector<double> > rPIV;
  std::vector<double> scaling,r00;
  std::vector<double> nl_skin;
  std::vector<double> fmass;
  std::vector<bool> dosort;
  std::vector<Vector> compos;
  std::vector<string> sw;
  std::vector<NeighborList *> nl;
  std::vector<NeighborList *> nlcom;
  std::vector<Vector> m_deriv;
  Tensor m_virial;
  bool Svol,cross,direct,doneigh,test,CompDer,com;
public:
  static void registerKeywords( Keywords& keys );
  explicit PIV(const ActionOptions&);
  ~PIV();
  // active methods:
  virtual void calculate();
  void checkFieldsAllowed() {}
};

PLUMED_REGISTER_ACTION(PIV,"PIV")

void PIV::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered","SWITCH","The switching functions parameter."
           "You should specify a Switching function for all PIV blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("compulsory","PRECISION","the precision for approximating reals with integers in sorting.");
  keys.add("compulsory","REF_FILE","PDB file name that contains the \\f$i\\f$th reference structure.");
  keys.add("compulsory","PIVATOMS","Number of atoms to use for PIV.");
  keys.add("compulsory","SORT","Whether to sort or not the PIV block.");
  keys.add("compulsory","ATOMTYPES","The atom types to use for PIV.");
  keys.add("optional","SFACTOR","Scale the PIV-distance by such block-specific factor");
  keys.add("optional","VOLUME","Scale atom-atom distances by the cubic root of the cell volume. The input volume is used to scale the R_0 value of the switching function. ");
  keys.add("optional","UPDATEPIV","Frequency (in steps) at which the PIV is updated.");
  keys.addFlag("TEST",false,"Print the actual and reference PIV and exit");
  keys.addFlag("COM",false,"Use centers of mass of groups of atoms instead of atoms as specified in the Pdb file");
  keys.addFlag("ONLYCROSS",false,"Use only cross-terms (A-B, A-C, B-C, ...) in PIV");
  keys.addFlag("ONLYDIRECT",false,"Use only direct-terms (A-A, B-B, C-C, ...) in PIV");
  keys.addFlag("DERIVATIVES",false,"Activate the calculation of the PIV for every class (needed for numerical derivatives).");
  keys.addFlag("NLIST",false,"Use a neighbor list for distance calculations.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("TIMER",false,"Perform timing analysis on heavy loops.");
  keys.add("optional","NL_CUTOFF","Neighbor lists cutoff.");
  keys.add("optional","NL_STRIDE","Update neighbor lists every NL_STRIDE steps.");
  keys.add("optional","NL_SKIN","The maximum atom displacement tolerated for the neighbor lists update.");
  keys.reset_style("SWITCH","compulsory");
}

PIV::PIV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  timer(false),
  updatePIV(1),
  Nprec(1000),
  Natm(1),
  Nlist(1),
  NLsize(1),
  Fvol(1.),
  Vol0(0.),
  m_PIVdistance(0.),
  rPIV(std:: vector<std:: vector<double> >(Nlist)),
  scaling(std:: vector<double>(Nlist)),
  r00(std:: vector<double>(Nlist)),
  nl_skin(std:: vector<double>(Nlist)),
  fmass(std:: vector<double>(Nlist)),
  dosort(std:: vector<bool>(Nlist)),
  compos(std:: vector<Vector>(NLsize)),
  sw(std:: vector<string>(Nlist)),
  nl(std:: vector<NeighborList *>(Nlist)),
  nlcom(std:: vector<NeighborList *>(NLsize)),
  m_deriv(std:: vector<Vector>(1)),
  Svol(false),
  cross(true),
  direct(true),
  doneigh(false),
  test(false),
  CompDer(false),
  com(false)
{
  log << "Starting PIV Constructor\n";

  // Precision on the real-to-integer transformation for the sorting
  parse("PRECISION",Nprec);
  if(Nprec<2) error("Precision must be => 2");

  // PBC
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }

  // SERIAL/PARALLEL
  parseFlag("SERIAL",serial);
  if(serial) {
    log << "Serial PIV construction\n";
  } else     {
    log << "Parallel PIV construction\n";
  }

  // Derivatives
  parseFlag("DERIVATIVES",CompDer);
  if(CompDer) log << "Computing Derivatives\n";

  // Timing
  parseFlag("TIMER",timer);
  if(timer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }

  // Test
  parseFlag("TEST",test);

  // UPDATEPIV
  if(keywords.exists("UPDATEPIV")) {
    parse("UPDATEPIV",updatePIV);
  }

  // Test
  parseFlag("COM",com);
  if(com) log << "Building PIV using COMs\n";

  // Volume Scaling
  parse("VOLUME",Vol0);
  if (Vol0>0) {
    Svol=true;
  }

  // PIV direct and cross blocks
  bool oc=false,od=false;
  parseFlag("ONLYCROSS",oc);
  parseFlag("ONLYDIRECT",od);
  if (oc&&od) {
    error("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  if(oc) {
    direct=false;
    log << "Using only CROSS-PIV blocks\n";
  }
  if(od) {
    cross=false;
    log << "Using only DIRECT-PIV blocks\n";
  }

  // Atoms for PIV
  parse("PIVATOMS",Natm);
  std:: vector<string> atype(Natm);
  parseVector("ATOMTYPES",atype);
  //if(atype.size()!=getNumberOfArguments() && atype.size()!=0) error("not enough values for ATOMTYPES");

  // Reference PDB file
  parse("REF_FILE",ref_file);
  PDB mypdb;
  FILE* fp=fopen(ref_file.c_str(),"r");
  if (fp!=NULL) {
    log<<"Opening PDB file with reference frame: "<<ref_file.c_str()<<"\n";
    mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    fclose (fp);
  } else {
    error("Error in reference PDB file");
  }

  // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
  // Atomlist or Plist used to build pair lists
  std:: vector<std:: vector<AtomNumber> > Plist(Natm);
  // Atomlist used to build list of atoms for each COM
  std:: vector<std:: vector<AtomNumber> > comatm(1);
  // NLsize is the number of atoms in the pdb cell
  NLsize=mypdb.getAtomNumbers().size();
  // In the following P stands for Point (either an Atom or a COM)
  unsigned resnum=0;
  // Presind (array size: number of residues) contains the contains the residue number
  //   this is because the residue numbers may not alwyas be ordered from 1 to resnum
  std:: vector<unsigned> Presind;
  // Build Presind
  for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
    unsigned rind=mypdb.getResidueNumber(mypdb.getAtomNumbers()[i]);
    bool oldres=false;
    for (unsigned j=0; j<Presind.size(); j++) {
      if(rind==Presind[j]) {
        oldres=true;
      }
    }
    if(!oldres) {
      Presind.push_back(rind);
    }
  }
  resnum=Presind.size();

  // Pind0 is the atom/COM used in Nlists (for COM Pind0 is the first atom in the pdb belonging to that COM)
  unsigned Pind0size;
  if(com) {
    Pind0size=resnum;
  } else {
    Pind0size=NLsize;
  }
  std:: vector<unsigned> Pind0(Pind0size);
  // If COM resize important arrays
  comatm.resize(NLsize);
  if(com) {
    nlcom.resize(NLsize);
    compos.resize(NLsize);
    fmass.resize(NLsize,0.);
  }
  log << "Total COM/Atoms: " << Natm*resnum << " \n";
  // Build lists of Atoms/COMs for NLists
  //   comatm filled also for non_COM calculation for analysis purposes
  for (unsigned j=0; j<Natm; j++) {
    unsigned oind;
    for (unsigned i=0; i<Pind0.size(); i++) {
      Pind0[i]=0;
    }
    for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
      // Residue/Atom AtomNumber: used to build NL for COMS/Atoms pairs.
      AtomNumber anum=mypdb.getAtomNumbers()[i];
      // ResidueName/Atomname associated to atom
      string rname=mypdb.getResidueName(anum);
      string aname=mypdb.getAtomName(anum);
      // Index associated to residue/atom: used to separate COM-lists
      unsigned rind=mypdb.getResidueNumber(anum);
      unsigned aind=anum.index();
      // This builds lists for NL
      string Pname;
      unsigned Pind;
      if(com) {
        Pname=rname;
        for(unsigned l=0; l<resnum; l++) {
          if(rind==Presind[l]) {
            Pind=l;
          }
        }
      } else {
        Pname=aname;
        Pind=aind;
      }
      if(Pname==atype[j]) {
        if(Pind0[Pind]==0) {
          // adding the atomnumber to the atom/COM list for pairs
          Plist[j].push_back(anum);
          Pind0[Pind]=aind+1;
          oind=Pind;
        }
        // adding the atomnumber to list of atoms for every COM/Atoms
        comatm[Pind0[Pind]-1].push_back(anum);
      }
    }
    // Output Lists
    log << "  Groups of type  " << j << ": " << Plist[j].size() << " \n";
    string gname;
    unsigned gsize;
    if(com) {
      gname=mypdb.getResidueName(comatm[Pind0[oind]-1][0]);
      gsize=comatm[Pind0[oind]-1].size();
    } else {
      gname=mypdb.getAtomName(comatm[Pind0[oind]-1][0]);
      gsize=1;
    }
    log.printf("    %6s %3s %13s %10i %6s\n", "type  ", gname.c_str(),"   containing ",gsize," atoms");
  }

  // This is to build the list with all the atoms
  std:: vector<AtomNumber> listall;
  for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
    listall.push_back(mypdb.getAtomNumbers()[i]);
  }

  // PIV blocks and Neighbour Lists
  Nlist=0;
  // Direct adds the A-A ad B-B blocks (N)
  if(direct) {
    Nlist=Nlist+unsigned(Natm);
  }
  // Cross adds the A-B blocks (N*(N-1)/2)
  if(cross) {
    Nlist=Nlist+unsigned(double(Natm*(Natm-1))/2.);
  }
  // Resize vectors according to Nlist
  rPIV.resize(Nlist);

  // PIV scaled option
  scaling.resize(Nlist);
  for(unsigned j=0; j<Nlist; j++) {
    scaling[j]=1.;
  }
  if(keywords.exists("SFACTOR")) {
    parseVector("SFACTOR",scaling);
    //if(scaling.size()!=getNumberOfArguments() && scaling.size()!=0) error("not enough values for SFACTOR");
  }
  // Neighbour Lists option
  parseFlag("NLIST",doneigh);
  nl.resize(Nlist);
  nl_skin.resize(Nlist);
  if(doneigh) {
    std:: vector<double> nl_cut(Nlist,0.);
    std:: vector<int> nl_st(Nlist,0);
    parseVector("NL_CUTOFF",nl_cut);
    //if(nl_cut.size()!=getNumberOfArguments() && nl_cut.size()!=0) error("not enough values for NL_CUTOFF");
    parseVector("NL_STRIDE",nl_st);
    //if(nl_st.size()!=getNumberOfArguments() && nl_st.size()!=0) error("not enough values for NL_STRIDE");
    parseVector("NL_SKIN",nl_skin);
    //if(nl_skin.size()!=getNumberOfArguments() && nl_skin.size()!=0) error("not enough values for NL_SKIN");
    for (unsigned j=0; j<Nlist; j++) {
      if(nl_cut[j]<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
      if(nl_st[j]<=0) error("NL_STRIDE should be explicitly specified and positive");
      if(nl_skin[j]<=0.) error("NL_SKIN should be explicitly specified and positive");
      nl_cut[j]=nl_cut[j]+nl_skin[j];
    }
    log << "Creating Neighbor Lists \n";
    // WARNING: is nl_cut meaningful here?
    nlall= new NeighborList(listall,pbc,getPbc(),nl_cut[0],nl_st[0]);
    if(com) {
      //Build lists of Atoms for every COM
      for (unsigned i=0; i<compos.size(); i++) {
        // WARNING: is nl_cut meaningful here?
        nlcom[i]= new NeighborList(comatm[i],pbc,getPbc(),nl_cut[0],nl_st[0]);
      }
    }
    unsigned ncnt=0;
    // Direct blocks AA, BB, CC, ...
    if(direct) {
      for (unsigned j=0; j<Natm; j++) {
        nl[ncnt]= new NeighborList(Plist[j],pbc,getPbc(),nl_cut[j],nl_st[j]);
        ncnt+=1;
      }
    }
    // Cross blocks AB, AC, BC, ...
    if(cross) {
      for (unsigned j=0; j<Natm; j++) {
        for (unsigned i=j+1; i<Natm; i++) {
          nl[ncnt]= new NeighborList(Plist[i],Plist[j],false,pbc,getPbc(),nl_cut[ncnt],nl_st[ncnt]);
          ncnt+=1;
        }
      }
    }
  } else {
    log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
    nlall= new NeighborList(listall,pbc,getPbc());
    for (unsigned j=0; j<Nlist; j++) {
      nl[j]= new NeighborList(Plist[j],Plist[j],true,pbc,getPbc());
    }
  }
  // Output Nlist
  log << "Total Nlists: " << Nlist << " \n";
  for (unsigned j=0; j<Nlist; j++) {
    log << "  list " << j+1 << "   size " << nl[j]->size() << " \n";
  }
  // Calculate COM masses once and for all from lists
  if(com) {
    for(unsigned j=0; j<compos.size(); j++) {
      double commass=0.;
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        commass+=mypdb.getOccupancy()[andx];
      }
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        if(commass>0.) {
          fmass[andx]=mypdb.getOccupancy()[andx]/commass;
        } else {
          fmass[andx]=1.;
        }
      }
    }
  }

  // Sorting
  dosort.resize(Nlist);
  std:: vector<int> ynsort(Nlist);
  parseVector("SORT",ynsort);
  for (unsigned i=0; i<Nlist; i++) {
    if(ynsort[i]==0||CompDer) {
      dosort[i]=false;
    } else {
      dosort[i]=true;
    }
  }

  //build box vectors and correct for pbc
  log << "Building the box from PDB data ... \n";
  Tensor Box=mypdb.getBoxVec();
  log << "  Done! A,B,C vectors in Cartesian space:  \n";
  log.printf("  A:  %12.6f%12.6f%12.6f\n", Box[0][0],Box[0][1],Box[0][2]);
  log.printf("  B:  %12.6f%12.6f%12.6f\n", Box[1][0],Box[1][1],Box[1][2]);
  log.printf("  C:  %12.6f%12.6f%12.6f\n", Box[2][0],Box[2][1],Box[2][2]);
  log << "Changing the PBC according to the new box \n";
  Pbc mypbc;
  mypbc.setBox(Box);
  log << "The box volume is " << mypbc.getBox().determinant() << " \n";

  //Compute scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/mypbc.getBox().determinant());
    log << "Scaling atom distances by  " << Fvol << " \n";
  } else {
    log << "Using unscaled atom distances \n";
  }

  r00.resize(Nlist);
  sw.resize(Nlist);
  for (unsigned j=0; j<Nlist; j++) {
    if( !parseNumbered( "SWITCH", j+1, sw[j] ) ) break;
  }
  if(CompDer) {
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
  }

  // build COMs from positions if requested
  if(com) {
    for(unsigned j=0; j<compos.size(); j++) {
      compos[j][0]=0.;
      compos[j][1]=0.;
      compos[j][2]=0.;
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        compos[j]+=fmass[andx]*mypdb.getPositions()[andx];
      }
    }
  }
  // build the rPIV distances (transformation and sorting is done afterwards)
  if(CompDer) {
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
  }
  for(unsigned j=0; j<Nlist; j++) {
    for(unsigned i=0; i<nl[j]->size(); i++) {
      unsigned i0=(nl[j]->getClosePairAtomNumber(i).first).index();
      unsigned i1=(nl[j]->getClosePairAtomNumber(i).second).index();
      //calculate/get COM position of centers i0 and i1
      Vector Pos0,Pos1;
      if(com) {
        //if(pbc) makeWhole();
        Pos0=compos[i0];
        Pos1=compos[i1];
      } else {
        Pos0=mypdb.getPositions()[i0];
        Pos1=mypdb.getPositions()[i1];
      }
      Vector ddist;
      if(pbc) {
        ddist=mypbc.distance(Pos0,Pos1);
      } else {
        ddist=delta(Pos0,Pos1);
      }
      double df=0.;
      // Transformation and sorting done at the first timestep to solve the r0 definition issue
      if(CompDer) {
        rPIV[j].push_back(sfs[j].calculate(ddist.modulo()*Fvol, df));
      } else {
        rPIV[j].push_back(ddist.modulo()*Fvol);
      }
    }
    if(CompDer) {
      if(dosort[j]) {
        std::sort(rPIV[j].begin(),rPIV[j].end());
      }
      int lmt0=0;
      int lmt1=0;
      for(unsigned i=0; i<rPIV[j].size(); i++) {
        if(int(rPIV[j][i]*double(Nprec-1))==0) {
          lmt0+=1;
        }
        if(int(rPIV[j][i]*double(Nprec-1))==1) {
          lmt1+=1;
        }
      }
      log.printf("       |%10i|%15zu|%15i|%15i|\n", j, rPIV[j].size(), lmt0, lmt1);
    }
  }

  checkRead();
  // From the plumed manual on how to build-up a new Colvar
  addValueWithDerivatives();
  requestAtoms(nlall->getFullAtomList());
  setNotPeriodic();
  // getValue()->setPeridodicity(false);
  // set size of derivative vector
  m_deriv.resize(getNumberOfAtoms());
}

// The following deallocates pointers
PIV::~PIV()
{
  for (unsigned j=0; j<Nlist; j++) {
    delete nl[j];
  }
  if(com) {
    for (unsigned j=0; j<NLsize; j++) {
      delete nlcom[j];
    }
  }
  delete nlall;
}

void PIV::calculate()
{

  // Local varaibles
  // The following are probably needed as static arrays
  static int prev_stp=-1;
  static int init_stp=1;
  static std:: vector<std:: vector<Vector> > prev_pos(Nlist);
  static std:: vector<std:: vector<double> > cPIV(Nlist);
  static std:: vector<std:: vector<int> > Atom0(Nlist);
  static std:: vector<std:: vector<int> > Atom1(Nlist);
  std:: vector<std:: vector<int> > A0(Nprec);
  std:: vector<std:: vector<int> > A1(Nprec);
  unsigned stride=1;
  unsigned rank=0;

  if(!serial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  } else {
    stride=1;
    rank=0;
  }

  // Transform (and sort) the rPIV before starting the dynamics
  if (((prev_stp==-1) || (init_stp==1)) &&!CompDer) {
    if(prev_stp!=-1) {init_stp=0;}
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    //Set switching function parameters
    log << "\n";
    log << "REFERENCE PDB # " << prev_stp+2 << " \n";
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
    //Transform and sort
    log << "Building Reference PIV Vector \n";
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    double df=0.;
    for (unsigned j=0; j<Nlist; j++) {
      for (unsigned i=0; i<rPIV[j].size(); i++) {
        rPIV[j][i]=sfs[j].calculate(rPIV[j][i], df);
      }
      if(dosort[j]) {
        std::sort(rPIV[j].begin(),rPIV[j].end());
      }
      int lmt0=0;
      int lmt1=0;
      for(unsigned i=0; i<rPIV[j].size(); i++) {
        if(int(rPIV[j][i]*double(Nprec-1))==0) {
          lmt0+=1;
        }
        if(int(rPIV[j][i]*double(Nprec-1))==1) {
          lmt1+=1;
        }
      }
      log.printf("       |%10i|%15zu|%15i|%15i|\n", j, rPIV[j].size(), lmt0, lmt1);
    }
    log << "\n";
  }
  // Do the sorting only once per timestep to avoid building the PIV N times for N rPIV PDB structures!
  if ((getStep()>prev_stp&&getStep()%updatePIV==0)||CompDer) {
    if (CompDer) log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV \n";
    //
    // build COMs from positions if requested
    if(com) {
      if(pbc) makeWhole();
      for(unsigned j=0; j<compos.size(); j++) {
        compos[j][0]=0.;
        compos[j][1]=0.;
        compos[j][2]=0.;
        for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
          unsigned andx=nlcom[j]->getFullAtomList()[i].index();
          compos[j]+=fmass[andx]*getPosition(andx);
        }
      }
    }
    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (doneigh) {
      bool doupdate=false;
      // For the first step build previous positions = actual positions
      if (prev_stp==-1) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            prev_pos[j].push_back(Pos);
          }
        }
        doupdate=true;
      }
      // Decide whether to update lists based on atom displacement, every stride
      std:: vector<std:: vector<Vector> > tmp_pos(Nlist);
      if (getStep() % nlall->getStride() ==0) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            tmp_pos[j].push_back(Pos);
            if (pbcDistance(tmp_pos[j][i],prev_pos[j][i]).modulo()>=nl_skin[j]) {
              doupdate=true;
            }
          }
        }
      }
      // Update Nlists if needed
      if (doupdate==true) {
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            prev_pos[j][i]=tmp_pos[j][i];
          }
          nl[j]->update(prev_pos[j]);
          log << " Step " << getStep() << "  Neighbour lists updated " << nl[j]->size() << " \n";
        }
      }
    }
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    Vector ddist;
    // Global to local variables
    bool doserial=serial;
    // Build "Nlist" PIV blocks
    for(unsigned j=0; j<Nlist; j++) {
      if(dosort[j]) {
        // from global to local variables to speedup the for loop with if statements
        bool docom=com;
        bool dopbc=pbc;
        // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
        std:: vector<int> OrdVec(Nprec,0);
        cPIV[j].resize(0);
        Atom0[j].resize(0);
        Atom1[j].resize(0);
        // Building distances for the PIV vector at time t
        if(timer) stopwatch.start("1 Build cPIV");
        for(unsigned i=rank; i<nl[j]->size(); i+=stride) {
          unsigned i0=(nl[j]->getClosePairAtomNumber(i).first).index();
          unsigned i1=(nl[j]->getClosePairAtomNumber(i).second).index();
          Vector Pos0,Pos1;
          if(docom) {
            Pos0=compos[i0];
            Pos1=compos[i1];
          } else {
            Pos0=getPosition(i0);
            Pos1=getPosition(i1);
          }
          if(dopbc) {
            ddist=pbcDistance(Pos0,Pos1);
          } else {
            ddist=delta(Pos0,Pos1);
          }
          double df=0.;
          //Integer sorting ... faster!
          //Transforming distances with the Switching function + real to integer transformation
          int Vint=int(sfs[j].calculate(ddist.modulo()*Fvol, df)*double(Nprec-1)+0.5);
          //Integer transformed distance values as index of the Ordering Vector OrdVec
          OrdVec[Vint]+=1;
          //Keeps track of atom indices for force and virial calculations
          A0[Vint].push_back(i0);
          A1[Vint].push_back(i1);
        }
        if(timer) stopwatch.stop("1 Build cPIV");
        if(timer) stopwatch.start("2 Sort cPIV");
        if(!doserial && comm.initialized()) {
          // Vectors keeping track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          std:: vector<int> Vdim(stride,0);
          std:: vector<int> Vpos(stride,0);
          // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
          std:: vector<int> OrdVecAll(stride*Nprec);
          // Big vectors contining all Atom indexes for every occupancy (Atom0O(Nprec,n) and Atom1O(Nprec,n) matrices in one vector)
          std:: vector<int> Atom0F;
          std:: vector<int> Atom1F;
          // Vector used to reconstruct arrays
          std:: vector<unsigned> k(stride,0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for(unsigned i=1; i<Nprec; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            // Can this be avoided ???
            Atom0F.insert(Atom0F.end(),A0[i].begin(),A0[i].end());
            Atom1F.insert(Atom1F.end(),A1[i].begin(),A1[i].end());
            A0[i].resize(0);
            A1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV block
          A0[0].resize(0);
          A1[0].resize(0);
          A0[Nprec-1].resize(0);
          A1[Nprec-1].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          OrdVec[0]=0;
          OrdVec[Nprec-1]=0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim=Atom0F.size();
          // Vdim and Vpos keep track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim,1,&Vdim[0],1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int Fdim=0;
          for(unsigned i=1; i<stride; i++) {
            Vpos[i]=Vpos[i-1]+Vdim[i-1];
            Fdim+=Vdim[i];
          }
          Fdim+=Vdim[0];
          // build big vectors for atom pairs on all ranks for all ranks
          std:: vector<int> Atom0FAll(Fdim);
          std:: vector<int> Atom1FAll(Fdim);
          // TO BE IMPROVED: Allgathers may be substituded by gathers by proc 0
          //   Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV
          comm.Allgather(&OrdVec[0],Nprec,&OrdVecAll[0],Nprec);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&Atom0F[0],Atom0F.size(),&Atom0FAll[0],&Vdim[0],&Vpos[0]);
          comm.Allgatherv(&Atom1F[0],Atom1F.size(),&Atom1FAll[0],&Vdim[0],&Vpos[0]);

          // Reconstruct the full vectors from collections of Allgathered parts (this is a serial step)
          // This is the tricky serial step, to assemble toghether PIV and atom-pair info from head-tail big vectors
          // Loop before on l and then on i would be better but the allgather should be modified
          // Loop on blocks
          //for(unsigned m=0;m<Nlist;m++) {
          // Loop on Ordering Vector size excluding zeros (i=1)
          if(timer) stopwatch.stop("2 Sort cPIV");
          if(timer) stopwatch.start("3 Reconstruct cPIV");
          for(unsigned i=1; i<Nprec; i++) {
            // Loop on the ranks
            for(unsigned l=0; l<stride; l++) {
              // Loop on the number of head-to-tail pieces
              for(unsigned m=0; m<OrdVecAll[i+l*Nprec]; m++) {
                // cPIV is the current PIV at time t
                cPIV[j].push_back(double(i)/double(Nprec-1));
                Atom0[j].push_back(Atom0FAll[k[l]+Vpos[l]]);
                Atom1[j].push_back(Atom1FAll[k[l]+Vpos[l]]);
                k[l]+=1;
              }
            }
          }
          if(timer) stopwatch.stop("3 Reconstruct cPIV");
        } else {
          for(unsigned i=1; i<Nprec; i++) {
            for(unsigned m=0; m<OrdVec[i]; m++) {
              cPIV[j].push_back(double(i)/double(Nprec-1));
              Atom0[j].push_back(A0[i][m]);
              Atom1[j].push_back(A1[i][m]);
            }
          }
        }
      }
    }
  }
  Vector distance;
  double dfunc=0.;
  // Calculate volume scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/getBox().determinant());
  }

  // This test may be run by specifying the TEST keyword as input, it pritnts rPIV and cPIV and quits
  if(test) {
    unsigned limit=0;
    for(unsigned j=0; j<Nlist; j++) {
      if(dosort[j]) {
        limit = cPIV[j].size();
      } else {
        limit = rPIV[j].size();
      }
      log.printf("PIV Block:  %6i %12s %6i \n", j, "      Size:", limit);
      log.printf("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV   ","    r-PIV   ","   i-j distance vector       ");
      for(unsigned i=0; i<limit; i++) {
        unsigned i0=0;
        unsigned i1=0;
        if(dosort[j]) {
          i0=Atom0[j][i];
          i1=Atom1[j][i];
        } else {
          i0=(nl[j]->getClosePairAtomNumber(i).first).index();
          i1=(nl[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector Pos0,Pos1;
        if(com) {
          Pos0=compos[i0];
          Pos1=compos[i1];
        } else {
          Pos0=getPosition(i0);
          Pos1=getPosition(i1);
        }
        if(pbc) {
          distance=pbcDistance(Pos0,Pos1);
        } else {
          distance=delta(Pos0,Pos1);
        }
        dfunc=0.;
        double cP,rP;
        if(dosort[j]) {
          cP = cPIV[j][i];
          rP = rPIV[j][rPIV[j].size()-cPIV[j].size()+i];
        } else {
          double dm=distance.modulo();
          cP = sfs[j].calculate(dm*Fvol, dfunc);
          rP = rPIV[j][i];
        }
        log.printf("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
      }
    }
    log.printf("This was a test, now exit \n");
    exit();
  }

  if(timer) stopwatch.start("4 Build For Derivatives");
  // non-global variables Nder and Scalevol defined to speedup if structures in cycles
  bool Nder=CompDer;
  bool Scalevol=Svol;
  if(getStep()%updatePIV==0) {
    // set to zero PIVdistance, derivatives and virial when they are calculated
    for(unsigned j=0; j<m_deriv.size(); j++) {
      for(unsigned k=0; k<3; k++) {m_deriv[j][k]=0.;}
    }
    for(unsigned j=0; j<3; j++) {
      for(unsigned k=0; k<3; k++) {
        m_virial[j][k]=0.;
      }
    }
    m_PIVdistance=0.;
    // Re-compute atomic distances for derivatives and compute PIV-PIV distance
    for(unsigned j=0; j<Nlist; j++) {
      unsigned limit=0;
      // dosorting definition is to speedup if structure in cycles with non-global variables
      bool dosorting=dosort[j];
      bool docom=com;
      bool dopbc=pbc;
      if(dosorting) {
        limit = cPIV[j].size();
      } else {
        limit = rPIV[j].size();
      }
      for(unsigned i=rank; i<limit; i+=stride) {
        unsigned i0=0;
        unsigned i1=0;
        if(dosorting) {
          i0=Atom0[j][i];
          i1=Atom1[j][i];
        } else {
          i0=(nl[j]->getClosePairAtomNumber(i).first).index();
          i1=(nl[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector Pos0,Pos1;
        if(docom) {
          Pos0=compos[i0];
          Pos1=compos[i1];
        } else {
          Pos0=getPosition(i0);
          Pos1=getPosition(i1);
        }
        if(dopbc) {
          distance=pbcDistance(Pos0,Pos1);
        } else {
          distance=delta(Pos0,Pos1);
        }
        dfunc=0.;
        // this is needed for dfunc and dervatives
        double dm=distance.modulo();
        double tPIV = sfs[j].calculate(dm*Fvol, dfunc);
        // PIV distance
        double coord=0.;
        if(!dosorting||Nder) {
          coord = tPIV - rPIV[j][i];
        } else {
          coord = cPIV[j][i] - rPIV[j][rPIV[j].size()-cPIV[j].size()+i];
        }
        // Calculate derivatives, virial, and variable=sum_j (scaling[j] *(cPIV-rPIV)_j^2)
        // WARNING: dfunc=dswf/(Fvol*dm)  (this may change in future Plumed versions)
        double tmp = 2.*scaling[j]*coord*Fvol*Fvol*dfunc;
        Vector tmpder = tmp*distance;
        // 0.5*(x_i-x_k)*f_ik         (force on atom k due to atom i)
        if(docom) {
          Vector dist;
          for(unsigned k=0; k<nlcom[i0]->getFullAtomList().size(); k++) {
            unsigned x0=nlcom[i0]->getFullAtomList()[k].index();
            m_deriv[x0] -= tmpder*fmass[x0];
            for(unsigned l=0; l<3; l++) {
              dist[l]=0.;
            }
            Vector P0=getPosition(x0);
            for(unsigned l=0; l<nlcom[i0]->getFullAtomList().size(); l++) {
              unsigned x1=nlcom[i0]->getFullAtomList()[l].index();
              Vector P1=getPosition(x1);
              if(dopbc) {
                dist+=pbcDistance(P0,P1);
              } else {
                dist+=delta(P0,P1);
              }
            }
            for(unsigned l=0; l<nlcom[i1]->getFullAtomList().size(); l++) {
              unsigned x1=nlcom[i1]->getFullAtomList()[l].index();
              Vector P1=getPosition(x1);
              if(dopbc) {
                dist+=pbcDistance(P0,P1);
              } else {
                dist+=delta(P0,P1);
              }
            }
            m_virial    -= 0.25*fmass[x0]*Tensor(dist,tmpder);
          }
          for(unsigned k=0; k<nlcom[i1]->getFullAtomList().size(); k++) {
            unsigned x1=nlcom[i1]->getFullAtomList()[k].index();
            m_deriv[x1] += tmpder*fmass[x1];
            for(unsigned l=0; l<3; l++) {
              dist[l]=0.;
            }
            Vector P1=getPosition(x1);
            for(unsigned l=0; l<nlcom[i1]->getFullAtomList().size(); l++) {
              unsigned x0=nlcom[i1]->getFullAtomList()[l].index();
              Vector P0=getPosition(x0);
              if(dopbc) {
                dist+=pbcDistance(P1,P0);
              } else {
                dist+=delta(P1,P0);
              }
            }
            for(unsigned l=0; l<nlcom[i0]->getFullAtomList().size(); l++) {
              unsigned x0=nlcom[i0]->getFullAtomList()[l].index();
              Vector P0=getPosition(x0);
              if(dopbc) {
                dist+=pbcDistance(P1,P0);
              } else {
                dist+=delta(P1,P0);
              }
            }
            m_virial    += 0.25*fmass[x1]*Tensor(dist,tmpder);
          }
        } else {
          m_deriv[i0] -= tmpder;
          m_deriv[i1] += tmpder;
          m_virial    -= tmp*Tensor(distance,distance);
        }
        if(Scalevol) {
          m_virial+=1./3.*tmp*dm*dm*Tensor::identity();
        }
        m_PIVdistance    += scaling[j]*coord*coord;
      }
    }

    if (!serial && comm.initialized()) {
      comm.Barrier();
      comm.Sum(&m_PIVdistance,1);
      if(!m_deriv.empty()) comm.Sum(&m_deriv[0][0],3*m_deriv.size());
      comm.Sum(&m_virial[0][0],9);
    }
  }
  prev_stp=getStep();

  //Timing
  if(timer) stopwatch.stop("4 Build For Derivatives");
  if(timer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }

  // Update derivatives, virial, and variable (PIV-distance^2)
  for(unsigned i=0; i<m_deriv.size(); ++i) setAtomsDerivatives(i,m_deriv[i]);
  setValue           (m_PIVdistance);
  setBoxDerivatives  (m_virial);
}
//Close Namespaces at the very beginning
}
}

