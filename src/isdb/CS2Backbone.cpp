/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

#define cutOffNB      0.60	// buffer distance for neighbour-lists 
#define cutOffDist    0.50  	// cut off distance for non-bonded pairwise forces
#define cutOnDist     0.32   	// cut off distance for non-bonded pairwise forces
#define cutOffNB2     cutOffNB*cutOffNB // squared buffer distance for neighbour-lists 
#define cutOffDist2   cutOffDist*cutOffDist
#define cutOnDist2    cutOnDist*cutOnDist
#define invswitch     1.0/((cutOffDist2-cutOnDist2)*(cutOffDist2-cutOnDist2)*(cutOffDist2-cutOnDist2))
#define cutOffDist4   cutOffDist2*cutOffDist2
#define cutMixed      cutOffDist2*cutOffDist2*cutOffDist2 -3.*cutOffDist2*cutOffDist2*cutOnDist2

#include <string>
#include <fstream>
#include <iterator>
#include <sstream>

#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "tools/Torsion.h"

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR CS2BACKBONE
/*
Calculates the backbone chemical shifts for a protein.

The functional form is that of CamShift \cite Kohlhoff:2009us. The chemical shift
of the selected nuclei can be saved as components. Alternatively one can calculate either
the CAMSHIFT score (useful as a collective variable \cite Granata:2013dk or as a scoring
function \cite Robustelli:2010dn) or a \ref METAINFERENCE score (using DOSCORE).
For these two latter cases experimental chemical shifts must be provided.

CS2BACKBONE calculation can be relatively heavy because it often uses a large number of atoms, it can
be run in parallel using MPI and \ref Openmp.

As a general rule, when using \ref CS2BACKBONE or other experimental restraints it may be better to
increase the accuracy of the constraint algorithm due to the increased strain on the bonded structure.
In the case of GROMACS it is safer to use lincs-iter=2 and lincs-order=6.

In general the system for which chemical shifts are calculated must be completely included in
ATOMS and a TEMPLATE pdb file for the same atoms should be provided as well in the folder DATADIR.
The system is made automatically whole unless NOPBC is used, in particular if the system is made
by multiple chains it is usually better to use NOPBC and make the molecule whole \ref WHOLEMOLECULES
selecting an appropriate order of the atoms. The pdb file is needed to the generate a simple topology of the protein.
For histidine residues in protonation states different from D the HIE/HSE HIP/HSP name should be used.
GLH and ASH can be used for the alternative protonation of GLU and ASP. Non-standard amino acids and other
molecules are not yet supported, but in principle they can be named UNK. If multiple chains are present
the chain identifier must be in the standard PDB format, together with the TER keyword at the end of each chain.
Termini groups like ACE or NME should be removed from the TEMPLATE pdb because they are not recognized by
CS2BACKBONE.

In addition to a pdb file one needs to provide a list of chemical shifts to be calculated using one
file per nucleus type (CAshifts.dat, CBshifts.dat, Cshifts.dat, Hshifts.dat, HAshifts.dat, Nshifts.dat),
add only the files for the nuclei you need, but each file should include all protein residues.
A chemical shift for a nucleus is calculated if a value greater than 0 is provided.
For practical purposes the value can correspond to the experimental value.
Residues numbers should match that used in the pdb file, but must be positive, so double check the pdb.
The first and last residue of each chain should be preceded by a # character.

\verbatim
CAshifts.dat:
#1 0.0
2 55.5
3 58.4
.
.
#last 0.0
#first of second chain
.
#last of second chain
\endverbatim

The default behavior is to store the values for the active nuclei in components (ca-#, cb-#,
co-#, ha-#, hn-#, nh-# and expca-#, expcb-#, expco-#, expha-#, exphn-#, exp-nh#) with NOEXP it is possible
to only store the back-calculated values, where # includes a chain and residue number.

One additional file is always needed in the folder DATADIR: camshift.db. This file includes all the parameters needed to
calculate the chemical shifts and can be found in regtest/isdb/rt-cs2backbone/data/ .

Additional material and examples can be also found in the tutorial \ref isdb-1 as well as in the cs2backbone regtests
in the isdb folder.

\par Examples

In this first example the chemical shifts are used to calculate a collective variable to be used
in NMR driven Metadynamics \cite Granata:2013dk :

\plumedfile
whole: GROUP ATOMS=2612-2514:-1,961-1:-1,2466-962:-1,2513-2467:-1
WHOLEMOLECULES ENTITY0=whole
cs: CS2BACKBONE ATOMS=1-2612 DATADIR=data/ TEMPLATE=template.pdb CAMSHIFT NOPBC
metad: METAD ARG=cs HEIGHT=0.5 SIGMA=0.1 PACE=200 BIASFACTOR=10
PRINT ARG=cs,metad.bias FILE=COLVAR STRIDE=100
\endplumedfile

In this second example the chemical shifts are used as replica-averaged restrained as in \cite Camilloni:2012je \cite Camilloni:2013hs.

\plumedfile
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/
encs: ENSEMBLE ARG=(cs\.hn-.*),(cs\.nh-.*)
stcs: STATS ARG=encs.* SQDEVSUM PARARG=(cs\.exphn-.*),(cs\.expnh-.*)
RESTRAINT ARG=stcs.sqdevsum AT=0 KAPPA=0 SLOPE=24

PRINT ARG=(cs\.hn-.*),(cs\.nh-.*) FILE=RESTRAINT STRIDE=100

\endplumedfile

This third example show how to use chemical shifts to calculate a \ref METAINFERENCE score .

\plumedfile
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/ SIGMA_MEAN0=1.0 DOSCORE
csbias: BIASVALUE ARG=cs.score

PRINT ARG=(cs\.hn-.*),(cs\.nh-.*) FILE=CS.dat STRIDE=1000
PRINT ARG=cs.score FILE=BIAS STRIDE=100
\endplumedfile

*/
//+ENDPLUMEDOC

class CS2BackboneDB {
  enum { STD, GLY, PRO};
  enum { HA_ATOM, H_ATOM, N_ATOM, CA_ATOM, CB_ATOM, C_ATOM };
  static const unsigned aa_kind = 3;
  static const unsigned atm_kind = 6;
  static const unsigned numXtraDists = 27;

  // ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
  double c_aa[aa_kind][atm_kind][20];
  double c_aa_prev[aa_kind][atm_kind][20];
  double c_aa_succ[aa_kind][atm_kind][20];
  double co_bb[aa_kind][atm_kind][16];
  double co_sc_[aa_kind][atm_kind][20][20];
  double co_xd[aa_kind][atm_kind][numXtraDists];
  double co_sphere[aa_kind][atm_kind][2][8];
  // for ring current effects
  // Phe, Tyr, Trp_1, Trp_2, His
  double co_ring[aa_kind][atm_kind][5];
  // for dihedral angles
  // co * (a * cos(3 * omega + c) + b * cos(omega + d))
  double co_da[aa_kind][atm_kind][3];
  double pars_da[aa_kind][atm_kind][3][5];

public:

  inline unsigned kind(const string &s) {
    if(s=="GLY") return GLY;
    else if(s=="PRO") return PRO;
    return STD;
  }

  inline unsigned atom_kind(const string &s) {
    if(s=="HA")return HA_ATOM;
    else if(s=="H") return H_ATOM;
    else if(s=="N") return N_ATOM;
    else if(s=="CA")return CA_ATOM;
    else if(s=="CB")return CB_ATOM;
    else if(s=="C") return C_ATOM;
    return -1;
  }

  unsigned get_numXtraDists() {return numXtraDists;}

  //PARAMETERS
  inline double * CONSTAACURR(const unsigned a_kind, const unsigned at_kind) {return c_aa[a_kind][at_kind];}
  inline double * CONSTAANEXT(const unsigned a_kind, const unsigned at_kind) {return c_aa_succ[a_kind][at_kind];}
  inline double * CONSTAAPREV(const unsigned a_kind, const unsigned at_kind) {return c_aa_prev[a_kind][at_kind];}
  inline double * CONST_BB2(const unsigned a_kind, const unsigned at_kind) {return co_bb[a_kind][at_kind];}
  inline double * CONST_SC2(const unsigned a_kind, const unsigned at_kind, unsigned res_type) { return co_sc_[a_kind][at_kind][res_type];}
  inline double * CONST_XD(const unsigned a_kind, const unsigned at_kind) { return co_xd[a_kind][at_kind];}
  inline double * CO_SPHERE(const unsigned a_kind, const unsigned at_kind, unsigned exp_type) { return co_sphere[a_kind][at_kind][exp_type];}
  inline double * CO_RING(const unsigned a_kind, const unsigned at_kind) { return co_ring[a_kind][at_kind];}
  inline double * CO_DA(const unsigned a_kind, const unsigned at_kind) { return co_da[a_kind][at_kind];}
  inline double * PARS_DA(const unsigned a_kind, const unsigned at_kind, const unsigned ang_kind) { return pars_da[a_kind][at_kind][ang_kind];}

  void parse(const string &file, const double dscale) {
    ifstream in;
    in.open(file.c_str());
    if(!in) plumed_merror("Unable to open DB file: " + file);

    unsigned c_kind = 0;
    unsigned c_atom = 0;
    unsigned nline = 0;

    for(unsigned i=0; i<3; i++) for(unsigned j=0; j<6; j++) {
        for(unsigned k=0; k<20; k++) {
          c_aa[i][j][k]=0.;
          c_aa_prev[i][j][k]=0.;
          c_aa_succ[i][j][k]=0.;
          for(unsigned m=0; m<20; m++) co_sc_[i][j][k][m]=0.;
        }
        for(unsigned k=0; k<16; k++) {co_bb[i][j][k]=0.; }
        for(unsigned k=0; k<8; k++) { co_sphere[i][j][0][k]=0.; co_sphere[i][j][1][k]=0.; }
        for(unsigned k=0; k<3; k++) {
          co_da[i][j][k]=0.;
          for(unsigned l=0; l<5; l++) pars_da[i][j][k][l]=0.;
        }
        for(unsigned k=0; k<5; k++) co_ring[i][j][k]=0.;
        for(unsigned k=0; k<numXtraDists; k++) co_xd[i][j][k]=0.;
      }

    while(!in.eof()) {
      string line;
      getline(in,line);
      ++nline;
      if(line.compare(0,1,"#")==0) continue;
      vector<string> tok;
      vector<string> tmp;
      tok = split(line,' ');
      for(unsigned q=0; q<tok.size(); q++)
        if(tok[q].size()) tmp.push_back(tok[q]);
      tok = tmp;
      if(tok.size()==0) continue;
      if(tok[0]=="PAR") {
        c_kind = kind(tok[2]);
        c_atom = atom_kind(tok[1]);
        continue;
      }
      else if(tok[0]=="WEIGHT") {
        continue;
      }
      else if(tok[0]=="FLATBTM") {
        continue;
      }
      else if (tok[0] == "SCALEHARM") {
        continue;
      }
      else if (tok[0] == "TANHAMPLI") {
        continue;
      }
      else if (tok[0] == "ENDHARMON") {
        continue;
      }
      else if (tok[0] == "MAXRCDEVI") {
        continue;
      }
      else if (tok[0] == "RANDCOIL") {
        continue;
      }
      else if (tok[0] == "CONST") {
        continue;
      }
      else if (tok[0] == "CONSTAA") {
        assign(c_aa[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "CONSTAA-1") {
        assign(c_aa_prev[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "CONSTAA+1") {
        assign(c_aa_succ[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "COBB1") {
        continue;
      }
      else if (tok[0] == "COBB2") {
        //angstrom to nm
        assign(co_bb[c_kind][c_atom],tok,dscale);
        continue;
      }
      else if (tok[0] == "SPHERE1") {
        // angstrom^-3 to nm^-3
        assign(co_sphere[c_kind][c_atom][0],tok,1./(dscale*dscale*dscale));
        continue;
      }
      else if (tok[0] == "SPHERE2") {
        //angstrom to nm
        assign(co_sphere[c_kind][c_atom][1],tok,dscale);
        continue;
      }
      else if (tok[0] == "DIHEDRALS") {
        assign(co_da[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "RINGS") {
        // angstrom^-3 to nm^-3
        assign(co_ring[c_kind][c_atom],tok,1./(dscale*dscale*dscale));
        for(unsigned i=1; i<tok.size(); i++)
          co_ring[c_kind][c_atom][i-1] *= 1000;
        continue;
      }
      else if (tok[0] == "HBONDS") {
        continue;
      }
      else if (tok[0] == "XTRADISTS") {
        //angstrom to nm
        assign(co_xd[c_kind][c_atom],tok,dscale);
        continue;
      }
      else if(tok[0]=="DIHEDPHI") {
        assign(pars_da[c_kind][c_atom][0],tok,1);
        continue;
      }
      else if(tok[0]=="DIHEDPSI") {
        assign(pars_da[c_kind][c_atom][1],tok,1);
        continue;
      }
      else if(tok[0]=="DIHEDCHI1") {
        assign(pars_da[c_kind][c_atom][2],tok,1);
        continue;
      }

      bool ok = false;
      string scIdent1 [] = {"COSCALA1", "COSCARG1", "COSCASN1", "COSCASP1", "COSCCYS1", "COSCGLN1", "COSCGLU1",
                            "COSCGLY1", "COSCHIS1", "COSCILE1", "COSCLEU1", "COSCLYS1", "COSCMET1", "COSCPHE1",
                            "COSCPRO1", "COSCSER1", "COSCTHR1", "COSCTRP1", "COSCTYR1", "COSCVAL1"
                           };

      for(unsigned scC = 0; scC < 20; scC++) {
        if(tok[0]==scIdent1[scC]) {
          ok = true;
          break;
        }
      }
      if(ok) continue;

      string scIdent2 [] = {"COSCALA2", "COSCARG2", "COSCASN2", "COSCASP2", "COSCCYS2", "COSCGLN2", "COSCGLU2",
                            "COSCGLY2", "COSCHIS2", "COSCILE2", "COSCLEU2", "COSCLYS2", "COSCMET2", "COSCPHE2",
                            "COSCPRO2", "COSCSER2", "COSCTHR2", "COSCTRP2", "COSCTYR2", "COSCVAL2"
                           };

      for(unsigned scC = 0; scC < 20; scC++) {
        if(tok[0]==scIdent2[scC]) {
          //angstrom to nm
          assign(co_sc_[c_kind][c_atom][scC],tok,dscale);
          ok = true; break;
        }
      }
      if(ok) continue;

      if(tok.size()) {
        string str_err = "DB WARNING: unrecognized token: " + tok[0];
        plumed_merror(str_err);
      }
    }
    in.close();
  }

private:

  vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
      elems.push_back(item);
    }
    return elems;
  }

  vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
  }

  void assign(double * f, const vector<string> & v, const double scale) {
    for(unsigned i=1; i<v.size(); i++) {
      f[i-1] = scale*(atof(v[i].c_str()));
      if(fabs(f[i-1])<0.000001) f[i-1]=0.;
    }
  }
};

class CS2Backbone : public MetainferenceBase {
  struct ChemicalShift {
    double exp_cs;              // a reference chemical shifts
    Value *comp;                // a pointer to the component
    unsigned res_kind;          // residue type (STD/GLY/PRO)
    unsigned atm_kind;          // nuclues (HA/CA/CB/CO/NH/HN)
    unsigned res_type_prev;     // previuos residue (ALA/VAL/..)
    unsigned res_type_curr;     // current residue (ALA/VAL/..)
    unsigned res_type_next;     // next residue (ALA/VAL/..)
    string res_name;            // residue name
    string nucleus;             // chemical shift
    bool has_chi1;              // does we have a chi1
    unsigned csatoms;           // fixed number of atoms used
    unsigned totcsatoms;        // number of atoms used
    unsigned res_num;           // residue number
    unsigned chain;             // chain number
    unsigned ipos;              // index of the atom for which we are calculating the chemical shifts
    vector<unsigned> bb;        // atoms for the previous, current and next backbone
    vector<unsigned> side_chain;// atoms for the current sidechain
    vector<int> xd1;            // additional couple of atoms
    vector<int> xd2;            // additional couple of atoms
    vector<unsigned> box_nb;    // non-bonded atoms

    ChemicalShift():
      exp_cs(0.),
      comp(NULL),
      res_kind(0),
      atm_kind(0),
      res_type_prev(0),
      res_type_curr(0),
      res_type_next(0),
      res_name(""),
      nucleus(""),
      has_chi1(true),
      csatoms(0),
      totcsatoms(0),
      res_num(0),
      chain(0),
      ipos(0)
    {
      xd1.reserve(26);
      xd2.reserve(26);
      box_nb.reserve(150);
    }
  };

  struct RingInfo {
    enum {R_PHE, R_TYR, R_TRP1, R_TRP2, R_HIS};
    unsigned rtype;    // one out of five different types
    unsigned atom[6];  // up to six member per ring
    unsigned numAtoms; // number of ring members (5 or 6)
    Vector position;   // center of ring coordinates
    Vector normVect;   // ring plane normal vector
    Vector g[6];       // vector of the vectors used for normVect
    double lengthN2;   // square of length of normVect
    double lengthNV;   // length of normVect
    RingInfo():
      rtype(0),
      numAtoms(0),
      lengthN2(NAN),
      lengthNV(NAN)
    {
      for(unsigned i=0; i<6; i++) atom[i]=0;
    }
  };

  enum aa_t {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK};
  enum sequence_t {Np, CAp, HAp, Cp, Op, Nc, Hc, CAc, HAc, Cc, Oc, Nn, Hn, CAn, HAn, Cn, CBc, CGc};

  CS2BackboneDB    db;
  vector<ChemicalShift> chemicalshifts;

  vector<RingInfo> ringInfo;
  vector<unsigned> type;
  vector<unsigned> res_num;
  unsigned         max_cs_atoms;
  unsigned         box_nupdate;
  unsigned         box_count;
  bool             camshift;
  bool             pbc;
  bool             serial;

  void init_cs(const string &file, const string &k, const PDB &pdb);
  void update_neighb();
  void compute_ring_parameters();
  void init_types(const PDB &pdb);
  void init_rings(const PDB &pdb);
  aa_t frag2enum(const string &aa);
  vector<string> side_chain_atoms(const string &s);
  bool isSP2(const string & resType, const string & atomName);
  bool is_chi1_cx(const string & frg, const string & atm);
  void xdist_name_map(string & name);

public:

  explicit CS2Backbone(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(CS2Backbone,"CS2BACKBONE")

void CS2Backbone::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","DATADIR","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file of the protein system.");
  keys.add("compulsory","NEIGH_FREQ","20","Period in step for neighbor list update.");
  keys.addFlag("CAMSHIFT",false,"Set to TRUE if you to calculate a single CamShift score.");
  keys.addFlag("NOEXP",false,"Set to TRUE if you don't want to have fixed components with the experimental values.");
  keys.addOutputComponent("ha","default","the calculated Ha hydrogen chemical shifts");
  keys.addOutputComponent("hn","default","the calculated H hydrogen chemical shifts");
  keys.addOutputComponent("nh","default","the calculated N nitrogen chemical shifts");
  keys.addOutputComponent("ca","default","the calculated Ca carbon chemical shifts");
  keys.addOutputComponent("cb","default","the calculated Cb carbon chemical shifts");
  keys.addOutputComponent("co","default","the calculated C' carbon chemical shifts");
  keys.addOutputComponent("expha","default","the experimental Ha hydrogen chemical shifts");
  keys.addOutputComponent("exphn","default","the experimental H hydrogen chemical shifts");
  keys.addOutputComponent("expnh","default","the experimental N nitrogen chemical shifts");
  keys.addOutputComponent("expca","default","the experimental Ca carbon chemical shifts");
  keys.addOutputComponent("expcb","default","the experimental Cb carbon chemical shifts");
  keys.addOutputComponent("expco","default","the experimental C' carbon chemical shifts");
}

CS2Backbone::CS2Backbone(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  max_cs_atoms(0),
  camshift(false),
  pbc(true),
  serial(false)
{
  vector<AtomNumber> used_atoms;
  parseAtomList("ATOMS",used_atoms);

  parseFlag("CAMSHIFT",camshift);
  if(camshift&&getDoScore()) plumed_merror("It is not possible to use CAMSHIFT and DOSCORE at the same time");

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parseFlag("SERIAL",serial);

  bool noexp=false;
  parseFlag("NOEXP",noexp);

  string stringa_data;
  parse("DATADIR",stringa_data);

  string stringa_template;
  parse("TEMPLATE",stringa_template);

  box_count=0;
  box_nupdate=20;
  parse("NEIGH_FREQ", box_nupdate);

  string stringadb  = stringa_data + string("/camshift.db");
  string stringapdb = stringa_data + string("/") + stringa_template;

  /* Lenght conversion (parameters are tuned for angstrom) */
  double scale=1.;
  if(!plumed.getAtoms().usingNaturalUnits()) {
    scale = 10.*atoms.getUnits().getLength();
  }

  log.printf("  Initialization of the predictor ...\n");
  db.parse(stringadb,scale);

  PDB pdb;
  if( !pdb.read(stringapdb,plumed.getAtoms().usingNaturalUnits(),1./scale) ) plumed_merror("missing input file " + stringapdb);

  // first of all we build the list of chemical shifts we want to predict
  log.printf("  Reading experimental data ...\n"); log.flush();
  stringadb = stringa_data + string("/CAshifts.dat");
  log.printf("  Initializing CA shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "CA", pdb);
  stringadb = stringa_data + string("/CBshifts.dat");
  log.printf("  Initializing CB shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "CB", pdb);
  stringadb = stringa_data + string("/Cshifts.dat");
  log.printf("  Initializing C' shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "C", pdb);
  stringadb = stringa_data + string("/HAshifts.dat");
  log.printf("  Initializing HA shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "HA", pdb);
  stringadb = stringa_data + string("/Hshifts.dat");
  log.printf("  Initializing H shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "H", pdb);
  stringadb = stringa_data + string("/Nshifts.dat");
  log.printf("  Initializing N shifts %s\n", stringadb.c_str());
  init_cs(stringadb, "N", pdb);

  if(chemicalshifts.size()==0) plumed_merror("There are no chemical shifts to calculate, there must be at least a not empty file (CA|CB|C|HA|H|N|shifts.dat)");

  init_types(pdb);
  init_rings(pdb);

  log<<"  Bibliography "
     <<plumed.cite("Kohlhoff K, Robustelli P, Cavalli A, Salvatella A, Vendruscolo M, J. Am. Chem. Soc. 131, 13894 (2009)");
  if(camshift) log<<plumed.cite("Granata D, Camilloni C, Vendruscolo M, Laio A, Proc. Natl. Acad. Sci. USA 110, 6817 (2013)");
  else log<<plumed.cite("Camilloni C, Robustelli P, De Simone A, Cavalli A, Vendruscolo M, J. Am. Chem. Soc. 134, 3968 (2012)");
  log<<plumed.cite("Bonomi M, Camilloni C, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  if(camshift) {
    noexp = true;
    addValueWithDerivatives();
    setNotPeriodic();
  } else {
    for(unsigned cs=0; cs<chemicalshifts.size(); cs++) {
      std::string num; Tools::convert(chemicalshifts[cs].res_num,num);
      std::string chain_num; Tools::convert(chemicalshifts[cs].chain,chain_num);
      if(getDoScore()) {
        addComponent(chemicalshifts[cs].nucleus+chain_num+"-"+num);
        componentIsNotPeriodic(chemicalshifts[cs].nucleus+chain_num+"-"+num);
        chemicalshifts[cs].comp = getPntrToComponent(chemicalshifts[cs].nucleus+chain_num+"-"+num);
        setParameter(chemicalshifts[cs].exp_cs);
      } else {
        addComponentWithDerivatives(chemicalshifts[cs].nucleus+chain_num+"-"+num);
        componentIsNotPeriodic(chemicalshifts[cs].nucleus+chain_num+"-"+num);
        chemicalshifts[cs].comp = getPntrToComponent(chemicalshifts[cs].nucleus+chain_num+"-"+num);
      }
    }
    if(getDoScore()) Initialise(chemicalshifts.size());
  }

  if(!noexp) {
    for(unsigned cs=0; cs<chemicalshifts.size(); cs++) {
      std::string num; Tools::convert(chemicalshifts[cs].res_num,num);
      std::string chain_num; Tools::convert(chemicalshifts[cs].chain,chain_num);
      addComponent("exp"+chemicalshifts[cs].nucleus+chain_num+"-"+num);
      componentIsNotPeriodic("exp"+chemicalshifts[cs].nucleus+chain_num+"-"+num);
      Value* comp=getPntrToComponent("exp"+chemicalshifts[cs].nucleus+chain_num+"-"+num);
      comp->set(chemicalshifts[cs].exp_cs);
    }
  }

  requestAtoms(used_atoms, false);
  setDerivatives();
  checkRead();
}

void CS2Backbone::init_cs(const string &file, const string &nucl, const PDB &pdb) {
  // number of chains
  vector<string> chains;
  pdb.getChainNames( chains );
  unsigned ichain=0;

  ifstream in;
  in.open(file.c_str());
  if(!in) return;
  istream_iterator<string> iter(in), end;
  unsigned begin=0;

  while(iter!=end) {
    string tok = *iter;
    ++iter;
    if(tok[0]=='#') {
      ++iter;
      if(begin==1) {
        begin=0;
        ichain++;
      } else begin=1;
      continue;
    }
    int ro = atoi(tok.c_str());
    if(ro<0) plumed_merror("Residue numbers should be positive\n");
    unsigned resnum = static_cast<unsigned> (ro);
    tok = *iter;
    ++iter;
    double cs = atof(tok.c_str());
    if(cs==0) continue;

    unsigned fres, lres;
    string errmsg;
    pdb.getResidueRange(chains[ichain], fres, lres, errmsg);
    if(resnum==fres||resnum==lres) plumed_merror("First and Last residue of each chain should be annotated as # in " + file + " Remember that residue numbers should match");

    // check in the PDB for the chain/residue/atom and enable the chemical shift
    string RES = pdb.getResidueName(resnum, chains[ichain]);
    if(RES=="HIE"||RES=="HIP"||RES=="HIS"||RES=="HSP"||RES=="HSE"||RES=="CYS"||RES=="GLH"||RES=="ASH"||RES=="UNK") continue;
    if(RES=="GLN"&&nucl=="CB") continue;
    if(RES=="ILE"&&nucl=="CB") continue;
    if(RES=="PRO"&&nucl=="N") continue;
    if(RES=="PRO"&&nucl=="H") continue;
    if(RES=="PRO"&&nucl=="CB") continue;
    if(RES=="GLY"&&nucl=="HA") continue;
    if(RES=="GLY"&&nucl=="CB") continue;

    ChemicalShift tmp_cs;

    tmp_cs.exp_cs = cs;
    if(nucl=="CA")      tmp_cs.nucleus = "ca-";
    else if(nucl=="CB") tmp_cs.nucleus = "cb-";
    else if(nucl=="C")  tmp_cs.nucleus = "co-";
    else if(nucl=="HA") tmp_cs.nucleus = "ha-";
    else if(nucl=="H")  tmp_cs.nucleus = "hn-";
    else if(nucl=="N")  tmp_cs.nucleus = "nh-";
    tmp_cs.chain = ichain;
    tmp_cs.res_num = resnum;
    tmp_cs.res_type_curr = frag2enum(RES);
    tmp_cs.res_type_prev = frag2enum(pdb.getResidueName(resnum-1, chains[ichain]));
    tmp_cs.res_type_next = frag2enum(pdb.getResidueName(resnum+1, chains[ichain]));
    tmp_cs.res_name = RES;
    tmp_cs.res_kind = db.kind(RES);
    tmp_cs.atm_kind = db.atom_kind(nucl);
    if(RES!="ALA"&&RES!="GLY") {tmp_cs.bb.resize(18); tmp_cs.has_chi1=true;}
    else {tmp_cs.bb.resize(16); tmp_cs.has_chi1=false;}

    vector<AtomNumber> res_atoms = pdb.getAtomsInResidue(resnum, chains[ichain]);
    // find the position of the nucleus and of the other backbone atoms as well as for phi/psi/chi
    for(unsigned a=0; a<res_atoms.size(); a++) {
      string AN = pdb.getAtomName(res_atoms[a]);
      if(nucl=="HA"&&(AN=="HA"||AN=="HA1"||AN=="HA3")) tmp_cs.ipos = res_atoms[a].index();
      else if(nucl=="H"&&(AN=="H"||AN=="HN"))          tmp_cs.ipos = res_atoms[a].index();
      else if(nucl=="N"&&AN=="N")                      tmp_cs.ipos = res_atoms[a].index();
      else if(nucl=="CA"&&AN=="CA")                    tmp_cs.ipos = res_atoms[a].index();
      else if(nucl=="CB"&&AN=="CB")                    tmp_cs.ipos = res_atoms[a].index();
      else if(nucl=="C"&&AN=="C" )                     tmp_cs.ipos = res_atoms[a].index();
    }

    vector<AtomNumber> prev_res_atoms = pdb.getAtomsInResidue(resnum-1, chains[ichain]);
    // find the position of the previous residues backbone atoms
    for(unsigned a=0; a<prev_res_atoms.size(); a++) {
      string AN = pdb.getAtomName(prev_res_atoms[a]);
      if(AN=="N")                             { tmp_cs.bb[Np]  = prev_res_atoms[a].index(); }
      else if(AN=="CA")                       { tmp_cs.bb[CAp] = prev_res_atoms[a].index(); }
      else if(AN=="HA"||AN=="HA1"||AN=="HA3") { tmp_cs.bb[HAp] = prev_res_atoms[a].index(); }
      else if(AN=="C" )                       { tmp_cs.bb[Cp]  = prev_res_atoms[a].index(); }
      else if(AN=="O" )                       { tmp_cs.bb[Op]  = prev_res_atoms[a].index(); }
    }

    for(unsigned a=0; a<res_atoms.size(); a++) {
      string AN = pdb.getAtomName(res_atoms[a]);
      if(AN=="N")                                         { tmp_cs.bb[Nc]  = res_atoms[a].index(); }
      else if(AN=="H" ||AN=="HN"||(AN=="CD"&&RES=="PRO")) { tmp_cs.bb[Hc]  = res_atoms[a].index(); }
      else if(AN=="CA")                                   { tmp_cs.bb[CAc] = res_atoms[a].index(); }
      else if(AN=="HA"||AN=="HA1"||AN=="HA3")             { tmp_cs.bb[HAc] = res_atoms[a].index(); }
      else if(AN=="C" )                                   { tmp_cs.bb[Cc]  = res_atoms[a].index(); }
      else if(AN=="O" )                                   { tmp_cs.bb[Oc]  = res_atoms[a].index(); }

      if(RES!="ALA"&&RES!="GLY") {
        if(AN=="CB") tmp_cs.bb[CBc] = res_atoms[a].index();
        if(is_chi1_cx(RES,AN)) tmp_cs.bb[CGc] = res_atoms[a].index();
      }
    }

    vector<AtomNumber> next_res_atoms = pdb.getAtomsInResidue(resnum+1, chains[ichain]);
    string NRES = pdb.getResidueName(resnum+1, chains[ichain]);
    // find the position of the previous residues backbone atoms
    for(unsigned a=0; a<next_res_atoms.size(); a++) {
      string AN = pdb.getAtomName(next_res_atoms[a]);
      if(AN=="N")                                          { tmp_cs.bb[Nn]  = next_res_atoms[a].index(); }
      else if(AN=="H" ||AN=="HN"||(AN=="CD"&&NRES=="PRO")) { tmp_cs.bb[Hn]  = next_res_atoms[a].index(); }
      else if(AN=="CA")                                    { tmp_cs.bb[CAn] = next_res_atoms[a].index(); }
      else if(AN=="HA"||AN=="HA1"||AN=="HA3")              { tmp_cs.bb[HAn] = next_res_atoms[a].index(); }
      else if(AN=="C" )                                    { tmp_cs.bb[Cn]  = next_res_atoms[a].index(); }
    }

    // set sidechain atoms
    vector<string> sc_atm = side_chain_atoms(RES);

    for(unsigned sc=0; sc<sc_atm.size(); sc++) {
      for(unsigned aa=0; aa<res_atoms.size(); aa++) {
        if(pdb.getAtomName(res_atoms[aa])==sc_atm[sc]) {
          tmp_cs.side_chain.push_back(res_atoms[aa].index());
        }
      }
    }

    // find atoms for extra distances
    const string atomsP1[] =  {"H", "H", "H", "C", "C", "C", "O", "O", "O", "N", "N", "N", "O", "O", "O", "N", "N", "N", "CG", "CG", "CG", "CG", "CG", "CG", "CG", "CA"};
    const int resOffsetP1[] = { 0,   0,   0,  -1,  -1,  -1,   0,   0,   0,   1,   1,   1,   -1,  -1,  -1,  0,   0,   0,   0,    0,    0,    0,    0,    -1,   1,    -1};

    const string atomsP2[] =  {"HA", "C", "CB", "HA", "C", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "CB", "HA", "N", "C", "C", "N", "CA", "CA", "CA"};
    const int resOffsetP2[] = { 0,    0,   0,    0,    0,   0,    0,    0,   0,    0,    0,   0,    -1,  -1,   -1,   -1,  -1,   -1,   0,    0,   0,   -1,  1,   0,    0,    1};

    for(unsigned q=0; q<db.get_numXtraDists()-1; q++) {
      vector<AtomNumber> at1;
      if(resOffsetP1[q]== 0) at1 = res_atoms;
      if(resOffsetP1[q]==-1) at1 = prev_res_atoms;
      if(resOffsetP1[q]==+1) at1 = next_res_atoms;

      vector<AtomNumber> at2;
      if(resOffsetP2[q]== 0) at2 = res_atoms;
      if(resOffsetP2[q]==-1) at2 = prev_res_atoms;
      if(resOffsetP2[q]==+1) at2 = next_res_atoms;

      int tmp1 = -1;
      for(unsigned a=0; a<at1.size(); a++) {
        string name = pdb.getAtomName(at1[a]);
        xdist_name_map(name);

        if(name==atomsP1[q]) {
          tmp1 = at1[a].index();
          break;
        }
      }

      int tmp2 = -1;
      for(unsigned a=0; a<at2.size(); a++) {
        string name = pdb.getAtomName(at2[a]);
        xdist_name_map(name);

        if(name==atomsP2[q]) {
          tmp2 = at2[a].index();
          break;
        }
      }

      tmp_cs.xd1.push_back(tmp1);
      tmp_cs.xd2.push_back(tmp2);
    }

    // ready to add a new chemical shifts
    tmp_cs.csatoms = 1 + 16 + tmp_cs.side_chain.size() + 2*tmp_cs.xd1.size();
    if(tmp_cs.res_name!="ALA"&&tmp_cs.res_name!="GLY") tmp_cs.csatoms += 2;
    chemicalshifts.push_back(tmp_cs);
  }

  in.close();
}

// this assigns an atom-type to each atom of the pdb
void CS2Backbone::init_types(const PDB &pdb) {
  enum atom_t {D_C, D_H, D_N, D_O, D_S, D_C2, D_N2, D_O2};
  vector<AtomNumber> aa = pdb.getAtomNumbers();
  for(unsigned i=0; i<aa.size(); i++) {
    unsigned frag = pdb.getResidueNumber(aa[i]);
    string fragName = pdb.getResidueName(aa[i]);
    string atom_name = pdb.getAtomName(aa[i]);
    char atom_type = atom_name[0];
    if(isdigit(atom_name[0])) atom_type = atom_name[1];
    res_num.push_back(frag);
    unsigned t = 0;
    if (!isSP2(fragName, atom_name)) {
      if (atom_type == 'C') t = D_C;
      else if (atom_type == 'O') t = D_O;
      else if (atom_type == 'H') t = D_H;
      else if (atom_type == 'N') t = D_N;
      else if (atom_type == 'S') t = D_S;
      else plumed_merror("Unknown atom type: " + atom_name);
    } else {
      if (atom_type == 'C') t = D_C2;
      else if (atom_type == 'O') t = D_O2;
      else if (atom_type == 'N') t = D_N2;
      else plumed_merror("Unknown atom type: " + atom_name);
    }
    type.push_back(t);
  }
}

void CS2Backbone::init_rings(const PDB &pdb)
{
  const string pheTyr_n[] = {"CG","CD1","CE1","CZ","CE2","CD2"};
  const string trp1_n[]   = {"CD2","CE2","CZ2","CH2","CZ3","CE3"};
  const string trp2_n[]   = {"CG","CD1","NE1","CE2","CD2"};
  const string his_n[]    = {"CG","ND1","CD2","CE1","NE2"};

  // number of chains
  vector<string> chains;
  pdb.getChainNames( chains );
  unsigned total_rings_atoms = 0;

  // cycle over chains
  for(unsigned i=0; i<chains.size(); i++) {
    unsigned start, end;
    string errmsg;
    pdb.getResidueRange( chains[i], start, end, errmsg );
    // cycle over residues
    for(unsigned res=start; res<end; res++) {
      string frg = pdb.getResidueName(res, chains[i]);
      if(!((frg=="PHE")||(frg=="TYR")||(frg=="TRP")||
           (frg=="HIS")||(frg=="HIP")||(frg=="HID")||
           (frg=="HIE")||(frg=="HSD")||(frg=="HSE")||
           (frg=="HSP"))) continue;

      vector<AtomNumber> frg_atoms = pdb.getAtomsInResidue(res,chains[i]);

      if(frg=="PHE"||frg=="TYR") {
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index();
          for(unsigned aa=0; aa<6; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==pheTyr_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 6;
        total_rings_atoms += 6;
        if(frg=="PHE") ri.rtype = RingInfo::R_PHE;
        if(frg=="TYR") ri.rtype = RingInfo::R_TYR;
        ringInfo.push_back(ri);

      } else if(frg=="TRP") {
        //First ring
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index();
          for(unsigned aa=0; aa<6; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==trp1_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 6;
        total_rings_atoms += 6;
        ri.rtype = RingInfo::R_TRP1;
        ringInfo.push_back(ri);
        //Second Ring
        RingInfo ri2;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index();
          for(unsigned aa=0; aa<5; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==trp2_n[aa]) {
              ri2.atom[aa] = atm;
              break;
            }
          }
        }
        ri2.numAtoms = 5;
        total_rings_atoms += 3;
        ri2.rtype = RingInfo::R_TRP2;
        ringInfo.push_back(ri2);

      } else if((frg=="HIS")||(frg=="HIP")||(frg=="HID")||
                (frg=="HIE")||(frg=="HSD")||(frg=="HSE")||
                (frg=="HSP")) {//HIS case
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index();
          for(unsigned aa=0; aa<5; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==his_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 5;
        total_rings_atoms += 3;
        ri.rtype = RingInfo::R_HIS;
        ringInfo.push_back(ri);
      } else {
        plumed_merror("Unknown Ring Fragment: " + frg);
      }
    }
  }

  for(unsigned cs=0; cs<chemicalshifts.size(); cs++) chemicalshifts[cs].csatoms += total_rings_atoms;
}

void CS2Backbone::calculate()
{
  if(pbc) makeWhole();
  if(getExchangeStep()) box_count=0;
  if(box_count==0) update_neighb();
  compute_ring_parameters();

  vector<double> camshift_sigma2(6);
  camshift_sigma2[0] = 0.08; // HA
  camshift_sigma2[1] = 0.30; // HN
  camshift_sigma2[2] = 9.00; // NH
  camshift_sigma2[3] = 1.30; // CA
  camshift_sigma2[4] = 1.56; // CB
  camshift_sigma2[5] = 1.70; // CO

  vector<Vector>   cs_derivs;
  vector<Vector>   aa_derivs;
  vector<unsigned> cs_atoms;
  vector<double>   all_shifts;

  cs_derivs.resize(chemicalshifts.size()*max_cs_atoms,Vector(0,0,0));
  cs_atoms.resize(chemicalshifts.size()*max_cs_atoms,0);
  all_shifts.resize(chemicalshifts.size(),0);
  if(camshift||getDoScore()) aa_derivs.resize(getNumberOfAtoms(),Vector(0,0,0));

  unsigned stride = comm.Get_size();
  unsigned rank   = comm.Get_rank();
  if(serial) {
    stride = 1;
    rank   = 0;
  }

  unsigned nt=OpenMP::getNumThreads();
  if(nt*stride*2>chemicalshifts.size()) nt=1;

  // a single loop over all chemical shifts
  #pragma omp parallel num_threads(nt)
  {
    #pragma omp for schedule(dynamic)
    for(unsigned cs=rank; cs<chemicalshifts.size(); cs+=stride) {
      const unsigned kdx=cs*max_cs_atoms;
      const ChemicalShift *myfrag = &chemicalshifts[cs];
      const unsigned aa_kind = myfrag->res_kind;
      const unsigned at_kind = myfrag->atm_kind;

      double shift = db.CONSTAAPREV(aa_kind,at_kind)[myfrag->res_type_prev] +
                     db.CONSTAACURR(aa_kind,at_kind)[myfrag->res_type_curr] +
                     db.CONSTAANEXT(aa_kind,at_kind)[myfrag->res_type_next];

      const unsigned ipos = myfrag->ipos;
      cs_atoms[kdx+0] = ipos;
      unsigned atom_counter = 1;

      //BACKBONE (PREV CURR NEXT)
      const double * CONST_BB2 = db.CONST_BB2(aa_kind,at_kind);
      const unsigned bbsize = 16;
      for(unsigned q=0; q<bbsize; q++) {
        const double cb2q = CONST_BB2[q];
        if(cb2q==0.) continue;
        const unsigned jpos = myfrag->bb[q];
        if(ipos==jpos) continue;
        const Vector distance = delta(getPosition(jpos),getPosition(ipos));
        const double d = distance.modulo();
        const double fact = cb2q/d;

        shift += cb2q*d;
        const Vector der = fact*distance;

        cs_derivs[kdx+0] += der;
        cs_derivs[kdx+q+atom_counter] = -der;
        cs_atoms[kdx+q+atom_counter] = jpos;
      }

      atom_counter += bbsize;

      //DIHEDRAL ANGLES
      const double *CO_DA = db.CO_DA(aa_kind,at_kind);
      //Phi
      {
        const Vector d0 = delta(getPosition(myfrag->bb[Nc]), getPosition(myfrag->bb[Cp]));
        const Vector d1 = delta(getPosition(myfrag->bb[CAc]), getPosition(myfrag->bb[Nc]));
        const Vector d2 = delta(getPosition(myfrag->bb[Cc]), getPosition(myfrag->bb[CAc]));
        Torsion t;
        Vector dd0, dd1, dd2;
        const double t_phi = t.compute(d0,d1,d2,dd0,dd1,dd2);
        const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,0);
        const double val1 = 3.*t_phi+PARS_DA[3];
        const double val2 = t_phi+PARS_DA[4];
        shift += CO_DA[0]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
        const double fact = -CO_DA[0]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

        cs_derivs[kdx+Cp+1] += fact*dd0;
        cs_derivs[kdx+Nc+1] += fact*(dd1-dd0);
        cs_derivs[kdx+CAc+1]+= fact*(dd2-dd1);
        cs_derivs[kdx+Cc+1] += -fact*dd2;
        cs_atoms[kdx+Cp+1] = myfrag->bb[Cp];
        cs_atoms[kdx+Nc+1] = myfrag->bb[Nc];
        cs_atoms[kdx+CAc+1]= myfrag->bb[CAc];
        cs_atoms[kdx+Cc+1] = myfrag->bb[Cc];
      }

      //Psi
      {
        const Vector d0 = delta(getPosition(myfrag->bb[CAc]), getPosition(myfrag->bb[Nc]));
        const Vector d1 = delta(getPosition(myfrag->bb[Cc]), getPosition(myfrag->bb[CAc]));
        const Vector d2 = delta(getPosition(myfrag->bb[Nn]), getPosition(myfrag->bb[Cc]));
        Torsion t;
        Vector dd0, dd1, dd2;
        const double t_psi = t.compute(d0,d1,d2,dd0,dd1,dd2);
        const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,1);
        const double val1 = 3.*t_psi+PARS_DA[3];
        const double val2 = t_psi+PARS_DA[4];
        shift += CO_DA[1]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
        const double fact = -CO_DA[1]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

        cs_derivs[kdx+Nc+1] += fact*dd0;
        cs_derivs[kdx+CAc+1] += fact*(dd1-dd0);
        cs_derivs[kdx+Cc+1] += fact*(dd2-dd1);
        cs_derivs[kdx+Nn+1] += -fact*dd2;
        cs_atoms[kdx+Nc+1] = myfrag->bb[Nc];
        cs_atoms[kdx+CAc+1]= myfrag->bb[CAc];
        cs_atoms[kdx+Cc+1] = myfrag->bb[Cc];
        cs_atoms[kdx+Nn+1] = myfrag->bb[Nn];
      }

      //Chi
      if(myfrag->has_chi1) {
        const Vector d0 = delta(getPosition(myfrag->bb[CAc]), getPosition(myfrag->bb[Nc]));
        const Vector d1 = delta(getPosition(myfrag->bb[CBc]), getPosition(myfrag->bb[CAc]));
        const Vector d2 = delta(getPosition(myfrag->bb[CGc]), getPosition(myfrag->bb[CBc]));
        Torsion t;
        Vector dd0, dd1, dd2;
        const double t_chi1 = t.compute(d0,d1,d2,dd0,dd1,dd2);
        const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,2);
        const double val1 = 3.*t_chi1+PARS_DA[3];
        const double val2 = t_chi1+PARS_DA[4];
        shift += CO_DA[2]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
        const double fact = -CO_DA[2]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

        cs_derivs[kdx+Nc+1] += fact*dd0;
        cs_derivs[kdx+CAc+1] += fact*(dd1-dd0);
        cs_derivs[kdx+CBc+1] += fact*(dd2-dd1);
        cs_derivs[kdx+CGc+1] += -fact*dd2;
        cs_atoms[kdx+Nc+1]  = myfrag->bb[Nc];
        cs_atoms[kdx+CAc+1] = myfrag->bb[CAc];
        cs_atoms[kdx+CBc+1] = myfrag->bb[CBc];
        cs_atoms[kdx+CGc+1] = myfrag->bb[CGc];

        atom_counter += 2;
      }
      //END OF DIHE

      //SIDE CHAIN
      const double * CONST_SC2 = db.CONST_SC2(aa_kind,at_kind,myfrag->res_type_curr);
      const unsigned sidsize = myfrag->side_chain.size();
      for(unsigned q=0; q<sidsize; q++) {
        const double cs2q = CONST_SC2[q];
        if(cs2q==0.) continue;
        const unsigned jpos = myfrag->side_chain[q];
        if(ipos==jpos) continue;
        const Vector distance = delta(getPosition(jpos),getPosition(ipos));
        const double d = distance.modulo();
        const double fact = cs2q/d;

        shift += cs2q*d;
        const Vector der = fact*distance;
        cs_derivs[kdx+0] += der;
        cs_derivs[kdx+q+atom_counter] = -der;
        cs_atoms[kdx+q+atom_counter] = jpos;
      }

      atom_counter += sidsize;

      //EXTRA DIST
      const double * CONST_XD  = db.CONST_XD(aa_kind,at_kind);
      const unsigned xdsize=myfrag->xd1.size();
      for(unsigned q=0; q<xdsize; q++) {
        const double cxdq = CONST_XD[q];
        if(cxdq==0.) continue;
        if(myfrag->xd1[q]==-1||myfrag->xd2[q]==-1) continue;
        const Vector distance = delta(getPosition(myfrag->xd1[q]),getPosition(myfrag->xd2[q]));
        const double d = distance.modulo();
        const double fact = cxdq/d;

        shift += cxdq*d;
        const Vector der = fact*distance;
        cs_derivs[kdx+2*q+atom_counter  ] = der;
        cs_derivs[kdx+2*q+atom_counter+1] = -der;
        cs_atoms[kdx+2*q+atom_counter] = myfrag->xd2[q];
        cs_atoms[kdx+2*q+atom_counter+1] = myfrag->xd1[q];
      }

      atom_counter += 2*xdsize;

      //RINGS
      const double *rc = db.CO_RING(aa_kind,at_kind);
      const unsigned rsize = ringInfo.size();
      // cycle over the list of rings
      for(unsigned q=0; q<rsize; q++) {
        // compute angle from ring middle point to current atom position
        // get distance vector from query atom to ring center and normal vector to ring plane
        const Vector n   = ringInfo[q].normVect;
        const double nL  = ringInfo[q].lengthNV;
        const double inL2 = ringInfo[q].lengthN2;

        const Vector d = delta(ringInfo[q].position, getPosition(ipos));
        const double dL2 = d.modulo2();
        double dL  = sqrt(dL2);
        const double idL3 = 1./(dL2*dL);

        const double dn    = dotProduct(d,n);
        const double dn2   = dn*dn;
        const double dLnL  = dL*nL;
        const double dL_nL = dL/nL;

        const double ang2 = dn2*inL2/dL2;
        const double u    = 1.-3.*ang2;
        const double cc   = rc[ringInfo[q].rtype];

        shift += cc*u*idL3;

        const double fUU    = -6.*dn*inL2;
        const double fUQ    = fUU/dL;
        const Vector gradUQ = fUQ*(dL2*n - dn*d);
        const Vector gradVQ = (3.*dL*u)*d;

        const double fact   = cc*idL3*idL3;
        cs_derivs[kdx+0] += fact*(gradUQ - gradVQ);

        const double fU       = fUU/nL;
        double OneOverN = 1./6.;
        if(ringInfo[q].numAtoms==5) OneOverN=1./3.;
        const Vector factor2  = OneOverN*n;
        const Vector factor4  = (OneOverN/dL_nL)*d;

        const Vector gradV    = -OneOverN*gradVQ;

        if(ringInfo[q].numAtoms==6) {
          // update forces on ring atoms
          for(unsigned at=0; at<6; at++) {
            const Vector ab = crossProduct(d,ringInfo[q].g[at]);
            const Vector c  = crossProduct(n,ringInfo[q].g[at]);
            const Vector factor3 = 0.5*dL_nL*c;
            const Vector factor1 = 0.5*ab;
            const Vector gradU   = fU*( dLnL*(factor1 - factor2) -dn*(factor3 - factor4) );
            cs_derivs[kdx+at+atom_counter] = fact*(gradU - gradV);
            cs_atoms[kdx+at+atom_counter] = ringInfo[q].atom[at];
          }
          atom_counter += 6;
        }  else {
          for(unsigned at=0; at<3; at++) {
            const Vector ab = crossProduct(d,ringInfo[q].g[at]);
            const Vector c  = crossProduct(n,ringInfo[q].g[at]);
            const Vector factor3 = dL_nL*c;
            const Vector factor1 = ab;
            const Vector gradU   = fU*( dLnL*(factor1 - factor2) -dn*(factor3 - factor4) );
            cs_derivs[kdx+at+atom_counter] = fact*(gradU - gradV);
          }
          cs_atoms[kdx+atom_counter] = ringInfo[q].atom[0];
          cs_atoms[kdx+atom_counter+1] = ringInfo[q].atom[2];
          cs_atoms[kdx+atom_counter+2] = ringInfo[q].atom[3];
          atom_counter += 3;
        }
      }
      //END OF RINGS

      //NON BOND
      const double * CONST_CO_SPHERE3 = db.CO_SPHERE(aa_kind,at_kind,0);
      const double * CONST_CO_SPHERE  = db.CO_SPHERE(aa_kind,at_kind,1);
      const unsigned boxsize = myfrag->box_nb.size();
      for(unsigned q=0; q<boxsize; q++) {
        const unsigned jpos = myfrag->box_nb[q];
        const Vector distance = delta(getPosition(jpos),getPosition(ipos));
        const double d2 = distance.modulo2();

        if(d2<cutOffDist2) {
          double factor1  = sqrt(d2);
          double dfactor1 = 1./factor1;
          double factor3  = dfactor1*dfactor1*dfactor1;
          double dfactor3 = -3.*factor3*dfactor1*dfactor1;

          if(d2>cutOnDist2) {
            const double af = cutOffDist2 - d2;
            const double bf = cutOffDist2 - 3.*cutOnDist2 + 2.*d2;
            const double cf = invswitch*af;
            const double df = cf*af*bf;
            factor1 *= df;
            factor3 *= df;

            const double d4  = d2*d2;
            const double af1 = 15.*cutOnDist2*d2;
            const double bf1 = -14.*d4;
            const double cf1 = -3.*cutOffDist2*cutOnDist2 + cutOffDist2*d2;
            const double df1 = af1+bf1+cf1;
            dfactor1 *= cf*(cutOffDist4+df1);

            const double af3 = +2.*cutOffDist2*cutOnDist2;
            const double bf3 = d2*(cutOffDist2+cutOnDist2);
            const double cf3 = -2.*d4;
            const double df3 = (af3+bf3+cf3)*d2;
            dfactor3 *= invswitch*(cutMixed+df3);
          }

          const unsigned t = type[jpos];
          shift += factor1*CONST_CO_SPHERE[t] + factor3*CONST_CO_SPHERE3[t] ;
          const double fact = dfactor1*CONST_CO_SPHERE[t]+dfactor3*CONST_CO_SPHERE3[t];
          const Vector der  = fact*distance;

          cs_derivs[kdx+0] += der;
          cs_derivs[kdx+q+atom_counter] = -der;
          cs_atoms[kdx+q+atom_counter] = jpos;
        }
      }
      //END NON BOND

      atom_counter += boxsize;
      all_shifts[cs] = shift;
    }
  }

  ++box_count;
  if(box_count == box_nupdate) box_count = 0;

  if(!camshift) {
    if(!serial) {
      if(!getDoScore()) {
        comm.Sum(&cs_derivs[0][0], 3*cs_derivs.size());
        comm.Sum(&cs_atoms[0], cs_atoms.size());
      }
      comm.Sum(&all_shifts[0], chemicalshifts.size());
    }
    for(unsigned cs=0; cs<chemicalshifts.size(); cs++) {
      Value *comp = chemicalshifts[cs].comp;
      comp->set(all_shifts[cs]);
      if(getDoScore()) setCalcData(cs, all_shifts[cs]);
      else {
        const unsigned kdx=cs*max_cs_atoms;
        Tensor csvirial;
        for(unsigned i=0; i<chemicalshifts[cs].totcsatoms; i++) {
          setAtomsDerivatives(comp,cs_atoms[kdx+i],cs_derivs[kdx+i]);
          csvirial-=Tensor(getPosition(cs_atoms[kdx+i]),cs_derivs[kdx+i]);
        }
        setBoxDerivatives(comp,csvirial);
      }
    }
    if(!getDoScore()) return;
  }

  double score = 0.;

  /* Metainference */
  if(getDoScore()) {
    score = getScore();
    for(unsigned cs=rank; cs<chemicalshifts.size(); cs+=stride) {
      const unsigned kdx=cs*max_cs_atoms;
      const double fact = getMetaDer(cs);
      for(unsigned i=0; i<chemicalshifts[cs].totcsatoms; i++) {
        aa_derivs[cs_atoms[kdx+i]] += cs_derivs[kdx+i]*fact;
      }
    }
  }

  /* camshift */
  if(camshift) {
    for(unsigned cs=rank; cs<chemicalshifts.size(); cs+=stride) {
      const unsigned kdx=cs*max_cs_atoms;
      score += (all_shifts[cs] - chemicalshifts[cs].exp_cs)*(all_shifts[cs] - chemicalshifts[cs].exp_cs)/camshift_sigma2[chemicalshifts[cs].atm_kind];
      double fact = 2.0*(all_shifts[cs] - chemicalshifts[cs].exp_cs)/camshift_sigma2[chemicalshifts[cs].atm_kind];
      for(unsigned i=0; i<chemicalshifts[cs].totcsatoms; i++) {
        aa_derivs[cs_atoms[kdx+i]] += cs_derivs[kdx+i]*fact;
      }
    }
  }

  if(!serial) {
    comm.Sum(&aa_derivs[0][0], 3*aa_derivs.size());
    if(camshift) comm.Sum(&score, 1);
  }

  Tensor virial;
  for(unsigned i=rank; i<getNumberOfAtoms(); i+=stride) {
    virial += Tensor(getPosition(i), aa_derivs[i]);
  }

  if(!serial) {
    comm.Sum(&virial[0][0], 9);
  }

  /* calculate final derivatives */
  Value* val;
  if(getDoScore()) {
    val=getPntrToComponent("score");
    setScore(score);
  } else {
    val=getPntrToValue();
    setValue(score);
  }

  /* at this point we cycle over all atoms */
  for(unsigned i=0; i<getNumberOfAtoms(); i++) setAtomsDerivatives(val, i,  aa_derivs[i]);
  setBoxDerivatives(val,-virial);
}

void CS2Backbone::update_neighb() {
  max_cs_atoms=0;
  // cycle over chemical shifts
  for(unsigned cs=0; cs<chemicalshifts.size(); cs++) {
    const unsigned boxsize = getNumberOfAtoms();
    chemicalshifts[cs].box_nb.clear();
    chemicalshifts[cs].box_nb.reserve(150);
    const unsigned res_curr = res_num[chemicalshifts[cs].ipos];
    for(unsigned bat=0; bat<boxsize; bat++) {
      const unsigned res_dist = abs(static_cast<int>(res_curr-res_num[bat]));
      if(res_dist<2) continue;
      const Vector distance = delta(getPosition(bat),getPosition(chemicalshifts[cs].ipos));
      const double d2=distance.modulo2();
      if(d2<cutOffNB2) chemicalshifts[cs].box_nb.push_back(bat);
    }
    chemicalshifts[cs].totcsatoms = chemicalshifts[cs].csatoms + chemicalshifts[cs].box_nb.size();
    if(chemicalshifts[cs].totcsatoms>max_cs_atoms) max_cs_atoms = chemicalshifts[cs].totcsatoms;
  }
}

void CS2Backbone::compute_ring_parameters() {
  for(unsigned i=0; i<ringInfo.size(); i++) {
    const unsigned size = ringInfo[i].numAtoms;
    if(size==6) {
      ringInfo[i].g[0] = delta(getPosition(ringInfo[i].atom[4]),getPosition(ringInfo[i].atom[2]));
      ringInfo[i].g[1] = delta(getPosition(ringInfo[i].atom[5]),getPosition(ringInfo[i].atom[3]));
      ringInfo[i].g[2] = delta(getPosition(ringInfo[i].atom[0]),getPosition(ringInfo[i].atom[4]));
      ringInfo[i].g[3] = delta(getPosition(ringInfo[i].atom[1]),getPosition(ringInfo[i].atom[5]));
      ringInfo[i].g[4] = delta(getPosition(ringInfo[i].atom[2]),getPosition(ringInfo[i].atom[0]));
      ringInfo[i].g[5] = delta(getPosition(ringInfo[i].atom[3]),getPosition(ringInfo[i].atom[1]));
      vector<Vector> a(6);
      a[0] = getPosition(ringInfo[i].atom[0]);
      // ring center
      Vector midP = a[0];
      for(unsigned j=1; j<size; j++) {
        a[j] = getPosition(ringInfo[i].atom[j]);
        midP += a[j];
      }
      ringInfo[i].position = midP/6.;
      // compute normal vector to plane
      Vector n1 = crossProduct(delta(a[0],a[4]), delta(a[0],a[2]));
      Vector n2 = crossProduct(delta(a[3],a[1]), delta(a[3],a[5]));
      ringInfo[i].normVect = 0.5*(n1 + n2);
    }  else {
      ringInfo[i].g[0] = delta(getPosition(ringInfo[i].atom[3]),getPosition(ringInfo[i].atom[2]));
      ringInfo[i].g[1] = delta(getPosition(ringInfo[i].atom[0]),getPosition(ringInfo[i].atom[3]));
      ringInfo[i].g[2] = delta(getPosition(ringInfo[i].atom[2]),getPosition(ringInfo[i].atom[0]));
      vector<Vector> a(size);
      for(unsigned j=0; j<size; j++) {
        a[j] = getPosition(ringInfo[i].atom[j]);
      }
      // ring center
      ringInfo[i].position = (a[0]+a[2]+a[3])/3.;
      // ring plane normal vector
      ringInfo[i].normVect = crossProduct(delta(a[0],a[3]), delta(a[0],a[2]));

    }
    // calculate squared length and length of normal vector
    ringInfo[i].lengthN2 = 1./ringInfo[i].normVect.modulo2();
    ringInfo[i].lengthNV = 1./sqrt(ringInfo[i].lengthN2);
  }
}

CS2Backbone::aa_t CS2Backbone::frag2enum(const string &aa) {
  aa_t type = ALA;
  if (aa == "ALA") type = ALA;
  else if (aa == "ARG") type = ARG;
  else if (aa == "ASN") type = ASN;
  else if (aa == "ASP") type = ASP;
  else if (aa == "ASH") type = ASP;
  else if (aa == "CYS") type = CYS;
  else if (aa == "CYM") type = CYS;
  else if (aa == "GLN") type = GLN;
  else if (aa == "GLU") type = GLU;
  else if (aa == "GLH") type = GLU;
  else if (aa == "GLY") type = GLY;
  else if (aa == "HIS") type = HIS;
  else if (aa == "HSE") type = HIS;
  else if (aa == "HIE") type = HIS;
  else if (aa == "HSP") type = HIS;
  else if (aa == "HIP") type = HIS;
  else if (aa == "HSD") type = HIS;
  else if (aa == "HID") type = HIS;
  else if (aa == "ILE") type = ILE;
  else if (aa == "LEU") type = LEU;
  else if (aa == "LYS") type = LYS;
  else if (aa == "MET") type = MET;
  else if (aa == "PHE") type = PHE;
  else if (aa == "PRO") type = PRO;
  else if (aa == "SER") type = SER;
  else if (aa == "THR") type = THR;
  else if (aa == "TRP") type = TRP;
  else if (aa == "TYR") type = TYR;
  else if (aa == "VAL") type = VAL;
  else if (aa == "UNK") type = UNK;
  else plumed_merror("Error converting string " + aa + " into amino acid index: not a valid 3-letter code");
  return type;
}

vector<string> CS2Backbone::side_chain_atoms(const string &s) {
  vector<string> sc;

  if(s=="ALA") {
    sc.push_back( "CB" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    return sc;
  } else if(s=="ARG") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD" );
    sc.push_back( "NE" );
    sc.push_back( "CZ" );
    sc.push_back( "NH1" );
    sc.push_back( "NH2" );
    sc.push_back( "NH3" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    sc.push_back( "HE" );
    sc.push_back( "HH11" );
    sc.push_back( "HH12" );
    sc.push_back( "HH21" );
    sc.push_back( "HH22" );
    sc.push_back( "1HH1" );
    sc.push_back( "2HH1" );
    sc.push_back( "1HH2" );
    sc.push_back( "2HH2" );
    return sc;
  } else if(s=="ASN") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "OD1" );
    sc.push_back( "ND2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HD21" );
    sc.push_back( "HD22" );
    sc.push_back( "1HD2" );
    sc.push_back( "2HD2" );
    return sc;
  } else if(s=="ASP"||s=="ASH") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "OD1" );
    sc.push_back( "OD2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    return sc;
  } else if(s=="CYS"||s=="CYM") {
    sc.push_back( "CB" );
    sc.push_back( "SG" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG" );
    return sc;
  } else if(s=="GLN") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD" );
    sc.push_back( "OE1" );
    sc.push_back( "NE2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    sc.push_back( "HE21" );
    sc.push_back( "HE22" );
    sc.push_back( "1HE2" );
    sc.push_back( "2HE2" );
    return sc;
  } else if(s=="GLU"||s=="GLH") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD" );
    sc.push_back( "OE1" );
    sc.push_back( "OE2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    return sc;
  } else if(s=="GLY") {
    sc.push_back( "HA2" );
    return sc;
  } else if(s=="HIS"||s=="HSE"||s=="HIE"||s=="HSD"||s=="HID"||s=="HIP"||s=="HSP") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "ND1" );
    sc.push_back( "CD2" );
    sc.push_back( "CE1" );
    sc.push_back( "NE2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HE1" );
    sc.push_back( "HE2" );
    return sc;
  } else if(s=="ILE") {
    sc.push_back( "CB" );
    sc.push_back( "CG1" );
    sc.push_back( "CG2" );
    sc.push_back( "CD" );
    sc.push_back( "HB" );
    sc.push_back( "HG11" );
    sc.push_back( "HG12" );
    sc.push_back( "HG21" );
    sc.push_back( "HG22" );
    sc.push_back( "HG23" );
    sc.push_back( "1HG1" );
    sc.push_back( "2HG1" );
    sc.push_back( "1HG2" );
    sc.push_back( "2HG2" );
    sc.push_back( "3HG2" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    return sc;
  } else if(s=="LEU") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD1" );
    sc.push_back( "CD2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG" );
    sc.push_back( "HD11" );
    sc.push_back( "HD12" );
    sc.push_back( "HD13" );
    sc.push_back( "HD21" );
    sc.push_back( "HD22" );
    sc.push_back( "HD23" );
    sc.push_back( "1HD1" );
    sc.push_back( "2HD1" );
    sc.push_back( "3HD1" );
    sc.push_back( "1HD2" );
    sc.push_back( "2HD2" );
    sc.push_back( "3HD2" );
    return sc;
  } else if(s=="LYS") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD" );
    sc.push_back( "CE" );
    sc.push_back( "NZ" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    sc.push_back( "HE1" );
    sc.push_back( "HE2" );
    sc.push_back( "HE3" );
    sc.push_back( "HZ1" );
    sc.push_back( "HZ2" );
    sc.push_back( "HZ3" );
    return sc;
  } else if(s=="MET") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "SD" );
    sc.push_back( "CE" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    sc.push_back( "HE1" );
    sc.push_back( "HE2" );
    sc.push_back( "HE3" );
    return sc;
  } else if(s=="PHE") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD1" );
    sc.push_back( "CD2" );
    sc.push_back( "CE1" );
    sc.push_back( "CE2" );
    sc.push_back( "CZ" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    sc.push_back( "HE1" );
    sc.push_back( "HE2" );
    sc.push_back( "HE3" );
    sc.push_back( "HZ" );
    return sc;
  } else if(s=="PRO") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG2" );
    sc.push_back( "HG3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    return sc;
  } else if(s=="SER") {
    sc.push_back( "CB" );
    sc.push_back( "OG" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG" );
    return sc;
  } else if(s=="THR") {
    sc.push_back( "CB" );
    sc.push_back( "OG1" );
    sc.push_back( "CG2" );
    sc.push_back( "HB" );
    sc.push_back( "HG1" );
    sc.push_back( "HG21" );
    sc.push_back( "HG22" );
    sc.push_back( "HG23" );
    sc.push_back( "1HG2" );
    sc.push_back( "2HG2" );
    sc.push_back( "3HG2" );
    return sc;
  } else if(s=="TRP") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD1" );
    sc.push_back( "CD2" );
    sc.push_back( "NE1" );
    sc.push_back( "CE2" );
    sc.push_back( "CE3" );
    sc.push_back( "CZ2" );
    sc.push_back( "CZ3" );
    sc.push_back( "CH2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HD1" );
    sc.push_back( "HE1" );
    sc.push_back( "HE3" );
    sc.push_back( "HZ2" );
    sc.push_back( "HZ3" );
    sc.push_back( "HH2" );
    return sc;
  } else if(s=="TYR") {
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "CD1" );
    sc.push_back( "CD2" );
    sc.push_back( "CE1" );
    sc.push_back( "CE2" );
    sc.push_back( "CZ" );
    sc.push_back( "OH" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HD1" );
    sc.push_back( "HD2" );
    sc.push_back( "HD3" );
    sc.push_back( "HE1" );
    sc.push_back( "HE2" );
    sc.push_back( "HE3" );
    sc.push_back( "HH" );
    return sc;
  } else if(s=="VAL") {
    sc.push_back( "CB" );
    sc.push_back( "CG1" );
    sc.push_back( "CG2" );
    sc.push_back( "HB" );
    sc.push_back( "HG11" );
    sc.push_back( "HG12" );
    sc.push_back( "HG13" );
    sc.push_back( "HG21" );
    sc.push_back( "HG22" );
    sc.push_back( "HG23" );
    sc.push_back( "1HG1" );
    sc.push_back( "2HG1" );
    sc.push_back( "3HG1" );
    sc.push_back( "1HG2" );
    sc.push_back( "2HG2" );
    sc.push_back( "3HG2" );
    return sc;
  } else plumed_merror("Sidechain atoms unknown: " + s);
}

bool CS2Backbone::isSP2(const string & resType, const string & atomName) {
  bool sp2 = false;
  if (atomName == "C") return true;
  if (atomName == "O") return true;

  if(resType == "TRP") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "CD1") sp2 = true;
    else if (atomName == "CD2") sp2 = true;
    else if (atomName == "CE2") sp2 = true;
    else if (atomName == "CE3") sp2 = true;
    else if (atomName == "CZ2") sp2 = true;
    else if (atomName == "CZ3") sp2 = true;
    else if (atomName == "CH2") sp2 = true;
  } else if (resType == "ASP") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "OD1") sp2 = true;
    else if (atomName == "OD2") sp2 = true;
  } else if (resType == "GLU") {
    if      (atomName == "CD")  sp2 = true;
    else if (atomName == "OE1") sp2 = true;
    else if (atomName == "OE2") sp2 = true;
  } else if (resType == "ARG") {
    if (atomName == "CZ") sp2 = true;
  } else if (resType == "HIS") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "ND1") sp2 = true;
    else if (atomName == "CD2") sp2 = true;
    else if (atomName == "CE1") sp2 = true;
    else if (atomName == "NE2") sp2 = true;
  } else if (resType == "PHE") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "CD1") sp2 = true;
    else if (atomName == "CD2") sp2 = true;
    else if (atomName == "CE1") sp2 = true;
    else if (atomName == "CE2") sp2 = true;
    else if (atomName == "CZ")  sp2 = true;
  } else if (resType == "TYR") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "CD1") sp2 = true;
    else if (atomName == "CD2") sp2 = true;
    else if (atomName == "CE1") sp2 = true;
    else if (atomName == "CE2") sp2 = true;
    else if (atomName == "CZ")  sp2 = true;
  } else if (resType == "ASN") {
    if      (atomName == "CG")  sp2 = true;
    else if (atomName == "OD1") sp2 = true;
  } else if (resType == "GLN") {
    if      (atomName == "CD")  sp2 = true;
    else if (atomName == "OE1") sp2 = true;
  }

  return sp2;
}

bool CS2Backbone::is_chi1_cx(const string & frg, const string & atm) {
  if(atm=="CG")                                        return true;
  if((frg == "CYS")&&(atm =="SG"))                     return true;
  if(((frg == "ILE")||(frg == "VAL"))&&(atm == "CG1")) return true;
  if((frg == "SER")&&(atm == "OG"))                    return true;
  if((frg == "THR")&&(atm == "OG1"))                   return true;

  return false;
}

void CS2Backbone::xdist_name_map(string & name) {
  if((name == "OT1")||(name == "OC1")) name = "O";
  else if ((name == "HN") || (name == "HT1") || (name == "H1")) name = "H";
  else if ((name == "CG1")|| (name == "OG")||
           (name == "SG") || (name == "OG1")) name = "CG";
  else if ((name == "HA1")|| (name == "HA3")) name = "HA";
}

void CS2Backbone::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
