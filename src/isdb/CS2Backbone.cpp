/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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

#define cutOffNB      0.70	// buffer distance for neighbour-lists 
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

The functional form is that of CamShift \cite Kohlhoff:2009us. The chemical shifts
of the selected nuclei/residues are saved as components. Reference experimental values
can also be stored as components. The two sets of components can then be used to calculate
either a scoring function as in \cite Robustelli:2010dn \cite Granata:2013dk, using
the keyword CAMSHIFT or to calculate ensemble averaged chemical shift as in \cite Camilloni:2012je
\cite Camilloni:2013hs (see \ref ENSEMBLE, \ref STATS and \ref RESTRAINT). Finally they can
also be used as input for \ref METAINFERENCE, \cite Bonomi:2016ip . In the current implementation there is
no need to pass the data to \ref METAINFERENCE because \ref CS2BACKBONE can internally enable Metainference
using the keywork DOSCORE.

CamShift calculation is relatively heavy because it often uses a large number of atoms, in order
to make it faster it is currently parallelised with \ref Openmp.

As a general rule, when using \ref CS2BACKBONE or other experimental restraints it is better to
increase the accuracy of the constraint algorithm due to the increased strain on the bonded structure.
In the case of GROMACS it is safer to use lincs-iter=2 and lincs-order=6.

In general the system for which chemical shifts are calculated must be completly included in
ATOMS and a TEMPLATE pdb file for the same atoms should be provided as well in the folder DATADIR.
The atoms are made automatically whole unless NOPBC is used, in particular if the system is made of
by multiple chains it is usually better to use NOPBC and make the molecule whole \ref WHOLEMOLECULES
selecting an appropriate order.

In addition to a pdb file one needs to provide a list of chemical shifts to be calculated using one
file per nucleus type (CAshifts.dat, CBshifts.dat, Cshifts.dat, Hshifts.dat, HAshifts.dat, Nshifts.dat),
all the six files should always be present. A chemical shift for a nucleus is calculated if a value
greater than 0 is provided. For practical purposes the value can correspond to the experimental value.
Residues numbers should go from 1 to N irrespectively of the numbers used in the pdb file. The first and
last residue of each chain should be preceeded by a # character. Termini groups like ACE or NME should
be removed from the PDB.

\verbatim
CAshifts.dat:
#1 0.0
2 55.5
3 58.4
.
.
#last 0.0
#last+1 (first) of second chain
.
#last of second chain
\endverbatim

The default behaviour is to store the values for the active nuclei in components (ca_#, cb_#,
co_#, ha_#, hn_#, nh_# and expca_#, expcb_#, expco_#, expha_#, exphn_#, exp_nh#) with NOEXP it is possible
to only store the backcalculated values.

A pdb file is needed to the generate a simple topology of the protein. For histidines in protonation
states different from D the HIE/HSE HIP/HSP name should be used. GLH and ASH can be used for the alternative
protonation of GLU and ASP. Non-standard amino acids and other molecules are not yet supported, but in principle
they can be named UNK. If multiple chains are present the chain identifier must be in the standard PDB format,
together with the TER keyword at the end of each chain.

One more standard file is also needed in the folder DATADIR: camshift.db. This file includes all the CamShift parameters
and can be found in regtest/isdb/rt-cs2backbone/data/ .

All the above files must be in a single folder that must be specified with the keyword DATADIR.

Additional material and examples can be also found in the tutorial \ref belfast-9

\par Examples

In this first example the chemical shifts are used to calculate a scoring function to be used
in NMR driven Metadynamics \cite Granata:2013dk :

\plumedfile
whole: GROUP ATOMS=2612-2514:-1,961-1:-1,2466-962:-1,2513-2467:-1
WHOLEMOLECULES ENTITY0=whole
cs: CS2BACKBONE ATOMS=1-2612 NRES=176 DATADIR=../data/ TEMPLATE=template.pdb CAMSHIFT NOPBC
metad: METAD ARG=cs HEIGHT=0.5 SIGMA=0.1 PACE=200 BIASFACTOR=10
PRINT ARG=cs,metad.bias FILE=COLVAR STRIDE=100
\endplumedfile

In this second example the chemical shifts are used as replica-averaged restrained as in \cite Camilloni:2012je \cite Camilloni:2013hs.

\plumedfile
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/ NRES=13
encs: ENSEMBLE ARG=(cs\.hn_.*),(cs\.nh_.*)
stcs: STATS ARG=encs.* SQDEVSUM PARARG=(cs\.exphn_.*),(cs\.expnh_.*)
RESTRAINT ARG=stcs.sqdevsum AT=0 KAPPA=0 SLOPE=24

PRINT ARG=(cs\.hn_.*),(cs\.nh_.*) FILE=RESTRAINT STRIDE=100

\endplumedfile

This third example show how to use chemical shifts to calculate a \ref METAINFERENCE score .

\plumedfile
cs: CS2BACKBONE ATOMS=1-174 DATADIR=data/ NRES=13 DOSCORE NDATA=24
csbias: BIASVALUE ARG=cs.score

PRINT ARG=(cs\.hn_.*),(cs\.nh_.*) FILE=CS.dat STRIDE=1000
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
  inline double * CONST_BB2_PREV(const unsigned a_kind, const unsigned at_kind) {return co_bb[a_kind][at_kind];}
  inline double * CONST_BB2_CURR(const unsigned a_kind, const unsigned at_kind) {return co_bb[a_kind][at_kind]+5;}
  inline double * CONST_BB2_NEXT(const unsigned a_kind, const unsigned at_kind) {return co_bb[a_kind][at_kind]+11;}
  inline double * CONST_SC2(const unsigned a_kind, const unsigned at_kind, unsigned res_type) { return co_sc_[a_kind][at_kind][res_type];}
  inline double * CONST_XD(const unsigned a_kind, const unsigned at_kind) { return co_xd[a_kind][at_kind];}
  inline double * CO_SPHERE(const unsigned a_kind, const unsigned at_kind, unsigned exp_type) { return co_sphere[a_kind][at_kind][exp_type];}
  inline double * CO_RING(const unsigned a_kind, const unsigned at_kind) { return co_ring[a_kind][at_kind];}
  inline double * CO_DA(const unsigned a_kind, const unsigned at_kind) { return co_da[a_kind][at_kind];}
  inline double * PARS_DA(const unsigned a_kind, const unsigned at_kind, const unsigned ang_kind) { return pars_da[a_kind][at_kind][ang_kind];}

  void parse(const string &file, const double dscale) {
    ifstream in;
    in.open(file.c_str());
    if(!in) plumed_merror("Unable to open CS2Backbone DB file " +file);

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
        string str_err = "CS2Backbone DB WARNING: unrecognized token: " + tok[0];
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
  struct Fragment {
    vector<Value*> comp;
    vector<double> exp_cs;
    unsigned res_type_prev;
    unsigned res_type_curr;
    unsigned res_type_next;
    unsigned res_kind;
    unsigned fd;
    string res_name;
    vector<int> pos;
    vector<int> prev;
    vector<int> curr;
    vector<int> next;
    vector<int> side_chain;
    vector<int> xd1;
    vector<int> xd2;
    vector<unsigned> box_nb;
    vector<int> phi;
    vector<int> psi;
    vector<int> chi1;
    double t_phi;
    double t_psi;
    double t_chi1;
    vector<Vector> dd0, dd10, dd21, dd2;

    Fragment() {
      comp.resize(6);
      exp_cs.resize(6,0);
      res_type_prev = res_type_curr = res_type_next = 0;
      res_kind = 0;
      fd = 0;
      res_name = "";
      pos.resize(6,-1);
      prev.reserve(5);
      curr.reserve(6);
      next.reserve(5);
      side_chain.reserve(20);
      xd1.reserve(27);
      xd2.reserve(27);
      box_nb.reserve(250);
      phi.reserve(4);
      psi.reserve(4);
      chi1.reserve(4);
      t_phi = t_psi = t_chi1 = 0;
      dd0.resize(3);
      dd10.resize(3);
      dd21.resize(3);
      dd2.resize(3);
    }

  };

  struct RingInfo {
    enum {R_PHE, R_TYR, R_TRP1, R_TRP2, R_HIS};
    unsigned rtype;    // one out of five different types
    unsigned atom[6];  // up to six member per ring
    unsigned numAtoms; // number of ring members (5 or 6)
    Vector position;   // center of ring coordinates
    Vector normVect;   // ring plane normal vector
    Vector n1, n2;     // two atom plane normal vectors used to compute ring plane normal
    Vector g[6];       // vector of the vectors used to calculate n1,n2
    double lengthN2;   // square of length of normVect
    double lengthNV;   // length of normVect
    RingInfo():
      rtype(0),numAtoms(0),
      lengthN2(NAN),lengthNV(NAN)
    {for(unsigned i=0; i<6; i++) atom[i]=0;}
  };

  enum aa_t {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK};
  enum atom_t {D_C, D_H, D_N, D_O, D_S, D_C2, D_N2, D_O2};

  CS2BackboneDB    db;
  vector<vector<Fragment> > atom;
  vector<vector<vector<unsigned> > > index_cs;
  vector<RingInfo> ringInfo;
  vector<unsigned> seg_last;
  vector<unsigned> type;
  vector<unsigned> res_num;
  unsigned         box_nupdate;
  unsigned         box_count;
  bool             camshift;
  bool             pbc;

  void remove_problematic(const string &res, const string &nucl);
  void read_cs(const string &file, const string &k);
  void update_neighb();
  void compute_ring_parameters();
  void compute_dihedrals();
  void init_backbone(const PDB &pdb);
  void init_sidechain(const PDB &pdb);
  void init_xdist(const PDB &pdb);
  void init_types(const PDB &pdb);
  void init_rings(const PDB &pdb);
  aa_t frag2enum(const string &aa);
  vector<string> side_chain_atoms(const string &s);
  bool isSP2(const string & resType, const string & atomName);
  bool is_chi1_cx(const string & frg, const string & atm);
  unsigned frag_segment(const unsigned p);
  unsigned frag_relitive_index(const unsigned p, const unsigned s);
  void debug_report();
  void xdist_name_map(string & name);

public:

  explicit CS2Backbone(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate();
  void update();
};

PLUMED_REGISTER_ACTION(CS2Backbone,"CS2BACKBONE")

void CS2Backbone::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  MetainferenceBase::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","DATADIR","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file of the protein system to initialise ALMOST.");
  keys.add("compulsory","NEIGH_FREQ","20","Period in step for neighbour list update.");
  keys.add("compulsory","NRES","Number of residues, corresponding to the number of chemical shifts.");
  keys.addFlag("CAMSHIFT",false,"Set to TRUE if you to calculate a single CamShift score.");
  keys.addFlag("NOEXP",false,"Set to TRUE if you don't want to have fixed components with the experimetnal values.");
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
  camshift(false),
  pbc(true)
{
  string stringadb;
  string stringapdb;

  parseFlag("CAMSHIFT",camshift);
  if(camshift&&getDoScore()) error("It is not possible to use CAMSHIFT together with DOSCORE");

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  bool noexp=false;
  parseFlag("NOEXP",noexp);

  string stringa_data;
  parse("DATADIR",stringa_data);

  string stringa_template;
  parse("TEMPLATE",stringa_template);

  box_count=0;
  box_nupdate=20;
  parse("NEIGH_FREQ", box_nupdate);

  unsigned numResidues;
  parse("NRES", numResidues);

  stringadb  = stringa_data + string("/camshift.db");
  stringapdb = stringa_data + string("/") + stringa_template;

  /* Lenght conversion (parameters are tuned for angstrom) */
  double scale=1.;
  if(!plumed.getAtoms().usingNaturalUnits()) {
    scale = 10.*atoms.getUnits().getLength();
  }

  log.printf("  Initialization of the predictor ...\n"); log.flush();
  db.parse(stringadb,scale);
  PDB pdb;
  if( !pdb.read(stringapdb,plumed.getAtoms().usingNaturalUnits(),1./scale) ) error("missing input file " + stringapdb);
  init_backbone(pdb);
  init_sidechain(pdb);
  init_xdist(pdb);
  init_types(pdb);
  init_rings(pdb);
#ifndef NDEBUG
  debug_report();
#endif
  log.printf("  Reading experimental data ...\n"); log.flush();
  stringadb = stringa_data + string("/CAshifts.dat");
  log.printf("  Initializing CA shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "CA");
  stringadb = stringa_data + string("/CBshifts.dat");
  log.printf("  Initializing CB shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "CB");
  stringadb = stringa_data + string("/Cshifts.dat");
  log.printf("  Initializing C' shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "C");
  stringadb = stringa_data + string("/HAshifts.dat");
  log.printf("  Initializing HA shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "HA");
  stringadb = stringa_data + string("/Hshifts.dat");
  log.printf("  Initializing H shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "H");
  stringadb = stringa_data + string("/Nshifts.dat");
  log.printf("  Initializing N shifts %s\n", stringadb.c_str()); log.flush();
  read_cs(stringadb, "N");

  /* this is a workaround for those chemical shifts that can result in too large forces */
  remove_problematic("GLN", "CB");
  remove_problematic("ILE", "CB");
  remove_problematic("PRO", "N");
  remove_problematic("PRO", "H");
  remove_problematic("PRO", "CB");
  remove_problematic("GLY", "HA");
  remove_problematic("GLY", "CB");
  /* this is a workaround for those chemical shifts that are not parameterized */
  remove_problematic("HIE", "HA"); remove_problematic("HIP", "HA"); remove_problematic("HSP", "HA");
  remove_problematic("HIE", "H");  remove_problematic("HIP", "H");  remove_problematic("HSP", "H");
  remove_problematic("HIE", "N");  remove_problematic("HIP", "N");  remove_problematic("HSP", "N");
  remove_problematic("HIE", "CA"); remove_problematic("HIP", "CA"); remove_problematic("HSP", "CA");
  remove_problematic("HIE", "CB"); remove_problematic("HIP", "CB"); remove_problematic("HSP", "CB");
  remove_problematic("HIE", "C");  remove_problematic("HIP", "C");  remove_problematic("HSP", "C");
  remove_problematic("GLH", "HA"); remove_problematic("ASH", "HA"); remove_problematic("HSE", "HA");
  remove_problematic("GLH", "H");  remove_problematic("ASH", "H");  remove_problematic("HSE", "H");
  remove_problematic("GLH", "N");  remove_problematic("ASH", "N");  remove_problematic("HSE", "N");
  remove_problematic("GLH", "CA"); remove_problematic("ASH", "CA"); remove_problematic("HSE", "CA");
  remove_problematic("GLH", "CB"); remove_problematic("ASH", "CB"); remove_problematic("HSE", "CB");
  remove_problematic("GLH", "C");  remove_problematic("ASH", "C");  remove_problematic("HSE", "C");
  remove_problematic("UNK", "HA");
  remove_problematic("UNK", "H");
  remove_problematic("UNK", "N");
  remove_problematic("UNK", "CA");
  remove_problematic("UNK", "CB");
  remove_problematic("UNK", "C");
  remove_problematic("CYS", "HA");
  remove_problematic("CYS", "H");
  remove_problematic("CYS", "N");
  remove_problematic("CYS", "CA");
  remove_problematic("CYS", "CB");
  remove_problematic("CYS", "C");
  remove_problematic("HIS", "HA");
  remove_problematic("HIS", "H");
  remove_problematic("HIS", "N");
  remove_problematic("HIS", "CA");
  remove_problematic("HIS", "CB");
  remove_problematic("HIS", "C");
  /* done */

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);

  log<<"  Bibliography "
     <<plumed.cite("Kohlhoff K, Robustelli P, Cavalli A, Salvatella A, Vendruscolo M, J. Am. Chem. Soc. 131, 13894 (2009)");
  if(camshift) log<<plumed.cite("Granata D, Camilloni C, Vendruscolo M, Laio A, Proc. Natl. Acad. Sci. USA 110, 6817 (2013)");
  else log<<plumed.cite("Camilloni C, Robustelli P, De Simone A, Cavalli A, Vendruscolo M, J. Am. Chem. Soc. 134, 3968 (2012)");
  log<<plumed.cite("Bonomi M, Camilloni C, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  const string str_cs[] = {"ha_","hn_","nh_","ca_","cb_","co_"};
  unsigned index=0;
  if(camshift) {
    noexp = true;
    addValueWithDerivatives();
    setNotPeriodic();
  } else {
    if(getDoScore()) {
      index_cs.resize(atom.size(), vector<vector<unsigned> >());
      for(unsigned i=0; i<atom.size(); i++) {
        index_cs[i].resize(atom[i].size(), vector<unsigned>(6));
      }
    }
    unsigned l=0;
    for(unsigned i=0; i<atom.size(); i++) {
      for(unsigned a=0; a<atom[i].size(); a++) {
        unsigned res=index+a;
        std::string num; Tools::convert(res,num);
        for(unsigned at_kind=0; at_kind<6; at_kind++) {
          if(atom[i][a].exp_cs[at_kind]!=0) {
            if(getDoScore()) {
              addComponent(str_cs[at_kind]+num);
              componentIsNotPeriodic(str_cs[at_kind]+num);
              atom[i][a].comp[at_kind] = getPntrToComponent(str_cs[at_kind]+num);
              index_cs[i][a][at_kind]=l;
              setParameter(atom[i][a].exp_cs[at_kind]);
              l++;
            } else {
              addComponentWithDerivatives(str_cs[at_kind]+num);
              componentIsNotPeriodic(str_cs[at_kind]+num);
              atom[i][a].comp[at_kind] = getPntrToComponent(str_cs[at_kind]+num);
            }
          }
        }
      }
      index += atom[i].size();
    }
    if(getDoScore()) Initialise(l);
  }

  if(!noexp) {
    index = 0;
    for(unsigned i=0; i<atom.size(); i++) {
      for(unsigned a=0; a<atom[i].size(); a++) {
        unsigned res=index+a;
        std::string num; Tools::convert(res,num);
        for(unsigned at_kind=0; at_kind<6; at_kind++) {
          if(atom[i][a].exp_cs[at_kind]!=0) {
            addComponent("exp"+str_cs[at_kind]+num);
            componentIsNotPeriodic("exp"+str_cs[at_kind]+num);
            Value* comp=getPntrToComponent("exp"+str_cs[at_kind]+num);
            comp->set(atom[i][a].exp_cs[at_kind]);
          }
        }
      }
      index += atom[i].size();
    }
  }

  /* temporary check, the idea is that I can remove NRES completely */
  index=0;
  for(unsigned i=0; i<atom.size(); i++) {
    index += atom[i].size();
  }
  if(index!=numResidues) error("NRES and the number of residues in the PDB do not match!");

  requestAtoms(atoms);
  setDerivatives();
  checkRead();
}

void CS2Backbone::remove_problematic(const string &res, const string &nucl) {
  unsigned n;
  if(nucl=="HA")     n=0;
  else if(nucl=="H") n=1;
  else if(nucl=="N") n=2;
  else if(nucl=="CA")n=3;
  else if(nucl=="CB")n=4;
  else if(nucl=="C") n=5;
  else return;

  for(unsigned i=0; i<atom.size(); i++) {
    for(unsigned a=0; a<atom[i].size(); a++) {
      if(atom[i][a].res_name.c_str()==res) {
        atom[i][a].exp_cs[n] = 0;
      }
    }
  }
}

void CS2Backbone::read_cs(const string &file, const string &nucl) {
  ifstream in;
  in.open(file.c_str());
  if(!in) error("CS2Backbone: Unable to open " + file);
  istream_iterator<string> iter(in), end;

  unsigned n;
  if(nucl=="HA")     n=0;
  else if(nucl=="H") n=1;
  else if(nucl=="N") n=2;
  else if(nucl=="CA")n=3;
  else if(nucl=="CB")n=4;
  else if(nucl=="C") n=5;
  else return;

  int oldseg = -1;
  int oldp = -1;
  while(iter!=end) {
    string tok;
    tok = *iter; ++iter;
    if(tok[0]=='#') { ++iter; continue;}
    int p = atoi(tok.c_str());
    p = p - 1;
    const int seg = frag_segment(p);
    p = frag_relitive_index(p,seg);
    if(oldp==-1) oldp=p;
    if(oldseg==-1) oldseg=seg;
    if(p<oldp&&seg==oldseg) {
      string errmsg = "Error while reading " + file + "! The same residue number has been used twice!";
      error(errmsg);
    }
    tok = *iter; ++iter;
    double cs = atof(tok.c_str());
    if(atom[seg][p].pos[n]<=0) cs=0;
    else atom[seg][p].exp_cs[n] = cs;
    oldseg = seg;
    oldp = p;
  }
  in.close();
}

void CS2Backbone::calculate()
{
  if(pbc) makeWhole();
  if(getExchangeStep()) box_count=0;
  if(box_count==0) update_neighb();

  compute_ring_parameters();
  compute_dihedrals();

  double score = 0.;

  vector<double> camshift_sigma2(6);
  camshift_sigma2[0] = 0.08; // HA
  camshift_sigma2[1] = 0.30; // HN
  camshift_sigma2[2] = 9.00; // NH
  camshift_sigma2[3] = 1.30; // CA
  camshift_sigma2[4] = 1.56; // CB
  camshift_sigma2[5] = 1.70; // CO

  const unsigned chainsize = atom.size();
  const unsigned atleastned = 72+ringInfo.size()*6;

  vector<vector<unsigned> > all_list;
  vector<vector<Vector> > all_ff;
  if(getDoScore()) {
    all_list.resize(getNarg(), vector<unsigned>());
    all_ff.resize(getNarg(), vector<Vector>());
  }

  // CYCLE OVER MULTIPLE CHAINS
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  for(unsigned s=0; s<chainsize; s++) {
    const unsigned psize = atom[s].size();
    vector<Vector> omp_deriv;
    if(camshift) omp_deriv.resize(getNumberOfAtoms(), Vector(0,0,0));
    #pragma omp for reduction(+:score)
    // SKIP FIRST AND LAST RESIDUE OF EACH CHAIN
    for(unsigned a=1; a<psize-1; a++) {

      const Fragment *myfrag = &atom[s][a];
      const unsigned aa_kind = myfrag->res_kind;
      const unsigned needed_atoms = atleastned+myfrag->box_nb.size();

      /* Extra Distances are the same for each residue */
      const unsigned xdsize=myfrag->xd1.size();
      vector<Vector> ext_distances(xdsize);
      vector<double> ext_d(xdsize);
      for(unsigned q=0; q<xdsize; q++) {
        if(myfrag->xd1[q]==-1||myfrag->xd2[q]==-1) continue;
        const Vector distance = delta(getPosition(myfrag->xd1[q]),getPosition(myfrag->xd2[q]));
        ext_d[q] = distance.modulo();
        ext_distances[q] = distance/ext_d[q];
      }

      // CYCLE OVER THE SIX BACKBONE CHEMICAL SHIFTS
      for(unsigned at_kind=0; at_kind<6; at_kind++) {
        if(atom[s][a].exp_cs[at_kind]!=0) {
          // Common constant and AATYPE
          const double * CONSTAACURR = db.CONSTAACURR(aa_kind,at_kind);
          const double * CONSTAANEXT = db.CONSTAANEXT(aa_kind,at_kind);
          const double * CONSTAAPREV = db.CONSTAAPREV(aa_kind,at_kind);

          double cs = CONSTAACURR[myfrag->res_type_curr] +
                      CONSTAANEXT[myfrag->res_type_next] +
                      CONSTAAPREV[myfrag->res_type_prev];
          // this is the atom for which we are calculating the chemical shift
          const unsigned ipos = myfrag->pos[at_kind];

          vector<unsigned> list;
          list.reserve(needed_atoms);
          list.push_back(ipos);
          vector<Vector> ff;
          ff.reserve(needed_atoms);
          ff.push_back(Vector(0,0,0));


          //PREV
          const double * CONST_BB2_PREV = db.CONST_BB2_PREV(aa_kind,at_kind);
          const unsigned presize = myfrag->prev.size();
          for(unsigned q=0; q<presize; q++) {
            const double cb2pq = CONST_BB2_PREV[q];
            if(cb2pq==0.) continue;
            const unsigned jpos = myfrag->prev[q];
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = cb2pq/d;

            cs += cb2pq*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //CURR
          const double * CONST_BB2_CURR = db.CONST_BB2_CURR(aa_kind,at_kind);
          const unsigned cursize = myfrag->curr.size();
          for(unsigned q=0; q<cursize; q++) {
            const double cb2cq = CONST_BB2_CURR[q];
            if(cb2cq==0.) continue;
            const unsigned jpos = myfrag->curr[q];
            if(ipos==jpos) continue;
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = cb2cq/d;

            cs += cb2cq*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //NEXT
          const double * CONST_BB2_NEXT = db.CONST_BB2_NEXT(aa_kind,at_kind);
          const unsigned nexsize = myfrag->next.size();
          for(unsigned q=0; q<nexsize; q++) {
            const double cb2nq = CONST_BB2_NEXT[q];
            if(cb2nq==0.) continue;
            const unsigned jpos = myfrag->next[q];
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = cb2nq/d;

            cs += cb2nq*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //SIDE CHAIN
          const double * CONST_SC2 = db.CONST_SC2(aa_kind,at_kind,myfrag->res_type_curr);
          const unsigned sidsize = myfrag->side_chain.size();
          for(unsigned q=0; q<sidsize; q++) {
            const double cs2q = CONST_SC2[q];
            if(cs2q==0.) continue;
            const unsigned jpos = myfrag->side_chain[q];
            if(ipos==jpos) continue;
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = cs2q/d;

            cs += cs2q*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //EXTRA DIST
          const double * CONST_XD  = db.CONST_XD(aa_kind,at_kind);
          for(unsigned q=0; q<xdsize; q++) {
            const double cxdq = CONST_XD[q];
            if(cxdq==0.) continue;
            if(myfrag->xd1[q]==-1||myfrag->xd2[q]==-1) continue;
            list.push_back(myfrag->xd2[q]);
            list.push_back(myfrag->xd1[q]);
            cs += cxdq*ext_d[q];
            const Vector der = cxdq*ext_distances[q];
            ff.push_back( der);
            ff.push_back(-der);
          }

          //NON BOND
          {
            const double * CONST_CO_SPHERE3 = db.CO_SPHERE(aa_kind,at_kind,0);
            const double * CONST_CO_SPHERE  = db.CO_SPHERE(aa_kind,at_kind,1);
            const unsigned boxsize = myfrag->box_nb.size();
            for(unsigned bat=0; bat<boxsize; bat++) {
              const unsigned jpos = myfrag->box_nb[bat];
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
                cs += factor1*CONST_CO_SPHERE[t] + factor3*CONST_CO_SPHERE3[t] ;
                const double fact = dfactor1*CONST_CO_SPHERE[t]+dfactor3*CONST_CO_SPHERE3[t];
                const Vector der  = fact*distance;

                list.push_back(jpos);
                ff[0] += der;
                ff.push_back(-der);
              }
            }
          }
          //END NON BOND

          //RINGS
          {
            const double *rc = db.CO_RING(aa_kind,at_kind);
            const unsigned rsize = ringInfo.size();
            for(unsigned i=0; i<rsize; i++) {
              // compute angle from ring middle point to current atom position
              // get distance vector from query atom to ring center and normal vector to ring plane
              const Vector n   = ringInfo[i].normVect;
              const double nL  = ringInfo[i].lengthNV;
              const double inL2 = ringInfo[i].lengthN2;

              const Vector d = delta(ringInfo[i].position, getPosition(ipos));
              const double dL2 = d.modulo2();
              double dL  = sqrt(dL2);
              const double idL3 = 1./(dL2*dL);

              const double dn    = dotProduct(d,n);
              const double dn2   = dn*dn;
              const double dLnL  = dL*nL;
              const double dL_nL = dL/nL;

              const double ang2 = dn2*inL2/dL2;
              const double u    = 1.-3.*ang2;
              const double cc   = rc[ringInfo[i].rtype];

              cs += cc*u*idL3;

              const double fUU    = -6*dn*inL2;
              const double fUQ    = fUU/dL;
              const Vector gradUQ = fUQ*(dL2*n - dn*d);
              const Vector gradVQ = (3*dL*u)*d;

              const double fact   = cc*idL3*idL3;
              ff[0] += fact*(gradUQ - gradVQ);

              const double fU       = fUU/nL;
              double OneOverN = 1./6.;
              if(ringInfo[i].numAtoms==5) OneOverN=1./3.;
              const Vector factor2  = OneOverN*n;
              const Vector factor4  = (OneOverN/dL_nL)*d;

              const Vector gradV    = -OneOverN*gradVQ;

              if(ringInfo[i].numAtoms==6) {
                // update forces on ring atoms
                for(unsigned at=0; at<6; at++) {
                  const Vector ab = crossProduct(d,ringInfo[i].g[at]);
                  const Vector c  = crossProduct(n,ringInfo[i].g[at]);
                  const Vector factor3 = 0.5*dL_nL*c;
                  const Vector factor1 = 0.5*ab;
                  const Vector gradU   = fU*( dLnL*(factor1 - factor2) -dn*(factor3 - factor4) );
                  ff.push_back(fact*(gradU - gradV));
                  list.push_back(ringInfo[i].atom[at]);
                }
              }  else {
                for(unsigned at=0; at<3; at++) {
                  const Vector ab = crossProduct(d,ringInfo[i].g[at]);
                  const Vector c  = crossProduct(n,ringInfo[i].g[at]);
                  const Vector factor3 = dL_nL*c;
                  const Vector factor1 = ab;
                  const Vector gradU   = fU*( dLnL*(factor1 - factor2) -dn*(factor3 - factor4) );
                  ff.push_back(fact*(gradU - gradV));
                }
                list.push_back(ringInfo[i].atom[0]);
                list.push_back(ringInfo[i].atom[2]);
                list.push_back(ringInfo[i].atom[3]);
              }
            }
          }
          //END OF RINGS

          //DIHEDRAL ANGLES
          {
            const double *CO_DA = db.CO_DA(aa_kind,at_kind);
            if(myfrag->phi.size()==4) {
              const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,0);
              const double val1 = 3.*myfrag->t_phi+PARS_DA[3];
              const double val2 = myfrag->t_phi+PARS_DA[4];
              cs += CO_DA[0]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = -CO_DA[0]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(fact*myfrag->dd0[0]);
              ff.push_back(fact*myfrag->dd10[0]);
              ff.push_back(fact*myfrag->dd21[0]);
              ff.push_back(-fact*myfrag->dd2[0]);
              list.push_back(myfrag->phi[0]);
              list.push_back(myfrag->phi[1]);
              list.push_back(myfrag->phi[2]);
              list.push_back(myfrag->phi[3]);
            }

            if(myfrag->psi.size()==4) {
              const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,1);
              const double val1 = 3.*myfrag->t_psi+PARS_DA[3];
              const double val2 = myfrag->t_psi+PARS_DA[4];
              cs += CO_DA[1]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = -CO_DA[1]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(fact*myfrag->dd0[1]);
              ff.push_back(fact*myfrag->dd10[1]);
              ff.push_back(fact*myfrag->dd21[1]);
              ff.push_back(-fact*myfrag->dd2[1]);
              list.push_back(myfrag->psi[0]);
              list.push_back(myfrag->psi[1]);
              list.push_back(myfrag->psi[2]);
              list.push_back(myfrag->psi[3]);
            }

            //Chi
            if(myfrag->chi1.size()==4) {
              const double *PARS_DA = db.PARS_DA(aa_kind,at_kind,2);
              const double val1 = 3.*myfrag->t_chi1+PARS_DA[3];
              const double val2 = myfrag->t_chi1+PARS_DA[4];
              cs += CO_DA[2]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = -CO_DA[2]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(fact*myfrag->dd0[2]);
              ff.push_back(fact*myfrag->dd10[2]);
              ff.push_back(fact*myfrag->dd21[2]);
              ff.push_back(-fact*myfrag->dd2[2]);
              list.push_back(myfrag->chi1[0]);
              list.push_back(myfrag->chi1[1]);
              list.push_back(myfrag->chi1[2]);
              list.push_back(myfrag->chi1[3]);
            }
          }
          //END OF DIHE

          if(getDoScore()) {
            setCalcData(index_cs[s][a][at_kind], cs);
            all_list[index_cs[s][a][at_kind]] = list;
            all_ff[index_cs[s][a][at_kind]] = ff;
            Value *comp = atom[s][a].comp[at_kind];
            comp->set(cs);
          } else {
            if(camshift) {
              score += (cs - atom[s][a].exp_cs[at_kind])*(cs - atom[s][a].exp_cs[at_kind])/camshift_sigma2[at_kind];
              double fact = 2.0*(cs - atom[s][a].exp_cs[at_kind])/camshift_sigma2[at_kind];
              for(unsigned i=0; i<list.size(); i++) omp_deriv[list[i]] += fact*ff[i];
            } else {
              Value *comp = atom[s][a].comp[at_kind];
              comp->set(cs);
              Tensor virial;
              for(unsigned i=0; i<list.size(); i++) {
                setAtomsDerivatives(comp,list[i],ff[i]);
                virial-=Tensor(getPosition(list[i]),ff[i]);
              }
              setBoxDerivatives(comp,virial);
            }
          }
        }
      }
    }
    #pragma omp critical
    if(camshift) for(unsigned i=0; i<getPositions().size(); i++) setAtomsDerivatives(getPntrToValue(),i,omp_deriv[i]);
  }

  if(getDoScore()) {
    /* Metainference */
    double score = getScore();
    setScore(score);

    /* calculate final derivatives */
    Value* val=getPntrToComponent("score");

    Tensor virial;

    for(unsigned i=0; i<all_list.size(); i++) {
      const double fact = getMetaDer(i);
      for(unsigned j=0; j<all_list[i].size(); j++) {
        setAtomsDerivatives(val, all_list[i][j],  all_ff[i][j]*fact);
        virial-=Tensor(getPosition(all_list[i][j]), all_ff[i][j]*fact);
      }
    }
    setBoxDerivatives(val,virial);
  }

  // in the case of camshift we calculate the virial at the end
  if(camshift) {
    Tensor virial;
    unsigned nat=getNumberOfAtoms();
    Value* val=getPntrToValue();
    for(unsigned i=0; i<nat; i++) virial-=Tensor(getPosition(i),
                                            Vector(val->getDerivative(3*i+0),
                                                val->getDerivative(3*i+1),
                                                val->getDerivative(3*i+2)));
    setBoxDerivatives(val,virial);
    setValue(score);
  }

  ++box_count;
  if(box_count == box_nupdate) box_count = 0;
}

void CS2Backbone::update_neighb() {
  const unsigned chainsize = atom.size();
  for(unsigned s=0; s<chainsize; s++) {
    const unsigned psize = atom[s].size();
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for(unsigned a=1; a<psize-1; a++) {
      const unsigned boxsize = getNumberOfAtoms();
      atom[s][a].box_nb.clear();
      atom[s][a].box_nb.reserve(300);
      const unsigned res_curr = res_num[atom[s][a].pos[0]];
      for(unsigned bat=0; bat<boxsize; bat++) {
        const unsigned res_dist = abs(static_cast<int>(res_curr-res_num[bat]));
        if(res_dist<2) continue;
        for(unsigned at_kind=0; at_kind<6; at_kind++) {
          if(atom[s][a].exp_cs[at_kind]==0.) continue;
          const unsigned ipos = atom[s][a].pos[at_kind];
          const Vector distance = delta(getPosition(bat),getPosition(ipos));
          const double d2=distance.modulo2();
          if(d2<cutOffNB2) {
            atom[s][a].box_nb.push_back(bat);
            break;
          }
        }
      }
    }
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
      vector<Vector> a(size);
      a[0] = getPosition(ringInfo[i].atom[0]);
      // ring center
      Vector midP = a[0];
      for(unsigned j=1; j<size; j++) {
        a[j] = getPosition(ringInfo[i].atom[j]);
        midP += a[j];
      }
      midP /= (double) size;
      ringInfo[i].position = midP;
      // compute normal vector to plane containing first three atoms in array
      ringInfo[i].n1 = crossProduct(delta(a[0],a[4]), delta(a[0],a[2]));
      // compute normal vector to plane containing last three atoms in array
      // NB: third atom of five-membered ring used for both computations above
      ringInfo[i].n2 = crossProduct(delta(a[3],a[1]), delta(a[3],a[5]));
      // ring plane normal vector is average of n1 and n2
      ringInfo[i].normVect = 0.5*(ringInfo[i].n1 + ringInfo[i].n2);
    }  else {
      ringInfo[i].g[0] = delta(getPosition(ringInfo[i].atom[3]),getPosition(ringInfo[i].atom[2]));
      ringInfo[i].g[1] = delta(getPosition(ringInfo[i].atom[0]),getPosition(ringInfo[i].atom[3]));
      ringInfo[i].g[2] = delta(getPosition(ringInfo[i].atom[2]),getPosition(ringInfo[i].atom[0]));
      vector<Vector> a(size);
      for(unsigned j=0; j<size; j++) {
        a[j] = getPosition(ringInfo[i].atom[j]);
      }
      // ring center
      Vector midP = (a[0]+a[2]+a[3])/3.;
      ringInfo[i].position = midP;
      // compute normal vector to plane containing first three atoms in array
      ringInfo[i].n1 = crossProduct(delta(a[0],a[3]), delta(a[0],a[2]));
      // ring plane normal vector is average of n1 and n2
      ringInfo[i].normVect = ringInfo[i].n1;

    }
    // calculate squared length and length of normal vector
    ringInfo[i].lengthN2 = 1./ringInfo[i].normVect.modulo2();
    ringInfo[i].lengthNV = 1./sqrt(ringInfo[i].lengthN2);
  }
}

void CS2Backbone::compute_dihedrals() {
  const unsigned chainsize = atom.size();
  for(unsigned s=0; s<chainsize; s++) {
    const unsigned psize = atom[s].size();
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for(unsigned a=1; a<psize-1; a++) {
      const Fragment *myfrag = &atom[s][a];
      if(myfrag->phi.size()==4) {
        const Vector d0 = delta(getPosition(myfrag->phi[1]), getPosition(myfrag->phi[0]));
        const Vector d1 = delta(getPosition(myfrag->phi[2]), getPosition(myfrag->phi[1]));
        const Vector d2 = delta(getPosition(myfrag->phi[3]), getPosition(myfrag->phi[2]));
        Torsion t;
        Vector dd0, dd1, dd2;
        atom[s][a].t_phi = t.compute(d0,d1,d2,dd0,dd1,dd2);
        atom[s][a].dd0[0]  = dd0;
        atom[s][a].dd10[0] = dd1-dd0;
        atom[s][a].dd21[0] = dd2-dd1;
        atom[s][a].dd2[0]  = dd2;
      }
      if(myfrag->psi.size()==4) {
        const Vector d0 = delta(getPosition(myfrag->psi[1]), getPosition(myfrag->psi[0]));
        const Vector d1 = delta(getPosition(myfrag->psi[2]), getPosition(myfrag->psi[1]));
        const Vector d2 = delta(getPosition(myfrag->psi[3]), getPosition(myfrag->psi[2]));
        Torsion t;
        Vector dd0, dd1, dd2;
        atom[s][a].t_psi = t.compute(d0,d1,d2,dd0,dd1,dd2);
        atom[s][a].dd0[1]  = dd0;
        atom[s][a].dd10[1] = dd1-dd0;
        atom[s][a].dd21[1] = dd2-dd1;
        atom[s][a].dd2[1]  = dd2;
      }
      if(myfrag->chi1.size()==4) {
        const Vector d0 = delta(getPosition(myfrag->chi1[1]), getPosition(myfrag->chi1[0]));
        const Vector d1 = delta(getPosition(myfrag->chi1[2]), getPosition(myfrag->chi1[1]));
        const Vector d2 = delta(getPosition(myfrag->chi1[3]), getPosition(myfrag->chi1[2]));
        Torsion t;
        Vector dd0, dd1, dd2;
        atom[s][a].t_chi1 = t.compute(d0,d1,d2,dd0,dd1,dd2);
        atom[s][a].dd0[2]  = dd0;
        atom[s][a].dd10[2] = dd1-dd0;
        atom[s][a].dd21[2] = dd2-dd1;
        atom[s][a].dd2[2]  = dd2;
      }
    }
  }
}

void CS2Backbone::init_backbone(const PDB &pdb) {
  // number of chains
  vector<string> chains;
  pdb.getChainNames( chains );
  seg_last.resize(chains.size());
  unsigned old_size=0;

  for(unsigned i=0; i<chains.size(); i++) {
    unsigned start, end;
    string errmsg;
    pdb.getResidueRange( chains[i], start, end, errmsg );

    unsigned res_offset = start;
    unsigned resrange = end-start+1;

    if(i==0) seg_last[i] = resrange;
    else seg_last[i] = seg_last[i-1]+resrange;

    vector<int> N_;
    vector<int> H_;
    vector<int> CA_;
    vector<int> CB_;
    vector<int> HA_;
    vector<int> C_;
    vector<int> O_;
    vector<int> CX_;
    N_.resize (resrange,-1);
    H_.resize (resrange,-1);
    CA_.resize(resrange,-1);
    CB_.resize(resrange,-1);
    HA_.resize(resrange,-1);
    C_.resize (resrange,-1);
    O_.resize (resrange,-1);
    CX_.resize(resrange,-1);

    vector<AtomNumber> allatoms = pdb.getAtomsInChain(chains[i]);
    // cycle over all the atoms in the chain
    for(unsigned a=0; a<allatoms.size(); a++) {
      unsigned atm_index=a+old_size;
      unsigned f = pdb.getResidueNumber(allatoms[a]);
      unsigned f_idx = f-res_offset;
      string AN = pdb.getAtomName(allatoms[a]);
      string RES = pdb.getResidueName(allatoms[a]);
      if(AN=="N")                  N_ [f_idx] = atm_index;
      else if(AN=="H" ||AN=="HN" ) H_ [f_idx] = atm_index;
      else if(AN=="HA"||AN=="HA1") HA_[f_idx] = atm_index;
      else if(AN=="CA"           ) CA_[f_idx] = atm_index;
      else if(AN=="CB"           ) CB_[f_idx] = atm_index;
      else if(AN=="C"            ) C_ [f_idx] = atm_index;
      else if(AN=="O"            ) O_ [f_idx] = atm_index;
      else if(AN=="CD"&&RES=="PRO") H_[f_idx] = atm_index;
      if(is_chi1_cx(RES,AN))       CX_[f_idx] = atm_index;
    }
    old_size+=allatoms.size();

    // vector of residues for a given chain
    vector<Fragment> atm_;
    // cycle over all residues in the chain
    for(unsigned a=start; a<=end; a++) {
      unsigned f_idx = a - res_offset;
      Fragment at;
      at.pos[0] = HA_[f_idx];
      at.pos[1] =  H_[f_idx];
      at.pos[2] =  N_[f_idx];
      at.pos[3] = CA_[f_idx];
      at.pos[4] = CB_[f_idx];
      at.pos[5] =  C_[f_idx];
      at.res_type_prev = at.res_type_curr = at.res_type_next = 0;
      at.res_name = pdb.getResidueName(a, chains[i]);
      at.res_kind = db.kind(at.res_name);
      at.fd = a;
      //REGISTER PREV CURR NEXT
      {
        if(a>start) {
          at.prev.push_back( N_[f_idx-1]);
          at.prev.push_back(CA_[f_idx-1]);
          at.prev.push_back(HA_[f_idx-1]);
          at.prev.push_back( C_[f_idx-1]);
          at.prev.push_back( O_[f_idx-1]);
          at.res_type_prev = frag2enum(pdb.getResidueName(a-1, chains[i]));
        }

        at.curr.push_back( N_[f_idx]);
        at.curr.push_back( H_[f_idx]);
        at.curr.push_back(CA_[f_idx]);
        at.curr.push_back(HA_[f_idx]);
        at.curr.push_back( C_[f_idx]);
        at.curr.push_back( O_[f_idx]);
        at.res_type_curr = frag2enum(pdb.getResidueName(a, chains[i]));

        if(a<end) {
          at.next.push_back (N_[f_idx+1]);
          at.next.push_back (H_[f_idx+1]);
          at.next.push_back(CA_[f_idx+1]);
          at.next.push_back(HA_[f_idx+1]);
          at.next.push_back( C_[f_idx+1]);
          at.res_type_next = frag2enum(pdb.getResidueName(a+1, chains[i]));
        }

        //PHI | PSI | CH1
        if(a>start) {
          at.phi.push_back( C_[f_idx-1]);
          at.phi.push_back( N_[f_idx]);
          at.phi.push_back(CA_[f_idx]);
          at.phi.push_back( C_[f_idx]);
        }

        if(a<end) {
          at.psi.push_back( N_[f_idx]);
          at.psi.push_back(CA_[f_idx]);
          at.psi.push_back( C_[f_idx]);
          at.psi.push_back( N_[f_idx+1]);
        }

        if(CX_[f_idx]!=-1&&CB_[f_idx]!=-1) {
          at.chi1.push_back( N_[f_idx]);
          at.chi1.push_back(CA_[f_idx]);
          at.chi1.push_back(CB_[f_idx]);
          at.chi1.push_back(CX_[f_idx]);
        }
      }
      atm_.push_back(at);
    }
    atom.push_back(atm_);
  }
}

void CS2Backbone::init_sidechain(const PDB &pdb) {
  vector<string> chains;
  pdb.getChainNames( chains );
  unsigned old_size=0;
  // cycle over chains
  for(unsigned s=0; s<atom.size(); s++) {
    AtomNumber astart, aend;
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();
    // cycle over residues
    for(unsigned a=0; a<atom[s].size(); a++) {
      if(atom[s][a].res_name=="UNK") continue;
      vector<AtomNumber> atm = pdb.getAtomsInResidue(atom[s][a].fd, chains[s]);
      vector<string> sc_atm = side_chain_atoms(atom[s][a].res_name);

      for(unsigned sc=0; sc<sc_atm.size(); sc++) {
        for(unsigned aa=0; aa<atm.size(); aa++) {
          if(pdb.getAtomName(atm[aa])==sc_atm[sc]) {
            atom[s][a].side_chain.push_back(atm[aa].index()-atom_offset+old_size);
          }
        }
      }

    }
    old_size += aend.index()+1;
  }
}

void CS2Backbone::init_xdist(const PDB &pdb) {
  const string atomsP1[] = {"H", "H", "H", "C", "C", "C",
                            "O", "O", "O", "N", "N", "N",
                            "O", "O", "O", "N", "N", "N",
                            "CG", "CG", "CG", "CG", "CG",
                            "CG", "CG", "CA"
                           };

  const int resOffsetP1 [] = {0, 0, 0, -1, -1, -1,
                              0, 0, 0, 1, 1, 1,
                              -1, -1, -1, 0, 0, 0,
                              0, 0, 0, 0, 0, -1, 1, -1
                             };

  const string atomsP2[] = {"HA", "C", "CB", "HA", "C", "CB",
                            "HA", "N", "CB", "HA", "N", "CB",
                            "HA", "N", "CB", "HA", "N", "CB",
                            "HA", "N", "C", "C", "N", "CA", "CA", "CA"
                           };

  const int resOffsetP2 [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 1, 0, 0, 1};

  vector<string> chains;
  pdb.getChainNames( chains );
  unsigned old_size=0;
  for(unsigned s=0; s<atom.size(); s++) {
    AtomNumber astart, aend;
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();

    for(unsigned a=1; a<atom[s].size()-1; a++) {
      vector<AtomNumber> atm_curr = pdb.getAtomsInResidue(atom[s][a].fd,chains[s]);
      vector<AtomNumber> atm_prev = pdb.getAtomsInResidue(atom[s][a].fd-1,chains[s]);
      vector<AtomNumber> atm_next = pdb.getAtomsInResidue(atom[s][a].fd+1,chains[s]);

      for(unsigned q=0; q<db.get_numXtraDists()-1; q++) {
        vector<AtomNumber>::iterator at1, at1_end;
        vector<AtomNumber>::iterator at2, at2_end;

        bool init_p1=false;
        AtomNumber p1;
        bool init_p2=false;
        AtomNumber p2;

        if(resOffsetP1[q]== 0) { at1 = atm_curr.begin(); at1_end = atm_curr.end();}
        if(resOffsetP1[q]==-1) { at1 = atm_prev.begin(); at1_end = atm_prev.end();}
        if(resOffsetP1[q]==+1) { at1 = atm_next.begin(); at1_end = atm_next.end();}
        while(at1!=at1_end) {
          AtomNumber aa = *at1;
          ++at1;
          string name = pdb.getAtomName(aa);

          xdist_name_map(name);

          if(name==atomsP1[q]) {
            p1 = aa;
            init_p1=true;
            break;
          }
        }

        if(resOffsetP2[q]== 0) { at2 = atm_curr.begin(); at2_end = atm_curr.end();}
        if(resOffsetP2[q]==-1) { at2 = atm_prev.begin(); at2_end = atm_prev.end();}
        if(resOffsetP2[q]==+1) { at2 = atm_next.begin(); at2_end = atm_next.end();}
        while(at2!=at2_end) {
          AtomNumber aa = *at2;
          ++at2;
          string name = pdb.getAtomName(aa);

          xdist_name_map(name);

          if(name==atomsP2[q]) {
            p2 = aa;
            init_p2=true;
            break;
          }
        }
        int add1 = -1;
        int add2 = -1;
        if(init_p1) add1=p1.index()-atom_offset+old_size;
        if(init_p2) add2=p2.index()-atom_offset+old_size;
        atom[s][a].xd1.push_back(add1);
        atom[s][a].xd2.push_back(add2);
      }
    }
    old_size += aend.index()+1;
  }
}

void CS2Backbone::init_types(const PDB &pdb) {
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
      else plumed_merror("CS2Backbone:init_type: unknown atom type!\n");
    } else {
      if (atom_type == 'C') t = D_C2;
      else if (atom_type == 'O') t = D_O2;
      else if (atom_type == 'N') t = D_N2;
      else plumed_merror("CS2Backbone:init_type: unknown atom type!\n");
    }
    type.push_back(t);
  }
}

void CS2Backbone::init_rings(const PDB &pdb) {

  const string pheTyr_n[] = {"CG","CD1","CE1","CZ","CE2","CD2"};
  const string trp1_n[]   = {"CD2","CE2","CZ2","CH2","CZ3","CE3"};
  const string trp2_n[]   = {"CG","CD1","NE1","CE2","CD2"};
  const string his_n[]    = {"CG","ND1","CD2","CE1","NE2"};

  vector<string> chains;
  pdb.getChainNames( chains );
  vector<AtomNumber> allatoms = pdb.getAtomNumbers();
  unsigned old_size=0;

  for(unsigned s=0; s<atom.size(); s++) {
    AtomNumber astart, aend;
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();
    for(unsigned r=0; r<atom[s].size(); r++) {
      string frg = pdb.getResidueName(atom[s][r].fd);
      if(!((frg=="PHE")||(frg=="TYR")||(frg=="TRP")||
           (frg=="HIS")||(frg=="HIP")||(frg=="HID")||
           (frg=="HIE")||(frg=="HSD")||(frg=="HSE")||
           (frg=="HSP"))) continue;

      vector<AtomNumber> frg_atoms = pdb.getAtomsInResidue(atom[s][r].fd,chains[s]);

      if(frg=="PHE"||frg=="TYR") {
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0; aa<6; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==pheTyr_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 6;
        if(frg=="PHE") ri.rtype = RingInfo::R_PHE;
        if(frg=="TYR") ri.rtype = RingInfo::R_TYR;
        ringInfo.push_back(ri);

      } else if(frg=="TRP") {
        //First ring
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0; aa<6; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==trp1_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 6;
        ri.rtype = RingInfo::R_TRP1;
        ringInfo.push_back(ri);
        //Second Ring
        RingInfo ri2;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0; aa<5; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==trp2_n[aa]) {
              ri2.atom[aa] = atm;
              break;
            }
          }
        }
        ri2.numAtoms = 5;
        ri2.rtype = RingInfo::R_TRP2;
        ringInfo.push_back(ri2);

      } else if((frg=="HIS")||(frg=="HIP")||(frg=="HID")||
                (frg=="HIE")||(frg=="HSD")||(frg=="HSE")||
                (frg=="HSP")) {//HIS case
        RingInfo ri;
        for(unsigned a=0; a<frg_atoms.size(); a++) {
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0; aa<5; aa++) {
            if(pdb.getAtomName(frg_atoms[a])==his_n[aa]) {
              ri.atom[aa] = atm;
              break;
            }
          }
        }
        ri.numAtoms = 5;
        ri.rtype = RingInfo::R_HIS;
        ringInfo.push_back(ri);
      } else {
        plumed_merror("Unkwown Ring Fragment");
      }
    }
    old_size += aend.index()+1;
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
  else plumed_merror("CS2Backbone: Error converting string " + aa + " into amino acid index: not a valid 3-letter code");
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
  } else plumed_merror("CS2Backbone: side_chain_atoms unknown");
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

unsigned CS2Backbone::frag_segment(const unsigned p) {
  unsigned s = 0;
  for(unsigned i=0; i<seg_last.size()-1; i++) {
    if(p>seg_last[i]) s  = i+1;
    else break;
  }
  return s;
}

unsigned CS2Backbone::frag_relitive_index(const unsigned p, const unsigned s) {
  if(s==0) return p;
  return p-seg_last[s-1];
}

void CS2Backbone::debug_report() {
  printf("\t CS2Backbone Initialization report: \n");
  printf("\t -------------------------------\n");
  printf("\t Number of segments: %u\n", static_cast<unsigned>(atom.size()));
  printf("\t Segments size:      ");
  for(unsigned i=0; i<atom.size(); i++) {printf("%u ", static_cast<unsigned>(atom[i].size()));} printf("\n");
  printf("\t%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n",
         "Seg","N","AA","Prev","Curr","Next","SC","XD1","XD2","Phi","Psi","Chi1");
  for(unsigned i=0; i<atom.size(); i++) {
    for(unsigned j=0; j<atom[i].size(); j++) {
      printf("\t%8u %8u %8s %8u %8u %8u %8u %8u %8u %8u %8u %8u \n",
             i+1, j+1,
             atom[i][j].res_name.c_str(),
             (unsigned)atom[i][j].prev.size(),
             (unsigned)atom[i][j].curr.size(),
             (unsigned)atom[i][j].next.size(),
             (unsigned)atom[i][j].side_chain.size(),
             (unsigned)atom[i][j].xd1.size(),
             (unsigned)atom[i][j].xd2.size(),
             (unsigned)atom[i][j].phi.size(),
             (unsigned)atom[i][j].psi.size(),
             (unsigned)atom[i][j].chi1.size());

      for(unsigned k=0; k<atom[i][j].prev.size(); k++) { printf("%8i ", atom[i][j].prev[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].curr.size(); k++) { printf("%8i ", atom[i][j].curr[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].next.size(); k++) { printf("%8i ", atom[i][j].next[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].side_chain.size(); k++) { printf("%8i ", atom[i][j].side_chain[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].xd1.size(); k++) { printf("%8i ", atom[i][j].xd1[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].xd2.size(); k++) { printf("%8i ", atom[i][j].xd2[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].phi.size(); k++) { printf("%8i ", atom[i][j].phi[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].psi.size(); k++) { printf("%8i ", atom[i][j].psi[k]);} printf("\n");
      for(unsigned k=0; k<atom[i][j].chi1.size(); k++) { printf("%8i ", atom[i][j].chi1[k]);} printf("\n");

    }
  }

  printf("\t Rings: \n");
  printf("\t ------ \n");
  printf("\t Number of rings: %u\n", static_cast<unsigned>(ringInfo.size()));
  printf("\t%8s %8s %8s %8s\n", "Num","Type","RType","N.atoms");
  for(unsigned i=0; i<ringInfo.size(); i++) {
    printf("\t%8u %8u %8u \n",i+1,ringInfo[i].rtype,ringInfo[i].numAtoms);
    for(unsigned j=0; j<ringInfo[i].numAtoms; j++) {printf("%8u ", ringInfo[i].atom[j]);} printf("\n");
  }
}

void CS2Backbone::xdist_name_map(string & name) {
  if((name == "OT1")||(name == "OC1")) name = "O";
  else if ((name == "HN") || (name == "HT1") || (name == "H1")) name = "H";
  else if ((name == "CG1")|| (name == "OG")||
           (name == "SG") || (name == "OG1")) name = "CG";
  else if ((name == "HA1"))                   name = "HA";
}

void CS2Backbone::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
