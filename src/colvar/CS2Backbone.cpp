/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

#define cutOffNB      0.43	// squared buffer distance for neighbour-lists 
#define cutOffDist    0.55  	// cut off distance for non-bonded pairwise forces
#define cutOffDist2   0.3025 	// square of cutOffDist
#define cutOnDist     0.45   	// cut off distance for non-bonded pairwise forces
#define cutOnDist2    0.2025 	// square of cutOffDist
#define invswitch     1.0/((cutOffDist2-cutOnDist2)*(cutOffDist2-cutOnDist2)*(cutOffDist2-cutOnDist2))

#include <string>
#include <fstream>
#include <iterator>
#include <sstream>

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/OpenMP.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "tools/Torsion.h"

using namespace std;

namespace PLMD {
namespace colvar{

//+PLUMEDOC COLVAR CS2BACKBONE 
/*
This collective variable calculates the backbone chemical shifts for a protein. 

The functional form is that of CamShift \cite Kohlhoff:2009us. The chemical shifts
of the selected nuclei/residues are saved as components. Reference experimental values
can also be stored as components. The two kind of components can then be used to calculate
either a scoring function as in \cite Robustelli:2010dn \cite Granata:2013dk or to calculate
ensemble averages as in \cite Camilloni:2012je \cite Camilloni:2013hs (see \ref STATS and
\ref ENSEMBLE).

CamShift calculation is relatively heavy because it often uses a large number of atoms, in order
to make it faster it is currently parallelised with \ref OpenMP.

In general the system for which chemical shifts are to be calculated must be completly included in
ATOMS and a TEMPLATE pdb file for the same atoms should be provided as well in the folder DATA. 
The atoms are made automatically whole unless NOPBC is used, in particular if the pdb is composed
by multiple chains it is usually better to use NOPBC and make the molecule whole \ref WHOLEMOLECULES
selecting an appropriate order.
 
In addition to a pdb file one needs to provide a list of chemical shifts to be calculated using one
file per nuclues type (CAshifts.dat, CBshifts.dat, Cshifts.dat, Hshifts.dat, HAshifts.dat, Nshifts.dat), 
all the six files should always be present. A chemical shifts for a nucleus is calculated if a value
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

The default behaviour is to store the values for the active nuclei are stored in components (ca_#, cb_#,
co_#, ha_#, hn_#, nh_# and expca_#, expcb_#, expco_#, expha_#, exphn_#, exp_nh#) with NOEXP it is possible
to only store the backcalculated values.
 
A pdb file is needed to the generate a simple topology of the protein. For histidines in protonation 
states different from D the HIE/HSE HIP/HSP name should be used. GLH and ASH can be used for the alternative 
protonation of GLU and ASP. Non-standard amino acids and other molecules are not yet supported, but in principle
they can be named UNK. If multiple chains are present the chain identifier must be in the standard PDB format, 
together with the TER keyword at the end of each chain. 

One more standard file is also needed in the folder DATA: camshift.db. This file includes all the CamShift parameters
and can be found in regtest/basic/rt45/data/ . 

All the above files must be in a single folder that must be specified with the keyword DATA. 

Additional material and examples can be also found in the tutorial \ref belfast-9 

\par Examples

In this first example the chemical shifts are used to calculate a scoring function to be used
in NMR driven Metadynamics \cite Granata:2013dk:

\verbatim
whole: GROUP ATOMS=2612-2514:-1,961-1:-1,2466-962:-1,2513-2467:-1
WHOLEMOLECULES ENTITY0=whole
cs: CS2BACKBONE ATOMS=1-2612 NRES=176 DATA=../data/ TEMPLATE=template.pdb NEIGH_FREQ=10
score: STATS ARG=(cs\.hn_.*),(cs\.nh_.*),(cs\.ca_.*),(cs\.cb_.*),(cs\.co_.*),(cs\.ha_.*) PARARG=(cs\.exphn_.*),(cs\.expnh_.*),(cs\.expca_.*),(cs\.expcb.*),(cs\.expco_.*),(cs\.expha_.*) SQDEVSUM  
metad: METAD ARG=score.sqdevsum ...
PRINT ARG=(cs\.hn_.*),(cs\.nh_.*),(cs\.ca_.*),(cs\.cb_.*),(cs\.co_.*),(cs\.ha_.*) FILE=CS.dat STRIDE=100 
PRINT ARG=score FILE=COLVAR STRIDE=100 
\endverbatim

In this second example the chemical shifts are used as replica-averaged restrained as in \cite Camilloni:2012je \cite Camilloni:2013hs. 
 
\verbatim
cs: CS2BACKBONE ATOMS=1-174 DATA=data/ NRES=13 
encs: ENSEMBLE ARG=(cs\.hn_.*),(cs\.nh_.*)
stcs: STATS ARG=encs.* SQDEVSUM PARARG=(cs\.exphn_.*),(cs\.expnh_.*)
RESTRAINT ARG=stcs.sqdevsum AT=0 KAPPA=0 SLOPE=24

PRINT ARG=(cs\.hn_.*),(cs\.nh_.*) FILE=RESTRAINT STRIDE=100

\endverbatim

(See also \ref WHOLEMOLECULES, \ref STATS, \ref METAD, \ref RESTRAINT and \ref PRINT)

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

  inline unsigned kind(const string &s){
    if(s=="GLY") return GLY;
    else if(s=="PRO") return PRO;      
    return STD;
  }

  inline unsigned atom_kind(const string &s){
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
  inline double * CONST_BB2_PREV(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind];}
  inline double * CONST_BB2_CURR(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind]+5;}
  inline double * CONST_BB2_NEXT(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind]+11;}
  inline double * CONST_SC2(const unsigned a_kind, const unsigned at_kind, unsigned res_type){ return co_sc_[a_kind][at_kind][res_type];}
  inline double * CONST_XD(const unsigned a_kind, const unsigned at_kind){ return co_xd[a_kind][at_kind];}
  inline double * CO_SPHERE(const unsigned a_kind, const unsigned at_kind, unsigned exp_type){ return co_sphere[a_kind][at_kind][exp_type];}
  inline double * CO_RING(const unsigned a_kind, const unsigned at_kind){ return co_ring[a_kind][at_kind];}
  inline double * CO_DA(const unsigned a_kind, const unsigned at_kind){ return co_da[a_kind][at_kind];}
  inline double * PARS_DA(const unsigned a_kind, const unsigned at_kind, const unsigned ang_kind){ return pars_da[a_kind][at_kind][ang_kind];}

  void parse(const string &file, const double dscale){
    ifstream in;
    in.open(file.c_str());
    if(!in) plumed_merror("Unable to open CS2Backbone DB file " +file);

    unsigned c_kind = 0;
    unsigned c_atom = 0;
    unsigned nline = 0;

    for(unsigned i=0;i<3;i++) for(unsigned j=0;j<6;j++) {
      for(unsigned k=0;k<20;k++) {
        c_aa[i][j][k]=0.;
        c_aa_prev[i][j][k]=0.;
        c_aa_succ[i][j][k]=0.;
        for(unsigned m=0;m<20;m++) co_sc_[i][j][k][m]=0.;
      }
      for(unsigned k=0;k<16;k++) {co_bb[i][j][k]=0.; }
      for(unsigned k=0;k<8;k++) { co_sphere[i][j][0][k]=0.; co_sphere[i][j][1][k]=0.; }
      for(unsigned k=0;k<3;k++) {
        co_da[i][j][k]=0.;
        for(unsigned l=0;l<5;l++) pars_da[i][j][k][l]=0.;
      }
      for(unsigned k=0;k<5;k++) co_ring[i][j][k]=0.;
      for(unsigned k=0;k<numXtraDists;k++) co_xd[i][j][k]=0.;
    }

    while(!in.eof()){
      string line;
      getline(in,line);
      ++nline;
      if(line.find("#")==0) continue;
      vector<string> tok;
      vector<string> tmp;
      tok = split(line,' ');
      for(unsigned q=0;q<tok.size();q++)
        if(tok[q].size()) tmp.push_back(tok[q]);
      tok = tmp;
      if(tok.size()==0) continue;
      if(tok[0]=="PAR"){
        c_kind = kind(tok[2]);
        c_atom = atom_kind(tok[1]);
        continue;
      }
      else if(tok[0]=="WEIGHT"){
        continue;
      }
      else if(tok[0]=="FLATBTM"){
        continue;
      }
      else if (tok[0] == "SCALEHARM"){
        continue;
      }
      else if (tok[0] == "TANHAMPLI"){
        continue;
      }
      else if (tok[0] == "ENDHARMON"){
        continue;
      }
      else if (tok[0] == "MAXRCDEVI"){
        continue;
      }
      else if (tok[0] == "RANDCOIL"){
        continue;
      }
      else if (tok[0] == "CONST"){
        continue;
      }
      else if (tok[0] == "CONSTAA"){
        assign(c_aa[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "CONSTAA-1"){
        assign(c_aa_prev[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "CONSTAA+1"){
        assign(c_aa_succ[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "COBB1"){
        continue;
      }
      else if (tok[0] == "COBB2"){
        //angstrom to nm
        assign(co_bb[c_kind][c_atom],tok,dscale);
        continue;
      }
      else if (tok[0] == "SPHERE1"){
        // angstrom^-3 to nm^-3
        assign(co_sphere[c_kind][c_atom][0],tok,1./(dscale*dscale*dscale));
        continue;
      }
      else if (tok[0] == "SPHERE2"){
        //angstrom to nm
        assign(co_sphere[c_kind][c_atom][1],tok,dscale);
        continue;
      }
      else if (tok[0] == "DIHEDRALS"){
        assign(co_da[c_kind][c_atom],tok,1);
        continue;
      }
      else if (tok[0] == "RINGS"){
        // angstrom^-3 to nm^-3
        assign(co_ring[c_kind][c_atom],tok,1./(dscale*dscale*dscale));
        for(unsigned i=1;i<tok.size();i++)
          co_ring[c_kind][c_atom][i-1] *= 1000;
        continue;
      }
      else if (tok[0] == "HBONDS"){
        continue;
      }
      else if (tok[0] == "XTRADISTS"){
        //angstrom to nm
        assign(co_xd[c_kind][c_atom],tok,dscale);
        continue;
      }
      else if(tok[0]=="DIHEDPHI"){
        assign(pars_da[c_kind][c_atom][0],tok,1);
        continue;
      }
      else if(tok[0]=="DIHEDPSI"){
        assign(pars_da[c_kind][c_atom][1],tok,1);
        continue;
      }
      else if(tok[0]=="DIHEDCHI1"){
        assign(pars_da[c_kind][c_atom][2],tok,1);
        continue;
      }

      bool ok = false;
      string scIdent1 [] = {"COSCALA1", "COSCARG1", "COSCASN1", "COSCASP1", "COSCCYS1", "COSCGLN1", "COSCGLU1", 
                            "COSCGLY1", "COSCHIS1", "COSCILE1", "COSCLEU1", "COSCLYS1", "COSCMET1", "COSCPHE1", 
                            "COSCPRO1", "COSCSER1", "COSCTHR1", "COSCTRP1", "COSCTYR1", "COSCVAL1"};

      for(unsigned scC = 0; scC < 20; scC++){
        if(tok[0]==scIdent1[scC]){
          ok = true; 
          break;
        }
      }
      if(ok) continue;

      string scIdent2 [] = {"COSCALA2", "COSCARG2", "COSCASN2", "COSCASP2", "COSCCYS2", "COSCGLN2", "COSCGLU2", 
                            "COSCGLY2", "COSCHIS2", "COSCILE2", "COSCLEU2", "COSCLYS2", "COSCMET2", "COSCPHE2", 
                            "COSCPRO2", "COSCSER2", "COSCTHR2", "COSCTRP2", "COSCTYR2", "COSCVAL2"};

      for(unsigned scC = 0; scC < 20; scC++){
        if(tok[0]==scIdent2[scC]){
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
  
  void assign(double * f, const vector<string> & v, const double scale){
    for(unsigned i=1;i<v.size();i++)
      f[i-1] = scale*(atof(v[i].c_str()));
  }
};

class CS2Backbone : public Colvar {
  struct Fragment {
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
    vector<vector<unsigned> > box_nb;
    vector<int> phi;
    vector<int> psi;
    vector<int> chi1;
    double t_phi;
    double t_psi;
    double t_chi1;
    vector<Vector> dd0, dd1, dd2;
  };

  struct RingInfo{
    enum {R_PHE, R_TYR, R_TRP1, R_TRP2, R_HIS};
    unsigned rtype;    // one out of five different types
    unsigned atom[6];  // up to six member per ring
    unsigned numAtoms; // number of ring members (5 or 6)
    Vector position;   // center of ring coordinates
    Vector normVect;   // ring plane normal vector
    Vector n1, n2;     // two atom plane normal vectors used to compute ring plane normal
    double lengthN2;   // square of length of normVect
    double lengthNV;   // length of normVect
  };

  enum aa_t {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK};
  enum atom_t {D_C, D_H, D_N, D_O, D_S, D_C2, D_N2, D_O2};

  CS2BackboneDB    db;
  vector<unsigned> seg_last;
  vector<unsigned> type;
  vector<unsigned> res_num;
  vector<RingInfo> ringInfo;
  unsigned         box_nupdate;
  unsigned         box_count;
  vector<vector<Fragment> > atom;
  unsigned         numResidues;
  unsigned         **c_sh;
  bool             pbc;

  void remove_problematic(const string &res, const string &nucl);
  void read_cs(const string &file, const string &k);
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
  ~CS2Backbone();
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(CS2Backbone,"CS2BACKBONE")

void CS2Backbone::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","DATA","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file of the protein system to initialise ALMOST.");
  keys.add("compulsory","NEIGH_FREQ","10","Period in step for neighbour list update.");
  keys.add("compulsory","NRES","Number of residues, corresponding to the number of chemical shifts.");
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
PLUMED_COLVAR_INIT(ao)
{
  string stringadb;
  string stringapdb;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  bool noexp=false;
  parseFlag("NOEXP",noexp);

  string stringa_data;
  parse("DATA",stringa_data);

  string stringa_template;
  parse("TEMPLATE",stringa_template);

  box_count=0;
  box_nupdate=10;
  parse("NEIGH_FREQ", box_nupdate);

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
  checkRead();

  c_sh = new unsigned*[numResidues];
  c_sh[0] = new unsigned[numResidues*6];
  for(unsigned i=1;i<numResidues;i++) c_sh[i]=c_sh[i-1]+6;

  log<<"  Bibliography "
     <<plumed.cite("Kohlhoff K, Robustelli P, Cavalli A, Salvatella A, Vendruscolo M, J. Am. Chem. Soc. 131, 13894 (2009)")
     <<plumed.cite("Camilloni C, Robustelli P, De Simone A, Cavalli A, Vendruscolo M, J. Am. Chem. Soc. 134, 3968 (2012)") <<"\n";

  const string str_cs[] = {"ha_","hn_","nh_","ca_","cb_","co_"};
  unsigned k=0;
  unsigned index=0;
  for(unsigned i=0;i<atom.size();i++) {
    for(unsigned a=0;a<atom[i].size();a++) {
      unsigned res=index+a;
      std::string num; Tools::convert(res,num);
      for(unsigned at_kind=0;at_kind<6;at_kind++){
        if(atom[i][a].exp_cs[at_kind]>0) {
          addComponentWithDerivatives(str_cs[at_kind]+num);
          componentIsNotPeriodic(str_cs[at_kind]+num);
          c_sh[res][at_kind]=k;
          k++; 
          if(!noexp) { 
            addComponent("exp"+str_cs[at_kind]+num); 
            componentIsNotPeriodic("exp"+str_cs[at_kind]+num);
            Value* comp=getPntrToComponent("exp"+str_cs[at_kind]+num);
            comp->set(atom[i][a].exp_cs[at_kind]);
            k++;
          }
        }
      }
    }
    index += atom[i].size();
  }

  /* temporary check, the idea is that I can remove NRES completely */
  if(index!=numResidues) error("NRES and the number of residues in the PDB do not match!");
 
  requestAtoms(atoms);
}

CS2Backbone::~CS2Backbone()
{
  delete[] c_sh[0];
  delete[] c_sh;
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

  for(unsigned i=0;i<atom.size();i++){
    for(unsigned a=0;a<atom[i].size();a++){
      if(atom[i][a].res_name.c_str()==res) atom[i][a].exp_cs[n] = 0.;
    }
  }
}

void CS2Backbone::read_cs(const string &file, const string &k){
  ifstream in;
  in.open(file.c_str());
  if(!in) error("CS2Backbone: Unable to open " + file);
  istream_iterator<string> iter(in), end;
  while(iter!=end){
    string tok;
    tok = *iter; ++iter;
    if(tok[0]=='#'){ ++iter; continue;}
    unsigned p = atoi(tok.c_str());
    p = p - 1;
    unsigned seg = frag_segment(p);
    p = frag_relitive_index(p,seg);
    tok = *iter; ++iter;
    double cs = atof(tok.c_str());
    if(k=="HA")     { if(atom[seg][p].pos[0]<=0) cs=0; atom[seg][p].exp_cs[0] = cs; }
    else if(k=="H") { if(atom[seg][p].pos[1]<=0) cs=0; atom[seg][p].exp_cs[1] = cs; }
    else if(k=="N") { if(atom[seg][p].pos[2]<=0) cs=0; atom[seg][p].exp_cs[2] = cs; }
    else if(k=="CA"){ if(atom[seg][p].pos[3]<=0) cs=0; atom[seg][p].exp_cs[3] = cs; }
    else if(k=="CB"){ if(atom[seg][p].pos[4]<=0) cs=0; atom[seg][p].exp_cs[4] = cs; }
    else if(k=="C") { if(atom[seg][p].pos[5]<=0) cs=0; atom[seg][p].exp_cs[5] = cs; }
  }
}

void CS2Backbone::calculate()
{
  if(pbc) makeWhole();

  if(getExchangeStep()) box_count=0;
  bool update = false;
  if(box_count==0) update = true;

  compute_ring_parameters();

  compute_dihedrals();

  unsigned index=0;
  const unsigned chainsize = atom.size();
  // CYCLE OVER MULTIPLE CHAINS
  for(unsigned s=0;s<chainsize;s++){
    const unsigned psize = atom[s].size();
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) 
    // SKIP FIRST AND LAST RESIDUE OF EACH CHAIN
    for(unsigned a=1;a<psize-1;a++){
      const Fragment *myfrag = &atom[s][a];
      // CYCLE OVER THE SIX BACKBONE CHEMICAL SHIFTS
      for(unsigned at_kind=0;at_kind<6;at_kind++){
        if(myfrag->exp_cs[at_kind]>0){
          const unsigned aa_kind = myfrag->res_kind;
          const unsigned res_type_curr = myfrag->res_type_curr;
          const unsigned res_type_prev = myfrag->res_type_prev;
          const unsigned res_type_next = myfrag->res_type_next;

          const double * CONSTAACURR = db.CONSTAACURR(aa_kind,at_kind);
          const double * CONSTAANEXT = db.CONSTAANEXT(aa_kind,at_kind);
          const double * CONSTAAPREV = db.CONSTAAPREV(aa_kind,at_kind);

          // Common constant and AATYPE
          double cs = CONSTAACURR[res_type_curr] + 
                      CONSTAANEXT[res_type_next] + 
                      CONSTAAPREV[res_type_prev];
          // this is the atom for which we are calculating the chemical shift 
          const unsigned ipos = myfrag->pos[at_kind];

          vector<unsigned> list;
          list.reserve(3*numResidues);
          list.push_back(ipos);
          vector<Vector> ff;
          ff.reserve(3*numResidues); 
          ff.push_back(Vector(0,0,0));

          //PREV
          const double * CONST_BB2_PREV = db.CONST_BB2_PREV(aa_kind,at_kind);
          const unsigned presize = myfrag->prev.size();
          for(unsigned q=0;q<presize;q++){
            if(CONST_BB2_PREV[q]==0.) continue;
            const unsigned jpos = myfrag->prev[q];
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = CONST_BB2_PREV[q]/d;

            cs += CONST_BB2_PREV[q]*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //CURR
          const double * CONST_BB2_CURR = db.CONST_BB2_CURR(aa_kind,at_kind);
          const unsigned cursize = myfrag->curr.size();
          for(unsigned q=0;q<cursize;q++){
            if(CONST_BB2_CURR[q]==0.) continue;
            const unsigned jpos = myfrag->curr[q];
            if(ipos==jpos) continue;
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = CONST_BB2_CURR[q]/d;

            cs += CONST_BB2_CURR[q]*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //NEXT
          const double * CONST_BB2_NEXT = db.CONST_BB2_NEXT(aa_kind,at_kind);
          const unsigned nexsize = myfrag->next.size();
          for(unsigned q=0;q<nexsize;q++){
            if(CONST_BB2_NEXT[q]==0.) continue;
            const unsigned jpos = myfrag->next[q];
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = CONST_BB2_NEXT[q]/d;

            cs += CONST_BB2_NEXT[q]*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }

          //SIDE CHAIN
          const double * CONST_SC2 = db.CONST_SC2(aa_kind,at_kind,res_type_curr);
          const unsigned sidsize = myfrag->side_chain.size();
          for(unsigned q=0;q<sidsize;q++){
            if(CONST_SC2[q]==0.) continue;
            const unsigned jpos = myfrag->side_chain[q];
            if(ipos==jpos) continue;
            list.push_back(jpos);
            const Vector distance = delta(getPosition(jpos),getPosition(ipos));
            const double d = distance.modulo();
            const double fact = CONST_SC2[q]/d;

            cs += CONST_SC2[q]*d;
            const Vector der = fact*distance;
            ff[0] += der;
            ff.push_back(-der);
          }
        
          //EXTRA DIST
          const double * CONST_XD  = db.CONST_XD(aa_kind,at_kind);
          const unsigned xdsize = myfrag->xd1.size();
          for(unsigned q=0;q<xdsize;q++){
            if(CONST_XD[q]==0.) continue;
            if(myfrag->xd1[q]==-1||myfrag->xd2[q]==-1) continue;
            const unsigned kpos = myfrag->xd1[q];
            const unsigned jpos = myfrag->xd2[q];
            list.push_back(jpos);
            list.push_back(kpos);
            const Vector distance = delta(getPosition(kpos),getPosition(jpos));
            const double d = distance.modulo();
            const double fact = CONST_XD[q]/d;

            cs += CONST_XD[q]*d;
            const Vector der = fact*distance;
            ff.push_back( der);
            ff.push_back(-der);
          }
           
          //NON BOND
          {
            const double * CONST_CO_SPHERE3 = db.CO_SPHERE(aa_kind,at_kind,0);
            const double * CONST_CO_SPHERE  = db.CO_SPHERE(aa_kind,at_kind,1);
            const double af1 = cutOffDist2*cutOffDist2;
            const double af3 = cutOffDist2*cutOffDist2*cutOffDist2 -3.*cutOffDist2*cutOffDist2*cutOnDist2;
            unsigned boxsize, res_curr;
            if(!update) boxsize = myfrag->box_nb[at_kind].size();
            else { 
              boxsize = getNumberOfAtoms();
              atom[s][a].box_nb[at_kind].clear();
              atom[s][a].box_nb[at_kind].reserve(2*numResidues);
              res_curr = res_num[ipos];
            }
            for(unsigned bat=0; bat<boxsize; bat++) {
              unsigned jpos;
              Vector distance;
              double d2;
              if(!update) {
                jpos = myfrag->box_nb[at_kind][bat];
                distance = delta(getPosition(jpos),getPosition(ipos));
                d2 = distance.modulo2();
              } else {
                const unsigned res_dist = abs(static_cast<int>(res_curr-res_num[bat]));
                if(res_dist<2) continue;
                jpos = bat;
                distance = delta(getPosition(jpos),getPosition(ipos));
                d2 = distance.modulo2();
                if(d2<cutOffNB) atom[s][a].box_nb[at_kind].push_back(bat);
                else continue;
              }
            
              if(d2<cutOffDist2) {
                list.push_back(jpos);
    	        const double d = sqrt(d2);
    	        const double d4 = d2*d2;
      	        const double dinv = 1.0/d;
                const double invd3 = dinv*dinv*dinv;

                double factor1;
    	        double factor3;
                double dfactor1;
                double dfactor3;

                if(d2>cutOnDist2) {
                  const double af = (cutOffDist2 - d2);
                  const double bf = (cutOffDist2 + 2.*d2 - 3.*cutOnDist2);
                  const double cf = invswitch*af*af*bf;
		  factor1 = d*cf;
                  factor3 = invd3*cf;

                  const double bf1 = invswitch*af;
                  const double cf1 = +15.*cutOnDist2*d2;
                  const double df1 = -14.*d4;
                  const double ef1 = cutOffDist2*(-3.*cutOnDist2+d2);
                  dfactor1 = bf1*(af1+cf1+df1+ef1);

                  const double bf3 =  -3.*invswitch/d4;
                  const double cf3 = +2.*cutOffDist2*cutOnDist2*d2;
                  const double df3 = d4*(cutOffDist2+cutOnDist2);
                  const double ef3 = -2.*d4*d2;
                  dfactor3 = bf3*(af3+cf3+df3+ef3);
                } else {
                  factor1 = d;
    	          factor3 = invd3;
                  dfactor1 = 1.0;
                  dfactor3 = -3./d4;
                }

    	        const unsigned t = type[jpos];
                cs += factor3*CONST_CO_SPHERE3[t] + factor1*CONST_CO_SPHERE[t];

    	        const double fact1 = dfactor1*dinv;
    	        const double fact2 = dfactor3*dinv;
    	        const double fact = fact1*CONST_CO_SPHERE[t]+fact2*CONST_CO_SPHERE3[t];
                const Vector der = fact*distance;
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
            for(unsigned i=0; i<rsize; i++){
    	      // compute angle from ring middle point to current atom position
    	      // get distance vector from query atom to ring center and normal vector to ring plane
              const Vector d = delta(ringInfo[i].position, getPosition(ipos));
    	      const Vector n = ringInfo[i].normVect;
    	      // compute square distance and distance from query atom to ring center
    	      const double dn  = dotProduct(d,n);
    	      const double dn2 = dn * dn;
    	      const double dL2 = d.modulo2();
    	      const double dL  = sqrt(dL2);
    	      const double dL3 = dL2 * dL;
    	      const double dL4 = dL2 * dL2;

    	      const double nL  = ringInfo[i].lengthNV;
    	      const double nL2 = ringInfo[i].lengthN2;

    	      const double dLnL = dL * nL;
    	      const double dL3nL3 = dL3 * nL2 * nL;

    	      const double sqrangle = dn2/(dL2*nL2);
    	      const double u = 1.-3.*sqrangle;
              cs += rc[ringInfo[i].rtype]*u/dL3;
              
    	      // calculate terms resulting from differentiating energy function with respect to query and ring atom coordinates
    	      double factor = -6 * dn / (dL4 * nL2);
              const Vector gradUQ = factor * (dL2 * n - dn * d);
    	
    	      factor = 3 * dL;
              const Vector gradVQ = factor * d;
    	      const double invdL6 = 1./(dL3 * dL3);
    	
      	      const double fact = -rc[ringInfo[i].rtype] * invdL6;
    	      ff[0] += -fact * (gradUQ * dL3 - u * gradVQ);
    	
    	      const Vector nSum = ringInfo[i].n1 + ringInfo[i].n2;

    	      // update forces on ring atoms
    	      Vector g, ab, c;
    	      const unsigned limit = ringInfo[i].numAtoms - 3; // 2 for a 5-membered ring, 3 for a 6-membered ring
    	      for(unsigned at=0; at<ringInfo[i].numAtoms; at++) {
                // atoms 0,1 (5 member) or 0,1,2 (6 member)
    	        if (at < limit) g = delta(getPosition(ringInfo[i].atom[(at+2)%3]),getPosition(ringInfo[i].atom[(at+1)%3]));
                // atoms 3,4 (5 member) or 3,4,5 (6 member)
    	        else if (at >= ringInfo[i].numAtoms - limit) {
                  // 2 for a 5-membered ring, 3 for a 6-membered ring
    	          unsigned offset = ringInfo[i].numAtoms - 3;
                  g = delta(getPosition(ringInfo[i].atom[((at+2-offset) % 3) + offset]),getPosition(ringInfo[i].atom[((at+1-offset) % 3) + offset])); 
                // atom 2 (5-membered rings)
    	        } else g = delta(getPosition(ringInfo[i].atom[4]) , getPosition(ringInfo[i].atom[3])) + 
                           delta(getPosition(ringInfo[i].atom[1]) , getPosition(ringInfo[i].atom[2]));

                ab = crossProduct(d,g);
                c  = crossProduct(nSum,g);
    	    
    	        factor = -6 * dn / dL3nL3;
    	        const double factor2 = 0.25 * dL / nL;
    	        const double OneOverN = 1 / ((double) ringInfo[i].numAtoms);
    	        const double factor3 = nL / dL * OneOverN;
                const Vector gradU = factor * ((0.5 * ab - n *  OneOverN) * dLnL - dn * (factor2 * c - factor3 * d));
    	    
    	        factor = -3 * dL * OneOverN;
                const Vector gradV = factor * d;
    	    
    	        ff.push_back(-fact*(gradU * dL3 - u * gradV));
                list.push_back(ringInfo[i].atom[at]);
              }
            }
          }
          //END OF RINGS

          //DIHEDRAL ANGLES
          {
            const double *CO_DA = db.CO_DA(aa_kind,at_kind);
    	    double *PARS_DA;
            if(myfrag->phi.size()==4){
      	      PARS_DA = db.PARS_DA(aa_kind,at_kind,0);
              const double val1 = 3.*myfrag->t_phi+PARS_DA[3];
              const double val2 = 1.*myfrag->t_phi+PARS_DA[4];
    	      cs += CO_DA[0]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = CO_DA[0]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(-fact*myfrag->dd0[0]);
              ff.push_back(-fact*(myfrag->dd1[0]-myfrag->dd0[0]));
              ff.push_back(-fact*(myfrag->dd2[0]-myfrag->dd1[0]));
              ff.push_back( fact*myfrag->dd2[0]);
              list.push_back(myfrag->phi[0]);
              list.push_back(myfrag->phi[1]);
              list.push_back(myfrag->phi[2]);
              list.push_back(myfrag->phi[3]);
            }

            if(myfrag->psi.size()==4){
    	      PARS_DA = db.PARS_DA(aa_kind,at_kind,1);
              const double val1 = 3.*myfrag->t_psi+PARS_DA[3];
              const double val2 = 1.*myfrag->t_psi+PARS_DA[4];
    	      cs += CO_DA[1]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = CO_DA[1]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(-fact*myfrag->dd0[1]);
              ff.push_back(-fact*(myfrag->dd1[1]-myfrag->dd0[1]));
              ff.push_back(-fact*(myfrag->dd2[1]-myfrag->dd1[1]));
              ff.push_back( fact*myfrag->dd2[1]);
              list.push_back(myfrag->psi[0]);
              list.push_back(myfrag->psi[1]);
              list.push_back(myfrag->psi[2]);
              list.push_back(myfrag->psi[3]);
            }

            //Chi
            if(myfrag->chi1.size()==4){
    	      PARS_DA = db.PARS_DA(aa_kind,at_kind,2);
              const double val1 = 3.*myfrag->t_chi1+PARS_DA[3];
              const double val2 = 1.*myfrag->t_chi1+PARS_DA[4];
    	      cs += CO_DA[2]*(PARS_DA[0]*cos(val1)+PARS_DA[1]*cos(val2)+PARS_DA[2]);
              const double fact = CO_DA[2]*(+3.*PARS_DA[0]*sin(val1)+PARS_DA[1]*sin(val2));

              ff.push_back(-fact*myfrag->dd0[2]);
              ff.push_back(-fact*(myfrag->dd1[2]-myfrag->dd0[2]));
              ff.push_back(-fact*(myfrag->dd2[2]-myfrag->dd1[2]));
              ff.push_back( fact*myfrag->dd2[2]);
              list.push_back(myfrag->chi1[0]);
              list.push_back(myfrag->chi1[1]);
              list.push_back(myfrag->chi1[2]);
              list.push_back(myfrag->chi1[3]);
            }
          }
          //END OF DIHE

          Value* comp=getPntrToComponent(c_sh[index+a][at_kind]);
          comp->set(cs);
          Tensor virial;
          const unsigned listsize = list.size();
          for(unsigned i=0;i<listsize;i++) {
            setAtomsDerivatives(comp,list[i],ff[i]);
            virial-=Tensor(getPosition(list[i]),ff[i]);
          }
          setBoxDerivatives(comp,virial);
        } 
      }
    }
    index += psize;
  }

  ++box_count;
  if(box_count == box_nupdate) box_count = 0;
}

void CS2Backbone::compute_ring_parameters(){
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) 
  for(unsigned i=0;i<ringInfo.size();i++){
    const unsigned size = ringInfo[i].numAtoms;
    vector<Vector> a;
    a.resize(size);
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
    ringInfo[i].n1 = crossProduct(delta(a[1],a[0]), delta(a[1], a[2]));
    // compute normal vector to plane containing last three atoms in array
    // NB: third atom of five-membered ring used for both computations above
    ringInfo[i].n2 = crossProduct(delta(a[size-2], a[size-3]), delta(a[size-2], a[size-1]));
    // ring plane normal vector is average of n1 and n2
    ringInfo[i].normVect = 0.5*(ringInfo[i].n1 + ringInfo[i].n2);
    // calculate squared length and length of normal vector
    ringInfo[i].lengthN2 = ringInfo[i].normVect.modulo2(); 
    ringInfo[i].lengthNV = sqrt(ringInfo[i].lengthN2);
  }
}

void CS2Backbone::compute_dihedrals(){
  const unsigned chainsize = atom.size();
  // CYCLE OVER MULTIPLE CHAINS
  for(unsigned s=0;s<chainsize;s++){
    const unsigned psize = atom[s].size();
    #pragma omp parallel for num_threads(OpenMP::getNumThreads()) 
    // SKIP FIRST AND LAST RESIDUE OF EACH CHAIN
    for(unsigned a=1;a<psize-1;a++){
      const Fragment *myfrag = &atom[s][a];
      // CYCLE OVER THE SIX BACKBONE CHEMICAL SHIFTS
      for(unsigned at_kind=0;at_kind<6;at_kind++){
        if(myfrag->exp_cs[at_kind]>0){
          if(myfrag->phi.size()==4){
            const Vector d0 = delta(getPosition(myfrag->phi[1]), getPosition(myfrag->phi[0]));
            const Vector d1 = delta(getPosition(myfrag->phi[2]), getPosition(myfrag->phi[1]));
            const Vector d2 = delta(getPosition(myfrag->phi[3]), getPosition(myfrag->phi[2]));
            Torsion t;
            atom[s][a].t_phi = t.compute(d0,d1,d2,atom[s][a].dd0[0],atom[s][a].dd1[0],atom[s][a].dd2[0]);
          }
          if(myfrag->psi.size()==4){
            const Vector d0 = delta(getPosition(myfrag->psi[1]), getPosition(myfrag->psi[0]));
            const Vector d1 = delta(getPosition(myfrag->psi[2]), getPosition(myfrag->psi[1]));
            const Vector d2 = delta(getPosition(myfrag->psi[3]), getPosition(myfrag->psi[2]));
            Torsion t;
            atom[s][a].t_psi = t.compute(d0,d1,d2,atom[s][a].dd0[1],atom[s][a].dd1[1],atom[s][a].dd2[1]);
          }
          if(myfrag->chi1.size()==4){
            const Vector d0 = delta(getPosition(myfrag->chi1[1]), getPosition(myfrag->chi1[0]));
            const Vector d1 = delta(getPosition(myfrag->chi1[2]), getPosition(myfrag->chi1[1]));
            const Vector d2 = delta(getPosition(myfrag->chi1[3]), getPosition(myfrag->chi1[2]));
            Torsion t;
            atom[s][a].t_chi1 = t.compute(d0,d1,d2,atom[s][a].dd0[2],atom[s][a].dd1[2],atom[s][a].dd2[2]);
          }
          break;
        }
      }
    }
  }
}

void CS2Backbone::init_backbone(const PDB &pdb){
  // number of chains
  vector<string> chains;
  pdb.getChainNames( chains );
  seg_last.resize(chains.size());
  unsigned old_size=0;

  for(unsigned i=0;i<chains.size();i++){
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
    N_.resize (resrange);
    H_.resize (resrange);
    CA_.resize(resrange);
    CB_.resize(resrange);
    HA_.resize(resrange);
    C_.resize (resrange);
    O_.resize (resrange);
    CX_.resize(resrange);

    for(unsigned a=0;a<resrange;a++){
      N_[a]  = -1;
      H_[a]  = -1;
      CA_[a] = -1;
      CB_[a] = -1;
      HA_[a] = -1;
      C_[a]  = -1;
      O_[a]  = -1;
      CX_[a] = -1;
    }

    vector<AtomNumber> allatoms = pdb.getAtomsInChain(chains[i]);
    // cycle over all the atoms in the chain
    for(unsigned a=0;a<allatoms.size();a++){
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
    for(unsigned a=start;a<=end;a++){
      unsigned f_idx = a - res_offset;
      Fragment at;
      at.box_nb.resize(6); 
      at.pos.resize(6); 
      at.exp_cs.resize(6); 
      {
        at.pos[0] = HA_[f_idx];
        at.pos[1] =  H_[f_idx];
        at.pos[2] =  N_[f_idx];
        at.pos[3] = CA_[f_idx];
        at.pos[4] = CB_[f_idx];
        at.pos[5] =  C_[f_idx];
      }
      at.res_type_prev = at.res_type_curr = at.res_type_next = 0;
      at.res_name = pdb.getResidueName(a, chains[i]); 
      at.res_kind = db.kind(at.res_name);
      at.fd = a;
      //REGISTER PREV CURR NEXT
      {
        if(a>start){
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

        if(a<end){
          at.next.push_back (N_[f_idx+1]);
          at.next.push_back (H_[f_idx+1]);
          at.next.push_back(CA_[f_idx+1]);
          at.next.push_back(HA_[f_idx+1]);
          at.next.push_back( C_[f_idx+1]);
          at.res_type_next = frag2enum(pdb.getResidueName(a+1, chains[i]));
        }

        //PHI | PSI | CH1
        if(a>start){
          at.phi.push_back( C_[f_idx-1]);
          at.phi.push_back( N_[f_idx]);
          at.phi.push_back(CA_[f_idx]);
          at.phi.push_back( C_[f_idx]);
        }
        
        if(a<end){
          at.psi.push_back( N_[f_idx]);
          at.psi.push_back(CA_[f_idx]);
          at.psi.push_back( C_[f_idx]);
          at.psi.push_back( N_[f_idx+1]);	
        }

        if(CX_[f_idx]!=-1&&CB_[f_idx]!=-1){
          at.chi1.push_back( N_[f_idx]);
          at.chi1.push_back(CA_[f_idx]);
          at.chi1.push_back(CB_[f_idx]);
          at.chi1.push_back(CX_[f_idx]);
        }
        at.dd0.resize(3);
        at.dd1.resize(3);
        at.dd2.resize(3);
      }
      atm_.push_back(at);
    }
    atom.push_back(atm_);
  }
}

void CS2Backbone::init_sidechain(const PDB &pdb){
  vector<string> chains; 
  pdb.getChainNames( chains );
  unsigned old_size=0;
  // cycle over chains
  for(unsigned s=0; s<atom.size(); s++){
    AtomNumber astart, aend; 
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();
    // cycle over residues  
    for(unsigned a=0; a<atom[s].size(); a++){
      if(atom[s][a].res_name=="UNK") continue;
      vector<AtomNumber> atm = pdb.getAtomsInResidue(atom[s][a].fd, chains[s]);
      vector<string> sc_atm = side_chain_atoms(atom[s][a].res_name);

      for(unsigned sc=0;sc<sc_atm.size();sc++){
        for(unsigned aa=0;aa<atm.size();aa++){
          if(pdb.getAtomName(atm[aa])==sc_atm[sc]){
    	    atom[s][a].side_chain.push_back(atm[aa].index()-atom_offset+old_size);
          }
        }
      }

    }
    old_size += aend.index()+1; 
  }
}

void CS2Backbone::init_xdist(const PDB &pdb){
  const string atomsP1[] = {"H", "H", "H", "C", "C", "C", 
                            "O", "O", "O", "N", "N", "N", 
                            "O", "O", "O", "N", "N", "N", 
                            "CG", "CG", "CG", "CG", "CG", 
                            "CG", "CG", "CA"};

  const int resOffsetP1 [] = {0, 0, 0, -1, -1, -1, 
                              0, 0, 0, 1, 1, 1, 
                              -1, -1, -1, 0, 0, 0, 
                              0, 0, 0, 0, 0, -1, 1, -1};

  const string atomsP2[] = {"HA", "C", "CB", "HA", "C", "CB", 
                            "HA", "N", "CB", "HA", "N", "CB", 
                            "HA", "N", "CB", "HA", "N", "CB", 
                            "HA", "N", "C", "C", "N", "CA", "CA", "CA"};

  const int resOffsetP2 [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 1, 0, 0, 1};

  vector<string> chains; 
  pdb.getChainNames( chains );
  unsigned old_size=0;
  for(unsigned s=0; s<atom.size(); s++){
    AtomNumber astart, aend; 
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();

    for(unsigned a=1; a<atom[s].size()-1; a++){
      vector<AtomNumber> atm_curr = pdb.getAtomsInResidue(atom[s][a].fd,chains[s]);
      vector<AtomNumber> atm_prev = pdb.getAtomsInResidue(atom[s][a].fd-1,chains[s]);
      vector<AtomNumber> atm_next = pdb.getAtomsInResidue(atom[s][a].fd+1,chains[s]);

      for(unsigned q=0;q<db.get_numXtraDists()-1; q++){
        vector<AtomNumber>::iterator at1, at1_end;
        vector<AtomNumber>::iterator at2, at2_end;

        bool init_p1=false;
        AtomNumber p1;
        bool init_p2=false;
        AtomNumber p2;

        if(resOffsetP1[q]== 0){ at1 = atm_curr.begin(); at1_end = atm_curr.end();}
        if(resOffsetP1[q]==-1){ at1 = atm_prev.begin(); at1_end = atm_prev.end();}
        if(resOffsetP1[q]==+1){ at1 = atm_next.begin(); at1_end = atm_next.end();}
        while(at1!=at1_end){
          AtomNumber aa = *at1;
          ++at1;
          string name = pdb.getAtomName(aa);

          xdist_name_map(name);

          if(name==atomsP1[q]){
    	    p1 = aa;
            init_p1=true;
    	    break;
          }
        }

        if(resOffsetP2[q]== 0){ at2 = atm_curr.begin(); at2_end = atm_curr.end();}
        if(resOffsetP2[q]==-1){ at2 = atm_prev.begin(); at2_end = atm_prev.end();}
        if(resOffsetP2[q]==+1){ at2 = atm_next.begin(); at2_end = atm_next.end();}
        while(at2!=at2_end){
          AtomNumber aa = *at2;
          ++at2;
          string name = pdb.getAtomName(aa);

          xdist_name_map(name);

          if(name==atomsP2[q]){
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

void CS2Backbone::init_types(const PDB &pdb){
  vector<AtomNumber> aa = pdb.getAtomNumbers();
  for(unsigned i=0;i<aa.size();i++){
    unsigned frag = pdb.getResidueNumber(aa[i]);
    string fragName = pdb.getResidueName(aa[i]);
    string atom_name = pdb.getAtomName(aa[i]);
    char atom_type = atom_name[0];
    if(isdigit(atom_name[0])) atom_type = atom_name[1];
    res_num.push_back(frag);
    unsigned t = 0;
    if (!isSP2(fragName, atom_name)){
      if (atom_type == 'C') t = D_C;
      else if (atom_type == 'O') t = D_O;
      else if (atom_type == 'H') t = D_H;
      else if (atom_type == 'N') t = D_N;
      else if (atom_type == 'S') t = D_S;
      else plumed_merror("Camshift:init_type: unknown atom type!\n");
    }else{
      if (atom_type == 'C') t = D_C2;
      else if (atom_type == 'O') t = D_O2;
      else if (atom_type == 'N') t = D_N2;
      else plumed_merror("Camshift:init_type: unknown atom type!\n");
    }
    type.push_back(t);
  }
}

void CS2Backbone::init_rings(const PDB &pdb){

  const string pheTyr_n[] = {"CG","CD1","CE1","CZ","CE2","CD2"};
  const string trp1_n[]   = {"CD2","CE2","CZ2","CH2","CZ3","CE3"};
  const string trp2_n[]   = {"CG","CD1","NE1","CE2","CD2"};
  const string his_n[]    = {"CG","ND1","CD2","CE1","NE2"};

  vector<string> chains; 
  pdb.getChainNames( chains );
  vector<AtomNumber> allatoms = pdb.getAtomNumbers();
  unsigned old_size=0;

  for(unsigned s=0; s<atom.size(); s++){
    AtomNumber astart, aend; 
    string errmsg;
    pdb.getAtomRange( chains[s], astart, aend, errmsg );
    unsigned atom_offset = astart.index();
    for(unsigned r=0; r<atom[s].size(); r++){
      string frg = pdb.getResidueName(atom[s][r].fd);
      if(!((frg=="PHE")||(frg=="TYR")||(frg=="TRP")||
           (frg=="HIS")||(frg=="HIP")||(frg=="HID")||
           (frg=="HIE")||(frg=="HSD")||(frg=="HSE")||
           (frg=="HSP"))) continue;

      vector<AtomNumber> frg_atoms = pdb.getAtomsInResidue(atom[s][r].fd,chains[s]);

      if(frg=="PHE"||frg=="TYR"){
        RingInfo ri;
        for(unsigned a=0;a<frg_atoms.size();a++){
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0;aa<6;aa++){
            if(pdb.getAtomName(frg_atoms[a])==pheTyr_n[aa]){
    	      ri.atom[aa] = atm;
    	      break;
            }
          }
        }
        ri.numAtoms = 6;
        if(frg=="PHE") ri.rtype = RingInfo::R_PHE;
        if(frg=="TYR") ri.rtype = RingInfo::R_TYR;
        ringInfo.push_back(ri);

      } else if(frg=="TRP"){
        //First ring
        RingInfo ri;
        for(unsigned a=0;a<frg_atoms.size();a++){
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0;aa<6;aa++){
            if(pdb.getAtomName(frg_atoms[a])==trp1_n[aa]){
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
        for(unsigned a=0;a<frg_atoms.size();a++){
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0;aa<5;aa++){
            if(pdb.getAtomName(frg_atoms[a])==trp2_n[aa]){
              ri2.atom[aa] = atm;
              break;
            }
          }
        }
        ri2.numAtoms = 5;
        ri2.rtype = RingInfo::R_TRP2;
        ringInfo.push_back(ri2);

      } else { //HIS case
        RingInfo ri;
        for(unsigned a=0;a<frg_atoms.size();a++){
          unsigned atm = frg_atoms[a].index()-atom_offset+old_size;
          for(unsigned aa=0;aa<5;aa++){
            if(pdb.getAtomName(frg_atoms[a])==his_n[aa]){
    	      ri.atom[aa] = atm;
    	      break;
            }
          }
        }
        ri.numAtoms = 5;
        ri.rtype = RingInfo::R_HIS;
        ringInfo.push_back(ri);
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
  else if (aa == "GME") type = GLY;
  else if (aa == "HIS") type = HIS;
  else if (aa == "HSE") type = HIS;
  else if (aa == "HIE") type = HIS;
  else if (aa == "HSP") type = HIS;
  else if (aa == "HIP") type = HIS;
  else if (aa == "HSD") type = HIS;
  else if (aa == "HID") type = HIS;
  else if (aa == "ILE") type = ILE;
  else if (aa == "IME") type = ILE;
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

vector<string> CS2Backbone::side_chain_atoms(const string &s){
  vector<string> sc;

  if(s=="ALA"){
    sc.push_back( "CB" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    return sc;
  } else if(s=="ARG"){
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
  } else if(s=="ASN"){
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
  } else if(s=="ASP"||s=="ASH"){
    sc.push_back( "CB" );
    sc.push_back( "CG" );
    sc.push_back( "OD1" );
    sc.push_back( "OD2" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    return sc;
  } else if(s=="CYS"||s=="CYM"){
    sc.push_back( "CB" );
    sc.push_back( "SG" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG" );
    return sc;
  } else if(s=="GLN"){
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
  } else if(s=="GLU"||s=="GLH"){
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
  } else if(s=="GLY"||s=="GME"){
    sc.push_back( "HA2" );
    return sc;
  } else if(s=="HIS"||s=="HSE"||s=="HIE"||s=="HSD"||s=="HID"||s=="HIP"||s=="HSP"){
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
  } else if(s=="ILE"||s=="IME"){
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
  } else if(s=="LEU"){
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
  } else if(s=="LYS"){
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
  } else if(s=="MET"){
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
  } else if(s=="PHE"){
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
  } else if(s=="PRO"){
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
  } else if(s=="SER"){
    sc.push_back( "CB" );
    sc.push_back( "OG" );
    sc.push_back( "HB1" );
    sc.push_back( "HB2" );
    sc.push_back( "HB3" );
    sc.push_back( "HG1" );
    sc.push_back( "HG" );
    return sc;
  } else if(s=="THR"){
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
  } else if(s=="TRP"){
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
  } else if(s=="TYR"){
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
  } else if(s=="VAL"){
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

bool CS2Backbone::is_chi1_cx(const string & frg, const string & atm){
  if(atm=="CG")                                        return true;
  if((frg == "CYS")&&(atm =="SG"))                     return true;
  if(((frg == "ILE")||(frg == "VAL"))&&(atm == "CG1")) return true;
  if((frg == "SER")&&(atm == "OG"))                    return true;
  if((frg == "THR")&&(atm == "OG1"))                   return true;

  return false;
}

unsigned CS2Backbone::frag_segment(const unsigned p){
  unsigned s = 0;
  for(unsigned i=0;i<seg_last.size()-1;i++){
    if(p>seg_last[i]) s  = i+1;
    else break;
  }
  return s;
}

unsigned CS2Backbone::frag_relitive_index(const unsigned p, const unsigned s){
  if(s==0) return p;
  return p-seg_last[s-1];
}

void CS2Backbone::debug_report(){
  printf("\t CS2Backbone Initialization report: \n");
  printf("\t -------------------------------\n");
  printf("\t Number of segments: %u\n", static_cast<unsigned>(atom.size()));
  printf("\t Segments size:      ");
  for(unsigned i=0;i<atom.size();i++) printf("%u ", static_cast<unsigned>(atom[i].size())); printf("\n");
  printf("\t%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s \n",
          "Seg","N","AA","Prev","Curr","Next","SC","XD1","XD2","Phi","Psi","Chi1");
  for(unsigned i=0;i<atom.size();i++){
    for(unsigned j=0;j<atom[i].size();j++){
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

      for(unsigned k=0;k<atom[i][j].prev.size();k++) printf("%8i ", atom[i][j].prev[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].curr.size();k++) printf("%8i ", atom[i][j].curr[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].next.size();k++) printf("%8i ", atom[i][j].next[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].side_chain.size();k++) printf("%8i ", atom[i][j].side_chain[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].xd1.size();k++) printf("%8i ", atom[i][j].xd1[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].xd2.size();k++) printf("%8i ", atom[i][j].xd2[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].phi.size();k++) printf("%8i ", atom[i][j].phi[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].psi.size();k++) printf("%8i ", atom[i][j].psi[k]); printf("\n");
      for(unsigned k=0;k<atom[i][j].chi1.size();k++) printf("%8i ", atom[i][j].chi1[k]); printf("\n");

    }
  }

  printf("\t Rings: \n");
  printf("\t ------ \n");
  printf("\t Number of rings: %u\n", static_cast<unsigned>(ringInfo.size()));
  printf("\t%8s %8s %8s %8s\n", "Num","Type","RType","N.atoms");
  for(unsigned i=0;i<ringInfo.size();i++){
    printf("\t%8u %8u %8u \n",i+1,ringInfo[i].rtype,ringInfo[i].numAtoms);
    for(unsigned j=0;j<ringInfo[i].numAtoms;j++) printf("%8u ", ringInfo[i].atom[j]); printf("\n");
  } 
}

void CS2Backbone::xdist_name_map(string & name){
  if((name == "OT1")||(name == "OC1")) name = "O";
  else if ((name == "HN") || (name == "HT1") || (name == "H1")) name = "H";
  else if ((name == "CG1")|| (name == "OG")|| 
           (name == "SG") || (name == "OG1")) name = "CG";
  else if ((name == "HA1"))                   name = "HA";
}

}
}
