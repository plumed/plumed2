#ifdef __PLUMED_HAS_ALMOST

#ifndef __CAMSHIFT3__
#define __CAMSHIFT3__

#define cutOffDist    0.55  	// cut off distance for non-bonded pairwise forces
#define cutOffDist2   0.3025 	// square of cutOffDist
#define cutOnDist     0.45   	// cut off distance for non-bonded pairwise forces
#define cutOnDist2    0.2025 	// square of cutOffDist
#define invswitch     1.0/((cutOffDist2 - cutOnDist2)*(cutOffDist2 - cutOnDist2)*(cutOffDist2 - cutOnDist2))
#define CSDIM         3

#include <string>
#include <sstream>
#include <vector>

#include <almost/mdb.h>
#include <almost/molecules/molecules.h>

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/OpenMP.h"
#include "tools/Torsion.h"

using namespace std;

namespace PLMD {

  class CamShiftDB {
    enum { STD, GLY, PRO};

    //6 ATOMS TYPES
    enum { HA_ATOM, H_ATOM, N_ATOM, CA_ATOM, CB_ATOM, C_ATOM };

    static const int aa_kind = 3;
    static const int atm_kind = 6;
    static const int numXtraDists = 27;

    double dscale;

    double c_all[aa_kind][atm_kind];
    // ALA, ARG, ASN, ASP, CYS, GLU, GLN, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
    double c_aa[aa_kind][atm_kind][20];
    double c_aa_prev[aa_kind][atm_kind][20];
    double c_aa_succ[aa_kind][atm_kind][20];
    double co_bb[aa_kind][atm_kind][2][16];
    double co_sc_[aa_kind][atm_kind][20][2][20];
    double co_sphere[aa_kind][atm_kind][2][8];

    // for dihedral angles
    // co * (a * cos(3 * omega + c) + b * cos(omega + d))
    double co_da[aa_kind][atm_kind][3];
    double pars_da[aa_kind][atm_kind][3][5];

    // for ring current effects
    // Phe, Tyr, Trp_1, Trp_2, His
    double co_ring[aa_kind][atm_kind][5];

    // for extra distances
    double co_xd[aa_kind][atm_kind][numXtraDists];

  public:
    CamShiftDB(string file, const double scale){ 
      dscale = scale;
      parse(file);
    }

    int kind(string s){
      if(s=="GLY") return GLY;
      if(s=="PRO") return PRO;      
      return STD;
    }

    int atom_kind(string s){
      if(s=="HA")return HA_ATOM;
      if(s=="H") return H_ATOM;
      if(s=="N") return N_ATOM;
      if(s=="CA")return CA_ATOM;
      if(s=="CB")return CB_ATOM;
      if(s=="C") return C_ATOM;
      return -1;
    }

    //PARAMETERS
    double   CONST(int a_kind, int at_kind){ return c_all[a_kind][at_kind];}
    double * CONSTAACURR(int a_kind, int at_kind){return c_aa[a_kind][at_kind];}
    double * CONSTAANEXT(int a_kind, int at_kind){return c_aa_succ[a_kind][at_kind];}
    double * CONSTAAPREV(int a_kind, int at_kind){return c_aa_prev[a_kind][at_kind];}

    double * CONST_BB2_PREV(int a_kind, int at_kind){return co_bb[a_kind][at_kind][1];}
    double * CONST_BB2_CURR(int a_kind, int at_kind){return co_bb[a_kind][at_kind][1]+5;}
    double * CONST_BB2_NEXT(int a_kind, int at_kind){return co_bb[a_kind][at_kind][1]+11;}

    double * CONST_SC2(int a_kind, int at_kind, int res_type){ return co_sc_[a_kind][at_kind][res_type][1];}
    double * CONST_XD(int a_kind, int at_kind){ return co_xd[a_kind][at_kind];}

    double * CO_SPHERE(int a_kind, int at_kind, int exp_type){ return co_sphere[a_kind][at_kind][exp_type];}
    
    double * CO_RING(int a_kind, int at_kind){ return co_ring[a_kind][at_kind];}
    
    double * CO_DA(int a_kind, int at_kind){ return co_da[a_kind][at_kind];}
    double * PARS_DA(int a_kind, int at_kind, int ang_kind){ return pars_da[a_kind][at_kind][ang_kind];}

  private:

    void parse(string file){
      ifstream in;
      in.open(file.c_str());
      if(!in) plumed_merror("Unable to open CamShiftDB file\n");

      int c_kind = -1;
      int c_atom = -1;
      int nline = 0;

      for(int i=0;i<3;i++) for(int j=0;j<6;j++) {
        c_all[i][j]=0.;
        for(int k=0;k<20;k++) {
          c_aa[i][j][k]=0.;
          c_aa_prev[i][j][k]=0.;
          c_aa_succ[i][j][k]=0.;
          for(int l=0;l<2;l++) for(int m=0;m<20;m++) co_sc_[i][j][k][l][m]=0.;
        }
        for(int k=0;k<16;k++) {co_bb[i][j][0][k]=0.; co_bb[i][j][1][k]=0.;}
        for(int k=0;k<8;k++) { co_sphere[i][j][0][k]=0.; co_sphere[i][j][1][k]=0.; }
        for(int k=0;k<3;k++) {
          co_da[i][j][k]=0.;
          for(int l=0;l<5;l++) pars_da[i][j][k][l]=0.;
        }
        for(int k=0;k<5;k++) co_ring[i][j][k]=0.;
        for(int k=0;k<numXtraDists;k++) co_xd[i][j][k]=0.;
      }

      while(!in.eof()){
	string line;
	getline(in,line);
	++nline;
	if(line.find("#")==0) continue;
	vector<string> tok;
	vector<string> tmp;
	tok = split(line,' ');
	for(unsigned int q=0;q<tok.size();q++)
	  if(tok[q].size()) tmp.push_back(tok[q]);
	tok = tmp;
	//PROCESS
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
	  c_all[c_kind][c_atom]= atof(tok[1].c_str());
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
          //angstrom to nm
	  assign(co_bb[c_kind][c_atom][0],tok,dscale);
	  continue;
	}
	else if (tok[0] == "COBB2"){
          //angstrom to nm
	  assign(co_bb[c_kind][c_atom][1],tok,dscale);
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
	  for(unsigned int i=1;i<tok.size();i++)
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

	for (int scC = 0; scC < 20; scC++){
	  if(tok[0]==scIdent1[scC]){
            //angstrom to nm
	    assign(co_sc_[c_kind][c_atom][scC][0],tok,dscale);
	    ok = true; 
            break;
	  }
	}
	if(ok) continue;

	string scIdent2 [] = {"COSCALA2", "COSCARG2", "COSCASN2", "COSCASP2", "COSCCYS2", "COSCGLN2", "COSCGLU2", 
                              "COSCGLY2", "COSCHIS2", "COSCILE2", "COSCLEU2", "COSCLYS2", "COSCMET2", "COSCPHE2", 
                              "COSCPRO2", "COSCSER2", "COSCTHR2", "COSCTRP2", "COSCTYR2", "COSCVAL2"};

	for (int scC = 0; scC < 20; scC++){
	  if(tok[0]==scIdent2[scC]){
            //angstrom to nm
	    assign(co_sc_[c_kind][c_atom][scC][1],tok,dscale);
	    ok = true; break;
	  }
	}
	if(ok) continue;
	
	if(tok.size()) {
          string str_err = "CAMSHIFTDB WARNING: unrecognized token: " + tok[0];
          plumed_merror(str_err);
        }
      }
    }

    vector<string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
        elems.push_back(item);
      }
      return elems;
    }

    vector<string> split(const std::string &s, char delim) {
      std::vector<std::string> elems;
      split(s, delim, elems);
      return elems;
    }
    
    void assign(double * f, const vector<string> & v, const double scale){
      for(unsigned int i=1;i<v.size();i++)
	f[i-1] = scale*(atof(v[i].c_str()));
      
    }
  };

  class CamShift3
  {
    struct Fragment {
      //Order HA H N CA CB C
      vector<int> pos;
      vector<int> prev;
      vector<int> curr;
      vector<int> next;

      int res_type_prev;
      int res_type_curr;
      int res_type_next;

      vector<int> side_chain;
      vector<int> xd1;
      vector<int> xd2;
      string res_name;
      int res_kind;
      int fd;

      vector<int> phi;
      vector<int> psi;
      vector<int> chi1;

      vector<double> exp_cs;

      vector<vector<int> > box_nb;
 
      Fragment() {
	res_type_prev = res_type_curr = res_type_next = -1;
        res_kind = fd = -1;
        box_nb.resize(6); 
        for(unsigned int i=0;i<6;i++) box_nb[i].resize(1000);
        pos.resize(6); 
        exp_cs.resize(6); 
        for(unsigned int i=0;i<6;i++) pos[i] = -1;
      }
    };

    vector<int> seg_last; // i<seg_last[0] is i n 0...

    struct RingInfo{
      // one out of five different types
      int type;
      enum ring_t {R_PHE, R_TYR, R_TRP1, R_TRP2, R_HIS};
      int rtype;
      // up to six member per ring
      int atom[6];
      // number of ring members (5 or 6)
      int numAtoms;
      // center of ring coordinates
      double position[3];
      // ring plane normal vector
      double normVect[3];
      // square of length of normVect
      double lengthN2;
      // length of normVect
      double lengthNV;
      // two atom plane normal vectors used to compute ring plane normal
      vector<double> n1;
      vector<double> n2;
      //Static
      string atomNames[4][6];
      string types[4];
      int    init_;
      RingInfo() {
        type = rtype = -1;
        for(int i = 0; i < 6; i++) atom[i]=-1;
        numAtoms = -1;
        position[0] = position[1] = position[2] = -1;
        normVect[0] = normVect[1] = normVect[2] = -1;
        lengthN2 = lengthNV = -1;

        const char* pheTyr_n[] = {"CG","CD1","CE1","CZ","CE2","CD2"};
        const char* trp1_n[]   = {"CD2","CE2","CZ2","CH2","CZ3","CE3"};
        const char* trp2_n[]   = {"CG","CD1","NE1","CE2","CD2"};
        const char* his_n[]    = {"CG","ND1","CD2","CE1","NE2"};

        for (int i = 0; i < 6; i++){
          atomNames[0][i] = pheTyr_n[i];
          atomNames[1][i] = trp1_n[i];
          if (i < 5){
            atomNames[2][i] = trp2_n[i];
            atomNames[3][i] = his_n[i];
          }
        }
 
        types[0] = "PHE";
        types[1] = "TYR";
        types[2] = "TRP";
        types[3] = "HIS";
        init_  = 1;      
      } 
    };

    vector<RingInfo> ringInfo;

    //TYPES
    vector<int> type;
    vector<int> res_num;
    enum aa_t {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK};
    enum atom_t {D_C, D_H, D_N, D_O, D_S, D_C2, D_N2, D_O2};
    CamShiftDB db;
    static const int numXtraDists = 27;
    int box_nupdate;
    int box_count;

  public:
    vector<vector<Fragment> > atom;
  
    CamShift3(const Almost::Molecules & molecules, string file, const double scale):db(file,scale){
      //Init
      init_backbone(molecules);
      init_sidechain(molecules);
      init_xdist(molecules);
      init_types(molecules);
      init_rings(molecules);
      //Non-Bonded neighbour lists 
      box_nupdate=10;
      box_count=0;
    }

    void set_box_nupdate(const int v){box_nupdate = v;}
    void set_box_count(const int v){box_count = v;}
   
    void remove_problematic(string res, string nucl) {
      unsigned n;
      if(nucl=="HA")     n=0;
      else if(nucl=="H") n=1;
      else if(nucl=="N") n=2;
      else if(nucl=="CA")n=3;
      else if(nucl=="CB")n=4;
      else if(nucl=="C") n=5;
      else return;

      for(unsigned int i=0;i<atom.size();i++){
	for(unsigned int a=0;a<atom[i].size();a++){
          if(atom[i][a].res_name.c_str()==res) atom[i][a].exp_cs[n] = 0.;
        }
      }
    }

    void read_cs(string file, string k){
      ifstream in;
      in.open(file.c_str());
      if(!in){
        string str_err = "CS2Backbone: Unable to open " + file; 
	plumed_merror(str_err);
      }
      istream_iterator<string> iter(in), end;
      while(iter!=end){
	string tok;
	tok = *iter; ++iter;
	if(tok[0]=='#'){ ++iter; continue;}
	int p = atoi(tok.c_str());
	p = p - 1;
	int seg = frag_segment(p);
	p = frag_relitive_index(p,seg);
	tok = *iter; ++iter;
	double cs = atof(tok.c_str());
	if(k=="HA")     atom[seg][p].exp_cs[0] = cs;
	else if(k=="H") atom[seg][p].exp_cs[1] = cs;
	else if(k=="N") atom[seg][p].exp_cs[2] = cs;
	else if(k=="CA")atom[seg][p].exp_cs[3] = cs;
	else if(k=="CB")atom[seg][p].exp_cs[4] = cs;
	else if(k=="C") atom[seg][p].exp_cs[5] = cs;
      }
    }

    void new_calc_cs(const vector<double> & coor, vector<double> & ff, const int N,  double **sh){

      compute_ring_parameters(coor);

      bool update = false;
      if(box_count==0) update = true;

      unsigned index=0;
      // CYCLE OVER MULTIPLE CHAINS
      for(unsigned int s=0;s<atom.size();s++){
	// SKIP FIRST AND LAST RESIDUE OF EACH CHAIN
#pragma omp parallel for num_threads(OpenMP::getNumThreads())
	for(unsigned int a=0;a<atom[s].size();a++){
          // CYCLE OVER THE SIX BACKBONE CHEMICAL SHIFTS
	  for(unsigned int at_kind=0;at_kind<6;at_kind++){
	    double cs;
	    double cs_deriv = -1.;
            
	    if(atom[s][a].pos[at_kind]>0&&atom[s][a].exp_cs[at_kind]>0){

              // this is the counter to find your place in ff
              // place is residue*6*CSDIM*N + at_kind*CSDIM*N
              // residue is atoms[s-1].size()+a
              //int place = (iamresidue+a)*6*CSDIM*N+at_kind*CSDIM*N;
              int place = (index+a)*6*CSDIM*N+at_kind*CSDIM*N;

	      int aa_kind = atom[s][a].res_kind;
	      int res_type_curr = atom[s][a].res_type_curr;
	      int res_type_prev = atom[s][a].res_type_prev;
	      int res_type_next = atom[s][a].res_type_next;

	      double   CONST = db.CONST(aa_kind,at_kind);
	      double * CONSTAACURR = db.CONSTAACURR(aa_kind,at_kind);
	      double * CONSTAANEXT = db.CONSTAANEXT(aa_kind,at_kind);
	      double * CONSTAAPREV = db.CONSTAAPREV(aa_kind,at_kind);

	      double * CONST_BB2_PREV = db.CONST_BB2_PREV(aa_kind,at_kind);
	      double * CONST_BB2_NEXT = db.CONST_BB2_NEXT(aa_kind,at_kind);
	      double * CONST_BB2_CURR = db.CONST_BB2_CURR(aa_kind,at_kind);

	      double * CONST_SC2 = db.CONST_SC2(aa_kind,at_kind, res_type_curr);
	      double * CONST_XD  = db.CONST_XD(aa_kind,at_kind);

	      //1. Common constant
	      cs = CONST;

	      //2. AATYPE	  
	      cs = cs + CONSTAACURR[res_type_curr];
	      cs = cs + CONSTAANEXT[res_type_next];
	      cs = cs + CONSTAAPREV[res_type_prev];

	      //3. distances const
	    
	      int ipos = CSDIM*atom[s][a].pos[at_kind];

	      //PREV
	      for(unsigned int q=0;q<atom[s][a].prev.size();q++){
	        double d2,dx,dy,dz;
	        int jpos = CSDIM*atom[s][a].prev[q];
	        dx = coor[ipos] - coor[jpos];
	        dy = coor[ipos+1] - coor[jpos+1];
	        dz = coor[ipos+2] - coor[jpos+2];
	        d2 = dx*dx+dy*dy+dz*dz;

	        double d  = sqrt(d2);
	        cs = cs + CONST_BB2_PREV[q]*d;

	        double dinv = 1.0/d;

  	        double fact = -cs_deriv*CONST_BB2_PREV[q]*dinv;
	        ff[place+ipos]   = ff[place+ipos]   + fact*dx;
	        ff[place+ipos+1] = ff[place+ipos+1] + fact*dy;
	        ff[place+ipos+2] = ff[place+ipos+2] + fact*dz;

	        ff[place+jpos]   = ff[place+jpos]   - fact*dx;
	        ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
	        ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;
	      }

	      //CURR
	      for(unsigned int q=0;q<atom[s][a].curr.size();q++){
	        int jpos = CSDIM*atom[s][a].curr[q];
	        if(ipos==jpos) continue;
	        double d2,dx,dy,dz;
	        dx = coor[ipos] - coor[jpos];
	        dy = coor[ipos+1] - coor[jpos+1];
	        dz = coor[ipos+2] - coor[jpos+2];
	        d2 = dx*dx+dy*dy+dz*dz;

	        double d  = sqrt(d2);
	        cs = cs + CONST_BB2_CURR[q]*d;

	        double dinv = 1.0/d;
	        double fact = -cs_deriv*CONST_BB2_CURR[q]*dinv;
	        ff[place+ipos]   = ff[place+ipos]   + fact*dx;
	        ff[place+ipos+1] = ff[place+ipos+1] + fact*dy;
	        ff[place+ipos+2] = ff[place+ipos+2] + fact*dz;

	        ff[place+jpos]   = ff[place+jpos]   - fact*dx;
	        ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
	        ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;
	      }

	      //NEXT
	      for(unsigned int q=0;q<atom[s][a].next.size();q++){
	        double d2,dx,dy,dz;
	        int jpos = CSDIM*atom[s][a].next[q];
	        dx = coor[ipos] - coor[jpos];
	        dy = coor[ipos+1] - coor[jpos+1];
	        dz = coor[ipos+2] - coor[jpos+2];
	        d2 = dx*dx+dy*dy+dz*dz;
	    
	        double d  = sqrt(d2);
	        cs = cs + CONST_BB2_NEXT[q]*d;

	        double dinv = 1.0/d;
	        double fact = -cs_deriv*CONST_BB2_NEXT[q]*dinv;

	        ff[place+ipos]   = ff[place+ipos]   + fact*dx;
	        ff[place+ipos+1] = ff[place+ipos+1] + fact*dy;
	        ff[place+ipos+2] = ff[place+ipos+2] + fact*dz;

	        ff[place+jpos]   = ff[place+jpos]   - fact*dx;
	        ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
	        ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;	    
	      }

	      //SIDE CHAIN
	      for(unsigned int q=0;q<atom[s][a].side_chain.size();q++){
	        int jpos = CSDIM*atom[s][a].side_chain[q];
	        if(ipos==jpos) continue;
	        double d2,dx,dy,dz;
	        dx = coor[ipos] - coor[jpos];
	        dy = coor[ipos+1] - coor[jpos+1];
	        dz = coor[ipos+2] - coor[jpos+2];
	        d2 = dx*dx+dy*dy+dz*dz;
	    
	        double d  = sqrt(d2);
	        cs = cs + CONST_SC2[q]*d;

	        double dinv = 1.0/d;
	        double fact = -cs_deriv*CONST_SC2[q]*dinv;
	        ff[place+ipos]   = ff[place+ipos]   + fact*dx;
	        ff[place+ipos+1] = ff[place+ipos+1] + fact*dy;
	        ff[place+ipos+2] = ff[place+ipos+2] + fact*dz;
	        ff[place+jpos]   = ff[place+jpos]   - fact*dx;
	        ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
	        ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;	    	    
	      }
	    
	      //EXTRA DIST
	      for(unsigned int q=0;q<atom[s][a].xd1.size();q++){
	        double d2,dx,dy,dz;
	        if(atom[s][a].xd1[q]==-1||atom[s][a].xd2[q]==-1) continue;
	        int kpos = CSDIM*atom[s][a].xd1[q];
	        int jpos = CSDIM*atom[s][a].xd2[q];
	        dx = coor[kpos] - coor[jpos];
	        dy = coor[kpos+1] - coor[jpos+1];
	        dz = coor[kpos+2] - coor[jpos+2];
	        d2 = dx*dx+dy*dy+dz*dz;

	        double d  = sqrt(d2);
	        cs = cs + CONST_XD[q]*d;

	        double dinv = 1.0/d;
	        double fact = -cs_deriv*CONST_XD[q]*dinv;

	        ff[place+kpos]   = ff[place+kpos]   + fact*dx;
	        ff[place+kpos+1] = ff[place+kpos+1] + fact*dy;
	        ff[place+kpos+2] = ff[place+kpos+2] + fact*dz;

	        ff[place+jpos]   = ff[place+jpos]   - fact*dx;
	        ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
	        ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;
	      }
	    
	      //NON BOND
	      {
	        double * CONST_CO_SPHERE3 = db.CO_SPHERE(aa_kind,at_kind,0);
	        double * CONST_CO_SPHERE  = db.CO_SPHERE(aa_kind,at_kind,1);
	        double dist_sum3[8] = {0,0,0,0,0,0,0,0};
	        double dist_sum [8] = {0,0,0,0,0,0,0,0};

	        int curr = atom[s][a].pos[at_kind];
                if(update) update_box(curr, atom[s][a].box_nb[at_kind], coor, 0.56);
  	        double ffi[3] = {0,0,0};
                 
                for(unsigned bat = 0; bat<atom[s][a].box_nb[at_kind].size(); bat++) {
                  int at = atom[s][a].box_nb[at_kind][bat];

		  int jpos = CSDIM*at;
		  double dx = coor[ipos] - coor[jpos];
		  double dy = coor[ipos+1] - coor[jpos+1];
		  double dz = coor[ipos+2] - coor[jpos+2];
		  double d2 = dx*dx+dy*dy+dz*dz;
                
                  if(d2<cutOffDist2) {
		    double d = sqrt(d2);
		    double d3inv = 1.0/(d*d2);
                    double factor1=d;
                    double factor3=d3inv;

  		    double dinv = 1.0/d;
		    double d4 = d2*d2;
		    double d4inv = 1.0/d4;
                    double dfactor1 = 1.0;
                    double dfactor3 = -3.*d4inv;

                    if(d2>cutOnDist2) {
                      double af = (cutOffDist2 - d2);
                      double bf = (cutOffDist2 + 2.*d2 - 3.*cutOnDist2);
                      double cf = invswitch*af*af*bf;
    		      factor1 = factor1*cf;
                      factor3 = factor3*cf;

                      double af1 = invswitch*af;
                      double bf1 = cutOffDist2*cutOffDist2;
                      double cf1 = +15.*cutOnDist2*d2;
                      double df1 = -14.*d4;
                      double ef1 = cutOffDist2*(-3.*cutOnDist2+d2);
                      dfactor1 = af1*(bf1+ef1+cf1+df1);

                      double af3 = cutOffDist2*cutOffDist2*cutOffDist2 -3.*cutOffDist2*cutOffDist2*cutOnDist2;
                      double cf3 = +2.*cutOffDist2*cutOnDist2*d2;
                      double df3 = d4*(cutOffDist2+cutOnDist2);
                      double ef3 = -2.*d4*d2;
                      dfactor3 = dfactor3*invswitch*(af3+cf3+df3+ef3);
                    }
		    int t = type[at];
		    dist_sum3[t] += factor3;
		    dist_sum [t] += factor1;

		    double fact1 = -cs_deriv*dfactor1*dinv;
		    double fact2 = -cs_deriv*dfactor3*dinv;
		    double fact = fact1*CONST_CO_SPHERE[t]+fact2*CONST_CO_SPHERE3[t];
		    ffi[0]   = ffi[0]   + fact*dx;
		    ffi[1]   = ffi[1]   + fact*dy;
		    ffi[2]   = ffi[2]   + fact*dz;
		    ff[place+jpos]   = ff[place+jpos]   - fact*dx;
		    ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
		    ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;
                  }
	        }
	        ff[place+ipos] = ff[place+ipos] + ffi[0];
	        ff[place+ipos+1] = ff[place+ipos+1] + ffi[1];
	        ff[place+ipos+2] = ff[place+ipos+2] + ffi[2];

	        for(unsigned int tt=0;tt<8;tt++){
		  cs += dist_sum3[tt]*CONST_CO_SPHERE3[tt];
		  cs += dist_sum [tt]*CONST_CO_SPHERE [tt];
	        }
	      }
	      //END NON BOND

	      //RINGS
	      {
	        double *rc = db.CO_RING(aa_kind,at_kind);
	        double contrib [] = {0.0, 0.0, 0.0, 0.0, 0.0};
	        double contribTot = 0.0;
	        for (unsigned int i = 0; i < ringInfo.size(); i++){
		  // compute angle from ring middle point to current atom position
		  // get distance vector from query atom to ring center and normal vector to ring plane
  		  double fact = cs_deriv*rc[ringInfo[i].rtype];
	          double d[3], n[3];
		  for (int j = 0; j < 3; j++){
		    d[j] = coor[ipos+j] - ringInfo[i].position[j];
		    n[j] = ringInfo[i].normVect[j];
		  }
		  // compute square distance and distance from query atom to ring center
		  double dL2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
		  double dL = sqrt(dL2);
		  double v = dL2 * dL;
		  double dL4 = dL2 * dL2;
		  double nL = ringInfo[i].lengthNV;
		  double nL2 = ringInfo[i].lengthN2;
		  double dLnL = dL * nL;
		  double dL3nL3 = v * nL2 * nL;
		
		  double dn = d[0]*n[0] + d[1]*n[1] + d[2]*n[2];
		  double dn2 = dn * dn;
		
  		  int aPos[6];
  		  for (int j = 0; j < ringInfo[i].numAtoms; j++) aPos[j] = CSDIM*ringInfo[i].atom[j];

		  // here we save the necessity of calculating the angle first and then the cos, 
                  // using that the acos and cos cancel each other out
		  double angle = fabs( dn/(dL*nL) );
		  contrib[ringInfo[i].rtype] += (1 - 3 * angle * angle) / v;
                  
		  // calculate terms resulting from differentiating energy function with respect to query and ring atom coordinates
		
		  double factor = -6 * dn / (dL4 * nL2);
		  double gradUQx = factor * (dL2 * n[0] - dn * d[0]);
		  double gradUQy = factor * (dL2 * n[1] - dn * d[1]);
		  double gradUQz = factor * (dL2 * n[2] - dn * d[2]);
		  double u = 1 - 3 * dn2 / (dL2 * nL2);
		
		  factor = 3 * dL;
		  double gradVQx = factor * d[0];
		  double gradVQy = factor * d[1];
		  double gradVQz = factor * d[2];
		  double v2 = v * v;
		
		  // update forces on query atom
		  ff[place+ipos  ] += -fact * (gradUQx * v - u * gradVQx) / v2;
		  ff[place+ipos+1] += -fact * (gradUQy * v - u * gradVQy) / v2;
		  ff[place+ipos+2] += -fact * (gradUQz * v - u * gradVQz) / v2;
		
		  double nSum[3] = {ringInfo[i].n1[0] + ringInfo[i].n2[0],
				    ringInfo[i].n1[1] + ringInfo[i].n2[1], 
				    ringInfo[i].n1[2] + ringInfo[i].n2[2]};
		  double g[3], ab[3], c[3];
		  int limit = ringInfo[i].numAtoms - 3; // 2 for a 5-membered ring, 3 for a 6-membered ring
		
		  // update forces on ring atoms
		  for (int at = 0; at < ringInfo[i].numAtoms; at++)
		  {
		    if (at < limit) // atoms 0,1 (5 member) or 0,1,2 (6 member)
		    {
		      g[0] = coor[aPos[(at+1)%3]  ] - coor[aPos[(at+2)%3]  ];
		      g[1] = coor[aPos[(at+1)%3]+1] - coor[aPos[(at+2)%3]+1];
		      g[2] = coor[aPos[(at+1)%3]+2] - coor[aPos[(at+2)%3]+2];
		    }
		    else if (at >= ringInfo[i].numAtoms - limit) // atoms 3,4 (5 member) or 3,4,5 (6 member)
		    {
		      int offset = ringInfo[i].numAtoms - 3; // 2 for a 5-membered ring, 3 for a 6-membered ring
		      g[0] = coor[aPos[((at + 1 - offset) % 3) + offset]  ] - coor[aPos[((at + 2 - offset) % 3) + offset]  ];
		      g[1] = coor[aPos[((at + 1 - offset) % 3) + offset]+1] - coor[aPos[((at + 2 - offset) % 3) + offset]+1];
		      g[2] = coor[aPos[((at + 1 - offset) % 3) + offset]+2] - coor[aPos[((at + 2 - offset) % 3) + offset]+2];
		    }
		    else // atom 2 (5-membered rings)
		    {
		      g[0] = coor[aPos[0]  ] - coor[aPos[1]  ] + coor[aPos[3]  ] - coor[aPos[4]  ];
		      g[1] = coor[aPos[0]+1] - coor[aPos[1]+1] + coor[aPos[3]+1] - coor[aPos[4]+1];
		      g[2] = coor[aPos[0]+2] - coor[aPos[1]+2] + coor[aPos[3]+2] - coor[aPos[4]+2];
		    }
		    ab[0] = d[1] * g[2] - d[2] * g[1];
		    ab[1] = d[2] * g[0] - d[0] * g[2];
		    ab[2] = d[0] * g[1] - d[1] * g[0];
		    c[0] = nSum[1] * g[2] - nSum[2] * g[1];
		    c[1] = nSum[2] * g[0] - nSum[0] * g[2];
		    c[2] = nSum[0] * g[1] - nSum[1] * g[0];
		    
		    factor = -6 * dn / dL3nL3;
		    double factor2 = 0.25 * dL / nL;
		    double OneOverN = 1 / ((double) ringInfo[i].numAtoms);
		    double factor3 = nL / dL * OneOverN;
		    double gradUx = factor * ((0.5 * ab[0] - n[0] * OneOverN) * dLnL - dn * (factor2 * c[0] - factor3 * d[0]));
		    double gradUy = factor * ((0.5 * ab[1] - n[1] * OneOverN) * dLnL - dn * (factor2 * c[1] - factor3 * d[1]));
		    double gradUz = factor * ((0.5 * ab[2] - n[2] * OneOverN) * dLnL - dn * (factor2 * c[2] - factor3 * d[2]));
		    
		    factor = -3 * dL * OneOverN;
		    double gradVx = factor * d[0];
		    double gradVy = factor * d[1];
		    double gradVz = factor * d[2];
		    
		    ff[place+aPos[at]  ] += -fact * (gradUx * v - u * gradVx) / v2;
		    ff[place+aPos[at]+1] += -fact * (gradUy * v - u * gradVy) / v2;
		    ff[place+aPos[at]+2] += -fact * (gradUz * v - u * gradVz) / v2;
	          }
                }
	        // then multiply with appropriate coefficient and sum up
	        for (int i = 0; i < 5; i++) contribTot += rc[i] * contrib[i];
	        cs += contribTot;
	      }
	      //END OF RINGS

	      //DIHEDRAL ANGLES
	      {
	        double *CO_DA = db.CO_DA(aa_kind,at_kind);
		double *PARS_DA;

	        if(atom[s][a].phi.size()==4){

                  Vector d0(coor[CSDIM*atom[s][a].phi[0]  ] - coor[CSDIM*atom[s][a].phi[1]],
		            coor[CSDIM*atom[s][a].phi[0]+1] - coor[CSDIM*atom[s][a].phi[1]+1],
		            coor[CSDIM*atom[s][a].phi[0]+2] - coor[CSDIM*atom[s][a].phi[1]+2]);
                                                                        
                  Vector d1(coor[CSDIM*atom[s][a].phi[1]  ] - coor[CSDIM*atom[s][a].phi[2]],
		            coor[CSDIM*atom[s][a].phi[1]+1] - coor[CSDIM*atom[s][a].phi[2]+1],
		            coor[CSDIM*atom[s][a].phi[1]+2] - coor[CSDIM*atom[s][a].phi[2]+2]);
                                                                        
                  Vector d2(coor[CSDIM*atom[s][a].phi[2]  ] - coor[CSDIM*atom[s][a].phi[3]],
		            coor[CSDIM*atom[s][a].phi[2]+1] - coor[CSDIM*atom[s][a].phi[3]+1],
		            coor[CSDIM*atom[s][a].phi[2]+2] - coor[CSDIM*atom[s][a].phi[3]+2]);

                  Vector dd0,dd1,dd2;
                  Torsion t;
                  double phi = t.compute(d0,d1,d2,dd0,dd1,dd2);

  		  PARS_DA = db.PARS_DA(aa_kind,at_kind,0);
		  cs += CO_DA[0] * ( PARS_DA[0] * cos(3 * phi + PARS_DA[3]) 
		      	            +PARS_DA[1] * cos(phi + PARS_DA[4])+ PARS_DA[2]);

                  double fact = cs_deriv * CO_DA[0] * ( -PARS_DA[0]*3*sin(3*phi+ PARS_DA[3]) 
		      	                                -PARS_DA[1] * sin(phi + PARS_DA[4]));

                  ff[place+CSDIM*atom[s][a].phi[0]  ] += -dd0[0]*fact;
                  ff[place+CSDIM*atom[s][a].phi[0]+1] += -dd0[1]*fact;
                  ff[place+CSDIM*atom[s][a].phi[0]+2] += -dd0[2]*fact;
                  ff[place+CSDIM*atom[s][a].phi[1]  ] += -(dd1[0]-dd0[0])*fact;
                  ff[place+CSDIM*atom[s][a].phi[1]+1] += -(dd1[1]-dd0[1])*fact;
                  ff[place+CSDIM*atom[s][a].phi[1]+2] += -(dd1[2]-dd0[2])*fact;
                  ff[place+CSDIM*atom[s][a].phi[2]  ] += -(dd2[0]-dd1[0])*fact;
                  ff[place+CSDIM*atom[s][a].phi[2]+1] += -(dd2[1]-dd1[1])*fact;
                  ff[place+CSDIM*atom[s][a].phi[2]+2] += -(dd2[2]-dd1[2])*fact;
                  ff[place+CSDIM*atom[s][a].phi[3]  ] += dd2[0]*fact;
                  ff[place+CSDIM*atom[s][a].phi[3]+1] += dd2[1]*fact;
                  ff[place+CSDIM*atom[s][a].phi[3]+2] += dd2[2]*fact;
                }

	        if(atom[s][a].psi.size()==4){
                  Vector d0(coor[CSDIM*atom[s][a].psi[0]  ] - coor[CSDIM*atom[s][a].psi[1]],
		            coor[CSDIM*atom[s][a].psi[0]+1] - coor[CSDIM*atom[s][a].psi[1]+1],
		            coor[CSDIM*atom[s][a].psi[0]+2] - coor[CSDIM*atom[s][a].psi[1]+2]);
                                                                        
                  Vector d1(coor[CSDIM*atom[s][a].psi[1]  ] - coor[CSDIM*atom[s][a].psi[2]],
		            coor[CSDIM*atom[s][a].psi[1]+1] - coor[CSDIM*atom[s][a].psi[2]+1],
		            coor[CSDIM*atom[s][a].psi[1]+2] - coor[CSDIM*atom[s][a].psi[2]+2]);
                                                                        
                  Vector d2(coor[CSDIM*atom[s][a].psi[2]  ] - coor[CSDIM*atom[s][a].psi[3]],
		            coor[CSDIM*atom[s][a].psi[2]+1] - coor[CSDIM*atom[s][a].psi[3]+1],
		            coor[CSDIM*atom[s][a].psi[2]+2] - coor[CSDIM*atom[s][a].psi[3]+2]);

                  Torsion t;
                  Vector dd0,dd1,dd2;
                  double psi = t.compute(d0,d1,d2,dd0,dd1,dd2);

		  PARS_DA = db.PARS_DA(aa_kind,at_kind,1);
		  cs += CO_DA[1] * ( PARS_DA[0] * cos(3 * psi + PARS_DA[3]) 
		      	            +PARS_DA[1] * cos(psi + PARS_DA[4])+PARS_DA[2]);

                  double fact = cs_deriv * CO_DA[1] * ( -PARS_DA[0] * 3 * sin(3 * psi + PARS_DA[3]) 
		      	            -PARS_DA[1] * sin(psi + PARS_DA[4]));

                  ff[place+CSDIM*atom[s][a].psi[0]  ] += -dd0[0]*fact;
                  ff[place+CSDIM*atom[s][a].psi[0]+1] += -dd0[1]*fact;
                  ff[place+CSDIM*atom[s][a].psi[0]+2] += -dd0[2]*fact;
                  ff[place+CSDIM*atom[s][a].psi[1]  ] += -(dd1[0]-dd0[0])*fact;
                  ff[place+CSDIM*atom[s][a].psi[1]+1] += -(dd1[1]-dd0[1])*fact;
                  ff[place+CSDIM*atom[s][a].psi[1]+2] += -(dd1[2]-dd0[2])*fact;
                  ff[place+CSDIM*atom[s][a].psi[2]  ] += -(dd2[0]-dd1[0])*fact;
                  ff[place+CSDIM*atom[s][a].psi[2]+1] += -(dd2[1]-dd1[1])*fact;
                  ff[place+CSDIM*atom[s][a].psi[2]+2] += -(dd2[2]-dd1[2])*fact;
                  ff[place+CSDIM*atom[s][a].psi[3]  ] += dd2[0]*fact;
                  ff[place+CSDIM*atom[s][a].psi[3]+1] += dd2[1]*fact;
                  ff[place+CSDIM*atom[s][a].psi[3]+2] += dd2[2]*fact;
                }

	        //Chi
	        if(atom[s][a].chi1.size()==4){
                  Vector d0(coor[CSDIM*atom[s][a].chi1[0]  ] - coor[CSDIM*atom[s][a].chi1[1]],
		            coor[CSDIM*atom[s][a].chi1[0]+1] - coor[CSDIM*atom[s][a].chi1[1]+1],
		            coor[CSDIM*atom[s][a].chi1[0]+2] - coor[CSDIM*atom[s][a].chi1[1]+2]);
                                                                         
                  Vector d1(coor[CSDIM*atom[s][a].chi1[1]  ] - coor[CSDIM*atom[s][a].chi1[2]],
		            coor[CSDIM*atom[s][a].chi1[1]+1] - coor[CSDIM*atom[s][a].chi1[2]+1],
		            coor[CSDIM*atom[s][a].chi1[1]+2] - coor[CSDIM*atom[s][a].chi1[2]+2]);
                                                                         
                  Vector d2(coor[CSDIM*atom[s][a].chi1[2]  ] - coor[CSDIM*atom[s][a].chi1[3]],
		            coor[CSDIM*atom[s][a].chi1[2]+1] - coor[CSDIM*atom[s][a].chi1[3]+1],
		            coor[CSDIM*atom[s][a].chi1[2]+2] - coor[CSDIM*atom[s][a].chi1[3]+2]);

                  Torsion t;
                  Vector dd0,dd1,dd2;
                  double chi = t.compute(d0,d1,d2,dd0,dd1,dd2);

		  PARS_DA = db.PARS_DA(aa_kind,at_kind,2);
		  cs += CO_DA[2] * ( PARS_DA[0] * cos(3 * chi + PARS_DA[3]) 
		      	            +PARS_DA[1] * cos(chi + PARS_DA[4])+PARS_DA[2]);

                  double fact = cs_deriv* CO_DA[2] * ( -PARS_DA[0] * 3 * sin(3 * chi + PARS_DA[3]) 
		      	            -PARS_DA[1] * sin(chi + PARS_DA[4]));

                  ff[place+CSDIM*atom[s][a].chi1[0]  ] += -dd0[0]*fact;
                  ff[place+CSDIM*atom[s][a].chi1[0]+1] += -dd0[1]*fact;
                  ff[place+CSDIM*atom[s][a].chi1[0]+2] += -dd0[2]*fact;
                  ff[place+CSDIM*atom[s][a].chi1[1]  ] += -(dd1[0]-dd0[0])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[1]+1] += -(dd1[1]-dd0[1])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[1]+2] += -(dd1[2]-dd0[2])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[2]  ] += -(dd2[0]-dd1[0])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[2]+1] += -(dd2[1]-dd1[1])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[2]+2] += -(dd2[2]-dd1[2])*fact;
                  ff[place+CSDIM*atom[s][a].chi1[3]  ] += dd2[0]*fact;
                  ff[place+CSDIM*atom[s][a].chi1[3]+1] += dd2[1]*fact;
                  ff[place+CSDIM*atom[s][a].chi1[3]+2] += dd2[2]*fact;
	        }
	      }
	      //END OF DIHE
              sh[index+a][at_kind] = cs;
	    } 
	  }
	}
        index += atom[s].size();
      }

      ++box_count;
      if(box_count == box_nupdate) box_count = 0;
    }

  private:

    inline void update_box(int curr, vector<int> & aa_box_i, const vector<double> & coor, const double cutnb2){
      aa_box_i.clear();
      int ipos = CSDIM*curr;
      int size = coor.size();
      for(int n=0;n<size;n++){
	int res_dist = abs(res_num[curr]-res_num[n]);
	if(res_dist<2) continue;

	int npos = CSDIM*n;
	double dx,dy,dz;
	dx = coor[ipos  ] - coor[npos];
	dy = coor[ipos+1] - coor[npos+1];
	dz = coor[ipos+2] - coor[npos+2];
	double d2 = dx*dx+dy*dy+dz*dz;
	if(d2<cutnb2) aa_box_i.push_back(n);
      }
    }

    void xdist_name_map(string & name){
      if((name == "OT1")) name = "O";
      else if ((name == "HN") || (name == "HT1")) name = "H";
      else if ((name == "CG1")|| (name == "OG")|| 
	       (name == "SG") || (name == "OG1")) name = "CG";
      else if ((name == "HA1"))                   name = "HA";
    }

    void init_backbone(const Almost::Molecules & molecules){

      seg_last.resize(molecules.protein_size());

      for(unsigned int i=0;i<molecules.protein_size();i++){
	if(i==0)
	  seg_last[i] = molecules.protein(i).fragment_size();
	else 
	  seg_last[i] = seg_last[i-1]+molecules.protein(i).fragment_size();

	const Almost::Protein & p = molecules.protein(i);
	int b = p.atom_offset();
	int e = b+p.atom_size();
	int fb = p.fragment_offset();

	vector<int> N_;
	vector<int> H_;
	vector<int> CA_;
	vector<int> CB_;
	vector<int> HA_;
	vector<int> C_;
	vector<int> O_;
	vector<int> CX_; //For chi1
	{	  
	  N_.resize(p.fragment_size());
	  H_.resize(p.fragment_size());
	  CA_.resize(p.fragment_size());
	  CB_.resize(p.fragment_size());
	  HA_.resize(p.fragment_size());
	  C_.resize(p.fragment_size());
	  O_.resize(p.fragment_size());
	  CX_.resize(p.fragment_size());
	}

	for(unsigned int a=0;a<N_.size();a++){
	  N_[a] = -1;
	  H_[a] = -1;
	  CA_[a] = -1;
	  CB_[a] = -1;
	  HA_[a] = -1;
	  C_[a] = -1;
	  O_[a] = -1;
	  CX_[a] = -1;
	}

	for(int a=b;a<e;a++){
	  if(molecules[a].name()=="N"){
	    int f = molecules.find_fragment(a);
	    N_[f-fb] = a;
	  } else if(molecules[a].name()=="H"||molecules[a].name()=="HN"){
	    int f = molecules.find_fragment(a);
	    H_[f-fb] = a;
	  } else if(molecules[a].name()=="HA"||molecules[a].name()=="HA1"){
	    int f = molecules.find_fragment(a);
	    HA_[f-fb] = a;
	  } else if(molecules[a].name()=="CA"){
	    int f = molecules.find_fragment(a);
	    CA_[f-fb] = a;
	  } else if(molecules[a].name()=="CB"){
	    int f = molecules.find_fragment(a);
	    CB_[f-fb] = a;
	  } else if(molecules[a].name()=="C"){
	    int f = molecules.find_fragment(a);
	    C_[f-fb] = a;
	  } else if(molecules[a].name()=="O"){
	    int f = molecules.find_fragment(a);
	    O_[f-fb] = a;
	  } else {
	    //SPECIAL CASE PRO
	    if(molecules[a].name()=="CD"){
	      int f = molecules.find_fragment(a);
	      if(molecules.fragment_name(f,Almost::Protein::SHORT)=="PRO")
		H_[f-fb] = a;
	    }
	  }

	  //CHI1 SIDE CHAIN
	  {
	    int f = molecules.find_fragment(a);
	    string frg = molecules.fragment_name(f,Almost::Protein::SHORT);
	    string atm = molecules[a].name();
	    if(is_chi1_cx(frg,atm)) CX_[f-fb] = a;
	  }
	}

	vector<Fragment> atm_;
	for(int a=b;a<e;a++){
	  if(molecules[a].name()=="CA"){
	    int fd = p.find_fragment(a);
	    int f_idx = fd - fb;
	    Fragment at;
	    {
	      at.pos[0] = HA_[f_idx];
	      at.pos[1] = H_[f_idx];
	      at.pos[2] = N_[f_idx];
	      at.pos[3] = CA_[f_idx];
	      at.pos[4] = CB_[f_idx];
	      at.pos[5] = C_[f_idx];
	    }
	    at.res_type_prev = at.res_type_curr = at.res_type_next = -1;
	    at.res_name = molecules.fragment_name(fd,Almost::Molecules::SHORT);
	    at.res_kind = db.kind(at.res_name);
	    at.fd = fd;
	    //REGISTER PREV CURR NEXT
	    {
	      //N H CA HA C O

	      //PREV
	      if(f_idx!=0){
		at.prev.push_back(N_[f_idx-1]);
		at.prev.push_back(CA_[f_idx-1]);
		at.prev.push_back(HA_[f_idx-1]);
		at.prev.push_back(C_[f_idx-1]);
		at.prev.push_back(O_[f_idx-1]);	
		at.res_type_prev = frag2enum(molecules.fragment_name(fd-1,Almost::Molecules::SHORT));
	      }

	      //CURR
	      {
		at.curr.push_back(N_[f_idx]);
		at.curr.push_back(H_[f_idx]);
		at.curr.push_back(CA_[f_idx]);
		at.curr.push_back(HA_[f_idx]);
		at.curr.push_back(C_[f_idx]);
		at.curr.push_back(O_[f_idx]);		
		at.res_type_curr = frag2enum(molecules.fragment_name(fd,Almost::Molecules::SHORT));
	      }

	      //NEXT
	      if(f_idx!=p.fragment_size()-1){
		at.next.push_back(N_[f_idx+1]);
		at.next.push_back(H_[f_idx+1]);
		at.next.push_back(CA_[f_idx+1]);
		at.next.push_back(HA_[f_idx+1]);
		at.next.push_back(C_[f_idx+1]);
		at.res_type_next = frag2enum(molecules.fragment_name(fd+1,Almost::Molecules::SHORT));
	      }

	      //PHI | PSI | CH1
	      if(f_idx!=0){
		at.phi.push_back(C_[f_idx-1]);
		at.phi.push_back(N_[f_idx]);
		at.phi.push_back(CA_[f_idx]);
		at.phi.push_back(C_[f_idx]);
	      }
	      
	      if(f_idx!=p.fragment_size()-1){
		at.psi.push_back(N_[f_idx]);
		at.psi.push_back(CA_[f_idx]);
		at.psi.push_back(C_[f_idx]);
		at.psi.push_back(N_[f_idx+1]);	
	      }

	      if(CX_[f_idx]!=-1&&CB_[f_idx]!=-1){
		at.chi1.push_back(N_[f_idx]);
		at.chi1.push_back(CA_[f_idx]);
		at.chi1.push_back(CB_[f_idx]);
		at.chi1.push_back(CX_[f_idx]);	
	      }
	    }
	    atm_.push_back(at);
	  }
	}
	atom.push_back(atm_);
      }
    }

    void init_sidechain(const Almost::Molecules & molecules){
      for(unsigned int s = 0; s<atom.size(); s++){
	for(unsigned int a = 0;a<atom[s].size();a++){
	  vector<int> atm = molecules.fragment_atoms(atom[s][a].fd);
          if(atom[s][a].res_name=="UNK") continue;
	  vector<string> sc_atm = side_chain_atoms(atom[s][a].res_name);

	  for(unsigned int sc=0;sc<sc_atm.size();sc++){
	    for(unsigned int aa=0;aa<atm.size();aa++){
	      if(molecules[atm[aa]].name()==sc_atm[sc])
		atom[s][a].side_chain.push_back(atm[aa]);
	    }
	  }
	}
      }
    }

    void init_xdist(const Almost::Molecules & molecules){
      string atomsP1[] = {"H", "H", "H", "C", "C", "C", 
                          "O", "O", "O", "N", "N", "N", 
                          "O", "O", "O", "N", "N", "N", 
                          "CG", "CG", "CG", "CG", "CG", "CG", "CG", "CA"};

      int resOffsetP1 [] = {0, 0, 0, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, -1};
      string atomsP2[] = {"HA", "C", "CB", "HA", "C", "CB", 
                          "HA", "N", "CB", "HA", "N", "CB", 
                          "HA", "N", "CB", "HA", "N", "CB", 
                          "HA", "N", "C", "C", "N", "CA", "CA", "CA"};

      int resOffsetP2 [] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 1, 0, 0, 1};

      for(unsigned int s = 0; s<atom.size(); s++){
	for(unsigned int a = 0;a<atom[s].size();a++){
	  int fd = atom[s][a].fd;
	  int fo = molecules.protein(s).fragment_offset();
	  int fs = molecules.protein(s).fragment_size();
	  int f_idx = fd-fo;
	  if((f_idx==0)||(f_idx==fs-1)) continue;
	  vector<int> atm_curr = molecules.fragment_atoms(atom[s][a].fd);
	  vector<int> atm_prev = molecules.fragment_atoms(atom[s][a].fd-1);
	  vector<int> atm_next = molecules.fragment_atoms(atom[s][a].fd+1);

	  for(unsigned int q=0;q<numXtraDists - 1; q++){
	    vector<int>::iterator at1,at1_end;
	    vector<int>::iterator at2,at2_end;

	    int p1 = -1;
	    int p2 = -1;

	    if(resOffsetP1[q]== 0){ at1 = atm_curr.begin(); at1_end = atm_curr.end();}
	    if(resOffsetP1[q]==-1){ at1 = atm_prev.begin(); at1_end = atm_prev.end();}
	    if(resOffsetP1[q]==+1){ at1 = atm_next.begin(); at1_end = atm_next.end();}
	    while(at1!=at1_end){
	      int aa = *at1; ++at1;
	      string name = molecules[aa].name();

	      xdist_name_map(name);

	      if(name==atomsP1[q]){
		p1 = aa;
		break;
	      }
	    }

	    if(resOffsetP2[q]== 0){ at2 = atm_curr.begin(); at2_end = atm_curr.end();}
	    if(resOffsetP2[q]==-1){ at2 = atm_prev.begin(); at2_end = atm_prev.end();}
	    if(resOffsetP2[q]==+1){ at2 = atm_next.begin(); at2_end = atm_next.end();}
	    while(at2!=at2_end){
	      int aa = *at2; ++at2;
	      string name = molecules[aa].name();

	      xdist_name_map(name);

	      if(name==atomsP2[q]){
		p2 = aa;
		break;
	      }
	    }
	    atom[s][a].xd1.push_back(p1);
	    atom[s][a].xd2.push_back(p2);
	  }
	}
      }
    }

    void init_types(const Almost::Molecules & molecules){
      for(unsigned int i=0;i<molecules.atom_size();i++){
	int frag = molecules.find_fragment(i);
	string fragName = molecules.fragment_name(frag,Almost::Molecules::SHORT);
	string atom_name = molecules[i].name();
	char atom_type = atom_name[0];
	res_num.push_back(frag);
	int t = -1;
	if (!isSP2(fragName, atom_name)){
	  if (atom_type == 'C') t = D_C;
	  else if (atom_type == 'O') t = D_O;
	  else if (atom_type == 'H') t = D_H;
	  else if (atom_type == 'N') t = D_N;
	  else if (atom_type == 'S') t = D_S;
	}else{
	  if (atom_type == 'C') t = D_C2;
	  else if (atom_type == 'O') t = D_O2;
	  else if (atom_type == 'N') t = D_N2;
	}
	type.push_back(t);
      }
    }

    void init_rings(const Almost::Molecules & m){
      for(unsigned int i=0;i<m.fragment_size();i++){
	string frg = m.fragment_name(i,Almost::Molecules::SHORT);
	if(!((frg=="PHE")||
	     (frg=="TYR")||
	     (frg=="TRP")||
	     (frg=="HIS")||
	     (frg=="HIP")||
	     (frg=="HID")||
	     (frg=="HIE")||
	     (frg=="HSD")||
	     (frg=="HSE")||
	     (frg=="HSP")))
	  continue;

	if(frg=="PHE"||frg=="TYR"){
	  RingInfo ri;
	  ri.type = 0;
	  vector<int> frg_atoms = m.fragment_atoms(i);
	  for(unsigned int a=0;a<frg_atoms.size();a++){
	    for(unsigned int aa=0;aa<6;aa++){
	      int atm = frg_atoms[a];
	      if(m[atm].name()==ri.atomNames[ri.type][aa]){
		ri.atom[aa] = atm;
		break;
	      }
	    }
	  }
	  ri.numAtoms = 6;
	  if(frg=="PHE") ri.rtype = RingInfo::R_PHE;
	  if(frg=="TYR") ri.rtype = RingInfo::R_TYR;
	  ringInfo.push_back(ri);
	}
	else if(frg=="TRP"){
	  //First ring
	  {
	    RingInfo ri;
	    ri.type = 1;
	    vector<int> frg_atoms = m.fragment_atoms(i);
	    for(unsigned int a=0;a<frg_atoms.size();a++){
	      for(unsigned int aa=0;aa<6;aa++){
		int atm = frg_atoms[a];
		if(m[atm].name()==ri.atomNames[ri.type][aa]){
		  ri.atom[aa] = atm;
		  break;
		}
	      }
	    }
	    ri.numAtoms = 6;
	    ri.rtype = RingInfo::R_TRP1;
	    ringInfo.push_back(ri);
	  }
	  //Second Ring
	  {
	    RingInfo ri;
	    ri.type = 2;
	    vector<int> frg_atoms = m.fragment_atoms(i);
	    for(unsigned int a=0;a<frg_atoms.size();a++){
	      for(unsigned int aa=0;aa<5;aa++){
		int atm = frg_atoms[a];
		if(m[atm].name()==ri.atomNames[ri.type][aa]){
		  ri.atom[aa] = atm;
		  break;
		}
	      }
	    }
	    ri.numAtoms = 5;
	    ri.rtype = RingInfo::R_TRP2;
	    ringInfo.push_back(ri);
	  }
	}
	else { //HIS case
	  RingInfo ri;
	  ri.type = 3;
	  vector<int> frg_atoms = m.fragment_atoms(i);
	  for(unsigned int a=0;a<frg_atoms.size();a++){
	    for(unsigned int aa=0;aa<5;aa++){
	      int atm = frg_atoms[a];
	      if(m[atm].name()==ri.atomNames[ri.type][aa]){
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
    }

    void compute_ring_parameters(const vector<double> & coor){
      for(unsigned int i=0;i<ringInfo.size();i++){
	int size = ringInfo[i].numAtoms;
	double a[6][3];
	for(int j = 0; j < size; j++)
	  for(unsigned int d=0;d<3;d++){
	    a[j][d] = coor[CSDIM*ringInfo[i].atom[j]+d];
	  }
	// calculate ring center
	double midP [3];
	for (int j = 0; j < 3; j++){
	  midP[j] = 0.0;
	  for (int k = 0; k < size; k++) midP[j] += a[k][j];
	  ringInfo[i].position[j] = midP[j] /= (double) size;
	}
	// compute normal vector to plane containing first three atoms in array
	ringInfo[i].n1 = xProduct(a[0][0] - a[1][0], a[0][1] - a[1][1], a[0][2] - a[1][2],
				  a[2][0] - a[1][0], a[2][1] - a[1][1], a[2][2] - a[1][2]);
	// compute normal vector to plane containing last three atoms in array
	// NB: third atom of five-membered ring used for both computations above
	ringInfo[i].n2 = xProduct(a[size-3][0] - a[size-2][0], a[size-3][1] - a[size-2][1], a[size-3][2] - a[size-2][2],
				  a[size-1][0] - a[size-2][0], a[size-1][1] - a[size-2][1], a[size-1][2] - a[size-2][2]);
	// ring plane normal vector is average of n1 and n2
	for (int j = 0; j < 3; j++) ringInfo[i].normVect[j] = 0.5*(ringInfo[i].n1[j] + ringInfo[i].n2[j]);
	// calculate squared length and length of normal vector
	ringInfo[i].lengthN2 = ringInfo[i].normVect[0]*ringInfo[i].normVect[0] + 
                               ringInfo[i].normVect[1]*ringInfo[i].normVect[1] + 
                               ringInfo[i].normVect[2]*ringInfo[i].normVect[2];
	ringInfo[i].lengthNV = sqrt(ringInfo[i].lengthN2);
	
      }
    }

    aa_t frag2enum(string aa)
    {
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
      else cerr << "Error converting string " << aa << " into amino acid index: not a valid 3-letter code" << endl;
      return type;
    }

    vector<string> side_chain_atoms(string s){
      vector<string> sc; sc.clear();

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
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
        sc.push_back( "HE" );
        sc.push_back( "HH11" );
        sc.push_back( "HH12" );
        sc.push_back( "HH21" );
        sc.push_back( "HH22" );
	return sc;
      } else if(s=="ASN"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "OD1" );
        sc.push_back( "ND2" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HD21" );
        sc.push_back( "HD22" );
	return sc;
      } else if(s=="ASP"||s=="ASH"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "OD1" );
        sc.push_back( "OD2" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
	return sc;
      } else if(s=="CYS"||s=="CYM"){
        sc.push_back( "CB" );
        sc.push_back( "SG" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
	return sc;
      } else if(s=="GLN"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "CD" );
        sc.push_back( "OE1" );
        sc.push_back( "NE2" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
        sc.push_back( "HE21" );
        sc.push_back( "HE22" );
	return sc;
      } else if(s=="GLU"||s=="GLH"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "CD" );
        sc.push_back( "OE1" );
        sc.push_back( "OE2" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
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
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
        sc.push_back( "HE1" );
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
        sc.push_back( "HG" );
        sc.push_back( "HD11" );
        sc.push_back( "HD12" );
        sc.push_back( "HD13" );
        sc.push_back( "HD21" );
        sc.push_back( "HD22" );
        sc.push_back( "HD23" );
	return sc;
      } else if(s=="LYS"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "CD" );
        sc.push_back( "CE" );
        sc.push_back( "NZ" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
        sc.push_back( "HE1" );
        sc.push_back( "HE2" );
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
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
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
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
        sc.push_back( "HE1" );
        sc.push_back( "HE2" );
        sc.push_back( "HZ" );
	return sc;
      } else if(s=="PRO"){
        sc.push_back( "CB" );
        sc.push_back( "CG" );
        sc.push_back( "CD" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
        sc.push_back( "HG2" );
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
	return sc;
      } else if(s=="SER"){
        sc.push_back( "CB" );
        sc.push_back( "OG" );
        sc.push_back( "HB1" );
        sc.push_back( "HB2" );
        sc.push_back( "HG1" );
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
        sc.push_back( "HD1" );
        sc.push_back( "HD2" );
        sc.push_back( "HE1" );
        sc.push_back( "HE2" );
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
	return sc;
      } else plumed_merror("CS2BackBone: side_chain_atoms unknown");
    }

    bool isSP2(const string & resType, const string & atomName)
    {
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

    bool is_chi1_cx(const string & frg, const string & atm){
      if(atm=="CG")                                        return true;
      if((frg == "CYS")&&(atm =="SG"))                     return true;
      if(((frg == "ILE")||(frg == "VAL"))&&(atm == "CG1")) return true;
      if((frg == "SER")&&(atm == "OG"))                    return true;
      if((frg == "THR")&&(atm == "OG1"))                   return true;

      return false;
    }

    static inline vector<double> xProduct(const double & x1, const double & y1, const double & z1, 
                                          const double & x2, const double & y2, const double & z2)
    {
      vector<double> result;
      result.reserve(3);
      result.push_back(y1 * z2 - z1 * y2);
      result.push_back(z1 * x2 - x1 * z2);
      result.push_back(x1 * y2 - y1 * x2);
      return result;
    }

    int frag_segment(int p){
      int s = 0;
      for(unsigned int i=0;i<seg_last.size()-1;i++){
	if(p>seg_last[i]) s  = i+1;
	else break;
      }
      return s;
    }

    int frag_relitive_index(int p, int s){
      if(s==0) return p;
      return p-seg_last[s-1];
    }
    
  };

}

#endif

#endif
