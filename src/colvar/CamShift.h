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

#ifndef __CAMSHIFT3__
#define __CAMSHIFT3__

#define cutOffDist    0.55  	// cut off distance for non-bonded pairwise forces
#define cutOffDist2   0.3025 	// square of cutOffDist
#define cutOnDist     0.45   	// cut off distance for non-bonded pairwise forces
#define cutOnDist2    0.2025 	// square of cutOffDist
#define invswitch     1.0/((cutOffDist2 - cutOnDist2)*(cutOffDist2 - cutOnDist2)*(cutOffDist2 - cutOnDist2))
#define CSDIM         3

#include <string>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>

#include "tools/OpenMP.h"
#include "tools/PDB.h"
#include "tools/Torsion.h"

using namespace std;

namespace PLMD {

  class CamShiftDB {
    enum { STD, GLY, PRO};
    enum { HA_ATOM, H_ATOM, N_ATOM, CA_ATOM, CB_ATOM, C_ATOM };

    static const unsigned aa_kind = 3;
    static const unsigned atm_kind = 6;
    static const unsigned numXtraDists = 27;

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
    CamShiftDB(const string &file, const double scale){ 
      dscale = scale;
      parse(file);
    }

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
    inline double   CONST(const unsigned a_kind, const unsigned at_kind){ return c_all[a_kind][at_kind];}
    inline double * CONSTAACURR(const unsigned a_kind, const unsigned at_kind){return c_aa[a_kind][at_kind];}
    inline double * CONSTAANEXT(const unsigned a_kind, const unsigned at_kind){return c_aa_succ[a_kind][at_kind];}
    inline double * CONSTAAPREV(const unsigned a_kind, const unsigned at_kind){return c_aa_prev[a_kind][at_kind];}
    inline double * CONST_BB2_PREV(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind][1];}
    inline double * CONST_BB2_CURR(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind][1]+5;}
    inline double * CONST_BB2_NEXT(const unsigned a_kind, const unsigned at_kind){return co_bb[a_kind][at_kind][1]+11;}
    inline double * CONST_SC2(const unsigned a_kind, const unsigned at_kind, unsigned res_type){ return co_sc_[a_kind][at_kind][res_type][1];}
    inline double * CONST_XD(const unsigned a_kind, const unsigned at_kind){ return co_xd[a_kind][at_kind];}
    inline double * CO_SPHERE(const unsigned a_kind, const unsigned at_kind, unsigned exp_type){ return co_sphere[a_kind][at_kind][exp_type];}
    inline double * CO_RING(const unsigned a_kind, const unsigned at_kind){ return co_ring[a_kind][at_kind];}
    inline double * CO_DA(const unsigned a_kind, const unsigned at_kind){ return co_da[a_kind][at_kind];}
    inline double * PARS_DA(const unsigned a_kind, const unsigned at_kind, const unsigned ang_kind){ return pars_da[a_kind][at_kind][ang_kind];}

  private:

    void parse(const string &file){
      ifstream in;
      in.open(file.c_str());
      if(!in) plumed_merror("Unable to open CamShiftDB file\n");

      unsigned c_kind = 0;
      unsigned c_atom = 0;
      unsigned nline = 0;

      for(unsigned i=0;i<3;i++) for(unsigned j=0;j<6;j++) {
        c_all[i][j]=0.;
        for(unsigned k=0;k<20;k++) {
          c_aa[i][j][k]=0.;
          c_aa_prev[i][j][k]=0.;
          c_aa_succ[i][j][k]=0.;
          for(unsigned l=0;l<2;l++) for(unsigned m=0;m<20;m++) co_sc_[i][j][k][l][m]=0.;
        }
        for(unsigned k=0;k<16;k++) {co_bb[i][j][0][k]=0.; co_bb[i][j][1][k]=0.;}
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

	for(unsigned scC = 0; scC < 20; scC++){
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

  class CamShift3
  {
    struct Fragment {
      //Order HA H N CA CB C
      vector<int> pos;
      vector<int> prev;
      vector<int> curr;
      vector<int> next;

      unsigned res_type_prev;
      unsigned res_type_curr;
      unsigned res_type_next;

      vector<int> side_chain;
      vector<int> xd1;
      vector<int> xd2;
      string res_name;
      unsigned res_kind;
      unsigned fd;

      vector<vector<unsigned> > box_nb;

      vector<int> phi;
      vector<int> psi;
      vector<int> chi1;

      vector<double> exp_cs;
 
      Fragment() {
	res_type_prev = res_type_curr = res_type_next = 0;
        res_kind = fd = 0;
        box_nb.resize(6); 
        pos.resize(6); 
        exp_cs.resize(6); 
        for(unsigned i=0;i<6;i++) {
           box_nb[i].resize(50);
           pos[i] = -1;
           exp_cs[i] = 0.;
        }
      }
    };

    struct RingInfo{
      enum {R_PHE, R_TYR, R_TRP1, R_TRP2, R_HIS};
      // one out of five different types
      unsigned rtype;
      // up to six member per ring
      unsigned atom[6];
      // number of ring members (5 or 6)
      unsigned numAtoms;
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
    };

    enum aa_t {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL, UNK};
    enum atom_t {D_C, D_H, D_N, D_O, D_S, D_C2, D_N2, D_O2};
    CamShiftDB db;
    vector<unsigned> seg_last; // i<seg_last[0] is i n 0...
    vector<unsigned> type;
    vector<unsigned> res_num;
    vector<RingInfo> ringInfo;
    unsigned box_nupdate;
    unsigned box_count;

  public:

    vector<vector<Fragment> > atom;
  
    CamShift3(const string &pdbfile, const string &dbfile, const bool NatUnits, const double scale):db(dbfile,scale){
      PDB pdb;
      if( !pdb.read(pdbfile,NatUnits,1./scale) ) plumed_merror("missing input file " + pdbfile );
      init_backbone(pdb);
      init_sidechain(pdb);
      init_xdist(pdb);
      init_types(pdb);
      init_rings(pdb);
#ifndef NDEBUG
      debug_report();
#endif
      box_nupdate=1;
      box_count=0;
    }

    void set_box_nupdate(const unsigned v){box_nupdate = v;}
    void set_box_count(const unsigned v){box_count = v;}
   
    void remove_problematic(const string &res, const string &nucl) {
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

    void read_cs(const string &file, const string &k){
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
	unsigned p = atoi(tok.c_str());
	p = p - 1;
	unsigned seg = frag_segment(p);
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

    void calc_cs(const vector<double> & coor, vector<double> & ff, const unsigned N,  double **sh){

      compute_ring_parameters(coor);

      bool update = false;
      if(box_count==0) update = true;

      unsigned index=0;
      // CYCLE OVER MULTIPLE CHAINS
      for(unsigned s=0;s<atom.size();s++){
	// SKIP FIRST AND LAST RESIDUE OF EACH CHAIN
#pragma omp parallel for num_threads(OpenMP::getNumThreads())
	for(unsigned a=0;a<atom[s].size();a++){
          // CYCLE OVER THE SIX BACKBONE CHEMICAL SHIFTS
	  for(unsigned at_kind=0;at_kind<6;at_kind++){
	    double cs = 0.;
	    if(atom[s][a].pos[at_kind]>0&&atom[s][a].exp_cs[at_kind]>0){

	      double cs_deriv = -1.;

              // this is the counter to find your place in ff
              // place is residue*6*CSDIM*N + at_kind*CSDIM*N
              // residue is atoms[s-1].size()+a
              //int place = (iamresidue+a)*6*CSDIM*N+at_kind*CSDIM*N;
              unsigned place = (index+a)*6*CSDIM*N+at_kind*CSDIM*N;

	      unsigned aa_kind = atom[s][a].res_kind;
	      unsigned res_type_curr = atom[s][a].res_type_curr;
	      unsigned res_type_prev = atom[s][a].res_type_prev;
	      unsigned res_type_next = atom[s][a].res_type_next;

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
	      unsigned ipos = CSDIM*atom[s][a].pos[at_kind];

	      //PREV
	      for(unsigned q=0;q<atom[s][a].prev.size();q++){
	        double d2,dx,dy,dz;
	        unsigned jpos = CSDIM*atom[s][a].prev[q];
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
	      for(unsigned q=0;q<atom[s][a].curr.size();q++){
	        unsigned jpos = CSDIM*atom[s][a].curr[q];
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
	      for(unsigned q=0;q<atom[s][a].next.size();q++){
	        double d2,dx,dy,dz;
	        unsigned jpos = CSDIM*atom[s][a].next[q];
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
	      for(unsigned q=0;q<atom[s][a].side_chain.size();q++){
	        unsigned jpos = CSDIM*atom[s][a].side_chain[q];
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
	      for(unsigned q=0;q<atom[s][a].xd1.size();q++){
	        double d2,dx,dy,dz;
	        if(atom[s][a].xd1[q]==-1||atom[s][a].xd2[q]==-1) continue;
	        unsigned kpos = CSDIM*atom[s][a].xd1[q];
	        unsigned jpos = CSDIM*atom[s][a].xd2[q];
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

	        unsigned curr = atom[s][a].pos[at_kind];
                if(update) update_box(atom[s][a].box_nb[at_kind], curr, coor, N, 0.56);
  	        double ffi[3] = {0,0,0};
                 
                for(unsigned bat = 0; bat<atom[s][a].box_nb[at_kind].size(); bat++) {
                  unsigned at = atom[s][a].box_nb[at_kind][bat];

		  unsigned jpos = CSDIM*at;
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
		    unsigned t = type[at];
		    dist_sum3[t] += factor3;
		    dist_sum [t] += factor1;

		    double fact1 = -cs_deriv*dfactor1*dinv;
		    double fact2 = -cs_deriv*dfactor3*dinv;
		    double fact = fact1*CONST_CO_SPHERE[t]+fact2*CONST_CO_SPHERE3[t];
		    ffi[0] = ffi[0] + fact*dx;
		    ffi[1] = ffi[1] + fact*dy;
		    ffi[2] = ffi[2] + fact*dz;
		    ff[place+jpos]   = ff[place+jpos]   - fact*dx;
		    ff[place+jpos+1] = ff[place+jpos+1] - fact*dy;
		    ff[place+jpos+2] = ff[place+jpos+2] - fact*dz;
                  }
	        }
	        ff[place+ipos] = ff[place+ipos] + ffi[0];
	        ff[place+ipos+1] = ff[place+ipos+1] + ffi[1];
	        ff[place+ipos+2] = ff[place+ipos+2] + ffi[2];

	        for(unsigned tt=0;tt<8;tt++){
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
	        for(unsigned i=0; i<ringInfo.size(); i++){
		  // compute angle from ring middle point to current atom position
		  // get distance vector from query atom to ring center and normal vector to ring plane
  		  double fact = cs_deriv*rc[ringInfo[i].rtype];
	          double d[3], n[3];
		  for(unsigned j=0; j<3; j++){
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
		
  		  unsigned aPos[6];
  		  for(unsigned j=0; j<ringInfo[i].numAtoms; j++) aPos[j] = CSDIM*ringInfo[i].atom[j];

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
		  unsigned limit = ringInfo[i].numAtoms - 3; // 2 for a 5-membered ring, 3 for a 6-membered ring
		
		  // update forces on ring atoms
		  for(unsigned at=0; at<ringInfo[i].numAtoms; at++)
		  {
		    if (at < limit) // atoms 0,1 (5 member) or 0,1,2 (6 member)
		    {
		      g[0] = coor[aPos[(at+1)%3]  ] - coor[aPos[(at+2)%3]  ];
		      g[1] = coor[aPos[(at+1)%3]+1] - coor[aPos[(at+2)%3]+1];
		      g[2] = coor[aPos[(at+1)%3]+2] - coor[aPos[(at+2)%3]+2];
		    }
		    else if (at >= ringInfo[i].numAtoms - limit) // atoms 3,4 (5 member) or 3,4,5 (6 member)
		    {
		      unsigned offset = ringInfo[i].numAtoms - 3; // 2 for a 5-membered ring, 3 for a 6-membered ring
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
	        for(unsigned i=0; i<5; i++) contribTot += rc[i]*contrib[i];
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
	    } 
            sh[index+a][at_kind] = cs;
	  }
	}
        index += atom[s].size();
      }

      ++box_count;
      if(box_count == box_nupdate) box_count = 0;
    }

  private:

    void init_backbone(const PDB &pdb){

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
	  }
	  atm_.push_back(at);
	}
	atom.push_back(atm_);
      }
    }

    void init_sidechain(const PDB &pdb){
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

    void init_xdist(const PDB &pdb){
      const string atomsP1[] = {"H", "H", "H", "C", "C", "C", 
                                "O", "O", "O", "N", "N", "N", 
                                "O", "O", "O", "N", "N", "N", 
                                "CG", "CG", "CG", "CG", "CG", 
                                "CG", "CG", "CA"};

      const int resOffsetP1 [] = {0, 0, 0, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, -1};

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

    void init_types(const PDB &pdb){
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

    void init_rings(const PDB &pdb){

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
            ri.n1.resize(3);
            ri.n2.resize(3);
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
            ri.n1.resize(3);
            ri.n2.resize(3);
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
            ri2.n1.resize(3);
            ri2.n2.resize(3);
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
            ri.n1.resize(3);
            ri.n2.resize(3);
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

    void compute_ring_parameters(const vector<double> & coor){
      for(unsigned i=0;i<ringInfo.size();i++){
	unsigned size = ringInfo[i].numAtoms;
	double a[6][3];
	for(unsigned j=0; j<size; j++) {
          a[j][0] = coor[CSDIM*ringInfo[i].atom[j]  ];
          a[j][1] = coor[CSDIM*ringInfo[i].atom[j]+1];
          a[j][2] = coor[CSDIM*ringInfo[i].atom[j]+2];
        }
	// calculate ring center
	double midP[3];
	for(unsigned j=0; j<3; j++){
	  midP[j] = a[0][j];
	  for(unsigned k=1; k<size; k++) midP[j] += a[k][j];
          midP[j] /= (double) size;
	  ringInfo[i].position[j] = midP[j];
	}
	// compute normal vector to plane containing first three atoms in array
	ringInfo[i].n1 = xProduct(a[0][0] - a[1][0], a[0][1] - a[1][1], a[0][2] - a[1][2],
				  a[2][0] - a[1][0], a[2][1] - a[1][1], a[2][2] - a[1][2]);
	// compute normal vector to plane containing last three atoms in array
	// NB: third atom of five-membered ring used for both computations above
	ringInfo[i].n2 = xProduct(a[size-3][0] - a[size-2][0], a[size-3][1] - a[size-2][1], a[size-3][2] - a[size-2][2],
				  a[size-1][0] - a[size-2][0], a[size-1][1] - a[size-2][1], a[size-1][2] - a[size-2][2]);
	// ring plane normal vector is average of n1 and n2
	ringInfo[i].normVect[0] = 0.5*(ringInfo[i].n1[0] + ringInfo[i].n2[0]);
	ringInfo[i].normVect[1] = 0.5*(ringInfo[i].n1[1] + ringInfo[i].n2[1]);
	ringInfo[i].normVect[2] = 0.5*(ringInfo[i].n1[2] + ringInfo[i].n2[2]);
	// calculate squared length and length of normal vector
	ringInfo[i].lengthN2 = ringInfo[i].normVect[0]*ringInfo[i].normVect[0] + 
                               ringInfo[i].normVect[1]*ringInfo[i].normVect[1] + 
                               ringInfo[i].normVect[2]*ringInfo[i].normVect[2];
	ringInfo[i].lengthNV = sqrt(ringInfo[i].lengthN2);
      }
    }

    aa_t frag2enum(const string &aa)
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
      else plumed_merror("CS2Backbone: Error converting string " + aa + " into amino acid index: not a valid 3-letter code");
      return type;
    }

    vector<string> side_chain_atoms(const string &s){
      vector<string> sc;
      sc.clear();

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
      result.push_back(y1 * z2 - z1 * y2);
      result.push_back(z1 * x2 - x1 * z2);
      result.push_back(x1 * y2 - y1 * x2);
      return result;
    }

    unsigned frag_segment(const unsigned p){
      unsigned s = 0;
      for(unsigned i=0;i<seg_last.size()-1;i++){
	if(p>seg_last[i]) s  = i+1;
	else break;
      }
      return s;
    }

    unsigned frag_relitive_index(const unsigned p, const unsigned s){
      if(s==0) return p;
      return p-seg_last[s-1];
    }

    void debug_report(){
      printf("\t CamShift3 Initialization report: \n");
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

    inline void update_box(vector<unsigned> & aa_box_i, const unsigned curr, const vector<double> & coor, const unsigned size, const double cutnb2){
      aa_box_i.clear();
      unsigned ipos = CSDIM*curr;
      unsigned res_curr = res_num[curr];
      for(unsigned n=0;n<size;n++){
	unsigned res_dist = abs(static_cast<int>(res_curr-res_num[n]));
	if(res_dist<2) continue;

	unsigned npos = CSDIM*n;
	double dx = coor[ipos  ] - coor[npos];
	double dy = coor[ipos+1] - coor[npos+1];
	double dz = coor[ipos+2] - coor[npos+2];
	double d2 = dx*dx+dy*dy+dz*dz;
	if(d2<cutnb2) aa_box_i.push_back(n);
      }
    }

    void xdist_name_map(string & name){
      if((name == "OT1")||(name == "OC1")) name = "O";
      else if ((name == "HN") || (name == "HT1") || (name == "H1")) name = "H";
      else if ((name == "CG1")|| (name == "OG")|| 
	       (name == "SG") || (name == "OG1")) name = "CG";
      else if ((name == "HA1"))                   name = "HA";
    }

  };

}

#endif
