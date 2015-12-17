/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#ifdef __PLUMED_HAS_ALMOST

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "CamShift.h"

#include <almost/mdb.h>
#include <almost/pdb.h>

using namespace std;
using namespace Almost;

string CamShift3::RingInfo::atomNames[4][6];
string CamShift3::RingInfo::types[4];
int    CamShift3::RingInfo::init_ = 0;

namespace PLMD{

//+PLUMEDOC COLVAR CS2BACKBONE 
/*
This collective variable calculates a scoring function based on the comparison of backcalculated and
experimental backbone chemical shifts for a protein (CA, CB, C', H, HA, N).

CamShift \cite Kohlhoff:2009us is employed to back calculate the chemical shifts that are then compared
with a set of experimental values to generate a score \cite Robustelli:2010dn \cite Granata:2013dk.

It is also possible to back-calculate the chemical shifts from multiple replicas and then average them
to perform Replica-Averaged Restrained MD simulations \cite Camilloni:2012je \cite Camilloni:2013hs.

In general the system for which chemical shifts are to be calculated must be completly included in
ATOMS. It should also be made whole \ref WHOLEMOLECULES before the the back-calculation. 

HOW TO COMPILE IT

\ref installingalmost on how to compile PLUMED with ALMOST.

HOW TO USE IT

To use CamShift a set of experimental data is needed. CamShift uses backbone and Cb chemical shifts 
that must be provided as text files:

\verbatim
CAshifts.dat:
CBshifts.dat:
Cshifts.dat:
Hshifts.dat:
HAshifts.dat:
Nshifts.dat:
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

All of them must always be there. If a chemical shift for an atom of a residue is not available 0.0 must be used. 
So if for example all the Cb are not available all the chemical shifts for all the residues should be set as 0.0.

A template.pdb file is needed to the generate a topology of the protein within ALMOST. For histidines in protonation 
states different from D the HIE/HIP name should be used in the template.pdb. GLH and ASH can be used for the alternative 
protonation of GLU and ASP. Non-standard amino acids and other molecules are not yet supported! If multiple chains are 
present the chain identifier must be in the standard PDB format, together with the TER keyword at the end of each chain.
Residues numbering should always go from 1 to N in both the chemical shifts files as well as in the template.pdb file.
Two more keywords can be used to setup the topology: CYS-DISU to tell ALMOST to look for disulphide bridges and TERMINI
to define the protonation state of the the chain's termini (currently only DEFAULT (NH3+, CO2-) and NONE (NH, CO)).

Two more standard files are also needed: an ALMOST force-field file, corresponding to the force-field family used in
the MD code, (a03_cs2_gromacs.mdb for the amber family or all22_gromacs.mdb for the charmm family) and camshift.db, 
both these files can be copied from almost/branches/almost-2.1/toppar.

All the above files must be in a single folder that must be specified with the keyword DATA. 

Additional material and examples can be also found in the tutorial \ref belfast-9 

\par Examples

case 1:

\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs: CS2BACKBONE ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=0.0 NRES=13 ENSEMBLE
cse: RESTRAINT ARG=cs SLOPE=24 KAPPA=0 AT=0.
PRINT ARG=cs,cse.bias
\endverbatim

case 2:

\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs: CS2BACKBONE ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=1.0 NRES=13 TERMINI=DEFAULT,NONE CYS-DISU WRITE_CS=1000
PRINT ARG=cs
\endverbatim

(See also \ref WHOLEMOLECULES, \ref RESTRAINT and \ref PRINT)

*/
//+ENDPLUMEDOC

class CS2Backbone : public Colvar {
  vector<CamShift3> cam_list;
  Molecules molecules;
  unsigned  numResidues;
  bool serial;
  bool pbc;
  bool noexp;
  double **sh;
  double len_pl2alm;
  double for_pl2alm;
  vector<double> coor;
  vector<double> csforces;
public:
  explicit CS2Backbone(const ActionOptions&);
  ~CS2Backbone();
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(CS2Backbone,"CS2BACKBONE")

void CS2Backbone::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  keys.addFlag("PARALLEL",false,"Perform the calculation in parallel - (experimental, bad scaling).");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","DATA","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","FF","a03_gromacs.mdb","The ALMOST force-field to map the atoms' names.");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file of the protein system to initialise ALMOST.");
  keys.add("compulsory","NEIGH_FREQ","10","Period in step for neighbour list update.");
  keys.add("compulsory","NRES","Number of residues, corresponding to the number of chemical shifts.");
  keys.add("optional","TERMINI","Defines the protonation states of the chain-termini.");
  keys.addFlag("CYS-DISU",false,"Set to TRUE if your system has disulphide bridges.");  
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
PLUMED_COLVAR_INIT(ao),
serial(true),
pbc(true),
noexp(false)
{
  string stringadb;
  string stringamdb;
  string stringapdb;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parseFlag("NOEXP",noexp);

  bool parallel=false;
  parseFlag("PARALLEL",parallel);
  serial=!parallel;

  string stringa_data;
  parse("DATA",stringa_data);

  string stringa_forcefield;
  parse("FF",stringa_forcefield);

  string stringa_template;
  parse("TEMPLATE",stringa_template);

  bool disu=false;
  parseFlag("CYS-DISU",disu);

  int neigh_f=10;
  parse("NEIGH_FREQ", neigh_f);

  parse("NRES", numResidues);

  stringadb  = stringa_data + string("/camshift.db");
  stringamdb = stringa_data + string("/") + stringa_forcefield;
  stringapdb = stringa_data + string("/") + stringa_template;

  log.printf("  loading force-field %s\n", stringamdb.c_str()); log.flush();
  Almost::MDB mdb((char*)stringamdb.c_str());
  log.printf("  loading template %s\n", stringapdb.c_str()); log.flush();
  Almost::PDB pdb((char*)stringapdb.c_str());

  vector<string> termini;
  string stringa_mol;
  parse("TERMINI",stringa_mol);
  if(stringa_mol.length()>0) {
    unsigned num_chains = pdb[0].size();
    vector<string> data=Tools::getWords(stringa_mol,",");
    if(data.size()!=2*num_chains) error("You have to define both the NTerm and the CTerm for each chain of your system!\n");
    for(unsigned i=0;i<data.size();i++) termini.push_back(data[i]);
  } else {
    unsigned num_chains = pdb[0].size();
    for(unsigned i=0;i<(2*num_chains);i++) termini.push_back("DEFAULT");
  }

  log.printf("  building molecule ..."); log.flush();
  for(unsigned i=0;i<pdb[0].size();i++){
    unsigned j=2*i;
    string str;
    str +='A'+i;
    Protein p(str);
    p.build_missing(pdb[0][i],mdb,termini[j],termini[j+1]);
    if(disu) p.auto_disu_bonds(2.9,mdb);
    molecules.add_protein(p);
  }
  log.printf(" done!\n"); log.flush(); 
  log.printf("  Writing converted template.pdb ...\n"); log.flush();
  mol2pdb(molecules,"converted-template.pdb");

  log.printf("  Initialization of the predictor ...\n"); log.flush();
  CamShift3 a = CamShift3(molecules, stringadb);

  log.printf("  Reading experimental data ...\n"); log.flush();
  stringadb = stringa_data + string("/CAshifts.dat");
  log.printf("  Initializing CA shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "CA");
  stringadb = stringa_data + string("/CBshifts.dat");
  log.printf("  Initializing CB shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "CB");
  stringadb = stringa_data + string("/Cshifts.dat");
  log.printf("  Initializing C' shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "C");
  stringadb = stringa_data + string("/HAshifts.dat");
  log.printf("  Initializing HA shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "HA");
  stringadb = stringa_data + string("/Hshifts.dat");
  log.printf("  Initializing H shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "H");
  stringadb = stringa_data + string("/Nshifts.dat");
  log.printf("  Initializing N shifts %s\n", stringadb.c_str()); log.flush();
  a.read_cs(stringadb, "N");
  /* this is a workaround for those chemical shifts that can result in too large forces */
  a.remove_problematic("GLN", "CB");
  a.remove_problematic("ILE", "CB");
  a.remove_problematic("PRO", "N");  
  a.remove_problematic("PRO", "H");
  a.remove_problematic("PRO", "CB");
  a.remove_problematic("GLY", "HA"); a.remove_problematic("GLY", "CB");
  /* this is a workaround for those chemical shifts that are not parameterized */
  a.remove_problematic("HIE", "HA"); a.remove_problematic("HIP", "HA"); a.remove_problematic("HSP", "HA");
  a.remove_problematic("HIE", "H");  a.remove_problematic("HIP", "H");  a.remove_problematic("HSP", "H"); 
  a.remove_problematic("HIE", "N");  a.remove_problematic("HIP", "N");  a.remove_problematic("HSP", "N"); 
  a.remove_problematic("HIE", "CA"); a.remove_problematic("HIP", "CA"); a.remove_problematic("HSP", "CA");
  a.remove_problematic("HIE", "CB"); a.remove_problematic("HIP", "CB"); a.remove_problematic("HSP", "CB");
  a.remove_problematic("HIE", "C");  a.remove_problematic("HIP", "C");  a.remove_problematic("HSP", "C"); 
  a.remove_problematic("GLH", "HA"); a.remove_problematic("ASH", "HA"); a.remove_problematic("HSE", "HA");
  a.remove_problematic("GLH", "H");  a.remove_problematic("ASH", "H");  a.remove_problematic("HSE", "H");
  a.remove_problematic("GLH", "N");  a.remove_problematic("ASH", "N");  a.remove_problematic("HSE", "N");
  a.remove_problematic("GLH", "CA"); a.remove_problematic("ASH", "CA"); a.remove_problematic("HSE", "CA");
  a.remove_problematic("GLH", "CB"); a.remove_problematic("ASH", "CB"); a.remove_problematic("HSE", "CB");
  a.remove_problematic("GLH", "C");  a.remove_problematic("ASH", "C");  a.remove_problematic("HSE", "C");
  a.remove_problematic("UNK", "HA");
  a.remove_problematic("UNK", "H");
  a.remove_problematic("UNK", "N");
  a.remove_problematic("UNK", "CA");
  a.remove_problematic("UNK", "CB");
  a.remove_problematic("UNK", "C");
  if(disu) { 
    a.remove_problematic("CYS", "HA");
    a.remove_problematic("CYS", "H");
    a.remove_problematic("CYS", "N");
    a.remove_problematic("CYS", "CA");
    a.remove_problematic("CYS", "CB");
    a.remove_problematic("CYS", "C");
  }
  /* done */

  log.printf("  Setting parameters ...\n"); log.flush();
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {stride=1; rank=0;}
  if(stride>1) log.printf("  Parallelized over %u processors\n", stride);
  a.set_mpi(stride, rank);
  
  a.set_box_nupdate(neigh_f);
  cam_list.push_back(a);

  sh = new double*[numResidues];
  sh[0] = new double[numResidues*6];
  for(unsigned i=1;i<numResidues;i++)  sh[i]=sh[i-1]+6;

  /* Energy and Lenght conversion */
  len_pl2alm = 10.00*plumed.getAtoms().getUnits().getLength();
  for_pl2alm = len_pl2alm;

  log.printf("  Conversion table from plumed to Almost:\n");
  log.printf("    Length %lf\n", len_pl2alm);
  log.printf("    Force %lf\n", for_pl2alm);

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite("Kohlhoff K, Robustelli P, Cavalli A, Salvatella A, Vendruscolo M, J. Am. Chem. Soc. 131, 13894 (2009)")
     <<plumed.cite("Camilloni C, Robustelli P, De Simone A, Cavalli A, Vendruscolo M, J. Am. Chem. Soc. 134, 3968 (2012)") <<"\n";


  unsigned k=0;
  for(unsigned i=0;i<cam_list[0].atom.size();i++) {
    for(unsigned a=0;a<cam_list[0].atom[i].size();a++) {
      std::string num; Tools::convert(k,num);
      if(cam_list[0].atom[i][a].exp_cs[0]>0) { 
        addComponentWithDerivatives("ha_"+num); componentIsNotPeriodic("ha_"+num); 
        if(!noexp) { 
          addComponent("expha_"+num); componentIsNotPeriodic("expha_"+num);
          Value* comp=getPntrToComponent("expha_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[0]);
        }
      }
      if(cam_list[0].atom[i][a].exp_cs[1]>0) { 
        addComponentWithDerivatives("hn_"+num); componentIsNotPeriodic("hn_"+num); 
        if(!noexp) { 
          addComponent("exphn_"+num); componentIsNotPeriodic("exphn_"+num);
          Value* comp=getPntrToComponent("exphn_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[1]);
        }
      }
      if(cam_list[0].atom[i][a].exp_cs[2]>0) { 
        addComponentWithDerivatives("nh_"+num); componentIsNotPeriodic("nh_"+num); 
        if(!noexp) { 
          addComponent("expnh_"+num); componentIsNotPeriodic("expnh_"+num);
          Value* comp=getPntrToComponent("expnh_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[2]);
        }
      }
      if(cam_list[0].atom[i][a].exp_cs[3]>0) { 
        addComponentWithDerivatives("ca_"+num); componentIsNotPeriodic("ca_"+num); 
        if(!noexp) { 
          addComponent("expca_"+num); componentIsNotPeriodic("expca_"+num);
          Value* comp=getPntrToComponent("expca_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[3]);
        }
      }
      if(cam_list[0].atom[i][a].exp_cs[4]>0) { 
        addComponentWithDerivatives("cb_"+num); componentIsNotPeriodic("cb_"+num); 
        if(!noexp) { 
          addComponent("expcb_"+num); componentIsNotPeriodic("expcb_"+num);
          Value* comp=getPntrToComponent("expcb_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[4]);
        }
      }
      if(cam_list[0].atom[i][a].exp_cs[5]>0) { 
        addComponentWithDerivatives("co_"+num); componentIsNotPeriodic("co_"+num);
        if(!noexp) { 
          addComponent("expco_"+num); componentIsNotPeriodic("expco_"+num);
          Value* comp=getPntrToComponent("expco_"+num); comp->set(cam_list[0].atom[i][a].exp_cs[5]);
        }
      }
      k++;
    }
  }

  coor.resize(CSDIM*atoms.size()); 
  csforces.resize(6*CSDIM*numResidues*atoms.size());

  requestAtoms(atoms);
  log.printf("  DONE!\n"); log.flush();
}

CS2Backbone::~CS2Backbone()
{
  delete[] sh[0];
  delete[] sh;
}

void CS2Backbone::calculate()
{
  unsigned N = getNumberOfAtoms();

  if(pbc) makeWhole();

  for(unsigned i=0;i<numResidues;i++) for(unsigned j=0;j<6;j++) sh[i][j]=0.;

  if(getExchangeStep()) cam_list[0].set_box_count(0);

  for (unsigned i=0;i<N;i++) {
     unsigned ipos = CSDIM*i;
     Vector Pos = getPosition(i);
     coor[ipos]   = len_pl2alm*Pos[0];
     coor[ipos+1] = len_pl2alm*Pos[1];
     coor[ipos+2] = len_pl2alm*Pos[2];
  }

  cam_list[0].new_calc_cs(coor, csforces, N, sh);

  if(!serial) {
    comm.Sum(&sh[0][0], numResidues*6);
    comm.Sum(csforces);
  }

  unsigned k=0;
  unsigned step=2;
  if(noexp) step=1;
 
  for(unsigned j=0;j<numResidues;j++) {
    unsigned placeres = CSDIM*N*6*j;
    for(unsigned cs=0;cs<6;cs++) {
      if(sh[j][cs]!=0.) {
        Value* comp=getPntrToComponent(k);
        comp->set(sh[j][cs]);
        unsigned place = placeres+cs*CSDIM*N;
        Tensor virial; 
        for(unsigned i=0;i<N;i++) {
          unsigned ipos = place+CSDIM*i;
          if(csforces[ipos]!=0||csforces[ipos+1]!=0||csforces[ipos+2]!=0) {
            Vector For(for_pl2alm*csforces[ipos], for_pl2alm*csforces[ipos+1], for_pl2alm*csforces[ipos+2]);
            csforces[ipos]=csforces[ipos+1]=csforces[ipos+2]=0.;
            setAtomsDerivatives(comp,i,For);
            virial-=Tensor(getPosition(i),For);
          }
        }
        setBoxDerivatives(comp,virial);
        k+=step;
      }
    }
  }

} 

}
#endif
