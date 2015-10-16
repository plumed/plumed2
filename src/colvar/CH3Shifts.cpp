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

#include <almost/mdb.h>
#include <almost/pdb.h>
#include <almost/camshift3/meth/methcs.h>
#include <almost/io/formostream.h>

using namespace std;
using namespace Almost;

namespace PLMD{

//+PLUMEDOC COLVAR CH3SHIFTS 
/*
This collective variable calculates a scoring function based on the comparison of calculated and
experimental methyl chemical shifts. 


CH3Shift \cite Sahakyan:2011bn is employed to back calculate the chemical shifts of methyl groups
(ALA:HB; ILE:HD,HG2; LEU:HD1,HD2; THR:HG2; VAL:HG1,HG2) that are then compared with a set of experimental 
values to generate a score \cite Robustelli:2010dn \cite Granata:2013dk.

It is also possible to backcalculate the chemical shifts from multiple replicas and then average them
to perform Replica-Averaged Restrained MD simulations \cite Camilloni:2012je \cite Camilloni:2013hs.

In general the system for which chemical shifts are to be calculated must be completly included in
ATOMS. It should also be made whole \ref WHOLEMOLECULES before the the back-calculation. 

HOW TO COMPILE IT

\ref installingalmost on how to compile PLUMED with ALMOST.

HOW TO USE IT

CH3Shift reads from a text file the experimental chemical shifts:

\verbatim
CH3shifts.dat:
1.596 28
0.956 46
0.576 3 HG2
0.536 3 HD1
0.836 13 HG2
0.666 13 HD1
0.716 23 HG2
0.506 23 HD1
\endverbatim

A template.pdb file is needed to the generate a topology of the protein within ALMOST. For histidines in protonation 
states different from D the HIE/HIP name should be used in the template.pdb. GLH and ASH can be used for the alternative 
protonation of GLU and ASP. Non-standard amino acids and other molecules are not yet supported! If multiple chains are 
present the chain identifier must be in the standard PDB format, together with the TER keyword at the end of each chain.
Residues numbering should always go from 1 to N in both the chemical shifts files as well as in the template.pdb file.
Two more keywords can be used to setup the topology: CYS-DISU to tell ALMOST to look for disulphide bridges and TERMINI
to define the protonation state of the the chain's termini (currently only DEFAULT (NH3+, CO2-) and NONE (NH, CO)).

One more standard file is also needed, that is an ALMOST force-field file, corresponding to the force-field family used in
the MD code, (a03_cs2_gromacs.mdb for the amber family or all22_gromacs.mdb for the charmm family).

All the above files must be in a single folder that must be specified with the keyword DATA (multiple definition of the 
CV can point to different folders). 

\par Examples

case 1:

\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs: CH3SHIFTS ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=0.0 NRES=13 ENSEMBLE
cse: RESTRAINT ARG=cs SLOPE=24 KAPPA=0 AT=0.
PRINT ARG=cs,cse.bias
\endverbatim

case 2:

\verbatim
WHOLEMOLECULES ENTITY0=1-174
cs: CH3SHIFTS ATOMS=1-174 DATA=data/ FF=a03_gromacs.mdb FLAT=1.0 NRES=13 TERMINI=DEFAULT,NONE CYS-DISU WRITE_CS=1000
PRINT ARG=cs
\endverbatim

(See also \ref WHOLEMOLECULES, \ref RESTRAINT and \ref PRINT)

*/
//+ENDPLUMEDOC

class CH3Shifts : public Colvar {
  vector<MethCS*> meth_list;
  Molecules molecules;
  int  numResidues;
  unsigned pperiod;
  unsigned ens_dim;
  bool ensemble;
  bool serial;
  double **sh;
  double ene_pl2alm;
  double len_pl2alm;
  double for_pl2alm;
public:
  explicit CH3Shifts(const ActionOptions&);
  ~CH3Shifts();
  static void registerKeywords( Keywords& keys );
  virtual void calculate();
};

PLUMED_REGISTER_ACTION(CH3Shifts,"CH3SHIFTS")

void CH3Shifts::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  //keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose.");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory","DATA","data/","The folder with the experimental chemical shifts.");
  keys.add("compulsory","FF","a03_gromacs.mdb","The ALMOST force-field to map the atoms' names.");
  keys.add("compulsory","FLAT","1.0","Flat region in the scoring function.");
  keys.add("compulsory","NEIGH_FREQ","10","Period in step for neighbour list update.");
  keys.add("compulsory","WRITE_CS","0","Write chemical shifts period.");
  keys.add("compulsory","NRES","Number of residues, corresponding to the number of chemical shifts.");
  keys.add("optional","TERMINI","Defines the protonation states of the chain-termini.");
  keys.addFlag("CYS-DISU",false,"Set to TRUE if your system has disulphide bridges.");  
  keys.addFlag("ENSEMBLE",false,"Set to TRUE if you want to average over multiple replicas.");  
  keys.remove("NOPBC");
}

CH3Shifts::CH3Shifts(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  string stringadb;
  string stringamdb;
  string stringapdb;

  serial=false;
  //parseFlag("SERIAL",serial);

  string stringa_data;
  parse("DATA",stringa_data);

  string stringa_forcefield;
  parse("FF",stringa_forcefield);

  bool disu=false;
  parseFlag("CYS-DISU",disu);

  double grains=1.0;
  parse("FLAT", grains);

  int neigh_f=10;
  parse("NEIGH_FREQ", neigh_f);

  unsigned w_period=0;
  parse("WRITE_CS", w_period);
  pperiod=w_period;

  parse("NRES", numResidues);

  ensemble=false;
  parseFlag("ENSEMBLE",ensemble);
  if(ensemble){
    if(comm.Get_rank()==0) {
      if(multi_sim_comm.Get_size()<2) error("You CANNOT run Replica-Averaged simulations without running multiple replicas!\n");
      ens_dim=multi_sim_comm.Get_size();
    } else ens_dim=0;
    comm.Sum(&ens_dim, 1);
  } else ens_dim=1;

  stringadb  = stringa_data + string("/CH3shifts.dat");
  stringamdb = stringa_data + string("/") + stringa_forcefield;
  stringapdb = stringa_data + string("/template.pdb");

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
  MethCS* a = new MethCS(molecules);

  log.printf("  Reading experimental data ...\n"); log.flush();
  a->read_cs(stringadb);
  log.printf("  Setting parameters ...\n"); log.flush();
  /*unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {stride=1; rank=0;}
  if(stride>1) log.printf("  Parallelized over %d processors\n", stride);
  a.set_mpi(stride, rank);*/
  
  if(ensemble) { log.printf("  ENSEMBLE averaging over %u replicas\n", ens_dim); }
  a->set_w_cs(1);
  a->set_flat_bottom_const(grains);
  a->set_box_nupdate(neigh_f);
  a->set_box_cutnb(11.); // cut-off for neigh-list
  meth_list.push_back(a);

  sh = new double*[numResidues];
  sh[0] = new double[numResidues*8];
  for (int i=1; i<numResidues; i++)  sh[i]=sh[i-1]+8; 

  /* Energy and Lenght conversion */
  ene_pl2alm = 4.186/plumed.getAtoms().getUnits().getEnergy();
  len_pl2alm = 10.00*plumed.getAtoms().getUnits().getLength();
  for_pl2alm = ene_pl2alm*len_pl2alm;
  log.printf("  Conversion table from plumed to Almost:\n");
  log.printf("    Energy %f\n", ene_pl2alm);
  log.printf("    Length %f\n", len_pl2alm);

  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite("Sahakyan AB, Vranken WF, Cavalli A, Vendruscolo M, J. Biomol. NMR 50, 331 (2011)")
     <<plumed.cite("Camilloni C, Robustelli P, De Simone A, Cavalli A, Vendruscolo M, J. Am. Chem. Soc. 134, 3968 (2012)") <<"\n";

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  log.printf("  DONE!\n"); log.flush();
}

CH3Shifts::~CH3Shifts()
{
  delete[] sh[0];
  delete[] sh;
}

void CH3Shifts::calculate()
{
  double energy=0.;
  Tensor virial;
  virial.zero();
  vector<Vector> deriv(getNumberOfAtoms());
  int N = getNumberOfAtoms();
  Coor<double> coor(N); 
  Coor<double> forces(N);

  forces.clear();
  for(int i=0; i<numResidues; i++) for(unsigned j=0; j<6; j++) sh[i][j]=0.;

  for (int i = 0; i < N; i++) {
     int ipos = 4 * i;
     Vector Pos = getPosition(i);
     coor.coor[ipos]   = len_pl2alm*Pos[0];
     coor.coor[ipos+1] = len_pl2alm*Pos[1];
     coor.coor[ipos+2] = len_pl2alm*Pos[2];
  }

  bool printout=false;
  string csfile;
  if(pperiod>0&&comm.Get_rank()==0) printout = (!(getStep()%pperiod));
  if(printout) {char tmp1[21]; sprintf(tmp1, "%ld", getStep()); csfile = string("cs")+getLabel()+"-"+tmp1+string(".dat");}

  double fact=1.0;
  if(!ensemble) { 
     energy = meth_list[0]->calc_cs_force(coor, forces);
     if(printout) meth_list[0]->write_cs(csfile.c_str());
  } else {
     meth_list[0]->calc_cs(coor);
     if(printout) meth_list[0]->write_cs(csfile.c_str());

     unsigned size = meth_list[0]->ala_calc_hb.size();
     for(unsigned j=0;j<size;j++) sh[0][j] = meth_list[0]->ala_calc_hb[j];
     size = meth_list[0]->ile_calc_hd.size();
     for(unsigned j=0;j<size;j++) sh[1][j] = meth_list[0]->ile_calc_hd[j];
     size = meth_list[0]->ile_calc_hg2.size();
     for(unsigned j=0;j<size;j++) sh[2][j] = meth_list[0]->ile_calc_hg2[j];
     size = meth_list[0]->leu_calc_hd1.size();
     for(unsigned j=0;j<size;j++) sh[3][j] = meth_list[0]->leu_calc_hd1[j];
     size = meth_list[0]->leu_calc_hd2.size();
     for(unsigned j=0;j<size;j++) sh[4][j] = meth_list[0]->leu_calc_hd2[j];
     size = meth_list[0]->thr_calc_hg2.size();
     for(unsigned j=0;j<size;j++) sh[5][j] = meth_list[0]->thr_calc_hg2[j];
     size = meth_list[0]->val_calc_hg1.size();
     for(unsigned j=0;j<size;j++) sh[6][j] = meth_list[0]->val_calc_hg1[j];
     size = meth_list[0]->val_calc_hg2.size();
     for(unsigned j=0;j<size;j++) sh[7][j] = meth_list[0]->val_calc_hg2[j];
     fact = 1./((double) ens_dim);
     if(comm.Get_rank()==0) { // I am the master of my replica
       // among replicas
       multi_sim_comm.Sum(&sh[0][0], numResidues*8);
       for(unsigned i=0;i<8;i++) for(int j=0;j<numResidues;j++) sh[j][i] *= fact; 
     } else for(unsigned i=0;i<8;i++) for(int j=0;j<numResidues;j++) sh[j][i] = 0.;
     // inside each replica
     comm.Sum(&sh[0][0], numResidues*8);
     // now send the averaged shifts back to almost
     size = meth_list[0]->ala_calc_hb.size();
     for(unsigned j=0;j<size;j++)  meth_list[0]->ala_calc_hb[j] = sh[0][j];
     size = meth_list[0]->ile_calc_hd.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->ile_calc_hd[j] = sh[1][j];
     size = meth_list[0]->ile_calc_hg2.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->ile_calc_hg2[j] = sh[2][j];
     size = meth_list[0]->leu_calc_hd1.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->leu_calc_hd1[j] = sh[3][j];
     size = meth_list[0]->leu_calc_hd2.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->leu_calc_hd2[j] = sh[4][j];
     size = meth_list[0]->thr_calc_hg2.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->thr_calc_hg2[j] = sh[5][j];
     size = meth_list[0]->val_calc_hg1.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->val_calc_hg1[j] = sh[6][j];
     size = meth_list[0]->val_calc_hg2.size();
     for(unsigned j=0;j<size;j++) meth_list[0]->val_calc_hg2[j] = sh[7][j];
     // calculate all the forces now
     energy = meth_list[0]->ens_calc_cs_force(coor, forces);
  }

  for (int i = 0; i < N; i++)
  {
    Vector For;
    int ipos = 4 * i;
    For[0] = forces.coor[ipos];
    For[1] = forces.coor[ipos+1];
    For[2] = forces.coor[ipos+2];
    deriv[i] = fact*for_pl2alm*For;
    virial=virial+(-1.*Tensor(getPosition(i),deriv[i]));
  }

  for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (ene_pl2alm*energy);
  setBoxDerivatives  (virial);
}

}
#endif
