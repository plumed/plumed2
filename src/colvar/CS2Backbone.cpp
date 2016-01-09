/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "CamShift.h"

using namespace std;

namespace PLMD{
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

class CS2Backbone : public Colvar {
  vector<CamShift3> cam_list;
  unsigned numResidues;
  bool pbc;
  bool noexp;
  double **sh;
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
PLUMED_COLVAR_INIT(ao),
pbc(true),
noexp(false)
{
  string stringadb;
  string stringapdb;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parseFlag("NOEXP",noexp);

  string stringa_data;
  parse("DATA",stringa_data);

  string stringa_template;
  parse("TEMPLATE",stringa_template);

  int neigh_f=10;
  parse("NEIGH_FREQ", neigh_f);

  parse("NRES", numResidues);

  stringadb  = stringa_data + string("/camshift.db");
  stringapdb = stringa_data + string("/") + stringa_template;

  /* Lenght conversion (parameters are tuned for angstrom) */
  double scale=1.;
  if(!plumed.getAtoms().usingNaturalUnits()) {
    scale = 10.*atoms.getUnits().getLength();
  }

  log.printf("  Initialization of the predictor ...\n"); log.flush();
  CamShift3 a = CamShift3(stringapdb, stringadb, plumed.getAtoms().usingNaturalUnits(), scale);

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
  a.remove_problematic("CYS", "HA");
  a.remove_problematic("CYS", "H");
  a.remove_problematic("CYS", "N");
  a.remove_problematic("CYS", "CA");
  a.remove_problematic("CYS", "CB");
  a.remove_problematic("CYS", "C");
  /* done */

  log.printf("  Setting parameters ...\n"); log.flush();
  
  a.set_box_nupdate(neigh_f);

  cam_list.push_back(a);

  sh = new double*[numResidues];
  sh[0] = new double[numResidues*6];
  for(unsigned i=1;i<numResidues;i++)  sh[i]=sh[i-1]+6;

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

  if(getExchangeStep()) cam_list[0].set_box_count(0);

#pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0;i<N;i++) {
     unsigned ipos = CSDIM*i;
     Vector Pos = getPosition(i);
     coor[ipos]   = Pos[0];
     coor[ipos+1] = Pos[1];
     coor[ipos+2] = Pos[2];
  }

  cam_list[0].calc_cs(coor, csforces, N, sh);

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
            Vector For(csforces[ipos], csforces[ipos+1], csforces[ipos+2]);
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
}
