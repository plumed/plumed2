/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2022 The plumed team
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
/*
 This class was originally written by Alexander Jussupow
 Arrayfire implementation by Alexander Jussupow and CC
 Extension for the middleman algorithm (now removed) by Max Muehlbauer
 Refactoring for hySAXS Martini beads structure factor for Nucleic Acids by Cristina Paissoni
*/

#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"

#include <map>
#include <iterator>
#include <iostream>

#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SAXS
/*
Calculates SAXS scattered intensity using either the Debye equation.

Intensities are calculated for a set of scattering length set using QVALUE keywords that are numbered starting from 0.
Structure factors can be either assigned using a polynomial expansion to any order using the PARAMETERS keywords;
automatically assigned to atoms using the ATOMISTIC flag reading a PDB file, a correction for the water density is
automatically added, with water density that by default is 0.334 but that can be set otherwise using SOLVDENS;
automatically assigned to Martini pseudo atoms using the MARTINI flag.
The calculated intensities can be scaled using the SCALEINT keywords. This is applied by rescaling the structure factors.
Experimental reference intensities can be added using the EXPINT keywords.
By default SAXS is calculated using Debye on CPU, by adding the GPU flag it is possible to solve the equation on a GPU
if the ARRAYFIRE libraries are installed and correctly linked.
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

\par Examples
in the following example the SAXS intensities are calculated using the single bead per residue approximation. structure factors
are obtained from the pdb file indicated in the MOLINFO.

\plumedfile
MOLINFO STRUCTURE=template.pdb

SAXS ...
LABEL=SAXS
ATOMS=1-355
ONEBEAD
QVALUE1=0.02 EXPINT1=1.0902
QVALUE2=0.05 EXPINT2=0.790632
QVALUE3=0.08 EXPINT3=0.453808
QVALUE4=0.11 EXPINT4=0.254737
QVALUE5=0.14 EXPINT5=0.154928
QVALUE6=0.17 EXPINT6=0.0921503
QVALUE7=0.2 EXPINT7=0.052633
QVALUE8=0.23 EXPINT8=0.0276557
QVALUE9=0.26 EXPINT9=0.0122775
QVALUE10=0.29 EXPINT10=0.00880634
QVALUE11=0.32 EXPINT11=0.0137301
QVALUE12=0.35 EXPINT12=0.0180036
QVALUE13=0.38 EXPINT13=0.0193374
QVALUE14=0.41 EXPINT14=0.0210131
QVALUE15=0.44 EXPINT15=0.0220506
... SAXS

PRINT ARG=(SAXS\.q-.*),(SAXS\.exp-.*) FILE=colvar STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC

class SAXS :
  public MetainferenceBase
{
private:
  enum { H, C, N, O, P, S, NTT };
  enum { ALA_BB, ARG_BB, ARG_SC1, ARG_SC2, ASN_BB, ASN_SC1, ASP_BB, ASP_SC1, CYS_BB, CYS_SC1,
         GLN_BB, GLN_SC1, GLU_BB, GLU_SC1, GLY_BB, HIS_BB, HIS_SC1, HIS_SC2, HIS_SC3, ILE_BB,
         ILE_SC1, LEU_BB, LEU_SC1, LYS_BB, LYS_SC1, LYS_SC2, MET_BB, MET_SC1, PHE_BB, PHE_SC1,
         PHE_SC2, PHE_SC3, PRO_BB, PRO_SC1, SER_BB, SER_SC1, THR_BB, THR_SC1, TRP_BB, TRP_SC1,
         TRP_SC2, TRP_SC3, TRP_SC4, TYR_BB, TYR_SC1, TYR_SC2, TYR_SC3, VAL_BB, VAL_SC1, A_BB1,
         A_BB2, A_BB3, A_SC1, A_SC2, A_SC3, A_SC4, A_3TE, A_5TE, A_TE3, A_TE5, C_BB1, C_BB2,
         C_BB3, C_SC1, C_SC2, C_SC3, C_3TE, C_5TE, C_TE3, C_TE5, G_BB1, G_BB2, G_BB3, G_SC1,
         G_SC2, G_SC3, G_SC4, G_3TE, G_5TE, G_TE3, G_TE5, U_BB1, U_BB2, U_BB3, U_SC1, U_SC2,
         U_SC3, U_3TE, U_5TE, U_TE3, U_TE5, DA_BB1, DA_BB2, DA_BB3, DA_SC1, DA_SC2, DA_SC3,
         DA_SC4, DA_3TE, DA_5TE, DA_TE3, DA_TE5, DC_BB1, DC_BB2, DC_BB3, DC_SC1, DC_SC2, DC_SC3,
         DC_3TE, DC_5TE, DC_TE3, DC_TE5, DG_BB1, DG_BB2, DG_BB3, DG_SC1, DG_SC2, DG_SC3, DG_SC4,
         DG_3TE, DG_5TE, DG_TE3, DG_TE5, DT_BB1, DT_BB2, DT_BB3, DT_SC1, DT_SC2, DT_SC3, DT_3TE,
         DT_5TE, DT_TE3, DT_TE5, NMARTINI
       };
  enum { TRP, TYR, PHE, HIS, HIP, HIE, ARG, LYS, CYS, ASP, GLU, ILE, LEU,
         MET, ASN, PRO, GLN, SER, THR, VAL, ALA, GLY, NONEBEAD
       };
  bool pbc;
  bool serial;
  bool gpu;
  bool onebead;
  bool isFirstStep;
  int  deviceid;
  unsigned nres;
  std::vector<unsigned> atoi;
  std::vector<unsigned> atoms_per_bead;
  std::vector<double>   atoms_masses;
  std::vector<double>   q_list;
  std::vector<double>   FF_rank;
  std::vector<std::vector<double> > FF_value_vacuum;
  std::vector<std::vector<double> > FF_value_water;
  std::vector<std::vector<double> > FF_value_mixed;
  std::vector<std::vector<double> > FF_value;
  std::vector<std::vector<float> >  FFf_value;

  std::vector<std::vector<double> > LCPOparam;
  std::vector<unsigned> residue_atom;

  double rho, rho_corr, sasa_cutoff;
  unsigned solv_stride;
  std::vector<double> Iq0_vac;
  std::vector<double> Iq0_wat;
  std::vector<double> Iq0_mix;
  double Iq0;

  void calculate_gpu(std::vector<Vector> &pos, std::vector<Vector> &deriv);
  void calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv);
  void getMartiniSFparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter);
  void getOnebeadparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac, std::vector<std::vector<long double> > &parameter_mix, std::vector<std::vector<long double> > &parameter_wat);
  void getOnebeadMapping(const std::vector<AtomNumber> &atoms);
  double calculateASF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double rho);
  std::map<std::string, std::vector<double> > setupLCPOparam();
  void readLCPOparam(const std::vector<std::vector<std::string> > &AtomResidueName, unsigned natoms);
  void calcNlist(std::vector<std::vector<int> > &Nlist);
  void sasa_calculate(std::vector<bool> &solv_res);

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXS(const ActionOptions&);
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(SAXS,"SAXS")

void SAXS::registerKeywords(Keywords& keys) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","DEVICEID","0","Identifier of the GPU to be used");
  keys.addFlag("GPU",false,"calculate SAXS using ARRAYFIRE on an accelerator device");
  keys.addFlag("ATOMISTIC",false,"calculate SAXS for an atomistic model");
  keys.addFlag("MARTINI",false,"calculate SAXS for a Martini model");
  keys.addFlag("ONEBEAD",false,"calculate SAXS for a single bead model");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("numbered","QVALUE","Selected scattering lengths in Angstrom are given as QVALUE1, QVALUE2, ... .");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the structure factor for the \\f$i\\f$th atom/bead.");
  keys.add("compulsory","SOLVDENS","0.334","Density of the water to be used for the correction of atomistic structure factors. (ONEBEAD only)");
  keys.add("compulsory","SOLVATION_CORRECTION","0.0","Hydration layer electron density correction.");
  keys.add("compulsory","SASA_CUTOFF","1.0","SASA value to consider a residue as exposed to the solvent.");
  keys.add("numbered","EXPINT","Add an experimental value for each q value.");
  keys.add("compulsory","SOLVATION_STRIDE","100","Number of steps between every new check of the residues solvation via LCPO estimate.");
  keys.addOutputComponent("q","default","the # SAXS of q");
  keys.addOutputComponent("exp","EXPINT","the # experimental intensity");
}

SAXS::SAXS(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  pbc(true),
  serial(false),
  gpu(false),
  onebead(false),
  isFirstStep(true),
  deviceid(0)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc)      log.printf("  using periodic boundary conditions\n");
  else         log.printf("  without periodic boundary conditions\n");

  parseFlag("GPU",gpu);
#ifndef  __PLUMED_HAS_ARRAYFIRE
  if(gpu) error("To use the GPU mode PLUMED must be compiled with ARRAYFIRE");
#endif

  parse("DEVICEID",deviceid);
#ifdef  __PLUMED_HAS_ARRAYFIRE
  if(gpu) {
    af::setDevice(deviceid);
    af::info();
  }
#endif

  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  if(atomistic) log.printf("  using atomistic structure factors\n");
  bool martini=false;
  parseFlag("MARTINI",martini);
  if(martini) log.printf("  using Martini structure factors\n");
  onebead=false;
  parseFlag("ONEBEAD",onebead);
  if(martini) log.printf("  using Single Bead structure factors\n");

  if(martini&&atomistic) error("You cannot use MARTINI and ATOMISTIC at the same time");
  if(martini&&onebead) error("You cannot use MARTINI and ONEBEAD at the same time");
  if(onebead&&atomistic) error("You cannot use ONEBEAD and ATOMISTIC at the same time");

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) break;
    if(t_list<=0.) error("QVALUE cannot be less or equal to zero!\n");
    if(onebead&&t_list>0.3) error("For ONEBEAD mapping QVALUE must be smaller or equal to 0.3.");
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;

  for(unsigned i=0; i<numq; ++i) {
    if(q_list[i]==0.) error("it is not possible to set q=0\n");
    if(i>0&&q_list[i]<q_list[i-1]) error("QVALUE must be in ascending order");
    log.printf("  my q: %lf \n",q_list[i]);
  }

  rho = 0.334;
  parse("SOLVDENS", rho);
  log.printf("  Solvent density of %lf\n", rho);

  solv_stride = 100;
  parse("SOLVATION_STRIDE", solv_stride);
  if(onebead) log.printf("  SASA calculated every %u steps\n", solv_stride);

  double correction = 0.00;
  parse("SOLVATION_CORRECTION", correction);
  if(correction>0&&!onebead) error("SOLVATION_CORRECTION can only be used with ONEBEAD");
  rho_corr=rho-correction;
  if(onebead) log.printf("  SASA Solvation density correction set to: %lf\n", rho_corr);

  sasa_cutoff = 1.0;
  parse("SASA_CUTOFF", sasa_cutoff);
  if(sasa_cutoff <= 0.) error("SASA_CUTOFF must be greater than 0.");

  // Here we perform the preliminary mapping for onebead representation
  if(onebead) {
    LCPOparam.resize(size);
    getOnebeadMapping(atoms);
    Iq0_vac.resize(nres);
    Iq0_wat.resize(nres);
    Iq0_mix.resize(nres);
    atoi.resize(nres);
  } else {
    atoi.resize(size);
  }

  Iq0=0;
  std::vector<std::vector<long double> > FF_tmp;
  std::vector<std::vector<long double> > FF_tmp_vac;
  std::vector<std::vector<long double> > FF_tmp_mix;
  std::vector<std::vector<long double> > FF_tmp_wat;
  std::vector<std::vector<long double> > parameter;
  if(!atomistic&&!martini&&!onebead) {
    //read in parameter std::vector
    parameter.resize(size);
    ntarget=0;
    for(unsigned i=0; i<size; ++i) {
      if( !parseNumberedVector( "PARAMETERS", i+1, parameter[i]) ) break;
      ntarget++;
    }
    if( ntarget!=size ) error("found wrong number of parameter std::vectors");
    FF_tmp.resize(numq,std::vector<long double>(size));
    for(unsigned i=0; i<size; ++i) {
      atoi[i]=i;
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[i][0];
    Iq0 *= Iq0;
  } else if(onebead) {
    //read in parameter std::vector
    FF_tmp_vac.resize(numq,std::vector<long double>(NONEBEAD));
    FF_tmp_mix.resize(numq,std::vector<long double>(NONEBEAD));
    FF_tmp_wat.resize(numq,std::vector<long double>(NONEBEAD));
    std::vector<std::vector<long double> > parameter_vac(NONEBEAD);
    std::vector<std::vector<long double> > parameter_mix(NONEBEAD);
    std::vector<std::vector<long double> > parameter_wat(NONEBEAD);
    getOnebeadparam(atoms, parameter_vac, parameter_mix, parameter_wat);
    for(unsigned i=0; i<NONEBEAD; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter_vac[i].size(); ++j) {
          FF_tmp_vac[k][i]+= parameter_vac[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
        for(unsigned j=0; j<parameter_mix[i].size(); ++j) {
          FF_tmp_mix[k][i]+= parameter_mix[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
        for(unsigned j=0; j<parameter_wat[i].size(); ++j) {
          FF_tmp_wat[k][i]+= parameter_wat[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<nres; ++i) {
      Iq0_vac[i]=parameter_vac[atoi[i]][0];
      Iq0_mix[i]=parameter_mix[atoi[i]][0];
      Iq0_wat[i]=parameter_wat[atoi[i]][0];
    }
  } else if(martini) {
    //read in parameter std::vector
    FF_tmp.resize(numq,std::vector<long double>(NMARTINI));
    parameter.resize(NMARTINI);
    getMartiniSFparam(atoms, parameter);
    for(unsigned i=0; i<NMARTINI; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) Iq0+=parameter[atoi[i]][0];
    Iq0 *= Iq0;
  } else if(atomistic) {
    FF_tmp.resize(numq,std::vector<long double>(NTT));
    Iq0=calculateASF(atoms, FF_tmp, rho);
    Iq0 *= Iq0;
  }

  std::vector<double> expint;
  expint.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) break;
    ntarget++;
  }
  bool exp=false;
  if(ntarget!=numq && ntarget!=0) error("found wrong number of EXPINT values");
  if(ntarget==numq) exp=true;
  if(getDoScore()&&!exp) error("with DOSCORE you need to set the EXPINT values");

  if(!gpu) {
    FF_rank.resize(numq);
    unsigned n_atom_types;
    if(onebead) {
      FF_value.resize(nres,std::vector<double>(numq));
      n_atom_types=NONEBEAD;
      FF_value_vacuum.resize(n_atom_types,std::vector<double>(numq));
      FF_value_water.resize(n_atom_types,std::vector<double>(numq));
      FF_value_mixed.resize(n_atom_types,std::vector<double>(numq));
    } else {
      FF_value.resize(size,std::vector<double>(numq));
    }
    for(unsigned k=0; k<numq; ++k) {
      if(!onebead) {
        for(unsigned i=0; i<size; ++i) FF_value[i][k] = static_cast<double>(FF_tmp[k][atoi[i]])/(std::sqrt(Iq0));
        for(unsigned i=0; i<size; ++i) FF_rank[k] += FF_value[i][k]*FF_value[i][k];
      } else {
        for(unsigned i=0; i<n_atom_types; ++i) {
          FF_value_vacuum[i][k] = static_cast<double>(FF_tmp_vac[k][i]);
          FF_value_mixed[i][k] = static_cast<double>(FF_tmp_mix[k][i]);
          FF_value_water[i][k] = static_cast<double>(FF_tmp_wat[k][i]);
        }
      }
    }
  } else {
    unsigned n_atom_types;
    if(onebead) {
      FFf_value.resize(numq,std::vector<float>(nres));
      n_atom_types=NONEBEAD;
      FF_value_vacuum.resize(n_atom_types,std::vector<double>(numq));
      FF_value_water.resize(n_atom_types,std::vector<double>(numq));
      FF_value_mixed.resize(n_atom_types,std::vector<double>(numq));
    } else {
      FFf_value.resize(numq,std::vector<float>(size));
    }
    for(unsigned k=0; k<numq; ++k) {
      if(!onebead) {
        for(unsigned i=0; i<size; ++i) {
          FFf_value[k][i] = static_cast<float>(FF_tmp[k][atoi[i]])/(std::sqrt(Iq0));
        }
      } else {
        for(unsigned i=0; i<n_atom_types; ++i) {
          FF_value_vacuum[i][k] = static_cast<double>(FF_tmp_vac[k][i]);
          FF_value_mixed[i][k] = static_cast<double>(FF_tmp_mix[k][i]);
          FF_value_water[i][k] = static_cast<double>(FF_tmp_wat[k][i]);
        }
      }
    }
  }

  if(!getDoScore()) {
    for(unsigned i=0; i<numq; ++i) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    if(exp) {
      for(unsigned i=0; i<numq; ++i) {
        std::string num; Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(expint[i]);
      }
    }
  } else {
    for(unsigned i=0; i<numq; ++i) {
      std::string num; Tools::convert(i,num);
      addComponent("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    for(unsigned i=0; i<numq; ++i) {
      std::string num; Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(expint[i]);
    }
  }

  // convert units to nm^-1
  for(unsigned i=0; i<numq; ++i) {
    q_list[i]=q_list[i]*10.0;    //factor 10 to convert from A^-1 to nm^-1
  }
  log<<"  Bibliography ";
  if(martini) {
    log<<plumed.cite("Niebling, Björling, Westenhoff, J Appl Crystallogr 47, 1190–1198 (2014).");
    log<<plumed.cite("Paissoni, Jussupow, Camilloni, J Appl Crystallogr 52, 394-402 (2019).");
  }
  if(atomistic) {
    log<<plumed.cite("Fraser, MacRae, Suzuki, J. Appl. Crystallogr., 11, 693–694 (1978).");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography C, 554–595 (International Union of Crystallography, 2006).");
  }

  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)");
  log<<"\n";

  requestAtoms(atoms, false);

  if(getDoScore()) {
    setParameters(expint);
    Initialise(numq);
  }
  setDerivatives();
  checkRead();
}

//calculates SASA neighbor list
void SAXS::calcNlist(std::vector<std::vector<int> > &Nlist)
{
  unsigned natoms = getNumberOfAtoms();
  for(unsigned i = 0; i < natoms; ++i) {
    if (LCPOparam[i].size()>0) {
      for (unsigned j = 0; j < i; ++j) {
        if (LCPOparam[j].size()>0) {
          double Delta_ij_mod = modulo(delta(getPosition(i), getPosition(j)))*10.;
          double overlapD = LCPOparam[i][0]+LCPOparam[j][0];
          if(Delta_ij_mod < overlapD) {
            Nlist.at(i).push_back(j);
            Nlist.at(j).push_back(i);
          }
        }
      }
    }
  }

}

//calculates SASA according to LCPO algorithm
void SAXS::sasa_calculate(std::vector<bool> &solv_res)
{
  unsigned natoms = getNumberOfAtoms();
  std::vector<std::vector<int> > Nlist(natoms);
  calcNlist(Nlist);
  std::vector<double> sasares(nres, 0.);

  for(unsigned i = 0; i < natoms; ++i) {
    if(LCPOparam[i].size()>1) {
      if(LCPOparam[i][1]>0.0) {
        double Aij = 0.0;
        double Aijk = 0.0;
        double Ajk = 0.0;
        double ri = LCPOparam[i][0];
        double S1 = 4.*M_PI*ri*ri;
        for (unsigned j = 0; j < Nlist[i].size(); ++j) {
          double d_ij = modulo(delta( getPosition(i), getPosition(Nlist[i][j]) ))*10.;
          double rj = LCPOparam[Nlist[i][j]][0];
          double Aijt = (2.*M_PI*ri*(ri-d_ij/2.-((ri*ri-rj*rj)/(2.*d_ij))));
          double Ajkt = 0.0;
          for (unsigned k = 0; k < Nlist[Nlist[i][j]].size(); ++k) {
            if (std::find (Nlist[i].begin(), Nlist[i].end(), Nlist[Nlist[i][j]][k]) !=  Nlist[i].end()) {
              double d_jk = modulo(delta( getPosition(Nlist[i][j]), getPosition(Nlist[Nlist[i][j]][k]) ))*10.;
              double rk = LCPOparam[Nlist[Nlist[i][j]][k]][0];
              double sjk =  (2.*M_PI*rj*(rj-d_jk/2.-((rj*rj-rk*rk)/(2.*d_jk))));
              Ajkt += sjk;
            }
          }
          Aijk += (Aijt * Ajkt);
          Aij += Aijt;
          Ajk += Ajkt;
        }
        double sasai = (LCPOparam[i][1]*S1+LCPOparam[i][2]*Aij+LCPOparam[i][3]*Ajk+LCPOparam[i][4]*Aijk);
        if (sasai > 0 ) {
          sasares[residue_atom[i]] += sasai/100.;
        }
      }
    }
  }
  for(unsigned i=0; i<nres; ++i) {
    if(sasares[i]>sasa_cutoff) solv_res[i] = 1;
    else solv_res[i] = 0;
  }
}

void SAXS::calculate_gpu(std::vector<Vector> &pos, std::vector<Vector> &deriv)
{
#ifdef __PLUMED_HAS_ARRAYFIRE
  unsigned size;
  if(onebead) size = nres;
  else size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  std::vector<float> sum;
  sum.resize(numq);

  std::vector<float> dd;
  dd.resize(size*3*numq);

  // on gpu only the master rank run the calculation
  if(comm.Get_rank()==0) {
    std::vector<float> posi;
    posi.resize(3*size);
    #pragma omp parallel for num_threads(OpenMP::getNumThreads())
    for (unsigned i=0; i<size; ++i) {
      const Vector tmp = pos[i];
      posi[3*i]   = static_cast<float>(tmp[0]);
      posi[3*i+1] = static_cast<float>(tmp[1]);
      posi[3*i+2] = static_cast<float>(tmp[2]);
    }

    // create array a and b containing atomic coordinates
    af::setDevice(deviceid);
    // 3,size,1,1
    af::array pos_a = af::array(3, size, &posi.front());
    // size,3,1,1
    pos_a = af::moddims(pos_a.T(), size, 3, 1);
    // size,3,1,1
    af::array pos_b = pos_a(af::span, af::span);
    // size,1,3,1
    pos_a = af::moddims(pos_a, size, 1, 3);
    // 1,size,3,1
    pos_b = af::moddims(pos_b, 1, size, 3);

    // size,size,3,1
    af::array pos_a_t = af::tile(pos_a, 1, size, 1);
    // size,size,3,1: for some reason we need this
    pos_a_t = af::moddims(pos_a_t, size, size, 3);
    // size,size,3,1
    af::array pos_b_t = af::tile(pos_b, size, 1, 1);
    // size,size,3,1: for some reason we need this
    pos_b_t = af::moddims(pos_b_t, size, size, 3);
    // size,size,3,1
    af::array xyz_dist = pos_a_t - pos_b_t;
    // size,size,1,1
    af::array square = af::sum(xyz_dist*xyz_dist,2);
    // size,size,1,1
    af::array dist_sqrt = af::sqrt(square);
    // replace the zero of square with one to avoid nan in the derivatives (the number does not matter because this are multiplied by zero)
    af::replace(square,!(af::iszero(square)),1.);
    // size,size,3,1
    xyz_dist = xyz_dist / af::tile(square, 1, 1, 3);
    // numq,1,1,1
    af::array sum_device   = af::constant(0, numq, f32);
    // numq,size,3,1
    af::array deriv_device = af::constant(0, numq, size, 3, f32);

    for (unsigned k=0; k<numq; ++k) {
      // calculate FF matrix
      // size,1,1,1
      af::array AFF_value(size, &FFf_value[k].front());
      // size,size,1,1
      af::array FFdist_mod = af::tile(AFF_value(af::span), 1, size)*af::transpose(af::tile(AFF_value(af::span), 1, size));

      // get q
      const float qvalue = static_cast<float>(q_list[k]);
      // size,size,1,1
      af::array dist_q = qvalue*dist_sqrt;
      // size,size,1
      af::array dist_sin = af::sin(dist_q)/dist_q;
      af::replace(dist_sin,!(af::isNaN(dist_sin)),1.);
      // 1,1,1,1
      sum_device(k) = af::sum(af::flat(dist_sin)*af::flat(FFdist_mod));

      // size,size,1,1
      af::array tmp = FFdist_mod*(dist_sin - af::cos(dist_q));
      // size,size,3,1
      af::array dd_all = af::tile(tmp, 1, 1, 3)*xyz_dist;
      // it should become 1,size,3
      deriv_device(k, af::span, af::span) = af::sum(dd_all,0);
    }

    // read out results
    sum_device.host(&sum.front());

    deriv_device = af::reorder(deriv_device, 2, 1, 0);
    deriv_device = af::flat(deriv_device);
    deriv_device.host(&dd.front());
  }

  comm.Bcast(dd, 0);
  comm.Bcast(sum, 0);

  for(unsigned k=0; k<numq; ++k) {
    std::string num; Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) setCalcData(k, sum[k]);
    for(unsigned i=0; i<size; ++i) {
      const unsigned di = k*size*3+i*3;
      deriv[k*size+i] = Vector(2.*dd[di+0],2.*dd[di+1],2.*dd[di+2]);
    }
  }
#endif
}

void SAXS::calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv)
{
  unsigned size;
  if(onebead) size = nres;
  else size = getNumberOfAtoms();
  const unsigned numq = q_list.size();

  unsigned stride = comm.Get_size();
  unsigned rank   = comm.Get_rank();
  if(serial) {
    stride = 1;
    rank   = 0;
  }
  std::vector<double> sum(numq,0);
  unsigned nt=OpenMP::getNumThreads();
  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(deriv.size());
    std::vector<double> omp_sum(numq,0);
    #pragma omp for nowait
    for (unsigned i=rank; i<size-1; i+=stride) {
      Vector posi = pos[i];
      for (unsigned j=i+1; j<size ; ++j) {
        Vector c_distances = delta(posi,pos[j]);
        double m_distances = c_distances.modulo();
        c_distances = c_distances/m_distances/m_distances;
        for (unsigned k=0; k<numq; ++k) {
          unsigned kdx=k*size;
          double qdist = q_list[k]*m_distances;
          double FFF = 2.*FF_value[i][k]*FF_value[j][k];
          double tsq = std::sin(qdist)/qdist;
          double tcq = std::cos(qdist);
          double tmp = FFF*(tcq-tsq);
          Vector dd  = c_distances*tmp;
          if(nt>1) {
            omp_deriv[kdx+i] -=dd;
            omp_deriv[kdx+j] +=dd;
            omp_sum[k] += FFF*tsq;
          } else {
            deriv[kdx+i] -= dd;
            deriv[kdx+j] += dd;
            sum[k] += FFF*tsq;
          }
        }
      }
    }
    #pragma omp critical
    if(nt>1) {
      for(unsigned i=0; i<deriv.size(); ++i) deriv[i]+=omp_deriv[i];
      for(unsigned k=0; k<numq; ++k) sum[k]+=omp_sum[k];
    }
  }

  if(!serial) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
    comm.Sum(&sum[0], numq);
  }

  for (unsigned k=0; k<numq; ++k) {
    sum[k]+=FF_rank[k];
    std::string num; Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) setCalcData(k, sum[k]);
  }
}

void SAXS::calculate()
{
  if(pbc) makeWhole();

  const size_t size = getNumberOfAtoms();
  const size_t numq = q_list.size();

  // these are the derivatives associated to the coarse graining
  std::vector<Vector> aa_deriv(size);

  size_t beads_size = size;
  if(onebead) beads_size = nres;
  // these are the derivatives particle,q
  std::vector<Vector> bd_deriv(numq*beads_size);

  std::vector<Vector> beads_pos(beads_size);
  if(onebead) {
    // mapping
    unsigned atom_id = 0;
    for(unsigned i=0; i<nres; ++i) {
      /* calculate center and derivatives */
      double sum_mass = 0.;
      Vector sum_pos = Vector(0,0,0);
      for(unsigned j=0; j<atoms_per_bead[i]; ++j) {
        aa_deriv[atom_id] = Vector(atoms_masses[atom_id],atoms_masses[atom_id],atoms_masses[atom_id]);
        sum_pos += atoms_masses[atom_id] * getPosition(atom_id); // getPosition(first_atom+atom_id)
        sum_mass += atoms_masses[atom_id];
        // Here I update the atom_id to stay sync'd with masses vector
        atom_id++;
      }
      beads_pos[i] = sum_pos/sum_mass;
      for(unsigned j=atom_id-atoms_per_bead[i]; j<atom_id; ++j) {
        aa_deriv[j] /= sum_mass;
      }
    }
    // SASA
    std::vector<bool> solv_res(nres, 0);
    if(getStep()%solv_stride == 0 || isFirstStep) {
      isFirstStep = 0;
      if(rho_corr!=rho) sasa_calculate(solv_res);
      Iq0=0.;
      for(unsigned i=0; i<nres; ++i) {
        if(solv_res[i] == 1 ) {
          Iq0 += std::sqrt((Iq0_vac[i]+(rho_corr*rho_corr)*Iq0_wat[i]-rho_corr*Iq0_mix[i]));
        } else {
          Iq0 += std::sqrt((Iq0_vac[i]+(rho*rho)*Iq0_wat[i]-rho*Iq0_mix[i]));
        }
      }
      // Form Factors
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned i=0; i<nres; ++i) {
          if(!gpu) {
            if(solv_res[i] == 0) { // buried
              FF_value[i][k] = std::sqrt(FF_value_vacuum[atoi[i]][k] + rho*rho*FF_value_water[atoi[i]][k] - rho*FF_value_mixed[atoi[i]][k])/Iq0;
            } else { // surface
              FF_value[i][k] = std::sqrt(FF_value_vacuum[atoi[i]][k] + rho_corr*rho_corr*FF_value_water[atoi[i]][k] - rho_corr*FF_value_mixed[atoi[i]][k])/Iq0;
            }
          } else {
            if(solv_res[i] == 0) { // buried
              FFf_value[k][i] = static_cast<float>(std::sqrt(FF_value_vacuum[atoi[i]][k] + rho*rho*FF_value_water[atoi[i]][k] - rho*FF_value_mixed[atoi[i]][k])/Iq0);
            } else { // surface
              FFf_value[k][i] = static_cast<float>(std::sqrt(FF_value_vacuum[atoi[i]][k] + rho_corr*rho_corr*FF_value_water[atoi[i]][k] - rho_corr*FF_value_mixed[atoi[i]][k])/Iq0);
            }
          }
        }
      }
      if(!gpu) {
        for(unsigned k=0; k<numq; ++k) {
          FF_rank[k]=0.;
          for(unsigned i=0; i<nres; ++i) {
            FF_rank[k]+=FF_value[i][k]*FF_value[i][k];
          }
        }
      }
    }
    // not ONEBEAD
  } else {
    for(unsigned i=0; i<size; ++i) {
      beads_pos[i] = getPosition(i);
    }
    aa_deriv = std::vector<Vector>(size,(Vector(1,1,1)));
  }

  if(gpu) calculate_gpu(beads_pos, bd_deriv);
  else calculate_cpu(beads_pos, bd_deriv);

  if(getDoScore()) {
    /* Metainference */
    double score = getScore();
    setScore(score);
  }

  for (unsigned k=0; k<numq; ++k) {
    const unsigned kdx=k*beads_size;
    Tensor deriv_box;
    Value* val;
    if(!getDoScore()) {
      std::string num; Tools::convert(k,num);
      val=getPntrToComponent("q-"+num);

      for(unsigned i=0; i<beads_size; ++i) {
        setAtomsDerivatives(val, i, Vector(aa_deriv[i][0]*bd_deriv[kdx+i][0], \
                                           aa_deriv[i][1]*bd_deriv[kdx+i][1], \
                                           aa_deriv[i][2]*bd_deriv[kdx+i][2]) );
        deriv_box += Tensor(getPosition(i),Vector(aa_deriv[i][0]*bd_deriv[kdx+i][0], \
                            aa_deriv[i][1]*bd_deriv[kdx+i][1], \
                            aa_deriv[i][2]*bd_deriv[kdx+i][2]) );
      }
    } else {
      val=getPntrToComponent("score");
      for(unsigned i=0; i<beads_size; ++i) {
        setAtomsDerivatives(val, i, Vector(aa_deriv[i][0]*bd_deriv[kdx+i][0]*getMetaDer(k),
                                           aa_deriv[i][1]*bd_deriv[kdx+i][1]*getMetaDer(k),
                                           aa_deriv[i][2]*bd_deriv[kdx+i][2]*getMetaDer(k)) );
        deriv_box += Tensor(getPosition(i),Vector(aa_deriv[i][0]*bd_deriv[kdx+i][0]*getMetaDer(k),
                            aa_deriv[i][1]*bd_deriv[kdx+i][1]*getMetaDer(k),
                            aa_deriv[i][2]*bd_deriv[kdx+i][2]*getMetaDer(k)) );
      }
    }
    setBoxDerivatives(val, -deriv_box);
  }
}

void SAXS::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

void SAXS::getOnebeadMapping(const std::vector<AtomNumber> &atoms) {
  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( moldat ) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    // RC: Here we assume that we are with a continuous sequence; maybe we can extend it to
    // discontinuous ones.
    unsigned first_res = moldat->getResidueNumber(atoms[0]);
    nres = moldat->getResidueNumber(atoms[atoms.size()-1]) - moldat->getResidueNumber(atoms[0]) + 1;
    atoms_per_bead.resize(nres);
    atoms_masses.resize(atoms.size());
    residue_atom.resize(atoms.size());       //@MOD: add vector resize
    std::vector<std::vector<std::string> > AtomResidueName;
    AtomResidueName.resize(2);

    log.printf("  Onebead residue mapping on %u residues\n", nres);

    for(unsigned i=0; i<atoms.size(); ++i) {
      atoms_per_bead[moldat->getResidueNumber(atoms[i])-first_res]++;
      std::string Aname = moldat->getAtomName(atoms[i]);
      std::string Rname = moldat->getResidueName(atoms[i]);
      AtomResidueName[0].push_back(Aname);
      AtomResidueName[1].push_back(Rname);
      residue_atom[i] = moldat->getResidueNumber(atoms[i])-first_res;
      char type;
      char first = Aname.at(0);

      // We assume that element symbol is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = Aname.at(1);
      }

      if (type == 'H') atoms_masses[i] = 1.008;
      else if (type == 'C') atoms_masses[i] = 12.011;
      else if (type == 'N') atoms_masses[i] = 14.007;
      else if (type == 'O') atoms_masses[i] = 15.999;
      else if (type == 'S') atoms_masses[i] = 32.065;
      else if (type == 'P') atoms_masses[i] = 30.974;
      else {
        error("Unknown element in mass extraction\n");
      }
    }
    readLCPOparam(AtomResidueName, atoms.size());
  } else {
    error("MOLINFO DATA not found\n");
  }
}

void SAXS::getMartiniSFparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter)
{
  parameter[ALA_BB].push_back(9.045);
  parameter[ALA_BB].push_back(-0.098114);
  parameter[ALA_BB].push_back(7.54281);
  parameter[ALA_BB].push_back(-1.97438);
  parameter[ALA_BB].push_back(-8.32689);
  parameter[ALA_BB].push_back(6.09318);
  parameter[ALA_BB].push_back(-1.18913);

  parameter[ARG_BB].push_back(10.729);
  parameter[ARG_BB].push_back(-0.0392574);
  parameter[ARG_BB].push_back(1.15382);
  parameter[ARG_BB].push_back(-0.155999);
  parameter[ARG_BB].push_back(-2.43619);
  parameter[ARG_BB].push_back(1.72922);
  parameter[ARG_BB].push_back(-0.33799);

  parameter[ARG_SC1].push_back(-2.796);
  parameter[ARG_SC1].push_back(0.472403);
  parameter[ARG_SC1].push_back(8.07424);
  parameter[ARG_SC1].push_back(4.37299);
  parameter[ARG_SC1].push_back(-10.7398);
  parameter[ARG_SC1].push_back(4.95677);
  parameter[ARG_SC1].push_back(-0.725797);

  parameter[ARG_SC2].push_back(15.396);
  parameter[ARG_SC2].push_back(0.0636736);
  parameter[ARG_SC2].push_back(-1.258);
  parameter[ARG_SC2].push_back(1.93135);
  parameter[ARG_SC2].push_back(-4.45031);
  parameter[ARG_SC2].push_back(2.49356);
  parameter[ARG_SC2].push_back(-0.410721);

  parameter[ASN_BB].push_back(10.738);
  parameter[ASN_BB].push_back(-0.0402162);
  parameter[ASN_BB].push_back(1.03007);
  parameter[ASN_BB].push_back(-0.254174);
  parameter[ASN_BB].push_back(-2.12015);
  parameter[ASN_BB].push_back(1.55535);
  parameter[ASN_BB].push_back(-0.30963);

  parameter[ASN_SC1].push_back(9.249);
  parameter[ASN_SC1].push_back(-0.0148678);
  parameter[ASN_SC1].push_back(5.52169);
  parameter[ASN_SC1].push_back(0.00853212);
  parameter[ASN_SC1].push_back(-6.71992);
  parameter[ASN_SC1].push_back(3.93622);
  parameter[ASN_SC1].push_back(-0.64973);

  parameter[ASP_BB].push_back(10.695);
  parameter[ASP_BB].push_back(-0.0410247);
  parameter[ASP_BB].push_back(1.03656);
  parameter[ASP_BB].push_back(-0.298558);
  parameter[ASP_BB].push_back(-2.06064);
  parameter[ASP_BB].push_back(1.53495);
  parameter[ASP_BB].push_back(-0.308365);

  parameter[ASP_SC1].push_back(9.476);
  parameter[ASP_SC1].push_back(-0.0254664);
  parameter[ASP_SC1].push_back(5.57899);
  parameter[ASP_SC1].push_back(-0.395027);
  parameter[ASP_SC1].push_back(-5.9407);
  parameter[ASP_SC1].push_back(3.48836);
  parameter[ASP_SC1].push_back(-0.569402);

  parameter[CYS_BB].push_back(10.698);
  parameter[CYS_BB].push_back(-0.0233493);
  parameter[CYS_BB].push_back(1.18257);
  parameter[CYS_BB].push_back(0.0684464);
  parameter[CYS_BB].push_back(-2.792);
  parameter[CYS_BB].push_back(1.88995);
  parameter[CYS_BB].push_back(-0.360229);

  parameter[CYS_SC1].push_back(8.199);
  parameter[CYS_SC1].push_back(-0.0261569);
  parameter[CYS_SC1].push_back(6.79677);
  parameter[CYS_SC1].push_back(-0.343845);
  parameter[CYS_SC1].push_back(-5.03578);
  parameter[CYS_SC1].push_back(2.7076);
  parameter[CYS_SC1].push_back(-0.420714);

  parameter[GLN_BB].push_back(10.728);
  parameter[GLN_BB].push_back(-0.0391984);
  parameter[GLN_BB].push_back(1.09264);
  parameter[GLN_BB].push_back(-0.261555);
  parameter[GLN_BB].push_back(-2.21245);
  parameter[GLN_BB].push_back(1.62071);
  parameter[GLN_BB].push_back(-0.322325);

  parameter[GLN_SC1].push_back(8.317);
  parameter[GLN_SC1].push_back(-0.229045);
  parameter[GLN_SC1].push_back(12.6338);
  parameter[GLN_SC1].push_back(-7.6719);
  parameter[GLN_SC1].push_back(-5.8376);
  parameter[GLN_SC1].push_back(5.53784);
  parameter[GLN_SC1].push_back(-1.12604);

  parameter[GLU_BB].push_back(10.694);
  parameter[GLU_BB].push_back(-0.0521961);
  parameter[GLU_BB].push_back(1.11153);
  parameter[GLU_BB].push_back(-0.491995);
  parameter[GLU_BB].push_back(-1.86236);
  parameter[GLU_BB].push_back(1.45332);
  parameter[GLU_BB].push_back(-0.29708);

  parameter[GLU_SC1].push_back(8.544);
  parameter[GLU_SC1].push_back(-0.249555);
  parameter[GLU_SC1].push_back(12.8031);
  parameter[GLU_SC1].push_back(-8.42696);
  parameter[GLU_SC1].push_back(-4.66486);
  parameter[GLU_SC1].push_back(4.90004);
  parameter[GLU_SC1].push_back(-1.01204);

  parameter[GLY_BB].push_back(9.977);
  parameter[GLY_BB].push_back(-0.0285799);
  parameter[GLY_BB].push_back(1.84236);
  parameter[GLY_BB].push_back(-0.0315192);
  parameter[GLY_BB].push_back(-2.88326);
  parameter[GLY_BB].push_back(1.87323);
  parameter[GLY_BB].push_back(-0.345773);

  parameter[HIS_BB].push_back(10.721);
  parameter[HIS_BB].push_back(-0.0379337);
  parameter[HIS_BB].push_back(1.06028);
  parameter[HIS_BB].push_back(-0.236143);
  parameter[HIS_BB].push_back(-2.17819);
  parameter[HIS_BB].push_back(1.58357);
  parameter[HIS_BB].push_back(-0.31345);

  parameter[HIS_SC1].push_back(-0.424);
  parameter[HIS_SC1].push_back(0.665176);
  parameter[HIS_SC1].push_back(3.4369);
  parameter[HIS_SC1].push_back(2.93795);
  parameter[HIS_SC1].push_back(-5.18288);
  parameter[HIS_SC1].push_back(2.12381);
  parameter[HIS_SC1].push_back(-0.284224);

  parameter[HIS_SC2].push_back(5.363);
  parameter[HIS_SC2].push_back(-0.0176945);
  parameter[HIS_SC2].push_back(2.9506);
  parameter[HIS_SC2].push_back(-0.387018);
  parameter[HIS_SC2].push_back(-1.83951);
  parameter[HIS_SC2].push_back(0.9703);
  parameter[HIS_SC2].push_back(-0.1458);

  parameter[HIS_SC3].push_back(5.784);
  parameter[HIS_SC3].push_back(-0.0293129);
  parameter[HIS_SC3].push_back(2.74167);
  parameter[HIS_SC3].push_back(-0.520875);
  parameter[HIS_SC3].push_back(-1.62949);
  parameter[HIS_SC3].push_back(0.902379);
  parameter[HIS_SC3].push_back(-0.139957);

  parameter[ILE_BB].push_back(10.699);
  parameter[ILE_BB].push_back(-0.0188962);
  parameter[ILE_BB].push_back(1.217);
  parameter[ILE_BB].push_back(0.242481);
  parameter[ILE_BB].push_back(-3.13898);
  parameter[ILE_BB].push_back(2.07916);
  parameter[ILE_BB].push_back(-0.392574);

  parameter[ILE_SC1].push_back(-4.448);
  parameter[ILE_SC1].push_back(1.20996);
  parameter[ILE_SC1].push_back(11.5141);
  parameter[ILE_SC1].push_back(6.98895);
  parameter[ILE_SC1].push_back(-19.1948);
  parameter[ILE_SC1].push_back(9.89207);
  parameter[ILE_SC1].push_back(-1.60877);

  parameter[LEU_BB].push_back(10.692);
  parameter[LEU_BB].push_back(-0.0414917);
  parameter[LEU_BB].push_back(1.1077);
  parameter[LEU_BB].push_back(-0.288062);
  parameter[LEU_BB].push_back(-2.17187);
  parameter[LEU_BB].push_back(1.59879);
  parameter[LEU_BB].push_back(-0.318545);

  parameter[LEU_SC1].push_back(-4.448);
  parameter[LEU_SC1].push_back(2.1063);
  parameter[LEU_SC1].push_back(6.72381);
  parameter[LEU_SC1].push_back(14.6954);
  parameter[LEU_SC1].push_back(-23.7197);
  parameter[LEU_SC1].push_back(10.7247);
  parameter[LEU_SC1].push_back(-1.59146);

  parameter[LYS_BB].push_back(10.706);
  parameter[LYS_BB].push_back(-0.0468629);
  parameter[LYS_BB].push_back(1.09477);
  parameter[LYS_BB].push_back(-0.432751);
  parameter[LYS_BB].push_back(-1.94335);
  parameter[LYS_BB].push_back(1.49109);
  parameter[LYS_BB].push_back(-0.302589);

  parameter[LYS_SC1].push_back(-2.796);
  parameter[LYS_SC1].push_back(0.508044);
  parameter[LYS_SC1].push_back(7.91436);
  parameter[LYS_SC1].push_back(4.54097);
  parameter[LYS_SC1].push_back(-10.8051);
  parameter[LYS_SC1].push_back(4.96204);
  parameter[LYS_SC1].push_back(-0.724414);

  parameter[LYS_SC2].push_back(3.070);
  parameter[LYS_SC2].push_back(-0.0101448);
  parameter[LYS_SC2].push_back(4.67994);
  parameter[LYS_SC2].push_back(-0.792529);
  parameter[LYS_SC2].push_back(-2.09142);
  parameter[LYS_SC2].push_back(1.02933);
  parameter[LYS_SC2].push_back(-0.137787);

  parameter[MET_BB].push_back(10.671);
  parameter[MET_BB].push_back(-0.0433724);
  parameter[MET_BB].push_back(1.13784);
  parameter[MET_BB].push_back(-0.40768);
  parameter[MET_BB].push_back(-2.00555);
  parameter[MET_BB].push_back(1.51673);
  parameter[MET_BB].push_back(-0.305547);

  parameter[MET_SC1].push_back(5.85);
  parameter[MET_SC1].push_back(-0.0485798);
  parameter[MET_SC1].push_back(17.0391);
  parameter[MET_SC1].push_back(-3.65327);
  parameter[MET_SC1].push_back(-13.174);
  parameter[MET_SC1].push_back(8.68286);
  parameter[MET_SC1].push_back(-1.56095);

  parameter[PHE_BB].push_back(10.741);
  parameter[PHE_BB].push_back(-0.0317275);
  parameter[PHE_BB].push_back(1.15599);
  parameter[PHE_BB].push_back(0.0276187);
  parameter[PHE_BB].push_back(-2.74757);
  parameter[PHE_BB].push_back(1.88783);
  parameter[PHE_BB].push_back(-0.363525);

  parameter[PHE_SC1].push_back(-0.636);
  parameter[PHE_SC1].push_back(0.527882);
  parameter[PHE_SC1].push_back(6.77612);
  parameter[PHE_SC1].push_back(3.18508);
  parameter[PHE_SC1].push_back(-8.92826);
  parameter[PHE_SC1].push_back(4.29752);
  parameter[PHE_SC1].push_back(-0.65187);

  parameter[PHE_SC2].push_back(-0.424);
  parameter[PHE_SC2].push_back(0.389174);
  parameter[PHE_SC2].push_back(4.11761);
  parameter[PHE_SC2].push_back(2.29527);
  parameter[PHE_SC2].push_back(-4.7652);
  parameter[PHE_SC2].push_back(1.97023);
  parameter[PHE_SC2].push_back(-0.262318);

  parameter[PHE_SC3].push_back(-0.424);
  parameter[PHE_SC3].push_back(0.38927);
  parameter[PHE_SC3].push_back(4.11708);
  parameter[PHE_SC3].push_back(2.29623);
  parameter[PHE_SC3].push_back(-4.76592);
  parameter[PHE_SC3].push_back(1.97055);
  parameter[PHE_SC3].push_back(-0.262381);

  parameter[PRO_BB].push_back(11.434);
  parameter[PRO_BB].push_back(-0.033323);
  parameter[PRO_BB].push_back(0.472014);
  parameter[PRO_BB].push_back(-0.290854);
  parameter[PRO_BB].push_back(-1.81409);
  parameter[PRO_BB].push_back(1.39751);
  parameter[PRO_BB].push_back(-0.280407);

  parameter[PRO_SC1].push_back(-2.796);
  parameter[PRO_SC1].push_back(0.95668);
  parameter[PRO_SC1].push_back(6.84197);
  parameter[PRO_SC1].push_back(6.43774);
  parameter[PRO_SC1].push_back(-12.5068);
  parameter[PRO_SC1].push_back(5.64597);
  parameter[PRO_SC1].push_back(-0.825206);

  parameter[SER_BB].push_back(10.699);
  parameter[SER_BB].push_back(-0.0325828);
  parameter[SER_BB].push_back(1.20329);
  parameter[SER_BB].push_back(-0.0674351);
  parameter[SER_BB].push_back(-2.60749);
  parameter[SER_BB].push_back(1.80318);
  parameter[SER_BB].push_back(-0.346803);

  parameter[SER_SC1].push_back(3.298);
  parameter[SER_SC1].push_back(-0.0366801);
  parameter[SER_SC1].push_back(5.11077);
  parameter[SER_SC1].push_back(-1.46774);
  parameter[SER_SC1].push_back(-1.48421);
  parameter[SER_SC1].push_back(0.800326);
  parameter[SER_SC1].push_back(-0.108314);

  parameter[THR_BB].push_back(10.697);
  parameter[THR_BB].push_back(-0.0242955);
  parameter[THR_BB].push_back(1.24671);
  parameter[THR_BB].push_back(0.146423);
  parameter[THR_BB].push_back(-2.97429);
  parameter[THR_BB].push_back(1.97513);
  parameter[THR_BB].push_back(-0.371479);

  parameter[THR_SC1].push_back(2.366);
  parameter[THR_SC1].push_back(0.0297604);
  parameter[THR_SC1].push_back(11.9216);
  parameter[THR_SC1].push_back(-9.32503);
  parameter[THR_SC1].push_back(1.9396);
  parameter[THR_SC1].push_back(0.0804861);
  parameter[THR_SC1].push_back(-0.0302721);

  parameter[TRP_BB].push_back(10.689);
  parameter[TRP_BB].push_back(-0.0265879);
  parameter[TRP_BB].push_back(1.17819);
  parameter[TRP_BB].push_back(0.0386457);
  parameter[TRP_BB].push_back(-2.75634);
  parameter[TRP_BB].push_back(1.88065);
  parameter[TRP_BB].push_back(-0.360217);

  parameter[TRP_SC1].push_back(0.084);
  parameter[TRP_SC1].push_back(0.752407);
  parameter[TRP_SC1].push_back(5.3802);
  parameter[TRP_SC1].push_back(4.09281);
  parameter[TRP_SC1].push_back(-9.28029);
  parameter[TRP_SC1].push_back(4.45923);
  parameter[TRP_SC1].push_back(-0.689008);

  parameter[TRP_SC2].push_back(5.739);
  parameter[TRP_SC2].push_back(0.0298492);
  parameter[TRP_SC2].push_back(4.60446);
  parameter[TRP_SC2].push_back(1.34463);
  parameter[TRP_SC2].push_back(-5.69968);
  parameter[TRP_SC2].push_back(2.84924);
  parameter[TRP_SC2].push_back(-0.433781);

  parameter[TRP_SC3].push_back(-0.424);
  parameter[TRP_SC3].push_back(0.388576);
  parameter[TRP_SC3].push_back(4.11859);
  parameter[TRP_SC3].push_back(2.29485);
  parameter[TRP_SC3].push_back(-4.76255);
  parameter[TRP_SC3].push_back(1.96849);
  parameter[TRP_SC3].push_back(-0.262015);

  parameter[TRP_SC4].push_back(-0.424);
  parameter[TRP_SC4].push_back(0.387685);
  parameter[TRP_SC4].push_back(4.12153);
  parameter[TRP_SC4].push_back(2.29144);
  parameter[TRP_SC4].push_back(-4.7589);
  parameter[TRP_SC4].push_back(1.96686);
  parameter[TRP_SC4].push_back(-0.261786);

  parameter[TYR_BB].push_back(10.689);
  parameter[TYR_BB].push_back(-0.0193526);
  parameter[TYR_BB].push_back(1.18241);
  parameter[TYR_BB].push_back(0.207318);
  parameter[TYR_BB].push_back(-3.0041);
  parameter[TYR_BB].push_back(1.99335);
  parameter[TYR_BB].push_back(-0.376482);

  parameter[TYR_SC1].push_back(-0.636);
  parameter[TYR_SC1].push_back(0.528902);
  parameter[TYR_SC1].push_back(6.78168);
  parameter[TYR_SC1].push_back(3.17769);
  parameter[TYR_SC1].push_back(-8.93667);
  parameter[TYR_SC1].push_back(4.30692);
  parameter[TYR_SC1].push_back(-0.653993);

  parameter[TYR_SC2].push_back(-0.424);
  parameter[TYR_SC2].push_back(0.388811);
  parameter[TYR_SC2].push_back(4.11851);
  parameter[TYR_SC2].push_back(2.29545);
  parameter[TYR_SC2].push_back(-4.7668);
  parameter[TYR_SC2].push_back(1.97131);
  parameter[TYR_SC2].push_back(-0.262534);

  parameter[TYR_SC3].push_back(4.526);
  parameter[TYR_SC3].push_back(-0.00381305);
  parameter[TYR_SC3].push_back(5.8567);
  parameter[TYR_SC3].push_back(-0.214086);
  parameter[TYR_SC3].push_back(-4.63649);
  parameter[TYR_SC3].push_back(2.52869);
  parameter[TYR_SC3].push_back(-0.39894);

  parameter[VAL_BB].push_back(10.691);
  parameter[VAL_BB].push_back(-0.0162929);
  parameter[VAL_BB].push_back(1.24446);
  parameter[VAL_BB].push_back(0.307914);
  parameter[VAL_BB].push_back(-3.27446);
  parameter[VAL_BB].push_back(2.14788);
  parameter[VAL_BB].push_back(-0.403259);

  parameter[VAL_SC1].push_back(-3.516);
  parameter[VAL_SC1].push_back(1.62307);
  parameter[VAL_SC1].push_back(5.43064);
  parameter[VAL_SC1].push_back(9.28809);
  parameter[VAL_SC1].push_back(-14.9927);
  parameter[VAL_SC1].push_back(6.6133);
  parameter[VAL_SC1].push_back(-0.964977);

  parameter[A_BB1].push_back(32.88500000);
  parameter[A_BB1].push_back(0.08339900);
  parameter[A_BB1].push_back(-7.36054400);
  parameter[A_BB1].push_back(2.19220300);
  parameter[A_BB1].push_back(-3.56523400);
  parameter[A_BB1].push_back(2.33326900);
  parameter[A_BB1].push_back(-0.39785500);

  parameter[A_BB2].push_back(3.80600000);
  parameter[A_BB2].push_back(-0.10727600);
  parameter[A_BB2].push_back(9.58854100);
  parameter[A_BB2].push_back(-6.23740500);
  parameter[A_BB2].push_back(-0.48267300);
  parameter[A_BB2].push_back(1.14119500);
  parameter[A_BB2].push_back(-0.21385600);

  parameter[A_BB3].push_back(3.59400000);
  parameter[A_BB3].push_back(0.04537300);
  parameter[A_BB3].push_back(9.59178900);
  parameter[A_BB3].push_back(-1.29202200);
  parameter[A_BB3].push_back(-7.10851000);
  parameter[A_BB3].push_back(4.05571200);
  parameter[A_BB3].push_back(-0.63372500);

  parameter[A_SC1].push_back(6.67100000);
  parameter[A_SC1].push_back(-0.00855300);
  parameter[A_SC1].push_back(1.63222400);
  parameter[A_SC1].push_back(-0.06466200);
  parameter[A_SC1].push_back(-1.48694200);
  parameter[A_SC1].push_back(0.78544600);
  parameter[A_SC1].push_back(-0.12083500);

  parameter[A_SC2].push_back(5.95100000);
  parameter[A_SC2].push_back(-0.02606600);
  parameter[A_SC2].push_back(2.54399900);
  parameter[A_SC2].push_back(-0.48436900);
  parameter[A_SC2].push_back(-1.55357400);
  parameter[A_SC2].push_back(0.86466900);
  parameter[A_SC2].push_back(-0.13509000);

  parameter[A_SC3].push_back(11.39400000);
  parameter[A_SC3].push_back(0.00871300);
  parameter[A_SC3].push_back(-0.23891300);
  parameter[A_SC3].push_back(0.48919400);
  parameter[A_SC3].push_back(-1.75289400);
  parameter[A_SC3].push_back(0.99267500);
  parameter[A_SC3].push_back(-0.16291300);

  parameter[A_SC4].push_back(6.45900000);
  parameter[A_SC4].push_back(0.01990600);
  parameter[A_SC4].push_back(4.17970400);
  parameter[A_SC4].push_back(0.97629900);
  parameter[A_SC4].push_back(-5.03297800);
  parameter[A_SC4].push_back(2.55576700);
  parameter[A_SC4].push_back(-0.39150500);

  parameter[A_3TE].push_back(4.23000000);
  parameter[A_3TE].push_back(0.00064800);
  parameter[A_3TE].push_back(0.92124600);
  parameter[A_3TE].push_back(0.08064300);
  parameter[A_3TE].push_back(-0.39054400);
  parameter[A_3TE].push_back(0.12429100);
  parameter[A_3TE].push_back(-0.01122700);

  parameter[A_5TE].push_back(4.23000000);
  parameter[A_5TE].push_back(0.00039300);
  parameter[A_5TE].push_back(0.92305100);
  parameter[A_5TE].push_back(0.07747500);
  parameter[A_5TE].push_back(-0.38792100);
  parameter[A_5TE].push_back(0.12323800);
  parameter[A_5TE].push_back(-0.01106600);

  parameter[A_TE3].push_back(7.82400000);
  parameter[A_TE3].push_back(-0.04881000);
  parameter[A_TE3].push_back(8.21557900);
  parameter[A_TE3].push_back(-0.89491400);
  parameter[A_TE3].push_back(-9.54293700);
  parameter[A_TE3].push_back(6.33122200);
  parameter[A_TE3].push_back(-1.16672900);

  parameter[A_TE5].push_back(8.03600000);
  parameter[A_TE5].push_back(0.01641200);
  parameter[A_TE5].push_back(5.14902200);
  parameter[A_TE5].push_back(0.83419700);
  parameter[A_TE5].push_back(-7.59068300);
  parameter[A_TE5].push_back(4.52063200);
  parameter[A_TE5].push_back(-0.78260800);

  parameter[C_BB1].push_back(32.88500000);
  parameter[C_BB1].push_back(0.08311100);
  parameter[C_BB1].push_back(-7.35432100);
  parameter[C_BB1].push_back(2.18610000);
  parameter[C_BB1].push_back(-3.55788300);
  parameter[C_BB1].push_back(2.32918700);
  parameter[C_BB1].push_back(-0.39720000);

  parameter[C_BB2].push_back(3.80600000);
  parameter[C_BB2].push_back(-0.10808100);
  parameter[C_BB2].push_back(9.61612600);
  parameter[C_BB2].push_back(-6.28595400);
  parameter[C_BB2].push_back(-0.45187000);
  parameter[C_BB2].push_back(1.13326000);
  parameter[C_BB2].push_back(-0.21320300);

  parameter[C_BB3].push_back(3.59400000);
  parameter[C_BB3].push_back(0.04484200);
  parameter[C_BB3].push_back(9.61919800);
  parameter[C_BB3].push_back(-1.33582800);
  parameter[C_BB3].push_back(-7.07200400);
  parameter[C_BB3].push_back(4.03952900);
  parameter[C_BB3].push_back(-0.63098200);

  parameter[C_SC1].push_back(5.95100000);
  parameter[C_SC1].push_back(-0.02911300);
  parameter[C_SC1].push_back(2.59700400);
  parameter[C_SC1].push_back(-0.55507700);
  parameter[C_SC1].push_back(-1.56344600);
  parameter[C_SC1].push_back(0.88956200);
  parameter[C_SC1].push_back(-0.14061300);

  parameter[C_SC2].push_back(11.62100000);
  parameter[C_SC2].push_back(0.01366100);
  parameter[C_SC2].push_back(-0.25959200);
  parameter[C_SC2].push_back(0.48918300);
  parameter[C_SC2].push_back(-1.52550500);
  parameter[C_SC2].push_back(0.83644100);
  parameter[C_SC2].push_back(-0.13407300);

  parameter[C_SC3].push_back(5.01900000);
  parameter[C_SC3].push_back(-0.03276100);
  parameter[C_SC3].push_back(5.53776900);
  parameter[C_SC3].push_back(-0.95105000);
  parameter[C_SC3].push_back(-3.71130800);
  parameter[C_SC3].push_back(2.16146000);
  parameter[C_SC3].push_back(-0.34918600);

  parameter[C_3TE].push_back(4.23000000);
  parameter[C_3TE].push_back(0.00057300);
  parameter[C_3TE].push_back(0.92174800);
  parameter[C_3TE].push_back(0.07964500);
  parameter[C_3TE].push_back(-0.38965700);
  parameter[C_3TE].push_back(0.12392500);
  parameter[C_3TE].push_back(-0.01117000);

  parameter[C_5TE].push_back(4.23000000);
  parameter[C_5TE].push_back(0.00071000);
  parameter[C_5TE].push_back(0.92082800);
  parameter[C_5TE].push_back(0.08150600);
  parameter[C_5TE].push_back(-0.39127000);
  parameter[C_5TE].push_back(0.12455900);
  parameter[C_5TE].push_back(-0.01126300);

  parameter[C_TE3].push_back(7.82400000);
  parameter[C_TE3].push_back(-0.05848300);
  parameter[C_TE3].push_back(8.29319900);
  parameter[C_TE3].push_back(-1.12563800);
  parameter[C_TE3].push_back(-9.42197600);
  parameter[C_TE3].push_back(6.35441700);
  parameter[C_TE3].push_back(-1.18356900);

  parameter[C_TE5].push_back(8.03600000);
  parameter[C_TE5].push_back(0.00493500);
  parameter[C_TE5].push_back(4.92622000);
  parameter[C_TE5].push_back(0.64810700);
  parameter[C_TE5].push_back(-7.05100000);
  parameter[C_TE5].push_back(4.26064400);
  parameter[C_TE5].push_back(-0.74819100);

  parameter[G_BB1].push_back(32.88500000);
  parameter[G_BB1].push_back(0.08325400);
  parameter[G_BB1].push_back(-7.35736000);
  parameter[G_BB1].push_back(2.18914800);
  parameter[G_BB1].push_back(-3.56154800);
  parameter[G_BB1].push_back(2.33120600);
  parameter[G_BB1].push_back(-0.39752300);

  parameter[G_BB2].push_back(3.80600000);
  parameter[G_BB2].push_back(-0.10788300);
  parameter[G_BB2].push_back(9.60930800);
  parameter[G_BB2].push_back(-6.27402500);
  parameter[G_BB2].push_back(-0.46192700);
  parameter[G_BB2].push_back(1.13737000);
  parameter[G_BB2].push_back(-0.21383100);

  parameter[G_BB3].push_back(3.59400000);
  parameter[G_BB3].push_back(0.04514500);
  parameter[G_BB3].push_back(9.61234700);
  parameter[G_BB3].push_back(-1.31542100);
  parameter[G_BB3].push_back(-7.09150500);
  parameter[G_BB3].push_back(4.04706200);
  parameter[G_BB3].push_back(-0.63201000);

  parameter[G_SC1].push_back(6.67100000);
  parameter[G_SC1].push_back(-0.00863200);
  parameter[G_SC1].push_back(1.63252300);
  parameter[G_SC1].push_back(-0.06567200);
  parameter[G_SC1].push_back(-1.48680500);
  parameter[G_SC1].push_back(0.78565600);
  parameter[G_SC1].push_back(-0.12088900);

  parameter[G_SC2].push_back(11.39400000);
  parameter[G_SC2].push_back(0.00912200);
  parameter[G_SC2].push_back(-0.22869000);
  parameter[G_SC2].push_back(0.49616400);
  parameter[G_SC2].push_back(-1.75039000);
  parameter[G_SC2].push_back(0.98649200);
  parameter[G_SC2].push_back(-0.16141600);

  parameter[G_SC3].push_back(10.90100000);
  parameter[G_SC3].push_back(0.02208700);
  parameter[G_SC3].push_back(0.17032800);
  parameter[G_SC3].push_back(0.73280800);
  parameter[G_SC3].push_back(-1.95292000);
  parameter[G_SC3].push_back(0.98357600);
  parameter[G_SC3].push_back(-0.14790900);

  parameter[G_SC4].push_back(6.45900000);
  parameter[G_SC4].push_back(0.02023700);
  parameter[G_SC4].push_back(4.17655400);
  parameter[G_SC4].push_back(0.98731800);
  parameter[G_SC4].push_back(-5.04352800);
  parameter[G_SC4].push_back(2.56059400);
  parameter[G_SC4].push_back(-0.39234300);

  parameter[G_3TE].push_back(4.23000000);
  parameter[G_3TE].push_back(0.00066300);
  parameter[G_3TE].push_back(0.92118800);
  parameter[G_3TE].push_back(0.08062700);
  parameter[G_3TE].push_back(-0.39041600);
  parameter[G_3TE].push_back(0.12419400);
  parameter[G_3TE].push_back(-0.01120500);

  parameter[G_5TE].push_back(4.23000000);
  parameter[G_5TE].push_back(0.00062800);
  parameter[G_5TE].push_back(0.92133500);
  parameter[G_5TE].push_back(0.08029900);
  parameter[G_5TE].push_back(-0.39015300);
  parameter[G_5TE].push_back(0.12411600);
  parameter[G_5TE].push_back(-0.01119900);

  parameter[G_TE3].push_back(7.82400000);
  parameter[G_TE3].push_back(-0.05177400);
  parameter[G_TE3].push_back(8.34606700);
  parameter[G_TE3].push_back(-1.02936300);
  parameter[G_TE3].push_back(-9.55211900);
  parameter[G_TE3].push_back(6.37776600);
  parameter[G_TE3].push_back(-1.17898000);

  parameter[G_TE5].push_back(8.03600000);
  parameter[G_TE5].push_back(0.00525100);
  parameter[G_TE5].push_back(4.71070600);
  parameter[G_TE5].push_back(0.66746900);
  parameter[G_TE5].push_back(-6.72538700);
  parameter[G_TE5].push_back(4.03644100);
  parameter[G_TE5].push_back(-0.70605700);

  parameter[U_BB1].push_back(32.88500000);
  parameter[U_BB1].push_back(0.08321400);
  parameter[U_BB1].push_back(-7.35634900);
  parameter[U_BB1].push_back(2.18826800);
  parameter[U_BB1].push_back(-3.56047400);
  parameter[U_BB1].push_back(2.33064700);
  parameter[U_BB1].push_back(-0.39744000);

  parameter[U_BB2].push_back(3.80600000);
  parameter[U_BB2].push_back(-0.10773100);
  parameter[U_BB2].push_back(9.60099900);
  parameter[U_BB2].push_back(-6.26131900);
  parameter[U_BB2].push_back(-0.46668300);
  parameter[U_BB2].push_back(1.13698100);
  parameter[U_BB2].push_back(-0.21351600);

  parameter[U_BB3].push_back(3.59400000);
  parameter[U_BB3].push_back(0.04544300);
  parameter[U_BB3].push_back(9.59625900);
  parameter[U_BB3].push_back(-1.29222200);
  parameter[U_BB3].push_back(-7.11143200);
  parameter[U_BB3].push_back(4.05687700);
  parameter[U_BB3].push_back(-0.63382800);

  parameter[U_SC1].push_back(5.95100000);
  parameter[U_SC1].push_back(-0.02924500);
  parameter[U_SC1].push_back(2.59668700);
  parameter[U_SC1].push_back(-0.56118700);
  parameter[U_SC1].push_back(-1.56477100);
  parameter[U_SC1].push_back(0.89265100);
  parameter[U_SC1].push_back(-0.14130800);

  parameter[U_SC2].push_back(10.90100000);
  parameter[U_SC2].push_back(0.02178900);
  parameter[U_SC2].push_back(0.18839000);
  parameter[U_SC2].push_back(0.72223100);
  parameter[U_SC2].push_back(-1.92581600);
  parameter[U_SC2].push_back(0.96654300);
  parameter[U_SC2].push_back(-0.14501300);

  parameter[U_SC3].push_back(5.24600000);
  parameter[U_SC3].push_back(-0.04586500);
  parameter[U_SC3].push_back(5.89978100);
  parameter[U_SC3].push_back(-1.50664700);
  parameter[U_SC3].push_back(-3.17054400);
  parameter[U_SC3].push_back(1.93717100);
  parameter[U_SC3].push_back(-0.31701000);

  parameter[U_3TE].push_back(4.23000000);
  parameter[U_3TE].push_back(0.00067500);
  parameter[U_3TE].push_back(0.92102300);
  parameter[U_3TE].push_back(0.08100800);
  parameter[U_3TE].push_back(-0.39084300);
  parameter[U_3TE].push_back(0.12441900);
  parameter[U_3TE].push_back(-0.01124900);

  parameter[U_5TE].push_back(4.23000000);
  parameter[U_5TE].push_back(0.00059000);
  parameter[U_5TE].push_back(0.92154600);
  parameter[U_5TE].push_back(0.07968200);
  parameter[U_5TE].push_back(-0.38950100);
  parameter[U_5TE].push_back(0.12382500);
  parameter[U_5TE].push_back(-0.01115100);

  parameter[U_TE3].push_back(7.82400000);
  parameter[U_TE3].push_back(-0.02968100);
  parameter[U_TE3].push_back(7.93783200);
  parameter[U_TE3].push_back(-0.33078100);
  parameter[U_TE3].push_back(-10.14120200);
  parameter[U_TE3].push_back(6.63334700);
  parameter[U_TE3].push_back(-1.22111200);

  parameter[U_TE5].push_back(8.03600000);
  parameter[U_TE5].push_back(-0.00909700);
  parameter[U_TE5].push_back(4.33193500);
  parameter[U_TE5].push_back(0.43416500);
  parameter[U_TE5].push_back(-5.80831400);
  parameter[U_TE5].push_back(3.52438800);
  parameter[U_TE5].push_back(-0.62382400);

  parameter[DA_BB1].push_back(32.88500000);
  parameter[DA_BB1].push_back(0.08179900);
  parameter[DA_BB1].push_back(-7.31735900);
  parameter[DA_BB1].push_back(2.15614500);
  parameter[DA_BB1].push_back(-3.52263200);
  parameter[DA_BB1].push_back(2.30604700);
  parameter[DA_BB1].push_back(-0.39270100);

  parameter[DA_BB2].push_back(3.80600000);
  parameter[DA_BB2].push_back(-0.10597700);
  parameter[DA_BB2].push_back(9.52537500);
  parameter[DA_BB2].push_back(-6.12991000);
  parameter[DA_BB2].push_back(-0.54092600);
  parameter[DA_BB2].push_back(1.15429100);
  parameter[DA_BB2].push_back(-0.21503500);

  parameter[DA_BB3].push_back(-1.35600000);
  parameter[DA_BB3].push_back(0.58928300);
  parameter[DA_BB3].push_back(6.71894100);
  parameter[DA_BB3].push_back(4.14050900);
  parameter[DA_BB3].push_back(-9.65859900);
  parameter[DA_BB3].push_back(4.43185000);
  parameter[DA_BB3].push_back(-0.64657300);

  parameter[DA_SC1].push_back(6.67100000);
  parameter[DA_SC1].push_back(-0.00871400);
  parameter[DA_SC1].push_back(1.63289100);
  parameter[DA_SC1].push_back(-0.06637700);
  parameter[DA_SC1].push_back(-1.48632900);
  parameter[DA_SC1].push_back(0.78551800);
  parameter[DA_SC1].push_back(-0.12087300);

  parameter[DA_SC2].push_back(5.95100000);
  parameter[DA_SC2].push_back(-0.02634300);
  parameter[DA_SC2].push_back(2.54864300);
  parameter[DA_SC2].push_back(-0.49015800);
  parameter[DA_SC2].push_back(-1.55386900);
  parameter[DA_SC2].push_back(0.86630200);
  parameter[DA_SC2].push_back(-0.13546200);

  parameter[DA_SC3].push_back(11.39400000);
  parameter[DA_SC3].push_back(0.00859500);
  parameter[DA_SC3].push_back(-0.25471400);
  parameter[DA_SC3].push_back(0.48718800);
  parameter[DA_SC3].push_back(-1.74520000);
  parameter[DA_SC3].push_back(0.99246200);
  parameter[DA_SC3].push_back(-0.16351900);

  parameter[DA_SC4].push_back(6.45900000);
  parameter[DA_SC4].push_back(0.01991800);
  parameter[DA_SC4].push_back(4.17962300);
  parameter[DA_SC4].push_back(0.97469100);
  parameter[DA_SC4].push_back(-5.02950400);
  parameter[DA_SC4].push_back(2.55371800);
  parameter[DA_SC4].push_back(-0.39113400);

  parameter[DA_3TE].push_back(4.23000000);
  parameter[DA_3TE].push_back(0.00062600);
  parameter[DA_3TE].push_back(0.92142000);
  parameter[DA_3TE].push_back(0.08016400);
  parameter[DA_3TE].push_back(-0.39000300);
  parameter[DA_3TE].push_back(0.12402500);
  parameter[DA_3TE].push_back(-0.01117900);

  parameter[DA_5TE].push_back(4.23000000);
  parameter[DA_5TE].push_back(0.00055500);
  parameter[DA_5TE].push_back(0.92183900);
  parameter[DA_5TE].push_back(0.07907600);
  parameter[DA_5TE].push_back(-0.38895100);
  parameter[DA_5TE].push_back(0.12359600);
  parameter[DA_5TE].push_back(-0.01111600);

  parameter[DA_TE3].push_back(2.87400000);
  parameter[DA_TE3].push_back(0.00112900);
  parameter[DA_TE3].push_back(12.51167200);
  parameter[DA_TE3].push_back(-7.67548000);
  parameter[DA_TE3].push_back(-2.02234000);
  parameter[DA_TE3].push_back(2.50837100);
  parameter[DA_TE3].push_back(-0.49458500);

  parameter[DA_TE5].push_back(8.03600000);
  parameter[DA_TE5].push_back(0.00473100);
  parameter[DA_TE5].push_back(4.65554400);
  parameter[DA_TE5].push_back(0.66424100);
  parameter[DA_TE5].push_back(-6.62131300);
  parameter[DA_TE5].push_back(3.96107400);
  parameter[DA_TE5].push_back(-0.69075800);

  parameter[DC_BB1].push_back(32.88500000);
  parameter[DC_BB1].push_back(0.08189900);
  parameter[DC_BB1].push_back(-7.32493500);
  parameter[DC_BB1].push_back(2.15976900);
  parameter[DC_BB1].push_back(-3.52612100);
  parameter[DC_BB1].push_back(2.31058600);
  parameter[DC_BB1].push_back(-0.39402700);

  parameter[DC_BB2].push_back(3.80600000);
  parameter[DC_BB2].push_back(-0.10559800);
  parameter[DC_BB2].push_back(9.52527700);
  parameter[DC_BB2].push_back(-6.12131700);
  parameter[DC_BB2].push_back(-0.54899400);
  parameter[DC_BB2].push_back(1.15592900);
  parameter[DC_BB2].push_back(-0.21494500);

  parameter[DC_BB3].push_back(-1.35600000);
  parameter[DC_BB3].push_back(0.55525700);
  parameter[DC_BB3].push_back(6.80305500);
  parameter[DC_BB3].push_back(4.05924700);
  parameter[DC_BB3].push_back(-9.61034700);
  parameter[DC_BB3].push_back(4.41253800);
  parameter[DC_BB3].push_back(-0.64315100);

  parameter[DC_SC1].push_back(5.95100000);
  parameter[DC_SC1].push_back(-0.02899900);
  parameter[DC_SC1].push_back(2.59587800);
  parameter[DC_SC1].push_back(-0.55388300);
  parameter[DC_SC1].push_back(-1.56395100);
  parameter[DC_SC1].push_back(0.88967400);
  parameter[DC_SC1].push_back(-0.14062500);

  parameter[DC_SC2].push_back(11.62100000);
  parameter[DC_SC2].push_back(0.01358100);
  parameter[DC_SC2].push_back(-0.24913000);
  parameter[DC_SC2].push_back(0.48787200);
  parameter[DC_SC2].push_back(-1.52867300);
  parameter[DC_SC2].push_back(0.83694900);
  parameter[DC_SC2].push_back(-0.13395300);

  parameter[DC_SC3].push_back(5.01900000);
  parameter[DC_SC3].push_back(-0.03298400);
  parameter[DC_SC3].push_back(5.54242800);
  parameter[DC_SC3].push_back(-0.96081500);
  parameter[DC_SC3].push_back(-3.71051600);
  parameter[DC_SC3].push_back(2.16500200);
  parameter[DC_SC3].push_back(-0.35023400);

  parameter[DC_3TE].push_back(4.23000000);
  parameter[DC_3TE].push_back(0.00055700);
  parameter[DC_3TE].push_back(0.92181400);
  parameter[DC_3TE].push_back(0.07924000);
  parameter[DC_3TE].push_back(-0.38916400);
  parameter[DC_3TE].push_back(0.12369900);
  parameter[DC_3TE].push_back(-0.01113300);

  parameter[DC_5TE].push_back(4.23000000);
  parameter[DC_5TE].push_back(0.00066500);
  parameter[DC_5TE].push_back(0.92103900);
  parameter[DC_5TE].push_back(0.08064600);
  parameter[DC_5TE].push_back(-0.39034900);
  parameter[DC_5TE].push_back(0.12417600);
  parameter[DC_5TE].push_back(-0.01120600);

  parameter[DC_TE3].push_back(2.87400000);
  parameter[DC_TE3].push_back(-0.05235500);
  parameter[DC_TE3].push_back(13.09201200);
  parameter[DC_TE3].push_back(-9.48128200);
  parameter[DC_TE3].push_back(-0.14958600);
  parameter[DC_TE3].push_back(1.75537200);
  parameter[DC_TE3].push_back(-0.39347500);

  parameter[DC_TE5].push_back(8.03600000);
  parameter[DC_TE5].push_back(-0.00513600);
  parameter[DC_TE5].push_back(4.67705700);
  parameter[DC_TE5].push_back(0.48333300);
  parameter[DC_TE5].push_back(-6.34511000);
  parameter[DC_TE5].push_back(3.83388500);
  parameter[DC_TE5].push_back(-0.67367800);

  parameter[DG_BB1].push_back(32.88500000);
  parameter[DG_BB1].push_back(0.08182900);
  parameter[DG_BB1].push_back(-7.32133900);
  parameter[DG_BB1].push_back(2.15767900);
  parameter[DG_BB1].push_back(-3.52369700);
  parameter[DG_BB1].push_back(2.30839600);
  parameter[DG_BB1].push_back(-0.39348300);

  parameter[DG_BB2].push_back(3.80600000);
  parameter[DG_BB2].push_back(-0.10618100);
  parameter[DG_BB2].push_back(9.54169000);
  parameter[DG_BB2].push_back(-6.15177600);
  parameter[DG_BB2].push_back(-0.53462400);
  parameter[DG_BB2].push_back(1.15581300);
  parameter[DG_BB2].push_back(-0.21567000);

  parameter[DG_BB3].push_back(-1.35600000);
  parameter[DG_BB3].push_back(0.57489100);
  parameter[DG_BB3].push_back(6.75164700);
  parameter[DG_BB3].push_back(4.11300900);
  parameter[DG_BB3].push_back(-9.63394600);
  parameter[DG_BB3].push_back(4.41675400);
  parameter[DG_BB3].push_back(-0.64339900);

  parameter[DG_SC1].push_back(6.67100000);
  parameter[DG_SC1].push_back(-0.00886600);
  parameter[DG_SC1].push_back(1.63333000);
  parameter[DG_SC1].push_back(-0.06892100);
  parameter[DG_SC1].push_back(-1.48683500);
  parameter[DG_SC1].push_back(0.78670800);
  parameter[DG_SC1].push_back(-0.12113900);

  parameter[DG_SC2].push_back(11.39400000);
  parameter[DG_SC2].push_back(0.00907900);
  parameter[DG_SC2].push_back(-0.22475500);
  parameter[DG_SC2].push_back(0.49535100);
  parameter[DG_SC2].push_back(-1.75324900);
  parameter[DG_SC2].push_back(0.98767400);
  parameter[DG_SC2].push_back(-0.16150800);

  parameter[DG_SC3].push_back(10.90100000);
  parameter[DG_SC3].push_back(0.02207600);
  parameter[DG_SC3].push_back(0.17932200);
  parameter[DG_SC3].push_back(0.73253200);
  parameter[DG_SC3].push_back(-1.95554900);
  parameter[DG_SC3].push_back(0.98339900);
  parameter[DG_SC3].push_back(-0.14763600);

  parameter[DG_SC4].push_back(6.45900000);
  parameter[DG_SC4].push_back(0.02018400);
  parameter[DG_SC4].push_back(4.17705400);
  parameter[DG_SC4].push_back(0.98531700);
  parameter[DG_SC4].push_back(-5.04354900);
  parameter[DG_SC4].push_back(2.56123700);
  parameter[DG_SC4].push_back(-0.39249300);

  parameter[DG_3TE].push_back(4.23000000);
  parameter[DG_3TE].push_back(0.00061700);
  parameter[DG_3TE].push_back(0.92140100);
  parameter[DG_3TE].push_back(0.08016400);
  parameter[DG_3TE].push_back(-0.39003500);
  parameter[DG_3TE].push_back(0.12406900);
  parameter[DG_3TE].push_back(-0.01119200);

  parameter[DG_5TE].push_back(4.23000000);
  parameter[DG_5TE].push_back(0.00064900);
  parameter[DG_5TE].push_back(0.92110500);
  parameter[DG_5TE].push_back(0.08031500);
  parameter[DG_5TE].push_back(-0.38997000);
  parameter[DG_5TE].push_back(0.12401200);
  parameter[DG_5TE].push_back(-0.01118100);

  parameter[DG_TE3].push_back(2.87400000);
  parameter[DG_TE3].push_back(0.00182000);
  parameter[DG_TE3].push_back(12.41507000);
  parameter[DG_TE3].push_back(-7.47384800);
  parameter[DG_TE3].push_back(-2.11864700);
  parameter[DG_TE3].push_back(2.50112600);
  parameter[DG_TE3].push_back(-0.48652200);

  parameter[DG_TE5].push_back(8.03600000);
  parameter[DG_TE5].push_back(0.00676400);
  parameter[DG_TE5].push_back(4.65989200);
  parameter[DG_TE5].push_back(0.78482500);
  parameter[DG_TE5].push_back(-6.86460600);
  parameter[DG_TE5].push_back(4.11675400);
  parameter[DG_TE5].push_back(-0.72249100);

  parameter[DT_BB1].push_back(32.88500000);
  parameter[DT_BB1].push_back(0.08220100);
  parameter[DT_BB1].push_back(-7.33006800);
  parameter[DT_BB1].push_back(2.16636500);
  parameter[DT_BB1].push_back(-3.53465700);
  parameter[DT_BB1].push_back(2.31447600);
  parameter[DT_BB1].push_back(-0.39445400);

  parameter[DT_BB2].push_back(3.80600000);
  parameter[DT_BB2].push_back(-0.10723000);
  parameter[DT_BB2].push_back(9.56675000);
  parameter[DT_BB2].push_back(-6.20236100);
  parameter[DT_BB2].push_back(-0.49550400);
  parameter[DT_BB2].push_back(1.14300600);
  parameter[DT_BB2].push_back(-0.21420000);

  parameter[DT_BB3].push_back(-1.35600000);
  parameter[DT_BB3].push_back(0.56737900);
  parameter[DT_BB3].push_back(6.76595400);
  parameter[DT_BB3].push_back(4.08976100);
  parameter[DT_BB3].push_back(-9.61512500);
  parameter[DT_BB3].push_back(4.40975100);
  parameter[DT_BB3].push_back(-0.64239800);

  parameter[DT_SC1].push_back(5.95100000);
  parameter[DT_SC1].push_back(-0.02926500);
  parameter[DT_SC1].push_back(2.59630300);
  parameter[DT_SC1].push_back(-0.56152200);
  parameter[DT_SC1].push_back(-1.56532600);
  parameter[DT_SC1].push_back(0.89322800);
  parameter[DT_SC1].push_back(-0.14142900);

  parameter[DT_SC2].push_back(10.90100000);
  parameter[DT_SC2].push_back(0.02183400);
  parameter[DT_SC2].push_back(0.19463000);
  parameter[DT_SC2].push_back(0.72393000);
  parameter[DT_SC2].push_back(-1.93199500);
  parameter[DT_SC2].push_back(0.96856300);
  parameter[DT_SC2].push_back(-0.14512600);

  parameter[DT_SC3].push_back(4.31400000);
  parameter[DT_SC3].push_back(-0.07745600);
  parameter[DT_SC3].push_back(12.49820300);
  parameter[DT_SC3].push_back(-7.64994200);
  parameter[DT_SC3].push_back(-3.00359600);
  parameter[DT_SC3].push_back(3.26263300);
  parameter[DT_SC3].push_back(-0.64498600);

  parameter[DT_3TE].push_back(4.23000000);
  parameter[DT_3TE].push_back(0.00062000);
  parameter[DT_3TE].push_back(0.92141100);
  parameter[DT_3TE].push_back(0.08030900);
  parameter[DT_3TE].push_back(-0.39021500);
  parameter[DT_3TE].push_back(0.12414000);
  parameter[DT_3TE].push_back(-0.01120100);

  parameter[DT_5TE].push_back(4.23000000);
  parameter[DT_5TE].push_back(0.00063700);
  parameter[DT_5TE].push_back(0.92130800);
  parameter[DT_5TE].push_back(0.08026900);
  parameter[DT_5TE].push_back(-0.39007500);
  parameter[DT_5TE].push_back(0.12406600);
  parameter[DT_5TE].push_back(-0.01118800);

  parameter[DT_TE3].push_back(2.87400000);
  parameter[DT_TE3].push_back(-0.00251200);
  parameter[DT_TE3].push_back(12.43576400);
  parameter[DT_TE3].push_back(-7.55343800);
  parameter[DT_TE3].push_back(-2.07363500);
  parameter[DT_TE3].push_back(2.51279300);
  parameter[DT_TE3].push_back(-0.49437100);

  parameter[DT_TE5].push_back(8.03600000);
  parameter[DT_TE5].push_back(0.00119900);
  parameter[DT_TE5].push_back(4.91762300);
  parameter[DT_TE5].push_back(0.65637000);
  parameter[DT_TE5].push_back(-7.23392500);
  parameter[DT_TE5].push_back(4.44636600);
  parameter[DT_TE5].push_back(-0.79467800);

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( moldat ) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {
      std::string Aname = moldat->getAtomName(atoms[i]);
      std::string Rname = moldat->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        if(Aname=="BB") {
          atoi[i]=ALA_BB;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ARG") {
        if(Aname=="BB") {
          atoi[i]=ARG_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ARG_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=ARG_SC2;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ASN") {
        if(Aname=="BB") {
          atoi[i]=ASN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASN_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ASP") {
        if(Aname=="BB") {
          atoi[i]=ASP_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASP_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="CYS") {
        if(Aname=="BB") {
          atoi[i]=CYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=CYS_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLN") {
        if(Aname=="BB") {
          atoi[i]=GLN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLN_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLU") {
        if(Aname=="BB") {
          atoi[i]=GLU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLU_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          atoi[i]=GLY_BB;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="HIS") {
        if(Aname=="BB") {
          atoi[i]=HIS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=HIS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=HIS_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=HIS_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="ILE") {
        if(Aname=="BB") {
          atoi[i]=ILE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ILE_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="LEU") {
        if(Aname=="BB") {
          atoi[i]=LEU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LEU_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="LYS") {
        if(Aname=="BB") {
          atoi[i]=LYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LYS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=LYS_SC2;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="MET") {
        if(Aname=="BB") {
          atoi[i]=MET_BB;
        } else if(Aname=="SC1") {
          atoi[i]=MET_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="PHE") {
        if(Aname=="BB") {
          atoi[i]=PHE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PHE_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=PHE_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=PHE_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="PRO") {
        if(Aname=="BB") {
          atoi[i]=PRO_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PRO_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="SER") {
        if(Aname=="BB") {
          atoi[i]=SER_BB;
        } else if(Aname=="SC1") {
          atoi[i]=SER_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="THR") {
        if(Aname=="BB") {
          atoi[i]=THR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=THR_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="TRP") {
        if(Aname=="BB") {
          atoi[i]=TRP_BB;
        } else if(Aname=="SC1") {
          atoi[i]=TRP_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=TRP_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=TRP_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=TRP_SC4;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="TYR") {
        if(Aname=="BB") {
          atoi[i]=TYR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=TYR_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=TYR_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=TYR_SC3;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="VAL") {
        if(Aname=="BB") {
          atoi[i]=VAL_BB;
        } else if(Aname=="SC1") {
          atoi[i]=VAL_SC1;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  A") {
        if(Aname=="BB1") {
          atoi[i]=A_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=A_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=A_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=A_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=A_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=A_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=A_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=A_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=A_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=A_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=A_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  C") {
        if(Aname=="BB1") {
          atoi[i]=C_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=C_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=C_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=C_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=C_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=C_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=C_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=C_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=C_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=C_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  G") {
        if(Aname=="BB1") {
          atoi[i]=G_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=G_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=G_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=G_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=G_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=G_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=G_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=G_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=G_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=G_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=G_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname=="  U") {
        if(Aname=="BB1") {
          atoi[i]=U_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=U_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=U_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=U_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=U_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=U_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=U_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=U_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=U_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=U_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DA") {
        if(Aname=="BB1") {
          atoi[i]=DA_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DA_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DA_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DA_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DA_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DA_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=DA_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=DA_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DA_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DA_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DA_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DC") {
        if(Aname=="BB1") {
          atoi[i]=DC_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DC_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DC_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DC_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DC_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DC_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=DC_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DC_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DC_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DC_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DG") {
        if(Aname=="BB1") {
          atoi[i]=DG_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DG_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DG_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DG_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DG_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DG_SC3;
        } else if(Aname=="SC4") {
          atoi[i]=DG_SC4;
        } else if(Aname=="3TE") {
          atoi[i]=DG_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DG_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DG_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DG_TE5;
        } else error("Atom name not known: "+Aname);
      } else if(Rname==" DT") {
        if(Aname=="BB1") {
          atoi[i]=DT_BB1;
        } else if(Aname=="BB2") {
          atoi[i]=DT_BB2;
        } else if(Aname=="BB3") {
          atoi[i]=DT_BB3;
        } else if(Aname=="SC1") {
          atoi[i]=DT_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=DT_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=DT_SC3;
        } else if(Aname=="3TE") {
          atoi[i]=DT_3TE;
        } else if(Aname=="5TE") {
          atoi[i]=DT_5TE;
        } else if(Aname=="TE3") {
          atoi[i]=DT_TE3;
        } else if(Aname=="TE5") {
          atoi[i]=DT_TE5;
        } else error("Atom name not known: "+Aname);
      } else error("Residue not known: "+Rname);
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

void SAXS::getOnebeadparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac, std::vector<std::vector<long double> > &parameter_mix, std::vector<std::vector<long double> > &parameter_wat)
{

  parameter_wat[TRP].push_back(60737.602499880035);
  parameter_wat[TRP].push_back(423.3538118566174);
  parameter_wat[TRP].push_back(-225729.57942426094);
  parameter_wat[TRP].push_back(-28863.96508764871);
  parameter_wat[TRP].push_back(703953.2853844669);
  parameter_wat[TRP].push_back(-750290.6360569759);
  parameter_wat[TRP].push_back(243441.79433046127);

  parameter_wat[TYR].push_back(46250.803600118066);
  parameter_wat[TYR].push_back(-257.310522603269);
  parameter_wat[TYR].push_back(-132302.7898977298);
  parameter_wat[TYR].push_back(-50352.95940852235);
  parameter_wat[TYR].push_back(406079.60996746755);
  parameter_wat[TYR].push_back(-375410.2955043933);
  parameter_wat[TYR].push_back(109181.47611497746);

  parameter_wat[PHE].push_back(42407.16489988089);
  parameter_wat[PHE].push_back(-101.43687190346151);
  parameter_wat[PHE].push_back(-127120.33522573594);
  parameter_wat[PHE].push_back(-38973.85574662809);
  parameter_wat[PHE].push_back(391846.2601350714);
  parameter_wat[PHE].push_back(-380079.9501473576);
  parameter_wat[PHE].push_back(115269.28702049167);

  parameter_wat[HIS].push_back(22888.66410011907);
  parameter_wat[HIS].push_back(-108.02829077408032);
  parameter_wat[HIS].push_back(-60158.91020245944);
  parameter_wat[HIS].push_back(-21412.911926812973);
  parameter_wat[HIS].push_back(175284.75094781365);
  parameter_wat[HIS].push_back(-160903.34874783095);
  parameter_wat[HIS].push_back(46623.44919862075);

  parameter_wat[HIP].push_back(24473.47360011933);
  parameter_wat[HIP].push_back(-113.35936100565031);
  parameter_wat[HIP].push_back(-65665.4454282035);
  parameter_wat[HIP].push_back(-23342.641495874075);
  parameter_wat[HIP].push_back(193991.90099311847);
  parameter_wat[HIP].push_back(-179359.8105131197);
  parameter_wat[HIP].push_back(52286.59943451184);

  parameter_wat[HIE].push_back(22888.66410011907);
  parameter_wat[HIE].push_back(-109.76722568752795);
  parameter_wat[HIE].push_back(-60653.77965364445);
  parameter_wat[HIE].push_back(-21845.181477770184);
  parameter_wat[HIE].push_back(178618.98076057335);
  parameter_wat[HIE].push_back(-164470.42658519786);
  parameter_wat[HIE].push_back(47778.241548568774);

  parameter_wat[ARG].push_back(34106.702399880836);
  parameter_wat[ARG].push_back(-31.968033811625222);
  parameter_wat[ARG].push_back(-108747.89030377955);
  parameter_wat[ARG].push_back(-30958.904014258027);
  parameter_wat[ARG].push_back(347163.4357312966);
  parameter_wat[ARG].push_back(-346276.2882554062);
  parameter_wat[ARG].push_back(107326.24946733958);

  parameter_wat[LYS].push_back(32292.089999880736);
  parameter_wat[LYS].push_back(-52.96029474662974);
  parameter_wat[LYS].push_back(-95909.23116199141);
  parameter_wat[LYS].push_back(-28184.175529381177);
  parameter_wat[LYS].push_back(297836.9203652572);
  parameter_wat[LYS].push_back(-293011.6716266663);
  parameter_wat[LYS].push_back(89967.65334564191);

  parameter_wat[CYS].push_back(11352.902500119093);
  parameter_wat[CYS].push_back(-42.586554888239874);
  parameter_wat[CYS].push_back(-20366.230938482615);
  parameter_wat[CYS].push_back(-5137.392457690132);
  parameter_wat[CYS].push_back(34975.323854506074);
  parameter_wat[CYS].push_back(-24243.22374505194);
  parameter_wat[CYS].push_back(5152.472178010549);

  parameter_wat[ASP].push_back(13511.73760011934);
  parameter_wat[ASP].push_back(-62.259809313175644);
  parameter_wat[ASP].push_back(-26335.159072036455);
  parameter_wat[ASP].push_back(-8002.9679940827655);
  parameter_wat[ASP].push_back(54030.09941724955);
  parameter_wat[ASP].push_back(-41082.18951195315);
  parameter_wat[ASP].push_back(9822.373386145666);

  parameter_wat[GLU].push_back(20443.28040011943);
  parameter_wat[GLU].push_back(-114.65211777672742);
  parameter_wat[GLU].push_back(-45705.214663013714);
  parameter_wat[GLU].push_back(-16305.612141085105);
  parameter_wat[GLU].push_back(113297.5758640635);
  parameter_wat[GLU].push_back(-93999.34286496713);
  parameter_wat[GLU].push_back(24699.32337589366);

  parameter_wat[ILE].push_back(27858.9481001196);
  parameter_wat[ILE].push_back(-160.09616278302047);
  parameter_wat[ILE].push_back(-61355.72128368112);
  parameter_wat[ILE].push_back(-21237.083952478268);
  parameter_wat[ILE].push_back(142633.447400751);
  parameter_wat[ILE].push_back(-113273.08967375412);
  parameter_wat[ILE].push_back(28319.143796630207);

  parameter_wat[LEU].push_back(27858.9481001196);
  parameter_wat[LEU].push_back(-164.2242755261456);
  parameter_wat[LEU].push_back(-62269.57074002018);
  parameter_wat[LEU].push_back(-22164.04215173996);
  parameter_wat[LEU].push_back(149446.70227598862);
  parameter_wat[LEU].push_back(-120581.31029609666);
  parameter_wat[LEU].push_back(30683.935787436156);

  parameter_wat[MET].push_back(25609.60090011981);
  parameter_wat[MET].push_back(-147.12983069832867);
  parameter_wat[MET].push_back(-65621.28306656817);
  parameter_wat[MET].push_back(-24719.313103059052);
  parameter_wat[MET].push_back(186430.47475376664);
  parameter_wat[MET].push_back(-165841.9203585998);
  parameter_wat[MET].push_back(46626.915939594095);

  parameter_wat[ASN].push_back(14376.010000119086);
  parameter_wat[ASN].push_back(-67.13662191717502);
  parameter_wat[ASN].push_back(-28168.8210696101);
  parameter_wat[ASN].push_back(-8448.880455615134);
  parameter_wat[ASN].push_back(56555.863212460674);
  parameter_wat[ASN].push_back(-42203.49388869736);
  parameter_wat[ASN].push_back(9840.827095275834);

  parameter_wat[PRO].push_back(16866.21690011945);
  parameter_wat[PRO].push_back(-77.79313454249083);
  parameter_wat[PRO].push_back(-32625.357427339673);
  parameter_wat[PRO].push_back(-9843.651473353346);
  parameter_wat[PRO].push_back(65875.36925812009);
  parameter_wat[PRO].push_back(-49499.88311947784);
  parameter_wat[PRO].push_back(11656.375656878794);

  parameter_wat[GLN].push_back(21503.289600119013);
  parameter_wat[GLN].push_back(-130.80754689004655);
  parameter_wat[GLN].push_back(-49582.982827439504);
  parameter_wat[GLN].push_back(-18648.86453544889);
  parameter_wat[GLN].push_back(128172.37310395157);
  parameter_wat[GLN].push_back(-107594.15890982642);
  parameter_wat[GLN].push_back(28581.982123379872);

  parameter_wat[SER].push_back(9181.472400119355);
  parameter_wat[SER].push_back(-29.299114175405876);
  parameter_wat[SER].push_back(-15273.60177234209);
  parameter_wat[SER].push_back(-3450.5500428922705);
  parameter_wat[SER].push_back(23773.483013633984);
  parameter_wat[SER].push_back(-15717.798552926966);
  parameter_wat[SER].push_back(3130.6799238799476);

  parameter_wat[THR].push_back(15020.9536001194);
  parameter_wat[THR].push_back(-59.886795759067695);
  parameter_wat[THR].push_back(-27450.363696477114);
  parameter_wat[THR].push_back(-7174.085244900661);
  parameter_wat[THR].push_back(48176.41855470643);
  parameter_wat[THR].push_back(-33484.54045901642);
  parameter_wat[THR].push_back(7115.493498450939);

  parameter_wat[VAL].push_back(19647.62890011937);
  parameter_wat[VAL].push_back(-90.37805699411572);
  parameter_wat[VAL].push_back(-38295.814997790214);
  parameter_wat[VAL].push_back(-11134.637603308249);
  parameter_wat[VAL].push_back(74170.4723318485);
  parameter_wat[VAL].push_back(-54051.48501366292);
  parameter_wat[VAL].push_back(12214.226839424213);

  parameter_wat[ALA].push_back(7515.156100119268);
  parameter_wat[ALA].push_back(-19.745812003651533);
  parameter_wat[ALA].push_back(-11691.3749145515);
  parameter_wat[ALA].push_back(-2278.6778311610356);
  parameter_wat[ALA].push_back(16179.890625578142);
  parameter_wat[ALA].push_back(-10070.196793821578);
  parameter_wat[ALA].push_back(1832.391603539423);

  parameter_wat[GLY].push_back(3594.002500119161);
  parameter_wat[GLY].push_back(-6.783782713939998);
  parameter_wat[GLY].push_back(-4918.1541249080265);
  parameter_wat[GLY].push_back(-769.4829100052555);
  parameter_wat[GLY].push_back(5761.026162132836);
  parameter_wat[GLY].push_back(-3312.0647071906815);
  parameter_wat[GLY].push_back(531.6151440835663);

  parameter_mix[TRP].push_back(48294.011756881504);
  parameter_mix[TRP].push_back(94.79679625134271);
  parameter_mix[TRP].push_back(-162545.8333320417);
  parameter_mix[TRP].push_back(-37693.98947041115);
  parameter_mix[TRP].push_back(526423.3675952561);
  parameter_mix[TRP].push_back(-543242.8087968872);
  parameter_mix[TRP].push_back(172770.31888199487);

  parameter_mix[TYR].push_back(36984.202403359115);
  parameter_mix[TYR].push_back(-222.0017930602204);
  parameter_mix[TYR].push_back(-99488.96547634256);
  parameter_mix[TYR].push_back(-39490.02998515562);
  parameter_mix[TYR].push_back(302058.8500035175);
  parameter_mix[TYR].push_back(-275632.27775690623);
  parameter_mix[TYR].push_back(79286.9486152365);

  parameter_mix[PHE].push_back(32119.469231338207);
  parameter_mix[PHE].push_back(-140.9766666035606);
  parameter_mix[PHE].push_back(-88429.06372052178);
  parameter_mix[PHE].push_back(-32440.111797288962);
  parameter_mix[PHE].push_back(276747.9437420169);
  parameter_mix[PHE].push_back(-262851.9318112437);
  parameter_mix[PHE].push_back(78470.28300750244);

  parameter_mix[HIS].push_back(21779.124723299235);
  parameter_mix[HIS].push_back(-125.81547101114079);
  parameter_mix[HIS].push_back(-51431.05749165743);
  parameter_mix[HIS].push_back(-20013.224366191025);
  parameter_mix[HIS].push_back(144694.92222746916);
  parameter_mix[HIS].push_back(-127819.49538087688);
  parameter_mix[HIS].push_back(35763.772645449855);

  parameter_mix[HIP].push_back(22833.36414923898);
  parameter_mix[HIP].push_back(-134.11759532585418);
  parameter_mix[HIP].push_back(-55189.4285078936);
  parameter_mix[HIP].push_back(-21799.88197554733);
  parameter_mix[HIP].push_back(158958.1762264112);
  parameter_mix[HIP].push_back(-141751.7818363874);
  parameter_mix[HIP].push_back(39997.324270833735);

  parameter_mix[HIE].push_back(21779.124723299235);
  parameter_mix[HIE].push_back(-129.6552781944993);
  parameter_mix[HIE].push_back(-51660.59291329289);
  parameter_mix[HIE].push_back(-20412.70089226015);
  parameter_mix[HIE].push_back(146341.19036764107);
  parameter_mix[HIE].push_back(-129182.83705926737);
  parameter_mix[HIE].push_back(36101.78353079046);

  parameter_mix[ARG].push_back(31385.401600920548);
  parameter_mix[ARG].push_back(-107.823749853272);
  parameter_mix[ARG].push_back(-96081.57737897056);
  parameter_mix[ARG].push_back(-34356.86448428557);
  parameter_mix[ARG].push_back(321119.58481540927);
  parameter_mix[ARG].push_back(-315674.66683333775);
  parameter_mix[ARG].push_back(96796.44214985141);

  parameter_mix[LYS].push_back(25511.358126718373);
  parameter_mix[LYS].push_back(-95.24657818372829);
  parameter_mix[LYS].push_back(-72472.87199806623);
  parameter_mix[LYS].push_back(-26437.446224277774);
  parameter_mix[LYS].push_back(237422.12523083566);
  parameter_mix[LYS].push_back(-231421.26921744435);
  parameter_mix[LYS].push_back(70582.85782972853);

  parameter_mix[CYS].push_back(11505.517261618916);
  parameter_mix[CYS].push_back(-28.00657295704436);
  parameter_mix[CYS].push_back(-17568.024889706485);
  parameter_mix[CYS].push_back(-3208.4239567270556);
  parameter_mix[CYS].push_back(23190.8240650359);
  parameter_mix[CYS].push_back(-14025.359913209422);
  parameter_mix[CYS].push_back(2432.996124564235);

  parameter_mix[ASP].push_back(13713.858501879387);
  parameter_mix[ASP].push_back(-56.16981725896633);
  parameter_mix[ASP].push_back(-24467.359478081828);
  parameter_mix[ASP].push_back(-6842.625685319477);
  parameter_mix[ASP].push_back(45409.799932497255);
  parameter_mix[ASP].push_back(-32675.494006878293);
  parameter_mix[ASP].push_back(7282.116459231734);

  parameter_mix[GLU].push_back(19156.03660739948);
  parameter_mix[GLU].push_back(-112.77021473479914);
  parameter_mix[GLU].push_back(-40559.88365261396);
  parameter_mix[GLU].push_back(-14896.089751789212);
  parameter_mix[GLU].push_back(97924.20440282047);
  parameter_mix[GLU].push_back(-78861.69840848455);
  parameter_mix[GLU].push_back(20003.69595126275);

  parameter_mix[ILE].push_back(20693.062159179222);
  parameter_mix[ILE].push_back(-101.97290141036149);
  parameter_mix[ILE].push_back(-40831.78703482502);
  parameter_mix[ILE].push_back(-12657.594458075528);
  parameter_mix[ILE].push_back(83296.23989690492);
  parameter_mix[ILE].push_back(-62017.47308548058);
  parameter_mix[ILE].push_back(14368.398225275843);

  parameter_mix[LEU].push_back(20693.062159179222);
  parameter_mix[LEU].push_back(-111.80857619890668);
  parameter_mix[LEU].push_back(-42131.62507792151);
  parameter_mix[LEU].push_back(-14009.106219647874);
  parameter_mix[LEU].push_back(91220.93610257414);
  parameter_mix[LEU].push_back(-69456.39325237977);
  parameter_mix[LEU].push_back(16496.74666112724);

  parameter_mix[MET].push_back(22400.800002738917);
  parameter_mix[MET].push_back(-133.76891528562035);
  parameter_mix[MET].push_back(-50562.14461060515);
  parameter_mix[MET].push_back(-19092.463028952076);
  parameter_mix[MET].push_back(131078.39403549655);
  parameter_mix[MET].push_back(-110255.40033841228);
  parameter_mix[MET].push_back(29347.828056519662);

  parameter_mix[ASN].push_back(14384.287416519466);
  parameter_mix[ASN].push_back(-54.368319093741725);
  parameter_mix[ASN].push_back(-25264.7008339227);
  parameter_mix[ASN].push_back(-6505.720384021238);
  parameter_mix[ASN].push_back(43726.56856322167);
  parameter_mix[ASN].push_back(-30312.120416038917);
  parameter_mix[ASN].push_back(6426.928446443609);

  parameter_mix[PRO].push_back(13503.79714565913);
  parameter_mix[PRO].push_back(-48.81691245595872);
  parameter_mix[PRO].push_back(-22749.678192064708);
  parameter_mix[PRO].push_back(-5863.718046448816);
  parameter_mix[PRO].push_back(39264.569464993685);
  parameter_mix[PRO].push_back(-27454.474446209842);
  parameter_mix[PRO].push_back(5909.564112806066);

  parameter_mix[GLN].push_back(19938.237246839024);
  parameter_mix[GLN].push_back(-125.86788493612829);
  parameter_mix[GLN].push_back(-43384.69315271189);
  parameter_mix[GLN].push_back(-16673.856129529086);
  parameter_mix[GLN].push_back(108716.8543788976);
  parameter_mix[GLN].push_back(-88501.68147616273);
  parameter_mix[GLN].push_back(22682.521341667303);

  parameter_mix[SER].push_back(8813.67020471935);
  parameter_mix[SER].push_back(-19.037976334608356);
  parameter_mix[SER].push_back(-12699.254690509915);
  parameter_mix[SER].push_back(-2154.45659945185);
  parameter_mix[SER].push_back(15784.762074798347);
  parameter_mix[SER].push_back(-9235.432880294426);
  parameter_mix[SER].push_back(1509.3565114265857);

  parameter_mix[THR].push_back(13233.997179639066);
  parameter_mix[THR].push_back(-37.56509776445388);
  parameter_mix[THR].push_back(-21203.56789307223);
  parameter_mix[THR].push_back(-4308.144424843769);
  parameter_mix[THR].push_back(30208.83306033092);
  parameter_mix[THR].push_back(-18848.577186687486);
  parameter_mix[THR].push_back(3411.913197678857);

  parameter_mix[VAL].push_back(15135.438016299164);
  parameter_mix[VAL].push_back(-52.635094911968245);
  parameter_mix[VAL].push_back(-26030.01153963802);
  parameter_mix[VAL].push_back(-6179.950964589099);
  parameter_mix[VAL].push_back(42022.56035826398);
  parameter_mix[VAL].push_back(-28010.86180138724);
  parameter_mix[VAL].push_back(5612.649276643284);

  parameter_mix[ALA].push_back(6586.942863819189);
  parameter_mix[ALA].push_back(-10.500985041159272);
  parameter_mix[ALA].push_back(-8682.564064794862);
  parameter_mix[ALA].push_back(-1167.7140089940317);
  parameter_mix[ALA].push_back(9157.822632690133);
  parameter_mix[ALA].push_back(-4869.150402902164);
  parameter_mix[ALA].push_back(660.7127436699247);

  parameter_mix[GLY].push_back(3596.0718542192735);
  parameter_mix[GLY].push_back(-3.872148741667527);
  parameter_mix[GLY].push_back(-4049.533333393464);
  parameter_mix[GLY].push_back(-426.8550732879424);
  parameter_mix[GLY].push_back(3595.6814286066124);
  parameter_mix[GLY].push_back(-1755.3404133215586);
  parameter_mix[GLY].push_back(199.02655960598142);

  parameter_vac[TRP].push_back(9599.949107129776);
  parameter_vac[TRP].push_back(-30.858609190810583);
  parameter_vac[TRP].push_back(-28551.048759771413);
  parameter_vac[TRP].push_back(-10231.240891996395);
  parameter_vac[TRP].push_back(96031.95748027436);
  parameter_vac[TRP].push_back(-94814.69691650089);
  parameter_vac[TRP].push_back(29200.454883200917);

  parameter_vac[TYR].push_back(7393.553846413058);
  parameter_vac[TYR].push_back(-47.23837907229961);
  parameter_vac[TYR].push_back(-18490.21435843763);
  parameter_vac[TYR].push_back(-7737.360235562757);
  parameter_vac[TYR].push_back(56063.21994893703);
  parameter_vac[TYR].push_back(-50556.859582519006);
  parameter_vac[TYR].push_back(14389.081569421436);

  parameter_vac[PHE].push_back(6081.874997705282);
  parameter_vac[PHE].push_back(-38.75296096284759);
  parameter_vac[PHE].push_back(-14795.200710070538);
  parameter_vac[PHE].push_back(-6406.629771873036);
  parameter_vac[PHE].push_back(46100.140106722465);
  parameter_vac[PHE].push_back(-42085.36076230338);
  parameter_vac[PHE].push_back(12122.907929159559);

  parameter_vac[HIS].push_back(5180.842705000208);
  parameter_vac[HIS].push_back(-32.89736234246629);
  parameter_vac[HIS].push_back(-10873.350997276972);
  parameter_vac[HIS].push_back(-4481.4769197362);
  parameter_vac[HIS].push_back(29212.72548906728);
  parameter_vac[HIS].push_back(-24529.701112742892);
  parameter_vac[HIS].push_back(6495.940243992503);

  parameter_vac[HIP].push_back(5325.791987063724);
  parameter_vac[HIP].push_back(-35.05625337809874);
  parameter_vac[HIP].push_back(-11441.518343889265);
  parameter_vac[HIP].push_back(-4859.835770461219);
  parameter_vac[HIP].push_back(31828.3405173923);
  parameter_vac[HIP].push_back(-27109.903007338115);
  parameter_vac[HIP].push_back(7283.467569951606);

  parameter_vac[HIE].push_back(5180.842705000208);
  parameter_vac[HIE].push_back(-33.29984675748154);
  parameter_vac[HIE].push_back(-10872.379635800535);
  parameter_vac[HIE].push_back(-4493.213445493111);
  parameter_vac[HIE].push_back(29107.17290177237);
  parameter_vac[HIE].push_back(-24328.802993028898);
  parameter_vac[HIE].push_back(6407.227555283756);

  parameter_vac[ARG].push_back(7220.306892486492);
  parameter_vac[ARG].push_back(-45.02676906795456);
  parameter_vac[ARG].push_back(-20987.6312443252);
  parameter_vac[ARG].push_back(-9443.110444762371);
  parameter_vac[ARG].push_back(74884.77524133732);
  parameter_vac[ARG].push_back(-72540.00466456635);
  parameter_vac[ARG].push_back(21959.97615859096);

  parameter_vac[LYS].push_back(5038.613120729024);
  parameter_vac[LYS].push_back(-30.37838672252398);
  parameter_vac[LYS].push_back(-13422.897999177683);
  parameter_vac[LYS].push_back(-6053.389263987596);
  parameter_vac[LYS].push_back(46835.63074102061);
  parameter_vac[LYS].push_back(-45089.04698246694);
  parameter_vac[LYS].push_back(13600.678970015228);

  parameter_vac[CYS].push_back(2915.0458981763995);
  parameter_vac[CYS].push_back(-3.595827832254589);
  parameter_vac[CYS].push_back(-3592.4507607863125);
  parameter_vac[CYS].push_back(-394.23885715070753);
  parameter_vac[CYS].push_back(3348.1185401489765);
  parameter_vac[CYS].push_back(-1606.9293130024218);
  parameter_vac[CYS].push_back(166.6775491017609);

  parameter_vac[ASP].push_back(3479.7507282248966);
  parameter_vac[ASP].push_back(-12.005347430085758);
  parameter_vac[ASP].push_back(-5600.883991352643);
  parameter_vac[ASP].push_back(-1389.3777545609155);
  parameter_vac[ASP].push_back(9245.650898756234);
  parameter_vac[ASP].push_back(-6156.207160613609);
  parameter_vac[ASP].push_back(1218.5833934642021);

  parameter_vac[GLU].push_back(4487.461543955492);
  parameter_vac[GLU].push_back(-27.373885533538896);
  parameter_vac[GLU].push_back(-8941.928757023);
  parameter_vac[GLU].push_back(-3371.688026217654);
  parameter_vac[GLU].push_back(21083.486308601627);
  parameter_vac[GLU].push_back(-16307.701537949622);
  parameter_vac[GLU].push_back(3909.6886537322757);

  parameter_vac[ILE].push_back(3842.5968201937767);
  parameter_vac[ILE].push_back(-13.4621492233332);
  parameter_vac[ILE].push_back(-6414.466724440508);
  parameter_vac[ILE].push_back(-1580.0124125764848);
  parameter_vac[ILE].push_back(10618.037884226278);
  parameter_vac[ILE].push_back(-7154.041016270302);
  parameter_vac[ILE].push_back(1452.4947569272492);

  parameter_vac[LEU].push_back(3842.5968201937767);
  parameter_vac[LEU].push_back(-15.681457814419176);
  parameter_vac[LEU].push_back(-6771.789004604622);
  parameter_vac[LEU].push_back(-1845.1900416951112);
  parameter_vac[LEU].push_back(12139.114106034092);
  parameter_vac[LEU].push_back(-8383.262864959117);
  parameter_vac[LEU].push_back(1747.4810835193168);

  parameter_vac[MET].push_back(4898.51289296739);
  parameter_vac[MET].push_back(-25.18325080293945);
  parameter_vac[MET].push_back(-9467.48983449123);
  parameter_vac[MET].push_back(-3202.202778951768);
  parameter_vac[MET].push_back(20961.877825984124);
  parameter_vac[MET].push_back(-16179.791857572387);
  parameter_vac[MET].push_back(3909.4934716887474);

  parameter_vac[ASN].push_back(3598.142399811549);
  parameter_vac[ASN].push_back(-10.193740412650593);
  parameter_vac[ASN].push_back(-5570.433196260574);
  parameter_vac[ASN].push_back(-1166.3919144288564);
  parameter_vac[ASN].push_back(8108.239221259332);
  parameter_vac[ASN].push_back(-5084.052169945779);
  parameter_vac[ASN].push_back(922.5084812557318);

  parameter_vac[PRO].push_back(2702.925890807491);
  parameter_vac[PRO].push_back(-6.344885973773661);
  parameter_vac[PRO].push_back(-3723.804313180885);
  parameter_vac[PRO].push_back(-725.0774829994464);
  parameter_vac[PRO].push_back(5083.307369521012);
  parameter_vac[PRO].push_back(-3153.6282584491305);
  parameter_vac[PRO].push_back(568.1054247577542);

  parameter_vac[GLN].push_back(4621.773132292558);
  parameter_vac[GLN].push_back(-29.611357037564705);
  parameter_vac[GLN].push_back(-9427.140597349113);
  parameter_vac[GLN].push_back(-3666.620711439403);
  parameter_vac[GLN].push_back(22875.745472246173);
  parameter_vac[GLN].push_back(-17863.872805261297);
  parameter_vac[GLN].push_back(4325.335756697777);

  parameter_vac[SER].push_back(2115.150465404398);
  parameter_vac[SER].push_back(-2.608961844911164);
  parameter_vac[SER].push_back(-2535.4783055939965);
  parameter_vac[SER].push_back(-285.24778630312363);
  parameter_vac[SER].push_back(2385.5343165459967);
  parameter_vac[SER].push_back(-1156.9360917654067);
  parameter_vac[SER].push_back(122.17917481954323);

  parameter_vac[THR].push_back(2914.906170715884);
  parameter_vac[THR].push_back(-4.775628956452242);
  parameter_vac[THR].push_back(-3871.750382748967);
  parameter_vac[THR].push_back(-527.0206564434785);
  parameter_vac[THR].push_back(4136.791177465494);
  parameter_vac[THR].push_back(-2169.9772449584666);
  parameter_vac[THR].push_back(279.42678356526335);

  parameter_vac[VAL].push_back(2914.8744247581994);
  parameter_vac[VAL].push_back(-6.0787188044112375);
  parameter_vac[VAL].push_back(-4122.984089016234);
  parameter_vac[VAL].push_back(-684.3239695392286);
  parameter_vac[VAL].push_back(5048.30676770002);
  parameter_vac[VAL].push_back(-2909.0035370830888);
  parameter_vac[VAL].push_back(460.72939731819946);

  parameter_vac[ALA].push_back(1443.3438146824442);
  parameter_vac[ALA].push_back(-1.079460116756403);
  parameter_vac[ALA].push_back(-1483.196871651843);
  parameter_vac[ALA].push_back(-116.4910603548198);
  parameter_vac[ALA].push_back(1109.6163998427003);
  parameter_vac[ALA].push_back(-462.4259514425249);
  parameter_vac[ALA].push_back(28.21768859752687);

  parameter_vac[GLY].push_back(899.5356000422943);
  parameter_vac[GLY].push_back(-0.4467373920268132);
  parameter_vac[GLY].push_back(-765.5891570102108);
  parameter_vac[GLY].push_back(-48.00532220136316);
  parameter_vac[GLY].push_back(493.4938562562617);
  parameter_vac[GLY].push_back(-189.2220341648864);
  parameter_vac[GLY].push_back(7.284562232361592);

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( moldat ) {
    unsigned init_resnum = moldat -> getResidueNumber(atoms[0]);
    for(unsigned i=0; i<atoms.size(); ++i) {
      std::string Rname = moldat->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=ALA;
      } else if(Rname=="ARG") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=ARG;
      } else if(Rname=="ASN") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=ASN;
      } else if(Rname=="ASP") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=ASP;
      } else if(Rname=="CYS") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=CYS;
      } else if(Rname=="GLN") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=GLN;
      } else if(Rname=="GLU") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=GLU;
      } else if(Rname=="GLY") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=GLY;
      } else if(Rname=="HIS") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIS;
      } else if(Rname=="HID") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIS;
      } else if(Rname=="HIE") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIE;
      } else if(Rname=="HIP") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIP;
        // CHARMM NAMING FOR PROTONATION STATES OF HISTIDINE
      } else if(Rname=="HSD") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIS;
      } else if(Rname=="HSE") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIE;
      } else if(Rname=="HSP") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=HIP;
      } else if(Rname=="ILE") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=ILE;
      } else if(Rname=="LEU") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=LEU;
      } else if(Rname=="LYS") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=LYS;
      } else if(Rname=="MET") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=MET;
      } else if(Rname=="PHE") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=PHE;
      } else if(Rname=="PRO") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=PRO;
      } else if(Rname=="SER") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=SER;
      } else if(Rname=="THR") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=THR;
      } else if(Rname=="TRP") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=TRP;
      } else if(Rname=="TYR") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=TYR;
      } else if(Rname=="VAL") {
        atoi[moldat->getResidueNumber(atoms[i])-init_resnum]=VAL;
      } else error("Residue not known: "+Rname);
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

double SAXS::calculateASF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double rho)
{
  std::map<std::string, unsigned> AA_map;
  AA_map["H"] = H;
  AA_map["C"] = C;
  AA_map["N"] = N;
  AA_map["O"] = O;
  AA_map["P"] = P;
  AA_map["S"] = S;

  std::vector<std::vector<double> > param_a;
  std::vector<std::vector<double> > param_b;
  std::vector<double> param_c;
  std::vector<double> param_v;

  param_a.resize(NTT, std::vector<double>(5));
  param_b.resize(NTT, std::vector<double>(5));
  param_c.resize(NTT);
  param_v.resize(NTT);

  param_a[H][0] = 0.493002; param_b[H][0] = 10.5109; param_c[H] = 0.003038;
  param_a[H][1] = 0.322912; param_b[H][1] = 26.1257; param_v[H] = 5.15;
  param_a[H][2] = 0.140191; param_b[H][2] = 3.14236;
  param_a[H][3] = 0.040810; param_b[H][3] = 57.7997;
  param_a[H][4] = 0.0;      param_b[H][4] = 1.0;

  param_a[C][0] = 2.31000; param_b[C][0] = 20.8439; param_c[C] = 0.215600;
  param_a[C][1] = 1.02000; param_b[C][1] = 10.2075; param_v[C] = 16.44;
  param_a[C][2] = 1.58860; param_b[C][2] = 0.56870;
  param_a[C][3] = 0.86500; param_b[C][3] = 51.6512;
  param_a[C][4] = 0.0;     param_b[C][4] = 1.0;

  param_a[N][0] = 12.2126; param_b[N][0] = 0.00570; param_c[N] = -11.529;
  param_a[N][1] = 3.13220; param_b[N][1] = 9.89330; param_v[N] = 2.49;
  param_a[N][2] = 2.01250; param_b[N][2] = 28.9975;
  param_a[N][3] = 1.16630; param_b[N][3] = 0.58260;
  param_a[N][4] = 0.0;     param_b[N][4] = 1.0;

  param_a[O][0] = 3.04850; param_b[O][0] = 13.2771; param_c[O] = 0.250800 ;
  param_a[O][1] = 2.28680; param_b[O][1] = 5.70110; param_v[O] = 9.13;
  param_a[O][2] = 1.54630; param_b[O][2] = 0.32390;
  param_a[O][3] = 0.86700; param_b[O][3] = 32.9089;
  param_a[O][4] = 0.0;     param_b[O][4] = 1.0;

  param_a[P][0] = 6.43450; param_b[P][0] = 1.90670; param_c[P] = 1.11490;
  param_a[P][1] = 4.17910; param_b[P][1] = 27.1570; param_v[P] = 5.73;
  param_a[P][2] = 1.78000; param_b[P][2] = 0.52600;
  param_a[P][3] = 1.49080; param_b[P][3] = 68.1645;
  param_a[P][4] = 0.0;     param_b[P][4] = 1.0;

  param_a[S][0] = 6.90530; param_b[S][0] = 1.46790; param_c[S] = 0.866900;
  param_a[S][1] = 5.20340; param_b[S][1] = 22.2151; param_v[S] = 19.86;
  param_a[S][2] = 1.43790; param_b[S][2] = 0.25360;
  param_a[S][3] = 1.58630; param_b[S][3] = 56.1720;
  param_a[S][4] = 0.0;     param_b[S][4] = 1.0;

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);

  double Iq0=0.;
  if( moldat ) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    // cycle over the atom types
    for(unsigned i=0; i<NTT; ++i) {
      const double volr = std::pow(param_v[i], (2.0/3.0)) /(4. * M_PI);
      // cycle over q
      for(unsigned k=0; k<q_list.size(); ++k) {
        const double q = q_list[k];
        const double s = q / (4. * M_PI);
        FF_tmp[k][i] = param_c[i];
        // SUM [a_i * EXP( - b_i * (q/4pi)^2 )] Waasmaier and Kirfel (1995)
        for(unsigned j=0; j<4; ++j) {
          FF_tmp[k][i] += param_a[i][j]*std::exp(-param_b[i][j]*s*s);
        }
        // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
        FF_tmp[k][i] -= rho*param_v[i]*std::exp(-volr*q*q);
      }
    }
    // cycle over the atoms to assign the atom type and calculate I0
    for(unsigned i=0; i<atoms.size(); ++i) {
      // get atom name
      std::string name = moldat->getAtomName(atoms[i]);
      char type;
      // get atom type
      char first = name.at(0);
      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = name.at(1);
      }
      std::string type_s = std::string(1,type);
      if(AA_map.find(type_s) != AA_map.end()) {
        const unsigned index=AA_map[type_s];
        atoi[i] = AA_map[type_s];
        for(unsigned j=0; j<4; ++j) Iq0 += param_a[index][j];
        Iq0 = Iq0 -rho*param_v[index] + param_c[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }

  return Iq0;
}

std::map<std::string, std::vector<double> > SAXS::setupLCPOparam() {
  std::map<std::string, std::vector<double> > lcpomap;

  //We arbitrarily set OC1 as the charged oxygen.

  lcpomap["ALA_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ALA_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ALA_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ALA_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_CB"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ALA_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ALA_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ALA_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["ASP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ASP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ASP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ASP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASP_OD1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["ASP_OD2"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["ASP_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ASP_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASP_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ASP_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["ASN_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ASN_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ASN_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASN_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASN_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ASN_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ASN_OD1"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASN_ND2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ASN_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ASN_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ASN_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ASN_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["ARG_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ARG_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ARG_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ARG_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ARG_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ARG_NE"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ARG_NH1"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ARG_NH2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["ARG_CZ"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ARG_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ARG_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ARG_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ARG_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["CYS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["CYS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["CYS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["CYS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["CYS_SG"] = { 1.9,  0.54581,  -0.19477,  -0.0012873,  0.00029247};
  lcpomap["CYS_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["CYS_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYS_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["CYS_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["GLU_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLU_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["GLU_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLU_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLU_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLU_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLU_CD"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLU_OE1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["GLU_OE2"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["GLU_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLU_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLU_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLU_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["GLN_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLN_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["GLN_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLN_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLN_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLN_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLN_CD"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLN_OE1"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLN_NE2"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["GLN_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLN_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLN_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLN_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["GLY_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLY_CA"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLY_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLY_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLY_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLY_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLY_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLY_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HIS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HIS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HIS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HIS_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIS_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIS_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIS_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIS_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HIS_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HIE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HIE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HIE_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HIE_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIE_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIE_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIE_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIE_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIE_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HIE_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HSE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HSE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HSE_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HSE_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSE_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSE_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSE_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSE_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSE_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HSE_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HID_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HID_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HID_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HID_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HID_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HID_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HID_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HID_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HID_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HID_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HID_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HID_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HSD_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSD_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HSD_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSD_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HSD_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HSD_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSD_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSD_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSD_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSD_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSD_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HSD_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HIP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HIP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HIP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HIP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HIP_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIP_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIP_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HIP_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HIP_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HIP_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["HSP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["HSP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["HSP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["HSP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["HSP_ND1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSP_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSP_NE2"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["HSP_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["HSP_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["HSP_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["ILE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ILE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ILE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ILE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ILE_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ILE_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ILE_CG1"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["ILE_CD1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ILE_CD"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ILE_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ILE_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ILE_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ILE_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["LEU_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["LEU_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LEU_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["LEU_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LEU_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LEU_CG"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LEU_CD1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["LEU_CD2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["LEU_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["LEU_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LEU_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["LEU_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["LYS_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["LYS_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["LYS_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["LYS_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LYS_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_CE"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["LYS_NZ"] = { 1.65,  0.73511,  -0.22116,  -0.00089148,  0.0002523};
  lcpomap["LYS_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["LYS_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["LYS_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["LYS_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["MET_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["MET_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["MET_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["MET_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["MET_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["MET_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["MET_SD"] = { 1.9,  0.54581,  -0.19477,  -0.0012873,  0.00029247};
  lcpomap["MET_CE"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["MET_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["MET_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["MET_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["MET_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["PHE_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["PHE_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["PHE_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PHE_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PHE_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PHE_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PHE_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CZ"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CE2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["PHE_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["PHE_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PHE_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["PHE_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["PRO_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["PRO_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["PRO_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["PRO_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PRO_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PRO_CG"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PRO_CD"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["PRO_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["PRO_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["PRO_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["PRO_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["SER_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["SER_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["SER_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["SER_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["SER_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["SER_OG"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["SER_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["SER_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["SER_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["SER_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["THR_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["THR_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["THR_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["THR_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["THR_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["THR_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["THR_OG1"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["THR_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["THR_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["THR_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["THR_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["TRP_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["TRP_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["TRP_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TRP_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["TRP_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_NE1"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["TRP_CE2"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_CZ2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CH2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CZ3"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CE3"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TRP_CD2"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TRP_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["TRP_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TRP_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["TRP_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["TYR_N"] = { 1.65,  0.062577,  -0.017874,  -8.312e-05,  1.9849e-05};
  lcpomap["TYR_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["TYR_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TYR_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["TYR_CG"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_CD1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CE1"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CZ"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["TYR_OH"] = { 1.6,  0.77914,  -0.25262,  -0.0016056,  0.00035071};
  lcpomap["TYR_CE2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_CD2"] = { 1.7,  0.51245,  -0.15966,  -0.00019781,  0.00016392};
  lcpomap["TYR_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["TYR_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["TYR_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["TYR_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  lcpomap["VAL_N"] = { 1.65,  0.062577,  -0.017874,  -8.312e-05,  1.9849e-05};
  lcpomap["VAL_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["VAL_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["VAL_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["VAL_CB"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["VAL_CG1"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["VAL_CG2"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["VAL_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["VAL_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["VAL_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["VAL_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};

  return lcpomap;
}

//assigns LCPO parameters to each atom reading from database
void SAXS::readLCPOparam(const std::vector<std::vector<std::string> > &AtomResidueName, unsigned natoms)
{
  std::map<std::string, std::vector<double> > lcpomap = setupLCPOparam();

  for(unsigned i=0; i<natoms; ++i)
  {
    if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
      std::string identifier = AtomResidueName[1][i]+"_"+AtomResidueName[0][i];
      std::vector<double> LCPOparamVector = lcpomap.at(identifier);
      double rs = 0.14;
      LCPOparam[i].push_back(LCPOparamVector[0]+rs*10.);
      LCPOparam[i].push_back(LCPOparamVector[1]);
      LCPOparam[i].push_back(LCPOparamVector[2]);
      LCPOparam[i].push_back(LCPOparamVector[3]);
      LCPOparam[i].push_back(LCPOparamVector[4]);
    }
  }

  for(unsigned i=0; i<natoms; ++i) {
    if (LCPOparam[i].size()==0 ) {
      if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S')) {
        std::cout << "Could not find LCPO paramaters for atom " << AtomResidueName[0][i] << " of residue " << AtomResidueName[1][i] << std::endl;
        error ("missing LCPO parameters\n");
      }
    }
  }

  if (AtomResidueName[0][0] == "N") {
    LCPOparam[0][1] = 7.3511e-01;
    LCPOparam[0][2] = -2.2116e-01;
    LCPOparam[0][3] = -8.9148e-04;
    LCPOparam[0][4] = 2.5230e-04;
  }

  if (AtomResidueName[0][natoms-1] == "O") {
    LCPOparam[natoms-1][1] = 8.8857e-01;
    LCPOparam[natoms-1][2] = -3.3421e-01;
    LCPOparam[natoms-1][3] = -1.8683e-03;
    LCPOparam[natoms-1][4] = 4.9372e-04;
  }
}


}//namespace isdb
}//namespace PLMD



