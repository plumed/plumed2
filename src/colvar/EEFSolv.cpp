/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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

/* This class was originally written by Thomas Loehr */

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/GenericMolInfo.h"
#include "tools/OpenMP.h"
#include <initializer_list>

#define INV_PI_SQRT_PI 0.179587122
#define KCAL_TO_KJ 4.184
#define ANG_TO_NM 0.1
#define ANG3_TO_NM3 0.001

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR EEFSOLV
/*
Calculates EEF1 solvation free energy for a group of atoms.

EEF1 is a solvent-accessible surface area based model, where the free energy of solvation is computed using a pairwise interaction term for non-hydrogen atoms:
\f[
    \Delta G^\mathrm{solv}_i = \Delta G^\mathrm{ref}_i - \sum_{j \neq i} f_i(r_{ij}) V_j
\f]
where \f$\Delta G^\mathrm{solv}_i\f$ is the free energy of solvation, \f$\Delta G^\mathrm{ref}_i\f$ is the reference solvation free energy, \f$V_j\f$ is the volume of atom \f$j\f$ and
\f[
    f_i(r) 4\pi r^2 = \frac{2}{\sqrt{\pi}} \frac{\Delta G^\mathrm{free}_i}{\lambda_i} \exp\left\{ - \frac{(r-R_i)^2}{\lambda^2_i}\right\}
\f]
where \f$\Delta G^\mathrm{free}_i\f$ is the solvation free energy of the isolated group, \f$\lambda_i\f$ is the correlation length equal to the width of the first solvation shell and \f$R_i\f$ is the van der Waals radius of atom \f$i\f$.

The output from this collective variable, the free energy of solvation, can be used with the \ref BIASVALUE keyword to provide implicit solvation to a system. All parameters are designed to be used with a modified CHARMM36 force field. It takes only non-hydrogen atoms as input, these can be conveniently specified using the \ref GROUP action with the NDX_GROUP parameter. To speed up the calculation, EEFSOLV internally uses a neighbor list with a cutoff dependent on the type of atom (maximum of 1.95 nm). This cutoff can be extended further by using the NL_BUFFER keyword.

\par Examples

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt77/peptide.pdb
MOLINFO MOLTYPE=protein STRUCTURE=peptide.pdb
WHOLEMOLECULES ENTITY0=1-111

# This allows us to select only non-hydrogen atoms
#SETTINGS AUXFILE=regtest/basic/rt77/index.ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# We extend the cutoff by 0.1 nm and update the neighbor list every 40 steps
solv: EEFSOLV ATOMS=protein-h

# Here we actually add our calculated energy back to the potential
bias: BIASVALUE ARG=solv

PRINT ARG=solv FILE=SOLV
\endplumedfile

*/
//+ENDPLUMEDOC

class EEFSolv : public Colvar {
private:
  bool pbc;
  bool serial;
  double delta_g_ref;
  double nl_buffer;
  unsigned nl_stride;
  unsigned nl_update;
  std::vector<std::vector<unsigned> > nl;
  std::vector<std::vector<bool> > nlexpo;
  std::vector<std::vector<double> > parameter;
  void setupConstants(const std::vector<AtomNumber> &atoms, std::vector<std::vector<double> > &parameter, bool tcorr);
  std::map<std::string, std::map<std::string, std::string> > setupTypeMap();
  std::map<std::string, std::vector<double> > setupValueMap();
  void update_neighb();

public:
  static void registerKeywords(Keywords& keys);
  explicit EEFSolv(const ActionOptions&);
  void calculate() override;
};

PLUMED_REGISTER_ACTION(EEFSolv,"EEFSOLV")

void EEFSolv::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "The atoms to be included in the calculation, e.g. the whole protein.");
  keys.add("compulsory", "NL_BUFFER", "0.1", "The buffer to the intrinsic cutoff used when calculating pairwise interactions.");
  keys.add("compulsory", "NL_STRIDE", "40", "The frequency with which the neighbor list is updated.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("TEMP_CORRECTION", false, "Correct free energy of solvation constants for temperatures different from 298.15 K");
}

EEFSolv::EEFSolv(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  delta_g_ref(0.),
  nl_buffer(0.1),
  nl_stride(40),
  nl_update(0)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  const unsigned size = atoms.size();
  bool tcorr = false;
  parseFlag("TEMP_CORRECTION", tcorr);
  parse("NL_BUFFER", nl_buffer);
  parse("NL_STRIDE", nl_stride);

  bool nopbc = !pbc;
  parseFlag("NOPBC", nopbc);
  pbc = !nopbc;

  parseFlag("SERIAL",serial);

  checkRead();

  log << "  Bibliography " << plumed.cite("Lazaridis T, Karplus M, Proteins Struct. Funct. Genet. 35, 133 (1999)"); log << "\n";

  nl.resize(size);
  nlexpo.resize(size);
  parameter.resize(size, std::vector<double>(4, 0));
  setupConstants(atoms, parameter, tcorr);

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
}

void EEFSolv::update_neighb() {
  const double lower_c2 = 0.24 * 0.24; // this is the cut-off for bonded atoms
  const unsigned size = getNumberOfAtoms();

  for (unsigned i=0; i<size; i++) {
    nl[i].clear();
    nlexpo[i].clear();
    const Vector posi = getPosition(i);
    // Loop through neighboring atoms, add the ones below cutoff
    for (unsigned j=i+1; j<size; j++) {
      if(parameter[i][1]==0&&parameter[j][1]==0) continue;
      const double d2 = delta(posi, getPosition(j)).modulo2();
      if (d2 < lower_c2 && j < i+14) {
        // crude approximation for i-i+1/2 interactions,
        // we want to exclude atoms separated by less than three bonds
        continue;
      }
      // We choose the maximum lambda value and use a more conservative cutoff
      double mlambda = 1./parameter[i][2];
      if (1./parameter[j][2] > mlambda) mlambda = 1./parameter[j][2];
      const double c2 = (2. * mlambda + nl_buffer) * (2. * mlambda + nl_buffer);
      if (d2 < c2 ) {
        nl[i].push_back(j);
        if(parameter[i][2] == parameter[j][2] && parameter[i][3] == parameter[j][3]) {
          nlexpo[i].push_back(true);
        } else nlexpo[i].push_back(false);
      }
    }
  }
}

void EEFSolv::calculate() {
  if(pbc) makeWhole();
  if(getExchangeStep()) nl_update = 0;
  if(nl_update==0) update_neighb();

  const unsigned size=getNumberOfAtoms();
  double bias = 0.0;
  std::vector<Vector> deriv(size, Vector(0,0,0));

  unsigned stride;
  unsigned rank;
  if(serial) {
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();
  if(nt*stride*10>size) nt=1;

  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> deriv_omp(size, Vector(0,0,0));
    #pragma omp for reduction(+:bias) nowait
    for (unsigned i=rank; i<size; i+=stride) {
      const Vector posi = getPosition(i);
      double fedensity = 0.0;
      Vector deriv_i;
      const double vdw_volume_i   = parameter[i][0];
      const double delta_g_free_i = parameter[i][1];
      const double inv_lambda_i   = parameter[i][2];
      const double vdw_radius_i   = parameter[i][3];

      // The pairwise interactions are unsymmetric, but we can get away with calculating the distance only once
      for (unsigned i_nl=0; i_nl<nl[i].size(); i_nl++) {
        const unsigned j = nl[i][i_nl];
        const double vdw_volume_j   = parameter[j][0];
        const double delta_g_free_j = parameter[j][1];
        const double inv_lambda_j   = parameter[j][2];
        const double vdw_radius_j   = parameter[j][3];

        const Vector dist     = delta(posi, getPosition(j));
        const double rij      = dist.modulo();
        const double inv_rij  = 1.0 / rij;
        const double inv_rij2 = inv_rij * inv_rij;
        const double fact_ij  = inv_rij2 * delta_g_free_i * vdw_volume_j * INV_PI_SQRT_PI * inv_lambda_i;
        const double fact_ji  = inv_rij2 * delta_g_free_j * vdw_volume_i * INV_PI_SQRT_PI * inv_lambda_j;

        // in this case we can calculate a single exponential
        if(!nlexpo[i][i_nl]) {
          // i-j interaction
          if(inv_rij > 0.5*inv_lambda_i && delta_g_free_i!=0.)
          {
            const double e_arg = (rij - vdw_radius_i)*inv_lambda_i;
            const double expo  = exp(-e_arg*e_arg);
            const double fact  = expo*fact_ij;
            const double e_deriv = inv_rij*fact*(inv_rij + e_arg*inv_lambda_i);
            const Vector dd    = e_deriv*dist;
            fedensity    += fact;
            deriv_i      += dd;
            if(nt>1) deriv_omp[j] -= dd;
            else deriv[j] -= dd;
          }

          // j-i interaction
          if(inv_rij > 0.5*inv_lambda_j && delta_g_free_j!=0.)
          {
            const double e_arg = (rij - vdw_radius_j)*inv_lambda_j;
            const double expo  = exp(-e_arg*e_arg);
            const double fact  = expo*fact_ji;
            const double e_deriv = inv_rij*fact*(inv_rij + e_arg*inv_lambda_j);
            const Vector dd    = e_deriv*dist;
            fedensity    += fact;
            deriv_i      += dd;
            if(nt>1) deriv_omp[j] -= dd;
            else deriv[j] -= dd;
          }
        } else {
          // i-j interaction
          if(inv_rij > 0.5*inv_lambda_i)
          {
            const double e_arg = (rij - vdw_radius_i)*inv_lambda_i;
            const double expo  = exp(-e_arg*e_arg);
            const double fact  = expo*(fact_ij + fact_ji);
            const double e_deriv = inv_rij*fact*(inv_rij + e_arg*inv_lambda_i);
            const Vector dd    = e_deriv*dist;
            fedensity    += fact;
            deriv_i      += dd;
            if(nt>1) deriv_omp[j] -= dd;
            else deriv[j] -= dd;
          }
        }

      }
      if(nt>1) deriv_omp[i] += deriv_i;
      else deriv[i] += deriv_i;
      bias += 0.5*fedensity;
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<size; i++) deriv[i]+=deriv_omp[i];
  }

  if(!serial) {
    comm.Sum(bias);
    if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
  }

  Tensor virial;
  for(unsigned i=0; i<size; i++) {
    setAtomsDerivatives(i, -deriv[i]);
    virial += Tensor(getPosition(i), -deriv[i]);
  }
  setBoxDerivatives(-virial);
  setValue(delta_g_ref - bias);

  // Keep track of the neighbourlist updates
  nl_update++;
  if (nl_update == nl_stride) {
    nl_update = 0;
  }
}

void EEFSolv::setupConstants(const std::vector<AtomNumber> &atoms, std::vector<std::vector<double> > &parameter, bool tcorr) {
  std::vector<std::vector<double> > parameter_temp;
  parameter_temp.resize(atoms.size(), std::vector<double>(7,0));
  std::map<std::string, std::vector<double> > valuemap;
  std::map<std::string, std::map<std::string, std::string> > typemap;
  valuemap = setupValueMap();
  typemap  = setupTypeMap();
  auto * moldat = plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  bool cter=false;
  if (moldat) {
    log<<"  MOLINFO DATA found with label " <<moldat->getLabel()<<", using proper atom names\n";
    for(unsigned i=0; i<atoms.size(); ++i) {

      // Get atom and residue names
      std::string Aname = moldat->getAtomName(atoms[i]);
      std::string Rname = moldat->getResidueName(atoms[i]);
      std::string Atype = typemap[Rname][Aname];

      // Check for terminal COOH or COO- (different atomtypes & parameters!)
      if (Aname == "OT1" || Aname == "OXT") {
        // We create a temporary AtomNumber object to access future atoms
        unsigned ai = atoms[i].index();
        AtomNumber tmp_an;
        tmp_an.setIndex(ai + 2);
        if (moldat->checkForAtom(tmp_an) && moldat->getAtomName(tmp_an) == "HT2") {
          // COOH
          Atype = "OB";
        } else {
          // COO-
          Atype = "OC";
        }
        cter = true;
      }
      if (Aname == "OT2" || (cter == true && Aname == "O")) {
        unsigned ai = atoms[i].index();
        AtomNumber tmp_an;
        tmp_an.setIndex(ai + 1);
        if (moldat->checkForAtom(tmp_an) && moldat->getAtomName(tmp_an) == "HT2") {
          // COOH
          Atype = "OH1";
        } else {
          // COO-
          Atype = "OC";
        }
      }

      // Check for H-atoms
      char type;
      char first = Aname.at(0);

      // GOLDEN RULE: type is first letter, if not a number
      if (!isdigit(first)) {
        type = first;
        // otherwise is the second
      } else {
        type = Aname.at(1);
      }

      if (type == 'H') {
        error("EEF1-SB does not allow the use of hydrogen atoms!\n");
      }

      // Lookup atomtype in table or throw exception if its not there
      try {
        parameter_temp[i] = valuemap.at(Atype);
      } catch (std::exception &e) {
        log << "Type: " << Atype << "  Name: " << Aname << "  Residue: " << Rname << "\n";
        error("Invalid atom type!\n");
      }

      // Temperature correction
      if (tcorr && parameter[i][1] > 0.0) {
        const double t0 = 298.15;
        const double delta_g_ref_t0 = parameter_temp[i][1];
        const double delta_h_ref_t0 = parameter_temp[i][3];
        const double delta_cp = parameter_temp[i][4];
        const double delta_s_ref_t0 = (delta_h_ref_t0 - delta_g_ref_t0) / t0;
        const double t = plumed.getAtoms().getKbT() / plumed.getAtoms().getKBoltzmann();
        parameter_temp[i][1] -= delta_s_ref_t0 * (t - t0) - delta_cp * t * std::log(t / t0) + delta_cp * (t - t0);
        parameter_temp[i][2] *= parameter_temp[i][1] / delta_g_ref_t0;
      }
      parameter[i][0] = parameter_temp[i][0];
      parameter[i][1] = parameter_temp[i][2];
      parameter[i][2] = parameter_temp[i][5];
      parameter[i][3] = parameter_temp[i][6];
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  for(unsigned i=0; i<atoms.size(); ++i) delta_g_ref += parameter_temp[i][1];
}

std::map<std::string, std::map<std::string, std::string> > EEFSolv::setupTypeMap()  {
  std::map<std::string, std::map<std::string, std::string> > typemap;
  typemap = {
    { "ACE", {
        {"CH3", "CT3"},
        {"HH31","HA3"},
        {"HH32","HA3"},
        {"HH33","HA3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "ALA", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT3"},
        {"HB1", "HA3"},
        {"HB2", "HA3"},
        {"HB3", "HA3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "ARG", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT2"},
        {"HB1",  "HA2"},
        {"HB2",  "HA2"},
        {"CG",   "CT2"},
        {"HG1",  "HA2"},
        {"HG2",  "HA2"},
        {"CD",   "CT2"},
        {"HD1",  "HA2"},
        {"HD2",  "HA2"},
        {"NE",   "NC2"},
        {"HE",   "HC" },
        {"CZ",   "C"  },
        {"NH1",  "NC2"},
        {"HH11", "HC" },
        {"HH12", "HC" },
        {"NH2",  "NC2"},
        {"HH21", "HC" },
        {"HH22", "HC" },
        {"C",    "C"  },
        {"O",    "O"  }
      }
    },
    { "ASN", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT2"},
        {"HB1",  "HA2"},
        {"HB2",  "HA2"},
        {"CG",   "CC" },
        {"OD1",  "O"  },
        {"ND2",  "NH2"},
        {"HD21", "H"  },
        {"HD22", "H"  },
        {"C",    "C"  },
        {"O",    "O"  }
      }
    },
    { "ASPP", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CD" },
        {"OD1", "OB" },
        {"OD2", "OH1"},
        {"HD2", "H"  },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "ASP", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CC" },
        {"OD1", "OC" },
        {"OD2", "OC" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "CYS", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"SG",  "S"  },
        {"HG1", "HS" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "GLN", {
        {"N",    "NH1" },
        {"HN",   "H"   },
        {"CA",   "CT1" },
        {"HA",   "HB1" },
        {"CB",   "CT2" },
        {"HB1",  "HA2" },
        {"HB2",  "HA2" },
        {"CG",   "CT2" },
        {"HG1",  "HA2" },
        {"HG2",  "HA2" },
        {"CD",   "CC"  },
        {"OE1",  "O"   },
        {"NE2",  "NH2" },
        {"HE21", "H"   },
        {"HE22", "H"   },
        {"C",    "C"   },
        {"O",    "O"   }
      }
    },
    { "GLUP", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CT2"},
        {"HG1", "HA2"},
        {"HG2", "HA2"},
        {"CD",  "CD" },
        {"OE1", "OB" },
        {"OE2", "OH1"},
        {"HE2", "H"  },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "GLU", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CT2"},
        {"HG1", "HA2"},
        {"HG2", "HA2"},
        {"CD",  "CC" },
        {"OE1", "OC" },
        {"OE2", "OC" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "GLY", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT2"},
        {"HA1", "HB2"},
        {"HA2", "HB2"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "HSD", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"ND1", "NR1"},
        {"HD1", "H"  },
        {"CG",  "CPH1"},
        {"CE1", "CPH2"},
        {"HE1", "HR1"},
        {"NE2", "NR2"},
        {"CD2", "CPH1"},
        {"HD2", "HR3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "HIS", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"ND1", "NR2"},
        {"CG",  "CPH1"},
        {"CE1", "CPH2"},
        {"HE1", "HR1"},
        {"NE2", "NR1"},
        {"HE2", "H"  },
        {"CD2", "CPH1"},
        {"HD2", "HR3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "HSE", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"ND1", "NR2"},
        {"CG",  "CPH1"},
        {"CE1", "CPH2"},
        {"HE1", "HR1"},
        {"NE2", "NR1"},
        {"HE2", "H"  },
        {"CD2", "CPH1"},
        {"HD2", "HR3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "HSP", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CD2", "CPH1"},
        {"HD2", "HR1"},
        {"CG",  "CPH1"},
        {"NE2", "NR3"},
        {"HE2", "H"  },
        {"ND1", "NR3"},
        {"HD1", "H"  },
        {"CE1", "CPH2"},
        {"HE1", "HR2"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "ILE", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT1"},
        {"HB",   "HA1"},
        {"CG2",  "CT3"},
        {"HG21", "HA3"},
        {"HG22", "HA3"},
        {"HG23", "HA3"},
        {"CG1",  "CT2"},
        {"HG11", "HA2"},
        {"HG12", "HA2"},
        {"CD",   "CT3"},
        {"HD1",  "HA3"},
        {"HD2",  "HA3"},
        {"HD3",  "HA3"},
        {"C",    "C"  },
        {"O",    "O"  }
      }
    },
    { "LEU", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT2"},
        {"HB1",  "HA2"},
        {"HB2",  "HA2"},
        {"CG",   "CT1"},
        {"HG",   "HA1"},
        {"CD1",  "CT3"},
        {"HD11", "HA3"},
        {"HD12", "HA3"},
        {"HD13", "HA3"},
        {"CD2",  "CT3"},
        {"HD21", "HA3"},
        {"HD22", "HA3"},
        {"HD23", "HA3"},
        {"C",    "C"  },
        {"O",    "O"  }
      }
    },
    { "LYS", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CT2"},
        {"HG1", "HA2"},
        {"HG2", "HA2"},
        {"CD",  "CT2"},
        {"HD1", "HA2"},
        {"HD2", "HA2"},
        {"CE",  "CT2"},
        {"HE1", "HA2"},
        {"HE2", "HA2"},
        {"NZ",  "NH3"},
        {"HZ1", "HC" },
        {"HZ2", "HC" },
        {"HZ3", "HC" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "MET", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CT2"},
        {"HG1", "HA2"},
        {"HG2", "HA2"},
        {"SD",  "S"  },
        {"CE",  "CT3"},
        {"HE1", "HA3"},
        {"HE2", "HA3"},
        {"HE3", "HA3"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "NMA", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CH3", "CT3"},
        {"HH31","HA3"},
        {"HH32","HA3"},
        {"HH33","HA3"},
      }
    },
    { "PHE", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CA" },
        {"CD1", "CA" },
        {"HD1", "HP" },
        {"CE1", "CA" },
        {"HE1", "HP" },
        {"CZ",  "CA" },
        {"HZ",  "HP" },
        {"CD2", "CA" },
        {"HD2", "HP" },
        {"CE2", "CA" },
        {"HE2", "HP" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "PRO", {
        {"N",   "N"  },
        {"CD",  "CP3"},
        {"HD1", "HA2"},
        {"HD2", "HA2"},
        {"CA",  "CP1"},
        {"HA",  "HB1"},
        {"CB",  "CP2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CP2"},
        {"HG1", "HA2"},
        {"HG2", "HA2"},
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "SER", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"OG",  "OH1"},
        {"HG1", "H"  },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "THR", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT1"},
        {"HB",   "HA1"},
        {"OG1",  "OH1"},
        {"HG1",  "H"  },
        {"CG2",  "CT3"},
        {"HG21", "HA3"},
        {"HG22", "HA3"},
        {"HG23", "HA3"},
        {"C",    "C"  },
        {"O",    "O"  }
      }
    },
    { "TRP", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CY" },
        {"CD1", "CA" },
        {"HD1", "HP" },
        {"NE1", "NY" },
        {"HE1", "H"  },
        {"CE2", "CPT"},
        {"CD2", "CPT"},
        {"CE3", "CAI"},
        {"HE3", "HP" },
        {"CZ3", "CA" },
        {"HZ3", "HP" },
        {"CZ2", "CAI"},
        {"HZ2", "HP" },
        {"CH2", "CA" },
        {"HH2", "HP" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "TYR", {
        {"N",   "NH1"},
        {"HN",  "H"  },
        {"CA",  "CT1"},
        {"HA",  "HB1"},
        {"CB",  "CT2"},
        {"HB1", "HA2"},
        {"HB2", "HA2"},
        {"CG",  "CA" },
        {"CD1", "CA" },
        {"HD1", "HP" },
        {"CE1", "CA" },
        {"HE1", "HP" },
        {"CZ",  "CA" },
        {"OH",  "OH1"},
        {"HH",  "H"  },
        {"CD2", "CA" },
        {"HD2", "HP" },
        {"CE2", "CA" },
        {"HE2", "HP" },
        {"C",   "C"  },
        {"O",   "O"  }
      }
    },
    { "VAL", {
        {"N",    "NH1"},
        {"HN",   "H"  },
        {"CA",   "CT1"},
        {"HA",   "HB1"},
        {"CB",   "CT1"},
        {"HB",   "HA1"},
        {"CG1",  "CT3"},
        {"HG11", "HA3"},
        {"HG12", "HA3"},
        {"HG13", "HA3"},
        {"CG2",  "CT3"},
        {"HG21", "HA3"},
        {"HG22", "HA3"},
        {"HG23", "HA3"},
        {"C",    "C"  },
        {"O",    "O"  }
      }
    }
  };
  return typemap;
}

std::map<std::string, std::vector<double> > EEFSolv::setupValueMap() {
  // Volume ∆Gref ∆Gfree ∆H ∆Cp λ vdw_radius
  std::map<std::string, std::vector<double> > valuemap;
  valuemap = {
    { "C", {
        ANG3_TO_NM3 * 14.720,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.0,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CD", {
        ANG3_TO_NM3 * 14.720,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.0,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CT1", {
        ANG3_TO_NM3 * 11.507,
        KCAL_TO_KJ * -0.187,
        KCAL_TO_KJ * -0.187,
        KCAL_TO_KJ * 0.876,
        KCAL_TO_KJ * 0.0,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CT2", {
        ANG3_TO_NM3 * 18.850,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * -0.610,
        KCAL_TO_KJ * 18.6,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CT2A", {
        ANG3_TO_NM3 * 18.666,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * -0.610,
        KCAL_TO_KJ * 18.6,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CT3", {
        ANG3_TO_NM3 * 27.941,
        KCAL_TO_KJ * 1.089,
        KCAL_TO_KJ * 1.089,
        KCAL_TO_KJ * -1.779,
        KCAL_TO_KJ * 35.6,
        1. / (ANG_TO_NM * 3.5),
        0.204,
      }
    },
    { "CPH1", {
        ANG3_TO_NM3 * 5.275,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * 0.080,
        KCAL_TO_KJ * -0.973,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.18,
      }
    },
    { "CPH2", {
        ANG3_TO_NM3 * 11.796,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * 0.080,
        KCAL_TO_KJ * -0.973,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.18,
      }
    },
    { "CPT", {
        ANG3_TO_NM3 * 4.669,
        KCAL_TO_KJ * -0.890,
        KCAL_TO_KJ * -0.890,
        KCAL_TO_KJ * 2.220,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.186,
      }
    },
    { "CY", {
        ANG3_TO_NM3 * 10.507,
        KCAL_TO_KJ * -0.890,
        KCAL_TO_KJ * -0.890,
        KCAL_TO_KJ * 2.220,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.199,
      }
    },
    { "CP1", {
        ANG3_TO_NM3 * 25.458,
        KCAL_TO_KJ * -0.187,
        KCAL_TO_KJ * -0.187,
        KCAL_TO_KJ * 0.876,
        KCAL_TO_KJ * 0.0,
        1. / (ANG_TO_NM * 3.5),
        0.227,
      }
    },
    { "CP2", {
        ANG3_TO_NM3 * 19.880,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * -0.610,
        KCAL_TO_KJ * 18.6,
        1. / (ANG_TO_NM * 3.5),
        0.217,
      }
    },
    { "CP3", {
        ANG3_TO_NM3 * 26.731,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * 0.372,
        KCAL_TO_KJ * -0.610,
        KCAL_TO_KJ * 18.6,
        1. / (ANG_TO_NM * 3.5),
        0.217,
      }
    },
    { "CC", {
        ANG3_TO_NM3 * 16.539,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.000,
        KCAL_TO_KJ * 0.0,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "CAI", {
        ANG3_TO_NM3 * 18.249,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * -0.973,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.199,
      }
    },
    { "CA", {
        ANG3_TO_NM3 * 18.249,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * 0.057,
        KCAL_TO_KJ * -0.973,
        KCAL_TO_KJ * 6.9,
        1. / (ANG_TO_NM * 3.5),
        0.199,
      }
    },
    { "N", {
        ANG3_TO_NM3 * 0.000,
        KCAL_TO_KJ * -1.000,
        KCAL_TO_KJ * -1.000,
        KCAL_TO_KJ * -1.250,
        KCAL_TO_KJ * 8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NR1", {
        ANG3_TO_NM3 * 15.273,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -9.059,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NR2", {
        ANG3_TO_NM3 * 15.111,
        KCAL_TO_KJ * -3.820,
        KCAL_TO_KJ * -3.820,
        KCAL_TO_KJ * -4.654,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NR3", {
        ANG3_TO_NM3 * 15.071,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -9.059,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NH1", {
        ANG3_TO_NM3 * 10.197,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -9.059,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NH2", {
        ANG3_TO_NM3 * 18.182,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -9.059,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NH3", {
        ANG3_TO_NM3 * 18.817,
        KCAL_TO_KJ * -20.000,
        KCAL_TO_KJ * -20.000,
        KCAL_TO_KJ * -25.000,
        KCAL_TO_KJ * -18.0,
        1. / (ANG_TO_NM * 6.0),
        0.185,
      }
    },
    { "NC2", {
        ANG3_TO_NM3 * 18.215,
        KCAL_TO_KJ * -10.000,
        KCAL_TO_KJ * -10.000,
        KCAL_TO_KJ * -12.000,
        KCAL_TO_KJ * -7.0,
        1. / (ANG_TO_NM * 6.0),
        0.185,
      }
    },
    { "NY", {
        ANG3_TO_NM3 * 12.001,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -5.950,
        KCAL_TO_KJ * -9.059,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.185,
      }
    },
    { "NP", {
        ANG3_TO_NM3 * 4.993,
        KCAL_TO_KJ * -20.000,
        KCAL_TO_KJ * -20.000,
        KCAL_TO_KJ * -25.000,
        KCAL_TO_KJ * -18.0,
        1. / (ANG_TO_NM * 6.0),
        0.185,
      }
    },
    { "O", {
        ANG3_TO_NM3 * 11.772,
        KCAL_TO_KJ * -5.330,
        KCAL_TO_KJ * -5.330,
        KCAL_TO_KJ * -5.787,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.170,
      }
    },
    { "OB", {
        ANG3_TO_NM3 * 11.694,
        KCAL_TO_KJ * -5.330,
        KCAL_TO_KJ * -5.330,
        KCAL_TO_KJ * -5.787,
        KCAL_TO_KJ * -8.8,
        1. / (ANG_TO_NM * 3.5),
        0.170,
      }
    },
    { "OC", {
        ANG3_TO_NM3 * 12.003,
        KCAL_TO_KJ * -10.000,
        KCAL_TO_KJ * -10.000,
        KCAL_TO_KJ * -12.000,
        KCAL_TO_KJ * -9.4,
        1. / (ANG_TO_NM * 6.0),
        0.170,
      }
    },
    { "OH1", {
        ANG3_TO_NM3 * 15.528,
        KCAL_TO_KJ * -5.920,
        KCAL_TO_KJ * -5.920,
        KCAL_TO_KJ * -9.264,
        KCAL_TO_KJ * -11.2,
        1. / (ANG_TO_NM * 3.5),
        0.177,
      }
    },
    { "OS", {
        ANG3_TO_NM3 * 6.774,
        KCAL_TO_KJ * -2.900,
        KCAL_TO_KJ * -2.900,
        KCAL_TO_KJ * -3.150,
        KCAL_TO_KJ * -4.8,
        1. / (ANG_TO_NM * 3.5),
        0.177,
      }
    },
    { "S", {
        ANG3_TO_NM3 * 20.703,
        KCAL_TO_KJ * -3.240,
        KCAL_TO_KJ * -3.240,
        KCAL_TO_KJ * -4.475,
        KCAL_TO_KJ * -39.9,
        1. / (ANG_TO_NM * 3.5),
        0.20,
      }
    },
    { "SM", {
        ANG3_TO_NM3 * 21.306,
        KCAL_TO_KJ * -3.240,
        KCAL_TO_KJ * -3.240,
        KCAL_TO_KJ * -4.475,
        KCAL_TO_KJ * -39.9,
        1. / (ANG_TO_NM * 3.5),
        0.197,
      }
    }
  };
  return valuemap;
}
}
}
