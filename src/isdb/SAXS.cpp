/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
 Refactoring for hySAXS Martini form factors for Nucleic Acids by Cristina Paissoni
 Refactoring for hySAS OneBead form factors with solvent correction by Federico Ballabio and Riccardo Capelli
 Resolution function by Henrique Musseli Cezar
*/

#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/GenericMolInfo.h"
#include "tools/MolDataClass.h"
#include "tools/Communicator.h"
#include "tools/Pbc.h"
#include "tools/PDB.h"
#include "tools/Tools.h"
#include "tools/IFile.h"

#include <map>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <cctype>

#ifdef __PLUMED_HAS_ARRAYFIRE
#include <arrayfire.h>
#include <af/util.h>
#ifdef __PLUMED_HAS_ARRAYFIRE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <af/cuda.h>
#elif __PLUMED_HAS_ARRAYFIRE_OCL
#include <af/opencl.h>
#endif
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR SAXS
/*
Calculates SAXS intensity.

SAXS intensities are calculated for a set of scattering vectors using QVALUE keywords numbered from 1.
Form factors can be assigned either by polynomial expansion of any order by using the PARAMETERS keywords, or
automatically matched to atoms using the ATOMISTIC flag by reading a PDB file. Alternatively to the atomistic
representation, two types of coarse-grained mapping are available:
- MARTINI.
- ONEBEAD.

Whether for PARAMETERS, ATOMISTIC, and ONEBEAD the user must provide an all-atom PDB file via MOLINFO before the
SAXS instruction. MARTINI requires a mapping scheme consisting of a PDB file that contains both the all-atom
and MARTINI representations, and a bead position file (e.g., bead1: CENTER ATOMS=1,5,7,11,12 WEIGHTS=14,12,12,
12,16).

ONEBEAD scheme consists in a single-bead per amino acid residue or three-bead for nucleic acid residue (one for
the phosphate group, one for the pentose sugar, one for the nucleobase). PLUMED creates a virtual bead on which
the SAXS calculations are performed, centred on the COM of all atoms belonging to the bead. It is possible to
account for the contribution of the solvation layer to the SAXS intensity by adding a correction term for the
solvent accessible beads only: the form factors of the amino acids / phosphate groups / pentose sugars /
nucleobases with a SASA (computed via LCPO algorithm) greater than a threshold are corrected according to an
electron density term. Both the surface cut-off threshold and the electron density term can be set by the user
with the SASA_CUTOFF and SOLVATION_CORRECTION keywords. Moreover, SASA stride calculation can be modified using
SOLVATION_STRIDE, which is set to 10 steps by default.
ONEBEAD requires an additional PDB file to perform mapping conversion, which must be provided via TEMPLATE
keyword. This PDB file should only include the atoms for which the SAXS intensity will be computed.
The AMBER OL3 (RNA) and OL15 (DNA) naming is required for nucleic acids.
Three additional bead types are available for DNA and RNA besides phosphate group, pentose sugar, and nucleobase:
- 5'-end pentose sugar capped with an hydroxyl moiety at C5' (the residue name in the PDB must be followed by
"5", e.g., DC5 or C5 for cytosine in DNA and RNA, respectively);
- 5'-phosphorylated end with an additional hydroxyl moiety at P (the residue name in the PDB must be followed by
"T", e.g., DCT or CT for cytosine in DNA and RNA, respectively);
- 3'-end pentose sugar capped with an hydroxyl moiety at C3' (the residue name in the PDB must be followed by
"3", e.g., DC3 or C3 for cytosine in DNA and RNA, respectively).
An additional bead type is available for proteins:
- Cysteine residue involved in disulfide bridge (the residue in the PDB must be named CYX).

Experimental reference intensities can be added using the EXPINT keywords. All these values must be normalised
to the SAXS intensity at q = 0. To facilitate this operation, the SCALE_EXPINT keyword can be used to provide
the intensity at q = 0. Each EXPINT is divided by SCALE_EXPINT.
The maximum QVALUE for ONEBEAD is set to 0.3 inverse angstroms.
The solvent density, that by default is set to 0.334 electrons per cubic angstrom (bulk water), can be modified
using the SOLVDENS keyword.

The ABSOLUTE flag can be used in order to calculate intensities in the absolute scale. It is only available for
the ATOMISTIC scheme and cannot be used with SCALE_EXPINT.

By default SAXS is calculated using Debye on CPU, by adding the GPU flag it is possible to solve the equation on
a GPU if the ARRAYFIRE libraries are installed and correctly linked.
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

\par Examples
in the following example the SAXS intensities are calculated using single-bead per residue approximation, with a
SASA threshold of 1 square nanometer and a solvation term of 0.04. Each experimental intensity is divided by
1.4002, which is the corresponding theoretical intensity value at q = 0. The form factors are selected according
to the PDB file specified by TEMPLATE keyword.

\plumedfile
MOLINFO STRUCTURE=template_AA.pdb

SAXS ...
LABEL=SAXS
ATOMS=1-355
ONEBEAD
TEMPLATE=template_AA.pdb
SOLVDENS=0.334
SOLVATION_CORRECTION=0.04
SOLVATION_STRIDE=1
SASA_CUTOFF=1.0
SCALE_EXPINT=1.4002
QVALUE1=0.03 EXPINT1=1.0902
QVALUE2=0.06 EXPINT2=0.790632
QVALUE3=0.09 EXPINT3=0.453808
QVALUE4=0.12 EXPINT4=0.254737
QVALUE5=0.15 EXPINT5=0.154928
QVALUE6=0.18 EXPINT6=0.0921503
QVALUE7=0.21 EXPINT7=0.052633
QVALUE8=0.24 EXPINT8=0.0276557
QVALUE9=0.27 EXPINT9=0.0122775
QVALUE10=0.30 EXPINT10=0.00880634
... SAXS

PRINT ARG=(SAXS\.q-.*),(SAXS\.exp-.*) FILE=saxsdata STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC ISDB_COLVAR SANS
/*
Calculates SANS intensity.

SANS intensities are calculated for a set of scattering vectors using QVALUE keywords numbered from 1.
Form factors are automatically assigned to atoms using the ATOMISTIC flag by reading a PDB file, by reading
the scattering lengths with the PARAMETERS keyword from input or with the PARAMETERSFILE keyword or, alternatively,
a ONEBEAD coarse-grained implementation is available.

Both for ATOMISTIC and ONEBEAD the user must provide an all-atom PDB file via MOLINFO before the SANS instruction.

ONEBEAD scheme consists in a single-bead per amino acid residue or three-bead for nucleic acid residue (one for
the phosphate group, one for the pentose sugar, one for the nucleobase). PLUMED creates a virtual bead on which
the SANS calculations are performed, centred on the COM of all atoms belonging to the bead. It is possible to
account for the contribution of the solvation layer to the SANS intensity by adding a correction term for the
solvent accessible beads only: the form factors of the amino acids / phosphate groups / pentose sugars /
nucleobases with a SASA (computed via LCPO algorithm) greater than a threshold are corrected according to an
electron density term. Both the surface cut-off threshold and the electron density term can be set by the user
with the SASA_CUTOFF and SOLVATION_CORRECTION keywords. Moreover, SASA stride calculation can be modified using
SOLVATION_STRIDE, which is set to 10 steps by default. The deuteration of the solvent-exposed residues is chosen
with a probability equal to the deuterium concentration in the buffer. The deuterated residues are updated with a
stride equal to SOLVATION_STRIDE. The fraction of deuterated water can be set with DEUTER_CONC, the default value
is 0.
ONEBEAD requires an additional PDB file to perform mapping conversion, which must be provided via TEMPLATE
keyword. This PDB file should only include the atoms for which the SANS intensity will be computed.
The AMBER OL3 (RNA) and OL15 (DNA) naming is required for nucleic acids.
Two additional bead types are available for DNA and RNA besides phosphate group, pentose sugar, and nucleobase:
- 5'-end pentose sugar capped with an hydroxyl moiety at C5' (the residue name in the PDB must be followed by "5",
e.g., DC5 or C5 for cytosine in DNA and RNA, respectively);
- 5'-phosphorylated end with an additional hydroxyl moiety at P (the residue name in the PDB must be followed by
"T", e.g., DCT or CT for cytosine in DNA and RNA, respectively);
- 3'-end pentose sugar capped with an hydroxyl moiety at C3' (the residue name in the PDB must be followed by "3",
e.g., DC3 or C3 for cytosine in DNA and RNA, respectively).
An additional bead type is available for proteins:
- Cysteine residue involved in disulfide bridge (the residue in the PDB must be named CYX).

PLEASE NOTE: at the moment, we DO NOT explicitly take into account deuterated residues in the ATOMISTIC
representation, but we correct the solvent contribution via the DEUTER_CONC keyword.

Experimental reference intensities can be added using the EXPINT keywords. All these values must be normalised
to the SANS intensity at q = 0. To facilitate this operation, the SCALE_EXPINT keyword can be used to provide
the intensity at q = 0. Each EXPINT is divided by SCALE_EXPINT.
The maximum QVALUE for ONEBEAD is set to 0.3 inverse angstroms.
The solvent density, that by default is set to 0.334 electrons per cubic angstrom (bulk water), can be modified
using the SOLVDENS keyword.

The ABSOLUTE flag can be used in order to calculate intensities in the absolute scale. It is only available for
the ATOMISTIC scheme and cannot be used with SCALE_EXPINT.

By default SANS is calculated using Debye on CPU, by adding the GPU flag it is possible to solve the equation on a
GPU if the ARRAYFIRE libraries are installed and correctly linked.
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

\par Examples
in the following example the SANS intensities are calculated at atomistic resolution. The form factors are assigned
according to the PDB file specified in the MOLINFO. Each experimental intensity is divided by 1.4002, which is the
corresponding theoretical intensity value at q = 0. The deuterated water fraction is set to 48%.

\plumedfile
MOLINFO STRUCTURE=template_AA.pdb

SANS ...
LABEL=SANS
ATOMS=1-355
ATOMISTIC
SCALE_EXPINT=1.4002
DEUTER_CONC=0.48
QVALUE1=0.03 EXPINT1=1.0902
QVALUE2=0.06 EXPINT2=0.790632
QVALUE3=0.09 EXPINT3=0.453808
QVALUE4=0.12 EXPINT4=0.254737
QVALUE5=0.15 EXPINT5=0.154928
QVALUE6=0.18 EXPINT6=0.0921503
QVALUE7=0.21 EXPINT7=0.052633
QVALUE8=0.24 EXPINT8=0.0276557
QVALUE9=0.27 EXPINT9=0.0122775
QVALUE10=0.30 EXPINT10=0.00880634
... SANS

PRINT ARG=(SANS\.q-.*),(SANS\.exp-.*) FILE=sansdata STRIDE=1

\endplumedfile

*/
//+ENDPLUMEDOC

class SAXS :
  public MetainferenceBase {
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
  enum { TRP,
         TYR,
         PHE,
         HIS,
         HIP,
         ARG,
         LYS,
         CYS,
         CYX,
         ASP,
         GLU,
         ILE,
         LEU,
         MET,
         ASN,
         PRO,
         GLN,
         SER,
         THR,
         VAL,
         ALA,
         GLY,
         BB_PO2,
         BB_PO3,
         BB_DNA,
         BB_DNA_5,
         BB_DNA_3,
         BB_RNA,
         BB_RNA_5,
         BB_RNA_3,
         BASE_A,
         BASE_C,
         BASE_T,
         BASE_G,
         BASE_U,
         NONEBEAD
       };
  struct SplineCoeffs {
    double a;
    double b;
    double c;
    double d;
    double x;
  };
  bool saxs;
  bool absolute;
  bool pbc;
  bool serial;
  bool gpu;
  bool onebead;
  bool resolution;
  bool isFirstStep;
  int  deviceid;
  unsigned nres;
  std::vector<unsigned> atoi;
  std::vector<unsigned> atoms_per_bead;
  std::vector<double>   atoms_masses;
  std::vector<double>   q_list;
  std::vector<double>   FF_rank;
  std::vector<std::vector<double> > FF_value_vacuum;
  std::vector<std::vector<double> > FF_value_solv;
  std::vector<std::vector<double> > FF_value_mixed;
  std::vector<std::vector<double> > FF_value;
  std::vector<std::vector<float> >  FFf_value;
  // SANS:
  std::vector<std::vector<double> > FF_value_vacuum_H;
  std::vector<std::vector<double> > FF_value_solv_H;
  std::vector<std::vector<double> > FF_value_mixed_H;
  std::vector<std::vector<double> > FF_value_vacuum_D;
  std::vector<std::vector<double> > FF_value_mixed_D;

  std::vector<std::vector<double> > LCPOparam;
  std::vector<unsigned> residue_atom;

  double rho, rho_corr, sasa_cutoff;
  double deuter_conc;
  unsigned solv_stride;
  std::vector<double> Iq0_vac;
  std::vector<double> Iq0_solv;
  std::vector<double> Iq0_mix;
  double Iq0;

  // SANS:
  std::vector<double> Iq0_vac_H;
  std::vector<double> Iq0_solv_H;
  std::vector<double> Iq0_mix_H;
  std::vector<double> Iq0_vac_D;
  std::vector<double> Iq0_mix_D;
  unsigned int Nj;
  std::vector<std::vector<double> > qj_list;
  std::vector<std::vector<double> > Rij;
  std::vector<double> sigma_res;

  // Chebyshev polynomial coefficients for i0e(x) used for resolution function
  // values taken from cephes library
  const std::vector<double> A = {
    -4.41534164647933937950E-18,
      3.33079451882223809783E-17,
      -2.43127984654795469359E-16,
      1.71539128555513303061E-15,
      -1.16853328779934516808E-14,
      7.67618549860493561688E-14,
      -4.85644678311192946090E-13,
      2.95505266312963983461E-12,
      -1.72682629144155570723E-11,
      9.67580903537323691224E-11,
      -5.18979560163526290666E-10,
      2.65982372468238665035E-9,
      -1.30002500998624804212E-8,
      6.04699502254191894932E-8,
      -2.67079385394061173391E-7,
      1.11738753912010371815E-6,
      -4.41673835845875056359E-6,
      1.64484480707288970893E-5,
      -5.75419501008210370398E-5,
      1.88502885095841655729E-4,
      -5.76375574538582365885E-4,
      1.63947561694133579842E-3,
      -4.32430999505057594430E-3,
      1.05464603945949983183E-2,
      -2.37374148058994688156E-2,
      4.93052842396707084878E-2,
      -9.49010970480476444210E-2,
      1.71620901522208775349E-1,
      -3.04682672343198398683E-1,
      6.76795274409476084995E-1
    };
  const std::vector<double> B = {
    -7.23318048787475395456E-18,
      -4.83050448594418207126E-18,
      4.46562142029675999901E-17,
      3.46122286769746109310E-17,
      -2.82762398051658348494E-16,
      -3.42548561967721913462E-16,
      1.77256013305652638360E-15,
      3.81168066935262242075E-15,
      -9.55484669882830764870E-15,
      -4.15056934728722208663E-14,
      1.54008621752140982691E-14,
      3.85277838274214270114E-13,
      7.18012445138366623367E-13,
      -1.79417853150680611778E-12,
      -1.32158118404477131188E-11,
      -3.14991652796324136454E-11,
      1.18891471078464383424E-11,
      4.94060238822496958910E-10,
      3.39623202570838634515E-9,
      2.26666899049817806459E-8,
      2.04891858946906374183E-7,
      2.89137052083475648297E-6,
      6.88975834691682398426E-5,
      3.36911647825569408990E-3,
      8.04490411014108831608E-1
    };

  void calculate_gpu(std::vector<Vector> &pos, std::vector<Vector> &deriv);
  void calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv);
  void getMartiniFFparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter);
  void getOnebeadparam(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac, std::vector<std::vector<long double> > &parameter_mix, std::vector<std::vector<long double> > &parameter_solv, const std::vector<unsigned> & residue_atom);
  unsigned getOnebeadMapping(const PDB &pdb, const std::vector<AtomNumber> &atoms);
  double calculateAFF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double rho);
  std::map<std::string, std::vector<double> > setupLCPOparam();
  void readLCPOparam(const std::vector<std::vector<std::string> > &AtomResidueName, unsigned natoms);
  void calcNlist(std::vector<std::vector<int> > &Nlist);
  void sasa_calculate(std::vector<bool> &solv_res);
  // SANS:
  void getOnebeadparam_sansH(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_H, std::vector<std::vector<long double> > &parameter_mix_H, std::vector<std::vector<long double> > &parameter_solv_H);
  void getOnebeadparam_sansD(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_D, std::vector<std::vector<long double> > &parameter_mix_D);
  double calculateAFFsans(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double deuter_conc);
  void resolution_function();
  std::vector<SplineCoeffs> spline_coeffs(std::vector<double> &x, std::vector<double> &y);
  inline double interpolation(std::vector<SplineCoeffs> &coeffs, double x);
  inline double i0e(double x);
  double chbevl(double x, const std::vector<double> &coeffs);

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXS(const ActionOptions&);
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(SAXS,"SAXS")
PLUMED_REGISTER_ACTION(SAXS,"SANS")

void SAXS::registerKeywords(Keywords& keys) {
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"Ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("compulsory","DEVICEID","-1","Identifier of the GPU to be used");
  keys.addFlag("GPU",false,"Calculate SAXS using ARRAYFIRE on an accelerator device");
  keys.addFlag("ABSOLUTE",false,"Absolute intensity: the intensities for each q-value are not normalised for the intensity at q=0.");
  keys.addFlag("ATOMISTIC",false,"Calculate SAXS for an atomistic model");
  keys.addFlag("MARTINI",false,"Calculate SAXS for a Martini model");
  keys.addFlag("ONEBEAD",false,"calculate SAXS for a single bead model");
  keys.add("compulsory","TEMPLATE","template.pdb","A PDB file is required for ONEBEAD mapping");
  keys.add("atoms","ATOMS","The atoms to be included in the calculation, e.g. the whole protein");
  keys.add("numbered","QVALUE","Selected scattering lengths in inverse angstroms are given as QVALUE1, QVALUE2, ...");
  keys.add("numbered","PARAMETERS","Used parameter Keywords like PARAMETERS1, PARAMETERS2. These are used to calculate the form factor for the \\f$i\\f$th atom/bead");
  keys.add("optional","PARAMETERSFILE","Read the PARAMETERS from a file");
  keys.add("compulsory","DEUTER_CONC","0.","Fraction of deuterated solvent");
  keys.add("compulsory","SOLVDENS","0.334","Density of the solvent to be used for the correction of atomistic form factors");
  keys.add("compulsory","SOLVATION_CORRECTION","0.0","Solvation layer electron density correction (ONEBEAD only)");
  keys.add("compulsory","SASA_CUTOFF","1.0","SASA value to consider a residue as exposed to the solvent (ONEBEAD only)");
  keys.add("numbered","EXPINT","Add an experimental value for each q value");
  keys.add("numbered","SIGMARES","Variance of Gaussian distribution describing the deviation in the scattering angle for each q value");
  keys.add("compulsory","N","10","Number of points in the resolution function integral");
  keys.add("compulsory","SOLVATION_STRIDE","10","Number of steps between every new residues solvation estimation via LCPO (ONEBEAD only)");
  keys.add("compulsory","SCALE_EXPINT","1.0","Scaling value for experimental data normalization");
  keys.addOutputComponent("q","default","scalar","The # SAXS of q");
  keys.addOutputComponent("exp","EXPINT","scalar","The # experimental intensity");
}

SAXS::SAXS(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  saxs(true),
  absolute(false),
  pbc(true),
  serial(false),
  gpu(false),
  onebead(false),
  isFirstStep(true),
  deviceid(-1) {
  if( getName().find("SAXS")!=std::string::npos) {
    saxs=true;
  } else if( getName().find("SANS")!=std::string::npos) {
    saxs=false;
  }

  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  unsigned size = atoms.size();

  parseFlag("SERIAL",serial);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  parseFlag("GPU",gpu);
#ifndef  __PLUMED_HAS_ARRAYFIRE
  if(gpu) {
    error("To use the GPU mode PLUMED must be compiled with ARRAYFIRE");
  }
#endif

  parse("DEVICEID",deviceid);
#ifdef  __PLUMED_HAS_ARRAYFIRE
  if(gpu&&comm.Get_rank()==0) {
    // if not set try to check the one set by the API
    if(deviceid==-1) {
      deviceid=plumed.getGpuDeviceId();
    }
    // if still not set use 0
    if(deviceid==-1) {
      deviceid=0;
    }
#ifdef  __PLUMED_HAS_ARRAYFIRE_CUDA
    af::setDevice(afcu::getNativeId(deviceid));
#elif   __PLUMED_HAS_ARRAYFIRE_OCL
    af::setDevice(afcl::getNativeId(deviceid));
#else
    af::setDevice(deviceid);
#endif
    af::info();
  }
#endif

  bool atomistic=false;
  parseFlag("ATOMISTIC",atomistic);
  if(atomistic) {
    log.printf("  using ATOMISTIC form factors\n");
  }
  bool martini=false;
  parseFlag("MARTINI",martini);
  if(martini) {
    log.printf("  using MARTINI form factors\n");
  }
  onebead=false;
  parseFlag("ONEBEAD",onebead);
  if(onebead) {
    log.printf("  using ONEBEAD form factors\n");
  }
  bool fromfile=false;
  std::string parametersfile;
  parse("PARAMETERSFILE",parametersfile);
  if (parametersfile.length() != 0) {
    fromfile=true;
  }
  if(fromfile) {
    log.printf("  will read form factors from file\n");
  }
  parseFlag("ABSOLUTE",absolute);

  if(martini&&atomistic) {
    error("You cannot use MARTINI and ATOMISTIC at the same time");
  }
  if(martini&&onebead) {
    error("You cannot use MARTINI and ONEBEAD at the same time");
  }
  if(onebead&&atomistic) {
    error("You cannot use ONEBEAD and ATOMISTIC at the same time");
  }
  if((martini)&&(!saxs)) {
    error("MARTINI cannot be used with SANS");
  }
  if((fromfile)&&((atomistic)||(martini)||(onebead))) {
    error("You cannot read parameters from file and use ATOMISTIC/MARTINI/ONEBEAD");
  }

  unsigned ntarget=0;
  for(unsigned i=0;; ++i) {
    double t_list;
    if( !parseNumbered( "QVALUE", i+1, t_list) ) {
      break;
    }
    if(t_list<=0.) {
      error("QVALUE cannot be less or equal to zero!\n");
    }
    if(onebead&&t_list>0.3) {
      error("ONEBEAD mapping QVALUE must be smaller or equal to 0.3");
    }
    q_list.push_back(t_list);
    ntarget++;
  }
  const unsigned numq = ntarget;

  for(unsigned i=0; i<numq; ++i) {
    if(q_list[i]==0.) {
      error("it is not possible to set q=0\n");
    }
    if(i>0&&q_list[i]<q_list[i-1]) {
      error("QVALUE must be in ascending order");
    }
    log.printf("  my q: %lf \n",q_list[i]);
  }

  rho = 0.334;
  parse("SOLVDENS", rho);
  log.printf("  Solvent density: %lf\n", rho);

  double scale_expint=1.;
  parse("SCALE_EXPINT",scale_expint);

  if((!atomistic&&absolute)||(absolute&&scale_expint!=1)) {
    error("ABSOLUTE can be used only combined with ATOMISTIC without SCALE_EXPINT");
  }
  if(atomistic) {
    log.printf("  Scale for intensities: %s\n", absolute ? "absolute" : "normalised");
  }

  double correction = 0.00;
  parse("SOLVATION_CORRECTION", correction);
  rho_corr=rho-correction;
  if(onebead) {
    log.printf("  Solvation density contribution: %lf\n", correction);
  }
  if((atomistic||martini||fromfile)&&(rho_corr!=rho)) {
    log.printf("  Solvation density contribution is taken into account in ONEBEAD only\n");
  }

  solv_stride = 10;
  parse("SOLVATION_STRIDE", solv_stride);
  if(solv_stride < 1.) {
    error("SOLVATION_STRIDE must be greater than 0");
  }
  if(onebead&&(rho_corr!=rho)) {
    log.printf("  SASA calculation stride: %u\n", solv_stride);
  }

  sasa_cutoff = 1.0;
  parse("SASA_CUTOFF", sasa_cutoff);
  if(sasa_cutoff <= 0.) {
    error("SASA_CUTOFF must be greater than 0");
  }

  deuter_conc = 0.;
  parse("DEUTER_CONC", deuter_conc);
  if ((deuter_conc)&&(fromfile)) {
    error("DEUTER_CONC cannot be used with PARAMETERSFILE");
  }
  if(deuter_conc < 0. || deuter_conc > 1.) {
    error("DEUTER_CONC must be in 0-1 range");
  }
  if ((atomistic||onebead)&&(!saxs)) {
    log.printf("  Solvent deuterium fraction: %lf/1.000000\n", deuter_conc);
  }

  PDB pdb;
  if(onebead) {
    std::string template_name;
    parse("TEMPLATE",template_name);
    log.printf("  Template for ONEBEAD mapping conversion: %s\n", template_name.c_str());
    if( !pdb.read(template_name,usingNaturalUnits(),1.) ) {
      plumed_merror("missing input file " + template_name);
    }
  }

  // preliminary mapping for onebead representation
  if(onebead) {
    LCPOparam.resize(size);
    nres = getOnebeadMapping(pdb, atoms);
    if(saxs) {
      Iq0_vac.resize(nres);
      Iq0_solv.resize(nres);
      Iq0_mix.resize(nres);
    } else { // SANS
      Iq0_vac_H.resize(nres);
      Iq0_solv_H.resize(nres);
      Iq0_mix_H.resize(nres);
      Iq0_vac_D.resize(nres);
      Iq0_mix_D.resize(nres);
    }
    atoi.resize(nres);
  } else {
    atoi.resize(size);
  }

  Iq0=0;
  std::vector<std::vector<long double> > FF_tmp;
  std::vector<std::vector<long double> > FF_tmp_vac;
  std::vector<std::vector<long double> > FF_tmp_mix;
  std::vector<std::vector<long double> > FF_tmp_solv;
  // SANS
  std::vector<std::vector<long double> > FF_tmp_vac_H;
  std::vector<std::vector<long double> > FF_tmp_mix_H;
  std::vector<std::vector<long double> > FF_tmp_solv_H;
  std::vector<std::vector<long double> > FF_tmp_vac_D;
  std::vector<std::vector<long double> > FF_tmp_mix_D;
  std::vector<std::vector<long double> > parameter_H;
  std::vector<std::vector<long double> > parameter_D;

  if(!atomistic&&!martini&&!onebead&&!fromfile) { // read PARAMETERS from PLUMED file
    if (saxs) {
      // read in parameter std::vector
      std::vector<std::vector<long double> > parameter;
      parameter.resize(size);
      ntarget=0;
      for(unsigned i=0; i<size; ++i) {
        if( !parseNumberedVector( "PARAMETERS", i+1, parameter[i]) ) {
          break;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of parameter std::vectors");
      }
      FF_tmp.resize(numq,std::vector<long double>(size));
      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter[i].size(); ++j) {
            FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i][0];
      }
      Iq0 *= Iq0;
    } else { // SANS
      std::vector<long double> parameter;
      parameter.resize(size);
      ntarget=0;
      for(unsigned i=0; i<size; ++i) {
        if( !parseNumbered( "PARAMETERS", i+1, parameter[i]) ) {
          break;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of parameter std::vectors");
      }
      FF_tmp.resize(numq,std::vector<long double>(size));
      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          FF_tmp[k][i]+= parameter[i];
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i];
      }
      Iq0 *= Iq0;
    }
  } else if (fromfile) { // read PARAMETERS from user-provided file
    log.printf("  Reading PARAMETERS from file: %s\n", parametersfile.c_str());
    if (saxs) {
      FF_tmp.resize(numq,std::vector<long double>(size));
      std::vector<std::vector<long double> > parameter;
      parameter.resize(size);

      IFile ifile;
      ifile.open(parametersfile);
      std::string line;

      ntarget=0;
      while(ifile.getline(line)) {
        Tools::ltrim(line);
        Tools::trimComments(line);
        if (line.empty()) {
          continue;
        }
        if (ntarget > size) {
          error("PARAMETERSFILE has more PARAMETERS than there are scattering centers");
        }
        std::string num;
        Tools::convert(ntarget+1,num);
        std::vector<std::string> lineread{line};
        if (!Tools::parseVector(lineread, "PARAMETERS"+num, parameter[ntarget], -1)) {
          error("Missing PARAMETERS or PARAMETERS not sorted");
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of PARAMETERS in file");
      }

      for(unsigned i=0; i<size; ++i) {
        atoi[i]=i;
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter[i].size(); ++j) {
            FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=parameter[i][0];
      }
      Iq0 *= Iq0;
    } else { // SANS
      FF_tmp.resize(numq,std::vector<long double>(size));

      IFile ifile;
      ifile.open(parametersfile);
      std::string line;

      ntarget=0;
      while(ifile.getline(line)) {
        Tools::ltrim(line);
        Tools::trimComments(line);
        if (line.empty()) {
          continue;
        }
        if (ntarget > size) {
          error("PARAMETERSFILE has more PARAMETERS than there are scattering centers");
        }
        std::string num;
        Tools::convert(ntarget+1,num);
        std::vector<std::string> lineread{line};
        long double scatlen;
        atoi[ntarget]=ntarget;
        if (!Tools::parse(lineread, "PARAMETERS"+num, scatlen, -1)) {
          error("Missing PARAMETERS or PARAMETERS not sorted");
        }
        for(unsigned k=0; k<numq; ++k) {
          FF_tmp[k][ntarget] = scatlen;
        }
        ntarget++;
      }
      if( ntarget!=size ) {
        error("found wrong number of PARAMETERS in file");
      }
      for(unsigned i=0; i<size; ++i) {
        Iq0+=FF_tmp[0][i];
      }
      Iq0 *= Iq0;
    }
  } else if(onebead) {
    if(saxs) {
      // read built-in ONEBEAD parameters
      FF_tmp_vac.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_solv.resize(numq,std::vector<long double>(NONEBEAD));
      std::vector<std::vector<long double> > parameter_vac(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix(NONEBEAD);
      std::vector<std::vector<long double> > parameter_solv(NONEBEAD);
      getOnebeadparam(pdb, atoms, parameter_vac, parameter_mix, parameter_solv, residue_atom);
      for(unsigned i=0; i<NONEBEAD; ++i) {
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter_vac[i].size(); ++j) {
            FF_tmp_vac[k][i]+= parameter_vac[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_mix[i].size(); ++j) {
            FF_tmp_mix[k][i]+= parameter_mix[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_solv[i].size(); ++j) {
            FF_tmp_solv[k][i]+= parameter_solv[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<nres; ++i) {
        Iq0_vac[i]=parameter_vac[atoi[i]][0];
        Iq0_mix[i]=parameter_mix[atoi[i]][0];
        Iq0_solv[i]=parameter_solv[atoi[i]][0];
      }
    } else { // SANS
      // read built-in ONEBEAD parameters
      FF_tmp_vac_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_solv_H.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_vac_D.resize(numq,std::vector<long double>(NONEBEAD));
      FF_tmp_mix_D.resize(numq,std::vector<long double>(NONEBEAD));
      std::vector<std::vector<long double> > parameter_vac_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_solv_H(NONEBEAD);
      std::vector<std::vector<long double> > parameter_vac_D(NONEBEAD);
      std::vector<std::vector<long double> > parameter_mix_D(NONEBEAD);
      getOnebeadparam_sansH(pdb, atoms, parameter_vac_H, parameter_mix_H, parameter_solv_H);
      getOnebeadparam_sansD(pdb, atoms, parameter_vac_D, parameter_mix_D);
      for(unsigned i=0; i<NONEBEAD; ++i) {
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned j=0; j<parameter_vac_H[i].size(); ++j) { // same number of parameters
            FF_tmp_vac_H[k][i]+= parameter_vac_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
            FF_tmp_vac_D[k][i]+= parameter_vac_D[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_mix_H[i].size(); ++j) {
            FF_tmp_mix_H[k][i]+= parameter_mix_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
            FF_tmp_mix_D[k][i]+= parameter_mix_D[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
          for(unsigned j=0; j<parameter_solv_H[i].size(); ++j) {
            FF_tmp_solv_H[k][i]+= parameter_solv_H[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
          }
        }
      }
      for(unsigned i=0; i<nres; ++i) {
        Iq0_vac_H[i]=parameter_vac_H[atoi[i]][0];
        Iq0_mix_H[i]=parameter_mix_H[atoi[i]][0];
        Iq0_solv_H[i]=parameter_solv_H[atoi[i]][0];
        Iq0_vac_D[i]=parameter_vac_D[atoi[i]][0];
        Iq0_mix_D[i]=parameter_mix_D[atoi[i]][0];
      }
    }
  } else if(martini) {
    // read built-in MARTINI parameters
    FF_tmp.resize(numq,std::vector<long double>(NMARTINI));
    std::vector<std::vector<long double> > parameter;
    parameter.resize(NMARTINI);
    getMartiniFFparam(atoms, parameter);
    for(unsigned i=0; i<NMARTINI; ++i) {
      for(unsigned k=0; k<numq; ++k) {
        for(unsigned j=0; j<parameter[i].size(); ++j) {
          FF_tmp[k][i]+= parameter[i][j]*std::pow(static_cast<long double>(q_list[k]),j);
        }
      }
    }
    for(unsigned i=0; i<size; ++i) {
      Iq0+=parameter[atoi[i]][0];
    }
    Iq0 *= Iq0;
  } else if(atomistic) {
    FF_tmp.resize(numq,std::vector<long double>(NTT));
    if(saxs) {
      Iq0=calculateAFF(atoms, FF_tmp, rho);
    } else {
      Iq0=calculateAFFsans(atoms, FF_tmp, deuter_conc);
    }
    Iq0 *= Iq0;
  }


  std::vector<double> expint;
  expint.resize( numq );
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "EXPINT", i+1, expint[i] ) ) {
      break;
    }
    ntarget++;
  }
  std::transform(expint.begin(), expint.begin() + ntarget, expint.begin(), [scale_expint](double x) {
    return x / scale_expint;
  });
  bool exp=false;
  if(ntarget!=numq && ntarget!=0) {
    error("found wrong number of EXPINT values");
  }
  if(ntarget==numq) {
    exp=true;
  }
  if(getDoScore()&&!exp) {
    error("with DOSCORE you need to set the EXPINT values");
  }

  sigma_res.resize( numq );
  resolution=false;
  ntarget=0;
  for(unsigned i=0; i<numq; ++i) {
    if( !parseNumbered( "SIGMARES", i+1, sigma_res[i] ) ) {
      break;
    }
    ntarget++;
  }
  if(ntarget!=numq && ntarget!=0) {
    error("found wrong number of SIGMARES values");
  }
  if(ntarget==numq) {
    resolution=true;
  }

  if(gpu && resolution) {
    error("Resolution function is not supported in GPUs");
  }

  Nj = 10;
  parse("N", Nj);
  if (Nj < 2) {
    error("N should be larger than 1");
  }
  if (resolution) {
    log.printf("  Resolution function with N: %d\n", Nj);
  }

  if(!gpu) {
    FF_rank.resize(numq);
    unsigned n_atom_types;
    if(onebead) {
      FF_value.resize(nres,std::vector<double>(numq));
      n_atom_types=NONEBEAD;
      if(saxs) {
        FF_value_vacuum.resize(n_atom_types,std::vector<double>(numq));
        FF_value_solv.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed.resize(n_atom_types,std::vector<double>(numq));
      } else {
        FF_value_vacuum_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_solv_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_vacuum_D.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed_D.resize(n_atom_types,std::vector<double>(numq));
      }
    } else {
      FF_value.resize(size,std::vector<double>(numq));
    }
    for(unsigned k=0; k<numq; ++k) {
      if(!onebead) {
        for(unsigned i=0; i<size; ++i) {
          FF_value[i][k] = static_cast<double>(FF_tmp[k][atoi[i]])/(std::sqrt(Iq0));
        }
        for(unsigned i=0; i<size; ++i) {
          FF_rank[k] += FF_value[i][k]*FF_value[i][k];
        }
      } else {
        if(saxs) {
          for(unsigned i=0; i<n_atom_types; ++i) {
            FF_value_vacuum[i][k] = static_cast<double>(FF_tmp_vac[k][i]);
            FF_value_mixed[i][k] = static_cast<double>(FF_tmp_mix[k][i]);
            FF_value_solv[i][k] = static_cast<double>(FF_tmp_solv[k][i]);
          }
        } else { // SANS
          for(unsigned i=0; i<n_atom_types; ++i) {
            FF_value_vacuum_H[i][k] = static_cast<double>(FF_tmp_vac_H[k][i]);
            FF_value_mixed_H[i][k] = static_cast<double>(FF_tmp_mix_H[k][i]);
            FF_value_solv_H[i][k] = static_cast<double>(FF_tmp_solv_H[k][i]);
            FF_value_vacuum_D[i][k] = static_cast<double>(FF_tmp_vac_D[k][i]);
            FF_value_mixed_D[i][k] = static_cast<double>(FF_tmp_mix_D[k][i]);
          }
        }
      }
    }
  } else {
    unsigned n_atom_types;
    if(onebead) {
      FFf_value.resize(numq,std::vector<float>(nres));
      n_atom_types=NONEBEAD;
      if(saxs) {
        FF_value_vacuum.resize(n_atom_types,std::vector<double>(numq));
        FF_value_solv.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed.resize(n_atom_types,std::vector<double>(numq));
      } else { // SANS
        FF_value_vacuum_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_solv_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed_H.resize(n_atom_types,std::vector<double>(numq));
        FF_value_vacuum_D.resize(n_atom_types,std::vector<double>(numq));
        FF_value_mixed_D.resize(n_atom_types,std::vector<double>(numq));
      }
    } else {
      FFf_value.resize(numq,std::vector<float>(size));
    }
    for(unsigned k=0; k<numq; ++k) {
      if(!onebead) {
        for(unsigned i=0; i<size; ++i) {
          FFf_value[k][i] = static_cast<float>(FF_tmp[k][atoi[i]])/(std::sqrt(Iq0));
        }
      } else {
        if(saxs) {
          for(unsigned i=0; i<n_atom_types; ++i) {
            FF_value_vacuum[i][k] = static_cast<double>(FF_tmp_vac[k][i]);
            FF_value_mixed[i][k] = static_cast<double>(FF_tmp_mix[k][i]);
            FF_value_solv[i][k] = static_cast<double>(FF_tmp_solv[k][i]);
          }
        } else { // SANS
          for(unsigned i=0; i<n_atom_types; ++i) {
            FF_value_vacuum_H[i][k] = static_cast<double>(FF_tmp_vac_H[k][i]);
            FF_value_mixed_H[i][k] = static_cast<double>(FF_tmp_mix_H[k][i]);
            FF_value_solv_H[i][k] = static_cast<double>(FF_tmp_solv_H[k][i]);
            FF_value_vacuum_D[i][k] = static_cast<double>(FF_tmp_vac_D[k][i]);
            FF_value_mixed_D[i][k] = static_cast<double>(FF_tmp_mix_D[k][i]);
          }
        }
      }
    }
  }

  if(!getDoScore()) {
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponentWithDerivatives("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    if(exp) {
      for(unsigned i=0; i<numq; ++i) {
        std::string num;
        Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(expint[i]);
      }
    }
  } else {
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponent("q-"+num);
      componentIsNotPeriodic("q-"+num);
    }
    for(unsigned i=0; i<numq; ++i) {
      std::string num;
      Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(expint[i]);
    }
  }

  // convert units to nm^-1
  for(unsigned i=0; i<numq; ++i) {
    q_list[i]=q_list[i]*10.0;    // factor 10 to convert from A^-1 to nm^-1
    if (resolution) {
      sigma_res[i]=sigma_res[i]*10.0;
    }
  }

  // compute resolution function after converting units
  if (resolution) {
    qj_list.resize(numq, std::vector<double>(Nj));
    Rij.resize(numq, std::vector<double>(Nj));
    // compute Rij and qj_list
    resolution_function();
  }

  log<<"  Bibliography ";
  if(onebead) {
    log<<plumed.cite("Ballabio, Paissoni, Bollati, de Rosa, Capelli, Camilloni, J. Chem. Theory Comput., 19, 22, 8401-8413 (2023)");
  }
  if(martini) {
    log<<plumed.cite("Niebling, Björling, Westenhoff, J. Appl. Crystallogr., 47, 1190–1198 (2014)");
    log<<plumed.cite("Paissoni, Jussupow, Camilloni, J. Appl. Crystallogr., 52, 394-402 (2019)");
  }
  if(atomistic) {
    log<<plumed.cite("Fraser, MacRae, Suzuki, J. Appl. Crystallogr., 11, 693–694 (1978)");
    log<<plumed.cite("Brown, Fox, Maslen, O'Keefe, Willis, International Tables for Crystallography, C, 554–595 (International Union of Crystallography, 2006)");
  }
  if(resolution) {
    log<<plumed.cite("Pedersen, Posselt, Mortensen, J. Appl. Crystallogr., 23, 321–333 (1990)");
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

// calculates SASA neighbor list
void SAXS::calcNlist(std::vector<std::vector<int> > &Nlist) {
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

// calculates SASA according to LCPO algorithm
void SAXS::sasa_calculate(std::vector<bool> &solv_res) {
  unsigned natoms = getNumberOfAtoms();
  std::vector<std::vector<int> > Nlist(natoms);
  calcNlist(Nlist);
  std::vector<double> sasares(nres, 0.);

  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    std::vector<double> private_sasares(nres, 0.);
    #pragma omp for
    for (unsigned i = 0; i < natoms; ++i) {
      if (LCPOparam[i].size() > 1 && LCPOparam[i][1] > 0.0) {
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
        if (sasai > 0) {
          private_sasares[residue_atom[i]] += sasai / 100.0;
        }
      }
    }
    #pragma omp critical
    {
      // combining private_sasares into sasares
      for (unsigned i = 0; i < nres; ++i) {
        sasares[i] += private_sasares[i];
      }
    }
  }
  for(unsigned i=0; i<nres; ++i) { // updating solv_res based on sasares
    if(sasares[i]>sasa_cutoff) {
      solv_res[i] = 1;
    } else {
      solv_res[i] = 0;
    }
  }
}

void SAXS::calculate_gpu(std::vector<Vector> &pos, std::vector<Vector> &deriv) {
#ifdef __PLUMED_HAS_ARRAYFIRE
  unsigned size;
  if(onebead) {
    size = nres;
  } else {
    size = getNumberOfAtoms();
  }
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
#ifdef  __PLUMED_HAS_ARRAYFIRE_CUDA
    af::setDevice(afcu::getNativeId(deviceid));
#elif   __PLUMED_HAS_ARRAYFIRE_OCL
    af::setDevice(afcl::getNativeId(deviceid));
#else
    af::setDevice(deviceid);
#endif
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
    std::string num;
    Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) {
      setCalcData(k, sum[k]);
    }
    for(unsigned i=0; i<size; ++i) {
      const unsigned di = k*size*3+i*3;
      deriv[k*size+i] = Vector(2.*dd[di+0],2.*dd[di+1],2.*dd[di+2]);
    }
  }
#endif
}

void SAXS::calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv) {
  unsigned size;
  if(onebead) {
    size = nres;
  } else {
    size = getNumberOfAtoms();
  }
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
      for(unsigned i=0; i<deriv.size(); ++i) {
        deriv[i]+=omp_deriv[i];
      }
      for(unsigned k=0; k<numq; ++k) {
        sum[k]+=omp_sum[k];
      }
    }
  }

  if(!serial) {
    comm.Sum(&deriv[0][0], 3*deriv.size());
    comm.Sum(&sum[0], numq);
  }

  if (resolution) {
    // get spline for scatering curve
    std::vector<SplineCoeffs> scatt_coeffs = spline_coeffs(q_list, sum);

    // get spline for the derivatives
    // copy the deriv to a new vector and zero deriv
    std::vector<Vector> old_deriv(deriv);
    memset(&deriv[0][0], 0.0, deriv.size() * sizeof deriv[0]);

    unsigned nt=OpenMP::getNumThreads();
    for (unsigned i=rank; i<size; i+=stride) {
      std::vector<double> deriv_i_x(numq);
      std::vector<double> deriv_i_y(numq);
      std::vector<double> deriv_i_z(numq);

      std::vector<SplineCoeffs> deriv_coeffs_x;
      std::vector<SplineCoeffs> deriv_coeffs_y;
      std::vector<SplineCoeffs> deriv_coeffs_z;
      for (unsigned k=0; k<numq; k++) {
        unsigned kdx = k*size;
        deriv_i_x[k] = old_deriv[kdx+i][0];
        deriv_i_y[k] = old_deriv[kdx+i][1];
        deriv_i_z[k] = old_deriv[kdx+i][2];
      }
      deriv_coeffs_x = spline_coeffs(q_list, deriv_i_x);
      deriv_coeffs_y = spline_coeffs(q_list, deriv_i_y);
      deriv_coeffs_z = spline_coeffs(q_list, deriv_i_z);

      // compute derivative with the smearing using the resolution function
      #pragma omp parallel for num_threads(nt)
      for (unsigned k=0; k<numq; k++) {
        unsigned kdx = k*size;
        double dq = qj_list[k][1] - qj_list[k][0];
        for (unsigned j=0; j<Nj; j++) {
          deriv[kdx+i][0] += Rij[k][j] * interpolation(deriv_coeffs_x, qj_list[k][j]) * dq;
          deriv[kdx+i][1] += Rij[k][j] * interpolation(deriv_coeffs_y, qj_list[k][j]) * dq;
          deriv[kdx+i][2] += Rij[k][j] * interpolation(deriv_coeffs_z, qj_list[k][j]) * dq;
        }
      }
    }

    if(!serial) {
      comm.Sum(&deriv[0][0], 3*deriv.size());
    }

    // compute the smeared spectra using the resolution function
    #pragma omp parallel for num_threads(nt)
    for (unsigned i=0; i<numq; i++) {
      sum[i] = 0.;
      double dq = qj_list[i][1] - qj_list[i][0];
      for (unsigned j=0; j<Nj; j++) {
        sum[i] += Rij[i][j] * interpolation(scatt_coeffs, qj_list[i][j]) * dq;
      }
    }
  }

  for (unsigned k=0; k<numq; ++k) {
    sum[k]+=FF_rank[k];
    std::string num;
    Tools::convert(k,num);
    Value* val=getPntrToComponent("q-"+num);
    val->set(sum[k]);
    if(getDoScore()) {
      setCalcData(k, sum[k]);
    }
  }
}

void SAXS::calculate() {
  if(pbc) {
    makeWhole();
  }

  const size_t size = getNumberOfAtoms();
  const size_t numq = q_list.size();

  // these are the derivatives associated to the coarse graining
  std::vector<Vector> aa_deriv(size);

  size_t beads_size = size;
  if(onebead) {
    beads_size = nres;
  }
  // these are the derivatives particle,q
  std::vector<Vector> bd_deriv(numq*beads_size);

  std::vector<Vector> beads_pos(beads_size);
  if(onebead) {
    for(unsigned resid=0; resid<nres; resid++) {
      double sum_mass = 0.;
      Vector sum_pos = Vector(0,0,0);
      for(unsigned atom_id=0; atom_id<size; atom_id++) {
        if(residue_atom[atom_id] == resid) {
          aa_deriv[atom_id] = Vector(atoms_masses[atom_id],atoms_masses[atom_id],atoms_masses[atom_id]);
          sum_pos += atoms_masses[atom_id] * getPosition(atom_id); // getPosition(first_atom+atom_id)
          sum_mass += atoms_masses[atom_id];
        }
      }
      beads_pos[resid] = sum_pos/sum_mass;
      for(unsigned atom_id=0; atom_id<size; atom_id++) {
        if(residue_atom[atom_id] == resid) {
          aa_deriv[atom_id] /= sum_mass;
        }
      }
    }
    // SASA
    std::vector<bool> solv_res(nres, 0);
    if(saxs) {
      if(getStep()%solv_stride == 0 || isFirstStep) {
        isFirstStep = 0;
        if(rho_corr!=rho) {
          sasa_calculate(solv_res);
        }
        Iq0=0.;
        for(unsigned i=0; i<nres; ++i) {
          if(solv_res[i] == 1 ) {
            Iq0 += std::sqrt((Iq0_vac[i]+(rho_corr*rho_corr)*Iq0_solv[i]-rho_corr*Iq0_mix[i]));
          } else {
            Iq0 += std::sqrt((Iq0_vac[i]+(rho*rho)*Iq0_solv[i]-rho*Iq0_mix[i]));
          }
        }
        // Form Factors
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned i=0; i<nres; ++i) {
            if(!gpu) {
              if(solv_res[i] == 0) { // buried
                FF_value[i][k] = std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho*rho*FF_value_solv[atoi[i]][k] - rho*FF_value_mixed[atoi[i]][k]))/Iq0;
              } else { // surface
                FF_value[i][k] = std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho_corr*rho_corr*FF_value_solv[atoi[i]][k] - rho_corr*FF_value_mixed[atoi[i]][k]))/Iq0;
              }
            } else {
              if(solv_res[i] == 0) { // buried
                FFf_value[k][i] = static_cast<float>(std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho*rho*FF_value_solv[atoi[i]][k] - rho*FF_value_mixed[atoi[i]][k]))/Iq0);
              } else { // surface
                FFf_value[k][i] = static_cast<float>(std::sqrt(std::fabs(FF_value_vacuum[atoi[i]][k] + rho_corr*rho_corr*FF_value_solv[atoi[i]][k] - rho_corr*FF_value_mixed[atoi[i]][k]))/Iq0);
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
    } else { // SANS
      std::vector<bool> deut_res(nres, 0);
      double solv_sc_length = 0.1*(0.580 + 2.*((1. - deuter_conc) * (-0.374) + deuter_conc * 0.667)); // per water electron (10 electrons)
      double rho_sans = rho * solv_sc_length;
      double rho_sans_corr = rho_corr * solv_sc_length;
      if(getStep()%solv_stride == 0 || isFirstStep) {
        isFirstStep = 0;
        if(deuter_conc!=0.||rho != rho_corr) {
          sasa_calculate(solv_res);
        }
        Iq0=0.;
        for(unsigned i=0; i<nres; ++i) {
          if(solv_res[i] == 1 ) {
            if(rand()/RAND_MAX<deuter_conc) {
              Iq0 += std::sqrt(std::fabs(Iq0_vac_D[i] + rho_sans_corr*rho_sans_corr*Iq0_solv_H[i] - rho_sans_corr*Iq0_mix_D[i]));
              deut_res[i] = 1;
            } else {
              Iq0 += std::sqrt(std::fabs(Iq0_vac_H[i] + rho_sans_corr*rho_sans_corr*Iq0_solv_H[i] - rho_sans_corr*Iq0_mix_H[i]));
            }
          } else {
            Iq0 += std::sqrt(std::fabs(Iq0_vac_H[i] + rho_sans*rho_sans*Iq0_solv_H[i] - rho_sans*Iq0_mix_H[i]));
          }
        }
        // Form Factors
        for(unsigned k=0; k<numq; ++k) {
          for(unsigned i=0; i<nres; ++i) {
            if(!gpu) {
              if(solv_res[i] == 0) { // hydrogen
                FF_value[i][k] = std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans*rho_sans*FF_value_solv_H[atoi[i]][k] - rho_sans*FF_value_mixed_H[atoi[i]][k]))/Iq0;
              } else {
                if(deut_res[i] == 0) {
                  FF_value[i][k] = std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_H[atoi[i]][k]))/Iq0;
                } else {
                  FF_value[i][k] = std::sqrt(std::fabs(FF_value_vacuum_D[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_D[atoi[i]][k]))/Iq0;
                }
              }
            } else {
              if(solv_res[i] == 0) { // hydrogen
                FFf_value[k][i] = static_cast<float>(std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans*rho_sans*FF_value_solv_H[atoi[i]][k] - rho_sans*FF_value_mixed_H[atoi[i]][k]))/Iq0);
              } else {
                if(deut_res[i] == 0) {
                  FFf_value[k][i] = static_cast<float>(std::sqrt(std::fabs(FF_value_vacuum_H[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_H[atoi[i]][k]))/Iq0);
                } else {
                  FFf_value[k][i] = static_cast<float>(std::sqrt(std::fabs(FF_value_vacuum_D[atoi[i]][k] + rho_sans_corr*rho_sans_corr*FF_value_solv_H[atoi[i]][k] - rho_sans_corr*FF_value_mixed_D[atoi[i]][k]))/Iq0);
                }
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
    }
    // not ONEBEAD
  } else {
    for(unsigned i=0; i<size; ++i) {
      beads_pos[i] = getPosition(i);
    }
    aa_deriv = std::vector<Vector>(size,(Vector(1,1,1)));
  }

  if(gpu) {
    calculate_gpu(beads_pos, bd_deriv);
  } else {
    calculate_cpu(beads_pos, bd_deriv);
  }

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
      std::string num;
      Tools::convert(k,num);
      val=getPntrToComponent("q-"+num);

      if(onebead) {
        unsigned atom_id=0;
        for(unsigned i=0; i<beads_size; ++i) {
          for(unsigned j=0; j<atoms_per_bead[i]; ++j) {
            setAtomsDerivatives(val, atom_id, Vector(aa_deriv[atom_id][0]*bd_deriv[kdx+i][0], \
                                aa_deriv[atom_id][1]*bd_deriv[kdx+i][1], \
                                aa_deriv[atom_id][2]*bd_deriv[kdx+i][2]) );
            deriv_box += Tensor(getPosition(atom_id),Vector(aa_deriv[atom_id][0]*bd_deriv[kdx+i][0], \
                                aa_deriv[atom_id][1]*bd_deriv[kdx+i][1], \
                                aa_deriv[atom_id][2]*bd_deriv[kdx+i][2]) );
            atom_id++;
          }
        }
      } else {
        for(unsigned i=0; i<beads_size; ++i) {
          setAtomsDerivatives(val, i, Vector(bd_deriv[kdx+i][0], \
                                             bd_deriv[kdx+i][1], \
                                             bd_deriv[kdx+i][2]) );
          deriv_box += Tensor(getPosition(i),Vector(bd_deriv[kdx+i][0], \
                              bd_deriv[kdx+i][1], \
                              bd_deriv[kdx+i][2]) );
        }
      }
    } else {
      val=getPntrToComponent("score");
      if(onebead) {
        unsigned atom_id=0;
        for(unsigned i=0; i<beads_size; ++i) {
          for(unsigned j=0; j<atoms_per_bead[i]; ++j) {
            setAtomsDerivatives(val, atom_id, Vector(aa_deriv[atom_id][0]*bd_deriv[kdx+i][0]*getMetaDer(k),
                                aa_deriv[atom_id][1]*bd_deriv[kdx+i][1]*getMetaDer(k),
                                aa_deriv[atom_id][2]*bd_deriv[kdx+i][2]*getMetaDer(k)) );
            deriv_box += Tensor(getPosition(atom_id),Vector(aa_deriv[atom_id][0]*bd_deriv[kdx+i][0]*getMetaDer(k),
                                aa_deriv[atom_id][1]*bd_deriv[kdx+i][1]*getMetaDer(k),
                                aa_deriv[atom_id][2]*bd_deriv[kdx+i][2]*getMetaDer(k)) );
            atom_id++;
          }
        }
      } else {
        for(unsigned i=0; i<beads_size; ++i) {
          setAtomsDerivatives(val, i, Vector(bd_deriv[kdx+i][0]*getMetaDer(k),
                                             bd_deriv[kdx+i][1]*getMetaDer(k),
                                             bd_deriv[kdx+i][2]*getMetaDer(k)) );
          deriv_box += Tensor(getPosition(i),Vector(bd_deriv[kdx+i][0]*getMetaDer(k),
                              bd_deriv[kdx+i][1]*getMetaDer(k),
                              bd_deriv[kdx+i][2]*getMetaDer(k)) );
        }
      }
    }
    setBoxDerivatives(val, -deriv_box);
  }
}

void SAXS::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) {
    writeStatus();
  }
}

unsigned SAXS::getOnebeadMapping(const PDB &pdb, const std::vector<AtomNumber> &atoms) {
  std::vector<std::string> chains;
  pdb.getChainNames( chains );
  std::vector<std::vector<std::string> > AtomResidueName;

  // cycle over chains
  for(unsigned i=0; i<chains.size(); ++i) {
    unsigned start, end;
    std::string errmsg;
    pdb.getResidueRange(chains[i], start, end, errmsg);
    AtomResidueName.resize(2);
    // cycle over residues
    for(unsigned res=start; res<=end; res++) {
      std::string Rname = pdb.getResidueName(res, chains[i]);
      Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
      std::vector<AtomNumber> res_atoms = pdb.getAtomsInResidue(res, chains[i]);
      std::vector<unsigned> tmp_residue_atom;
      tmp_residue_atom.resize(3,0);
      // cycle over atoms
      for(unsigned a=0; a<res_atoms.size(); a++) {
        // operations shared among all beads
        std::string Aname=pdb.getAtomName(res_atoms[a]);
        AtomResidueName[0].push_back(Aname);
        AtomResidueName[1].push_back(Rname);
        char type;
        char first = Aname.at(0);
        // We assume that element symbol is first letter, if not a number
        if (!isdigit(first)) {
          type = first;
          // otherwise is the second
        } else {
          type = Aname.at(1);
        }
        if (type == 'H') {
          atoms_masses.push_back(1.008);
        } else if(type == 'C') {
          atoms_masses.push_back(12.011);
        } else if(type == 'N') {
          atoms_masses.push_back(14.007);
        } else if(type == 'O') {
          atoms_masses.push_back(15.999);
        } else if(type == 'S') {
          atoms_masses.push_back(32.065);
        } else if(type == 'P') {
          atoms_masses.push_back(30.974);
        } else {
          error("Unknown element in mass extraction\n");
        }
        if(pdb.allowedResidue("protein",Rname)) {
          residue_atom.push_back(atoms_per_bead.size());
        } else {
          // check for nucleic acids
          // Pentose bead
          if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  || Aname=="O3'"  ||
              Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  || Aname=="C1'"  || Aname=="H5'"  ||
              Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  || Aname=="H2'"  || Aname=="H2''" ||
              Aname=="H2'2" || Aname=="H1'"  || Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" ||
              Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T"  ||
              Aname=="H3T" ) {
            residue_atom.push_back(atoms_per_bead.size()+0);
            tmp_residue_atom[0]++;
          }
          // Nucleobase bead
          else if(Aname=="N1"  || Aname=="N2"  || Aname=="N3"  || Aname=="N4"  || Aname=="N6"  ||
                  Aname=="N7"  || Aname=="N9"  || Aname=="C2"  || Aname=="C4"  || Aname=="C5"  ||
                  Aname=="C6"  || Aname=="C7"  || Aname=="C8"  || Aname=="O2"  || Aname=="O4"  ||
                  Aname=="O6"  || Aname=="H1"  || Aname=="H2"  || Aname=="H3"  || Aname=="H5"  ||
                  Aname=="H6"  || Aname=="H8"  || Aname=="H21" || Aname=="H22" || Aname=="H41" ||
                  Aname=="H42" || Aname=="H61" || Aname=="H62" || Aname=="H71" || Aname=="H72" ||
                  Aname=="H73" ) {
            residue_atom.push_back(atoms_per_bead.size()+1);
            tmp_residue_atom[1]++;
          }
          // PO bead
          else if(Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" || Aname=="O1P" ||
                  Aname=="O2P" || Aname=="O3P" || Aname=="HP"  || Aname=="HOP3" ) {
            residue_atom.push_back(atoms_per_bead.size()+2);
            tmp_residue_atom[2]++;
          }
          // error
          else {
            error("Atom name "+Aname+" cannot be indexed to any bead. Check the PDB.");
          }
        }
      }
      if(pdb.allowedResidue("protein",Rname)) {
        atoms_per_bead.push_back(res_atoms.size());
      } else {
        atoms_per_bead.push_back(tmp_residue_atom[0]);
        atoms_per_bead.push_back(tmp_residue_atom[1]);
        if(tmp_residue_atom[2]>0) {
          atoms_per_bead.push_back(tmp_residue_atom[2]);
        }
      }
    }
  }
  readLCPOparam(AtomResidueName, atoms.size());
  return atoms_per_bead.size();
}

void SAXS::getMartiniFFparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter) {
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
    for(unsigned i=0; i<atoms.size(); ++i) {
      std::string Aname = moldat->getAtomName(atoms[i]);
      std::string Rname = moldat->getResidueName(atoms[i]);
      if(Rname=="ALA") {
        if(Aname=="BB") {
          atoi[i]=ALA_BB;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="ARG") {
        if(Aname=="BB") {
          atoi[i]=ARG_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ARG_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=ARG_SC2;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="ASN") {
        if(Aname=="BB") {
          atoi[i]=ASN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASN_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="ASP") {
        if(Aname=="BB") {
          atoi[i]=ASP_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ASP_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="CYS") {
        if(Aname=="BB") {
          atoi[i]=CYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=CYS_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="GLN") {
        if(Aname=="BB") {
          atoi[i]=GLN_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLN_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="GLU") {
        if(Aname=="BB") {
          atoi[i]=GLU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=GLU_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="GLY") {
        if(Aname=="BB") {
          atoi[i]=GLY_BB;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="HIS") {
        if(Aname=="BB") {
          atoi[i]=HIS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=HIS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=HIS_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=HIS_SC3;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="ILE") {
        if(Aname=="BB") {
          atoi[i]=ILE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=ILE_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="LEU") {
        if(Aname=="BB") {
          atoi[i]=LEU_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LEU_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="LYS") {
        if(Aname=="BB") {
          atoi[i]=LYS_BB;
        } else if(Aname=="SC1") {
          atoi[i]=LYS_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=LYS_SC2;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="MET") {
        if(Aname=="BB") {
          atoi[i]=MET_BB;
        } else if(Aname=="SC1") {
          atoi[i]=MET_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="PHE") {
        if(Aname=="BB") {
          atoi[i]=PHE_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PHE_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=PHE_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=PHE_SC3;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="PRO") {
        if(Aname=="BB") {
          atoi[i]=PRO_BB;
        } else if(Aname=="SC1") {
          atoi[i]=PRO_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="SER") {
        if(Aname=="BB") {
          atoi[i]=SER_BB;
        } else if(Aname=="SC1") {
          atoi[i]=SER_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="THR") {
        if(Aname=="BB") {
          atoi[i]=THR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=THR_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="TYR") {
        if(Aname=="BB") {
          atoi[i]=TYR_BB;
        } else if(Aname=="SC1") {
          atoi[i]=TYR_SC1;
        } else if(Aname=="SC2") {
          atoi[i]=TYR_SC2;
        } else if(Aname=="SC3") {
          atoi[i]=TYR_SC3;
        } else {
          error("Atom name not known: "+Aname);
        }
      } else if(Rname=="VAL") {
        if(Aname=="BB") {
          atoi[i]=VAL_BB;
        } else if(Aname=="SC1") {
          atoi[i]=VAL_SC1;
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
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
        } else {
          error("Atom name not known: "+Aname);
        }
      } else {
        error("Residue not known: "+Rname);
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
}

void SAXS::getOnebeadparam(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac, std::vector<std::vector<long double> > &parameter_mix, std::vector<std::vector<long double> > &parameter_solv, const std::vector<unsigned> & residue_atom) {

  parameter_solv[TRP].push_back(60737.60249988003);
  parameter_solv[TRP].push_back(-77.75716755173752);
  parameter_solv[TRP].push_back(-205962.98557711052);
  parameter_solv[TRP].push_back(-62013.46984155453);
  parameter_solv[TRP].push_back(680710.7592231638);
  parameter_solv[TRP].push_back(-681336.8777362367);
  parameter_solv[TRP].push_back(211473.65530642506);

  parameter_solv[TYR].push_back(46250.80359987982);
  parameter_solv[TYR].push_back(-45.8287864681578);
  parameter_solv[TYR].push_back(-143872.91752817619);
  parameter_solv[TYR].push_back(-39049.68736409533);
  parameter_solv[TYR].push_back(441321.71874090104);
  parameter_solv[TYR].push_back(-434478.0972346327);
  parameter_solv[TYR].push_back(133179.3694641212);

  parameter_solv[PHE].push_back(42407.164900118914);
  parameter_solv[PHE].push_back(-159.1980754191431);
  parameter_solv[PHE].push_back(-123847.86192757386);
  parameter_solv[PHE].push_back(-41797.69041575073);
  parameter_solv[PHE].push_back(380283.7035277073);
  parameter_solv[PHE].push_back(-361432.67247521743);
  parameter_solv[PHE].push_back(107750.64978068044);

  parameter_solv[HIP].push_back(24473.47360011923);
  parameter_solv[HIP].push_back(-111.64156672747428);
  parameter_solv[HIP].push_back(-65826.16993707925);
  parameter_solv[HIP].push_back(-23305.91329798928);
  parameter_solv[HIP].push_back(194795.11911635034);
  parameter_solv[HIP].push_back(-180454.49458095312);
  parameter_solv[HIP].push_back(52699.374196745615);

  parameter_solv[ARG].push_back(34106.70239988039);
  parameter_solv[ARG].push_back(152.7472727640246);
  parameter_solv[ARG].push_back(-117086.49392248681);
  parameter_solv[ARG].push_back(-19664.229479267167);
  parameter_solv[ARG].push_back(364454.0909203641);
  parameter_solv[ARG].push_back(-382075.8018312776);
  parameter_solv[ARG].push_back(122775.75036605193);

  parameter_solv[LYS].push_back(32292.090000118922);
  parameter_solv[LYS].push_back(-111.97371180593888);
  parameter_solv[LYS].push_back(-91953.10997619898);
  parameter_solv[LYS].push_back(-30690.807047993283);
  parameter_solv[LYS].push_back(282092.40760143084);
  parameter_solv[LYS].push_back(-269503.2592457489);
  parameter_solv[LYS].push_back(80777.81552915688);

  parameter_solv[CYS].push_back(11352.902500119093);
  parameter_solv[CYS].push_back(-45.5226331859686);
  parameter_solv[CYS].push_back(-20925.085562607524);
  parameter_solv[CYS].push_back(-5662.685408989286);
  parameter_solv[CYS].push_back(38559.10376731146);
  parameter_solv[CYS].push_back(-27885.23426006181);
  parameter_solv[CYS].push_back(6280.15058191397);

  parameter_solv[CYX].push_back(10281.960000119348);
  parameter_solv[CYX].push_back(-42.315998754511);
  parameter_solv[CYX].push_back(-19328.174487327480);
  parameter_solv[CYX].push_back(-5523.191775626829);
  parameter_solv[CYX].push_back(38185.463172673335);
  parameter_solv[CYX].push_back(-28940.174693042034);
  parameter_solv[CYX].push_back(6925.390014187676);

  parameter_solv[ASP].push_back(13511.73760011933);
  parameter_solv[ASP].push_back(-59.929111107656595);
  parameter_solv[ASP].push_back(-25849.869639655575);
  parameter_solv[ASP].push_back(-7541.669448872824);
  parameter_solv[ASP].push_back(50760.92045144903);
  parameter_solv[ASP].push_back(-37677.87583269734);
  parameter_solv[ASP].push_back(8745.7056219399);

  parameter_solv[GLU].push_back(20443.280400119456);
  parameter_solv[GLU].push_back(-113.77561814283207);
  parameter_solv[GLU].push_back(-45587.79314626863);
  parameter_solv[GLU].push_back(-16187.556837331254);
  parameter_solv[GLU].push_back(112609.65830609271);
  parameter_solv[GLU].push_back(-93362.05323205091);
  parameter_solv[GLU].push_back(24519.557866124724);

  parameter_solv[ILE].push_back(27858.948100119596);
  parameter_solv[ILE].push_back(-159.27355145839834);
  parameter_solv[ILE].push_back(-61571.43463039565);
  parameter_solv[ILE].push_back(-21324.879474559468);
  parameter_solv[ILE].push_back(144070.7572894681);
  parameter_solv[ILE].push_back(-115021.81959095894);
  parameter_solv[ILE].push_back(28939.085108838968);

  parameter_solv[LEU].push_back(27858.948100119596);
  parameter_solv[LEU].push_back(-165.61892007509647);
  parameter_solv[LEU].push_back(-62564.568746500125);
  parameter_solv[LEU].push_back(-22465.332149768525);
  parameter_solv[LEU].push_back(151616.79489291538);
  parameter_solv[LEU].push_back(-122905.6119395393);
  parameter_solv[LEU].push_back(31436.664377885514);

  parameter_solv[MET].push_back(25609.60090011981);
  parameter_solv[MET].push_back(-135.38857843066708);
  parameter_solv[MET].push_back(-67771.01108177133);
  parameter_solv[MET].push_back(-25228.934337676077);
  parameter_solv[MET].push_back(199649.95030712147);
  parameter_solv[MET].push_back(-182251.94895101967);
  parameter_solv[MET].push_back(52502.88444247481);

  parameter_solv[ASN].push_back(14376.010000119095);
  parameter_solv[ASN].push_back(-67.65579048748472);
  parameter_solv[ASN].push_back(-28302.87809850141);
  parameter_solv[ASN].push_back(-8577.439830985548);
  parameter_solv[ASN].push_back(57532.879075695324);
  parameter_solv[ASN].push_back(-43261.79286366774);
  parameter_solv[ASN].push_back(10186.448634149085);

  parameter_solv[PRO].push_back(16866.21690011944);
  parameter_solv[PRO].push_back(-70.84327801054884);
  parameter_solv[PRO].push_back(-31465.84064925844);
  parameter_solv[PRO].push_back(-8653.3693368317);
  parameter_solv[PRO].push_back(58032.28250733714);
  parameter_solv[PRO].push_back(-41521.01146771431);
  parameter_solv[PRO].push_back(9184.530596102064);

  parameter_solv[GLN].push_back(21503.289600119);
  parameter_solv[GLN].push_back(-121.30164008960246);
  parameter_solv[GLN].push_back(-50468.580981118175);
  parameter_solv[GLN].push_back(-18462.49098408308);
  parameter_solv[GLN].push_back(132718.44904081387);
  parameter_solv[GLN].push_back(-113787.22666510186);
  parameter_solv[GLN].push_back(30920.348610969988);

  parameter_solv[SER].push_back(9181.472400119354);
  parameter_solv[SER].push_back(-28.77519915767741);
  parameter_solv[SER].push_back(-15205.543144104717);
  parameter_solv[SER].push_back(-3377.782176346411);
  parameter_solv[SER].push_back(23345.555771001076);
  parameter_solv[SER].push_back(-15312.694356014094);
  parameter_solv[SER].push_back(3013.8428466148);

  parameter_solv[THR].push_back(15020.953600119403);
  parameter_solv[THR].push_back(-61.91004832631006);
  parameter_solv[THR].push_back(-27814.537889259853);
  parameter_solv[THR].push_back(-7532.227289701552);
  parameter_solv[THR].push_back(50586.30566118166);
  parameter_solv[THR].push_back(-35943.866131120165);
  parameter_solv[THR].push_back(7880.093558764326);

  parameter_solv[VAL].push_back(19647.628900119355);
  parameter_solv[VAL].push_back(-89.04983250107853);
  parameter_solv[VAL].push_back(-38050.09958470928);
  parameter_solv[VAL].push_back(-10921.427112288537);
  parameter_solv[VAL].push_back(72774.32322962297);
  parameter_solv[VAL].push_back(-52689.060152305225);
  parameter_solv[VAL].push_back(11806.492503632868);

  parameter_solv[ALA].push_back(7515.156100119276);
  parameter_solv[ALA].push_back(-20.226381685697746);
  parameter_solv[ALA].push_back(-11761.841094237716);
  parameter_solv[ALA].push_back(-2341.4929468980367);
  parameter_solv[ALA].push_back(16545.385777961936);
  parameter_solv[ALA].push_back(-10397.175253025776);
  parameter_solv[ALA].push_back(1921.5264606725107);

  parameter_solv[GLY].push_back(3594.002500119159);
  parameter_solv[GLY].push_back(-6.910836154887606);
  parameter_solv[GLY].push_back(-4937.354220666574);
  parameter_solv[GLY].push_back(-785.4549468992149);
  parameter_solv[GLY].push_back(5852.854429532936);
  parameter_solv[GLY].push_back(-3391.2927115487832);
  parameter_solv[GLY].push_back(552.3280571490722);

  parameter_solv[HIS].push_back(22888.664100119073);
  parameter_solv[HIS].push_back(-133.86265270962434);
  parameter_solv[HIS].push_back(-57533.51591635819);
  parameter_solv[HIS].push_back(-21767.293192014684);
  parameter_solv[HIS].push_back(161255.14120001195);
  parameter_solv[HIS].push_back(-142176.64081149307);
  parameter_solv[HIS].push_back(39642.61185646193);

  parameter_mix[TRP].push_back(48294.0117571196);
  parameter_mix[TRP].push_back(-205.45879626487798);
  parameter_mix[TRP].push_back(-148816.1858118254);
  parameter_mix[TRP].push_back(-54968.030079609875);
  parameter_mix[TRP].push_back(491793.79967057955);
  parameter_mix[TRP].push_back(-476312.9117969879);
  parameter_mix[TRP].push_back(144159.96165644142);

  parameter_mix[TYR].push_back(36984.20240312081);
  parameter_mix[TYR].push_back(-83.86380083812203);
  parameter_mix[TYR].push_back(-108820.52211887162);
  parameter_mix[TYR].push_back(-33934.69818901515);
  parameter_mix[TYR].push_back(341504.736372253);
  parameter_mix[TYR].push_back(-334008.1748614056);
  parameter_mix[TYR].push_back(102033.08077851454);

  parameter_mix[PHE].push_back(32119.469231338233);
  parameter_mix[PHE].push_back(-172.96940450568917);
  parameter_mix[PHE].push_back(-85831.4326887122);
  parameter_mix[PHE].push_back(-33193.32405438845);
  parameter_mix[PHE].push_back(262940.64471909316);
  parameter_mix[PHE].push_back(-243540.06898907054);
  parameter_mix[PHE].push_back(71084.54387480798);

  parameter_mix[HIP].push_back(22833.36414923898);
  parameter_mix[HIP].push_back(-134.0493955562186);
  parameter_mix[HIP].push_back(-55325.55607328898);
  parameter_mix[HIP].push_back(-21898.314938881984);
  parameter_mix[HIP].push_back(159995.6912885654);
  parameter_mix[HIP].push_back(-142968.19796084083);
  parameter_mix[HIP].push_back(40417.44581470003);

  parameter_mix[ARG].push_back(31385.401600920715);
  parameter_mix[ARG].push_back(36.114094042884254);
  parameter_mix[ARG].push_back(-103730.44467490204);
  parameter_mix[ARG].push_back(-27036.249157905615);
  parameter_mix[ARG].push_back(347011.0339314942);
  parameter_mix[ARG].push_back(-358879.9736802336);
  parameter_mix[ARG].push_back(114432.18361399164);

  parameter_mix[LYS].push_back(25511.35812671878);
  parameter_mix[LYS].push_back(-130.4381491986372);
  parameter_mix[LYS].push_back(-69258.61236879184);
  parameter_mix[LYS].push_back(-27066.36783798707);
  parameter_mix[LYS].push_back(220092.65231165203);
  parameter_mix[LYS].push_back(-207794.5056092443);
  parameter_mix[LYS].push_back(61665.57004630315);

  parameter_mix[CYS].push_back(11505.517261618916);
  parameter_mix[CYS].push_back(-33.60468076978334);
  parameter_mix[CYS].push_back(-18328.882710004465);
  parameter_mix[CYS].push_back(-3956.9113649567626);
  parameter_mix[CYS].push_back(27546.35146501212);
  parameter_mix[CYS].push_back(-18024.826330595406);
  parameter_mix[CYS].push_back(3551.2207387570024);

  parameter_mix[CYX].push_back(10746.617793719070);
  parameter_mix[CYX].push_back(-37.082746200650);
  parameter_mix[CYX].push_back(-17871.552278655203);
  parameter_mix[CYX].push_back(-4512.203184574789);
  parameter_mix[CYX].push_back(30605.726711712588);
  parameter_mix[CYX].push_back(-21530.684072275839);
  parameter_mix[CYX].push_back(4694.601090758420);

  parameter_mix[ASP].push_back(13713.858501879382);
  parameter_mix[ASP].push_back(-51.33286241257164);
  parameter_mix[ASP].push_back(-23807.8549764091);
  parameter_mix[ASP].push_back(-6153.667315935503);
  parameter_mix[ASP].push_back(41296.118377286424);
  parameter_mix[ASP].push_back(-28740.28391184026);
  parameter_mix[ASP].push_back(6132.671533319127);

  parameter_mix[GLU].push_back(19156.03660739947);
  parameter_mix[GLU].push_back(-110.90600703589246);
  parameter_mix[GLU].push_back(-40319.3351514524);
  parameter_mix[GLU].push_back(-14679.813393816446);
  parameter_mix[GLU].push_back(96769.28565573556);
  parameter_mix[GLU].push_back(-77909.09315520026);
  parameter_mix[GLU].push_back(19770.047062759568);

  parameter_mix[ILE].push_back(20693.06215917923);
  parameter_mix[ILE].push_back(-102.87208880594848);
  parameter_mix[ILE].push_back(-41080.44036311675);
  parameter_mix[ILE].push_back(-12874.439649378206);
  parameter_mix[ILE].push_back(84947.33147117581);
  parameter_mix[ILE].push_back(-63779.07871450237);
  parameter_mix[ILE].push_back(14938.919981690511);

  parameter_mix[LEU].push_back(20693.062159179233);
  parameter_mix[LEU].push_back(-114.09539845409269);
  parameter_mix[LEU].push_back(-42417.3431074524);
  parameter_mix[LEU].push_back(-14393.801090829746);
  parameter_mix[LEU].push_back(93640.48403643962);
  parameter_mix[LEU].push_back(-71990.10354816525);
  parameter_mix[LEU].push_back(17299.01082057651);

  parameter_mix[MET].push_back(22400.800002738917);
  parameter_mix[MET].push_back(-138.14469221559762);
  parameter_mix[MET].push_back(-53013.97694299946);
  parameter_mix[MET].push_back(-21079.899452619244);
  parameter_mix[MET].push_back(148607.1089339919);
  parameter_mix[MET].push_back(-129827.63962878387);
  parameter_mix[MET].push_back(35882.3297822684);

  parameter_mix[ASN].push_back(14384.287416519475);
  parameter_mix[ASN].push_back(-55.24976731179147);
  parameter_mix[ASN].push_back(-25372.978199926372);
  parameter_mix[ASN].push_back(-6646.452004616925);
  parameter_mix[ASN].push_back(44594.5027556148);
  parameter_mix[ASN].push_back(-31202.511764907107);
  parameter_mix[ASN].push_back(6703.764135873442);

  parameter_mix[PRO].push_back(13503.797145659117);
  parameter_mix[PRO].push_back(-38.58316011847087);
  parameter_mix[PRO].push_back(-21446.17847324053);
  parameter_mix[PRO].push_back(-4480.55896170459);
  parameter_mix[PRO].push_back(31274.287350083254);
  parameter_mix[PRO].push_back(-19984.249229169505);
  parameter_mix[PRO].push_back(3782.272312712745);

  parameter_mix[GLN].push_back(19938.23724683901);
  parameter_mix[GLN].push_back(-121.24884503048865);
  parameter_mix[GLN].push_back(-43928.589472297834);
  parameter_mix[GLN].push_back(-16805.069757865473);
  parameter_mix[GLN].push_back(112831.61348476357);
  parameter_mix[GLN].push_back(-93979.08819184235);
  parameter_mix[GLN].push_back(24741.563493163732);

  parameter_mix[SER].push_back(8813.67020471935);
  parameter_mix[SER].push_back(-18.291615317790175);
  parameter_mix[SER].push_back(-12585.074732466266);
  parameter_mix[SER].push_back(-2064.454891600786);
  parameter_mix[SER].push_back(15273.905065790364);
  parameter_mix[SER].push_back(-8813.056005263466);
  parameter_mix[SER].push_back(1404.9812302289881);

  parameter_mix[THR].push_back(13233.997179639062);
  parameter_mix[THR].push_back(-39.40454157416847);
  parameter_mix[THR].push_back(-21430.58717233547);
  parameter_mix[THR].push_back(-4566.332853710876);
  parameter_mix[THR].push_back(31717.497780073558);
  parameter_mix[THR].push_back(-20299.614304281313);
  parameter_mix[THR].push_back(3837.207224537505);

  parameter_mix[VAL].push_back(15135.438016299158);
  parameter_mix[VAL].push_back(-51.415141550353205);
  parameter_mix[VAL].push_back(-25859.078442379723);
  parameter_mix[VAL].push_back(-6007.697291593915);
  parameter_mix[VAL].push_back(40997.969600345634);
  parameter_mix[VAL].push_back(-27036.257386814148);
  parameter_mix[VAL].push_back(5328.922363811635);

  parameter_mix[ALA].push_back(6586.942863819189);
  parameter_mix[ALA].push_back(-10.96713559950907);
  parameter_mix[ALA].push_back(-8758.836131731925);
  parameter_mix[ALA].push_back(-1223.1723720922605);
  parameter_mix[ALA].push_back(9475.182453543037);
  parameter_mix[ALA].push_back(-5124.611191433804);
  parameter_mix[ALA].push_back(721.7625962949869);

  parameter_mix[GLY].push_back(3596.0718542192762);
  parameter_mix[GLY].push_back(-4.079285964028122);
  parameter_mix[GLY].push_back(-4089.4217504382686);
  parameter_mix[GLY].push_back(-450.9650932154967);
  parameter_mix[GLY].push_back(3737.026778223427);
  parameter_mix[GLY].push_back(-1862.9856575810572);
  parameter_mix[GLY].push_back(222.97288276257262);

  parameter_mix[HIS].push_back(21779.124723299232);
  parameter_mix[HIS].push_back(-131.4603421188538);
  parameter_mix[HIS].push_back(-49068.74667421681);
  parameter_mix[HIS].push_back(-18685.909496392127);
  parameter_mix[HIS].push_back(127724.60792384286);
  parameter_mix[HIS].push_back(-107419.22159440348);
  parameter_mix[HIS].push_back(28577.228634530744);

  parameter_vac[TRP].push_back(9599.949107368187);
  parameter_vac[TRP].push_back(-66.35331786175249);
  parameter_vac[TRP].push_back(-26311.640290970638);
  parameter_vac[TRP].push_back(-11577.314600529338);
  parameter_vac[TRP].push_back(85847.52554160352);
  parameter_vac[TRP].push_back(-79417.17065742958);
  parameter_vac[TRP].push_back(23090.348430572863);

  parameter_vac[TYR].push_back(7393.553846412945);
  parameter_vac[TYR].push_back(-27.51954035778316);
  parameter_vac[TYR].push_back(-20329.10485615286);
  parameter_vac[TYR].push_back(-7444.276340508767);
  parameter_vac[TYR].push_back(66343.22315132803);
  parameter_vac[TYR].push_back(-64470.58721639446);
  parameter_vac[TYR].push_back(19614.63563898146);

  parameter_vac[PHE].push_back(6081.874997705279);
  parameter_vac[PHE].push_back(-40.474695969500104);
  parameter_vac[PHE].push_back(-14354.627390498901);
  parameter_vac[PHE].push_back(-6156.69750315959);
  parameter_vac[PHE].push_back(42580.84239395237);
  parameter_vac[PHE].push_back(-37704.09749809582);
  parameter_vac[PHE].push_back(10543.005717478625);

  parameter_vac[HIP].push_back(5325.791987063724);
  parameter_vac[HIP].push_back(-35.512112257530156);
  parameter_vac[HIP].push_back(-11488.443296477566);
  parameter_vac[HIP].push_back(-4916.724935318093);
  parameter_vac[HIP].push_back(32134.338675979947);
  parameter_vac[HIP].push_back(-27388.387595464188);
  parameter_vac[HIP].push_back(7359.899986748838);

  parameter_vac[ARG].push_back(7220.306892248294);
  parameter_vac[ARG].push_back(-20.65912886190997);
  parameter_vac[ARG].push_back(-22700.70129646048);
  parameter_vac[ARG].push_back(-8696.901551172636);
  parameter_vac[ARG].push_back(83641.36257312517);
  parameter_vac[ARG].push_back(-85237.33676336925);
  parameter_vac[ARG].push_back(26899.162688310953);

  parameter_vac[LYS].push_back(5038.613120729022);
  parameter_vac[LYS].push_back(-34.08366887546492);
  parameter_vac[LYS].push_back(-12812.921120433106);
  parameter_vac[LYS].push_back(-5843.761329082788);
  parameter_vac[LYS].push_back(42419.08427856609);
  parameter_vac[LYS].push_back(-39460.49038159249);
  parameter_vac[LYS].push_back(11542.320830663035);

  parameter_vac[CYS].push_back(2915.0458981763995);
  parameter_vac[CYS].push_back(-5.380571839804511);
  parameter_vac[CYS].push_back(-3865.366285883624);
  parameter_vac[CYS].push_back(-602.3275271136284);
  parameter_vac[CYS].push_back(4524.133086072617);
  parameter_vac[CYS].push_back(-2537.87137720241);
  parameter_vac[CYS].push_back(381.52870758240016);

  parameter_vac[CYX].push_back(2808.068549348085);
  parameter_vac[CYX].push_back(-7.372260063948);
  parameter_vac[CYX].push_back(-4017.492317531980);
  parameter_vac[CYX].push_back(-840.151815748375);
  parameter_vac[CYX].push_back(5800.074437790741);
  parameter_vac[CYX].push_back(-3637.868824045027);
  parameter_vac[CYX].push_back(659.528934122407);

  parameter_vac[ASP].push_back(3479.750728224898);
  parameter_vac[ASP].push_back(-10.33897891836596);
  parameter_vac[ASP].push_back(-5382.628188436401);
  parameter_vac[ASP].push_back(-1183.8060939576694);
  parameter_vac[ASP].push_back(8100.082378727997);
  parameter_vac[ASP].push_back(-5162.630696148773);
  parameter_vac[ASP].push_back(958.993022379732);

  parameter_vac[GLU].push_back(4487.461543955491);
  parameter_vac[GLU].push_back(-26.671865751817936);
  parameter_vac[GLU].push_back(-8829.738168429001);
  parameter_vac[GLU].push_back(-3297.668395415257);
  parameter_vac[GLU].push_back(20686.457747123466);
  parameter_vac[GLU].push_back(-16030.814134196151);
  parameter_vac[GLU].push_back(3858.4457728083275);

  parameter_vac[ILE].push_back(3842.5968201937776);
  parameter_vac[ILE].push_back(-13.848165050578492);
  parameter_vac[ILE].push_back(-6480.062699452774);
  parameter_vac[ILE].push_back(-1636.3888925440413);
  parameter_vac[ILE].push_back(10967.333210698738);
  parameter_vac[ILE].push_back(-7483.704914714421);
  parameter_vac[ILE].push_back(1548.5696047594895);

  parameter_vac[LEU].push_back(3842.5968201937785);
  parameter_vac[LEU].push_back(-16.2745108270949);
  parameter_vac[LEU].push_back(-6807.110269770606);
  parameter_vac[LEU].push_back(-1926.6303434106014);
  parameter_vac[LEU].push_back(12577.952756390941);
  parameter_vac[LEU].push_back(-8829.40489330961);
  parameter_vac[LEU].push_back(1882.919316016889);

  parameter_vac[MET].push_back(4898.512892967389);
  parameter_vac[MET].push_back(-30.588244886468207);
  parameter_vac[MET].push_back(-10159.093665859045);
  parameter_vac[MET].push_back(-4025.0261508449653);
  parameter_vac[MET].push_back(26007.394369425827);
  parameter_vac[MET].push_back(-21199.218680206573);
  parameter_vac[MET].push_back(5423.004225853842);

  parameter_vac[ASN].push_back(3598.1423998115492);
  parameter_vac[ASN].push_back(-10.357995638888545);
  parameter_vac[ASN].push_back(-5565.603011562138);
  parameter_vac[ASN].push_back(-1190.3294930971967);
  parameter_vac[ASN].push_back(8227.920711951123);
  parameter_vac[ASN].push_back(-5222.61541118056);
  parameter_vac[ASN].push_back(968.593406702772);

  parameter_vac[PRO].push_back(2702.925890807494);
  parameter_vac[PRO].push_back(-4.11690159421177);
  parameter_vac[PRO].push_back(-3395.325331307625);
  parameter_vac[PRO].push_back(-458.95242128002894);
  parameter_vac[PRO].push_back(3584.3640448715823);
  parameter_vac[PRO].push_back(-1921.4140769384692);
  parameter_vac[PRO].push_back(267.08577848319516);

  parameter_vac[GLN].push_back(4621.773132292556);
  parameter_vac[GLN].push_back(-29.511778489038818);
  parameter_vac[GLN].push_back(-9486.077450010192);
  parameter_vac[GLN].push_back(-3768.5756897489828);
  parameter_vac[GLN].push_back(23828.07111260487);
  parameter_vac[GLN].push_back(-19110.205836780202);
  parameter_vac[GLN].push_back(4791.718204894083);

  parameter_vac[SER].push_back(2115.1504654043965);
  parameter_vac[SER].push_back(-2.4158378234251234);
  parameter_vac[SER].push_back(-2488.1131972903822);
  parameter_vac[SER].push_back(-263.64072945387693);
  parameter_vac[SER].push_back(2251.376687850687);
  parameter_vac[SER].push_back(-1066.0790768852626);
  parameter_vac[SER].push_back(105.5155397911316);

  parameter_vac[THR].push_back(2914.9061707158835);
  parameter_vac[THR].push_back(-5.032844592364407);
  parameter_vac[THR].push_back(-3903.2546253886653);
  parameter_vac[THR].push_back(-559.4734271244915);
  parameter_vac[THR].push_back(4315.044828297787);
  parameter_vac[THR].push_back(-2331.211908177365);
  parameter_vac[THR].push_back(323.76849758109853);

  parameter_vac[VAL].push_back(2914.8744247581953);
  parameter_vac[VAL].push_back(-5.847217106105881);
  parameter_vac[VAL].push_back(-4096.370479502377);
  parameter_vac[VAL].push_back(-655.2917606620404);
  parameter_vac[VAL].push_back(4888.77261250007);
  parameter_vac[VAL].push_back(-2765.7552774385167);
  parameter_vac[VAL].push_back(421.9081598033515);

  parameter_vac[ALA].push_back(1443.3438146824446);
  parameter_vac[ALA].push_back(-1.1234573178567506);
  parameter_vac[ALA].push_back(-1492.4547663363514);
  parameter_vac[ALA].push_back(-121.47935619968672);
  parameter_vac[ALA].push_back(1139.689871538201);
  parameter_vac[ALA].push_back(-483.8336547914466);
  parameter_vac[ALA].push_back(32.48231950404626);

  parameter_vac[GLY].push_back(899.5356000422925);
  parameter_vac[GLY].push_back(-0.5200880084066986);
  parameter_vac[GLY].push_back(-787.5892053280859);
  parameter_vac[GLY].push_back(-56.07596224884467);
  parameter_vac[GLY].push_back(546.4212287680981);
  parameter_vac[GLY].push_back(-222.2667666932616);
  parameter_vac[GLY].push_back(12.474587265791476);

  parameter_vac[HIS].push_back(5180.842705000207);
  parameter_vac[HIS].push_back(-29.578973475252766);
  parameter_vac[HIS].push_back(-10323.417251934066);
  parameter_vac[HIS].push_back(-3788.977215582307);
  parameter_vac[HIS].push_back(24427.720792289427);
  parameter_vac[HIS].push_back(-19307.35836837878);
  parameter_vac[HIS].push_back(4780.831414992477);

  // NUCLEIC ACIDS
  parameter_solv[BB_PO3].push_back(1464.5929001192262);
  parameter_solv[BB_PO3].push_back(-2.316908934494931);
  parameter_solv[BB_PO3].push_back(-1882.4844584696532);
  parameter_solv[BB_PO3].push_back(-258.8660305554736);
  parameter_solv[BB_PO3].push_back(2007.5216385943972);
  parameter_solv[BB_PO3].push_back(-1087.6423562424877);
  parameter_solv[BB_PO3].push_back(154.89236486681165);

  parameter_solv[BB_PO2].push_back(575.5201001192197);
  parameter_solv[BB_PO2].push_back(-0.6126595489733868);
  parameter_solv[BB_PO2].push_back(-623.3371092254897);
  parameter_solv[BB_PO2].push_back(-68.05795957022156);
  parameter_solv[BB_PO2].push_back(561.8052621243662);
  parameter_solv[BB_PO2].push_back(-283.39573309540344);
  parameter_solv[BB_PO2].push_back(35.55001698010027);

  parameter_solv[BB_DNA].push_back(21211.009600118316);
  parameter_solv[BB_DNA].push_back(-90.18805990529991);
  parameter_solv[BB_DNA].push_back(-39731.1337351215);
  parameter_solv[BB_DNA].push_back(-10920.373563712878);
  parameter_solv[BB_DNA].push_back(72882.21702424977);
  parameter_solv[BB_DNA].push_back(-51747.487078112754);
  parameter_solv[BB_DNA].push_back(11308.67842901876);

  parameter_solv[BB_DNA_5].push_back(22737.624100119025);
  parameter_solv[BB_DNA_5].push_back(-102.72714886664163);
  parameter_solv[BB_DNA_5].push_back(-43685.329677789705);
  parameter_solv[BB_DNA_5].push_back(-12564.259374093454);
  parameter_solv[BB_DNA_5].push_back(83454.87540484876);
  parameter_solv[BB_DNA_5].push_back(-60367.15652138888);
  parameter_solv[BB_DNA_5].push_back(13507.33372986899);

  parameter_solv[BB_DNA_3].push_back(22737.62410011902);
  parameter_solv[BB_DNA_3].push_back(-101.57816684452263);
  parameter_solv[BB_DNA_3].push_back(-43488.53670557616);
  parameter_solv[BB_DNA_3].push_back(-12345.056184958417);
  parameter_solv[BB_DNA_3].push_back(81963.5236411489);
  parameter_solv[BB_DNA_3].push_back(-58791.59443618196);
  parameter_solv[BB_DNA_3].push_back(13003.199362335576);

  parameter_solv[BB_RNA].push_back(23953.752900120977);
  parameter_solv[BB_RNA].push_back(-117.35779348824401);
  parameter_solv[BB_RNA].push_back(-47644.41735332837);
  parameter_solv[BB_RNA].push_back(-14641.556643789863);
  parameter_solv[BB_RNA].push_back(96893.48627050382);
  parameter_solv[BB_RNA].push_back(-72249.62534169314);
  parameter_solv[BB_RNA].push_back(16792.05552105538);

  parameter_solv[BB_RNA_5].push_back(25574.406400119024);
  parameter_solv[BB_RNA_5].push_back(-131.99642772933734);
  parameter_solv[BB_RNA_5].push_back(-52136.51404531249);
  parameter_solv[BB_RNA_5].push_back(-16682.14273917604);
  parameter_solv[BB_RNA_5].push_back(110278.019216394);
  parameter_solv[BB_RNA_5].push_back(-83715.92027818545);
  parameter_solv[BB_RNA_5].push_back(19875.891337706045);

  parameter_solv[BB_RNA_3].push_back(25574.406400119024);
  parameter_solv[BB_RNA_3].push_back(-127.96875237036166);
  parameter_solv[BB_RNA_3].push_back(-51407.183917584385);
  parameter_solv[BB_RNA_3].push_back(-15922.900669975606);
  parameter_solv[BB_RNA_3].push_back(105078.58889106264);
  parameter_solv[BB_RNA_3].push_back(-78289.16276190648);
  parameter_solv[BB_RNA_3].push_back(18156.83214344118);

  parameter_solv[BASE_A].push_back(13282.562500119211);
  parameter_solv[BASE_A].push_back(-76.45124168404048);
  parameter_solv[BASE_A].push_back(-28376.06994108963);
  parameter_solv[BASE_A].push_back(-9972.910773722022);
  parameter_solv[BASE_A].push_back(65873.86341939073);
  parameter_solv[BASE_A].push_back(-52064.33492910885);
  parameter_solv[BASE_A].push_back(12931.608989412513);

  parameter_solv[BASE_C].push_back(10600.76160011891);
  parameter_solv[BASE_C].push_back(-49.1670871249108);
  parameter_solv[BASE_C].push_back(-20239.818742072875);
  parameter_solv[BASE_C].push_back(-6020.278780090207);
  parameter_solv[BASE_C].push_back(39632.13288981881);
  parameter_solv[BASE_C].push_back(-28954.779736165576);
  parameter_solv[BASE_C].push_back(6551.541109526305);

  parameter_solv[BASE_G].push_back(15470.384400119934);
  parameter_solv[BASE_G].push_back(-93.8013620200972);
  parameter_solv[BASE_G].push_back(-36188.29687013545);
  parameter_solv[BASE_G].push_back(-13717.685098209471);
  parameter_solv[BASE_G].push_back(95658.18473657136);
  parameter_solv[BASE_G].push_back(-81262.37811451119);
  parameter_solv[BASE_G].push_back(21841.903930943085);

  parameter_solv[BASE_T].push_back(17210.81610011936);
  parameter_solv[BASE_T].push_back(-93.10189802920208);
  parameter_solv[BASE_T].push_back(-36466.51927689957);
  parameter_solv[BASE_T].push_back(-12425.55615716932);
  parameter_solv[BASE_T].push_back(83847.427808925);
  parameter_solv[BASE_T].push_back(-66735.64997846584);
  parameter_solv[BASE_T].push_back(16757.3463987507);

  parameter_solv[BASE_U].push_back(10909.802500119395);
  parameter_solv[BASE_U].push_back(-46.17712672768298);
  parameter_solv[BASE_U].push_back(-20149.67695512526);
  parameter_solv[BASE_U].push_back(-5590.242961204435);
  parameter_solv[BASE_U].push_back(37169.2740983132);
  parameter_solv[BASE_U].push_back(-26475.631627167604);
  parameter_solv[BASE_U].push_back(5808.201015156168);

  parameter_mix[BB_PO3].push_back(3061.4050527391964);
  parameter_mix[BB_PO3].push_back(-2.022452986623699);
  parameter_mix[BB_PO3].push_back(-2998.2603666141363);
  parameter_mix[BB_PO3].push_back(-218.44256405859076);
  parameter_mix[BB_PO3].push_back(2113.3633950628328);
  parameter_mix[BB_PO3].push_back(-868.4021757095805);
  parameter_mix[BB_PO3].push_back(52.19052064107954);

  parameter_mix[BB_PO2].push_back(1487.2888381188868);
  parameter_mix[BB_PO2].push_back(-0.6155376265599789);
  parameter_mix[BB_PO2].push_back(-1181.5076027691764);
  parameter_mix[BB_PO2].push_back(-66.25027450710594);
  parameter_mix[BB_PO2].push_back(697.0421991965113);
  parameter_mix[BB_PO2].push_back(-261.8559466354847);
  parameter_mix[BB_PO2].push_back(9.974337082362194);

  parameter_mix[BB_DNA].push_back(17766.29474499878);
  parameter_mix[BB_DNA].push_back(-48.97330188566253);
  parameter_mix[BB_DNA].push_back(-28199.563596327207);
  parameter_mix[BB_DNA].push_back(-5623.82085602494);
  parameter_mix[BB_DNA].push_back(39646.22954828498);
  parameter_mix[BB_DNA].push_back(-24658.81157651943);
  parameter_mix[BB_DNA].push_back(4453.73906293146);

  parameter_mix[BB_DNA_5].push_back(18696.09744203927);
  parameter_mix[BB_DNA_5].push_back(-56.29408880833802);
  parameter_mix[BB_DNA_5].push_back(-30486.108599707608);
  parameter_mix[BB_DNA_5].push_back(-6524.195857141158);
  parameter_mix[BB_DNA_5].push_back(45280.80142686446);
  parameter_mix[BB_DNA_5].push_back(-29007.98616567993);
  parameter_mix[BB_DNA_5].push_back(5488.566965501818);

  parameter_mix[BB_DNA_3].push_back(18696.097442039274);
  parameter_mix[BB_DNA_3].push_back(-55.5645003501971);
  parameter_mix[BB_DNA_3].push_back(-30422.262113654506);
  parameter_mix[BB_DNA_3].push_back(-6409.659659309089);
  parameter_mix[BB_DNA_3].push_back(44605.76043515699);
  parameter_mix[BB_DNA_3].push_back(-28295.62152988411);
  parameter_mix[BB_DNA_3].push_back(5262.5765863484);

  parameter_mix[BB_RNA].push_back(21356.177105457366);
  parameter_mix[BB_RNA].push_back(-76.73490647754872);
  parameter_mix[BB_RNA].push_back(-36845.234814782816);
  parameter_mix[BB_RNA].push_back(-9066.559625582728);
  parameter_mix[BB_RNA].push_back(61167.998793390485);
  parameter_mix[BB_RNA].push_back(-41467.23384423218);
  parameter_mix[BB_RNA].push_back(8518.937793863257);

  parameter_mix[BB_RNA_5].push_back(22386.63276427916);
  parameter_mix[BB_RNA_5].push_back(-85.70426456933487);
  parameter_mix[BB_RNA_5].push_back(-39490.50298502025);
  parameter_mix[BB_RNA_5].push_back(-10223.702594972712);
  parameter_mix[BB_RNA_5].push_back(68450.60459618448);
  parameter_mix[BB_RNA_5].push_back(-47409.91098159006);
  parameter_mix[BB_RNA_5].push_back(10031.136138513202);

  parameter_mix[BB_RNA_3].push_back(22386.63276427916);
  parameter_mix[BB_RNA_3].push_back(-81.93760812351479);
  parameter_mix[BB_RNA_3].push_back(-39031.70571520093);
  parameter_mix[BB_RNA_3].push_back(-9666.316086142708);
  parameter_mix[BB_RNA_3].push_back(65120.07128126262);
  parameter_mix[BB_RNA_3].push_back(-44110.13603681317);
  parameter_mix[BB_RNA_3].push_back(9036.76498256983);

  parameter_mix[BASE_A].push_back(15897.31116611889);
  parameter_mix[BASE_A].push_back(-67.86385832953485);
  parameter_mix[BASE_A].push_back(-28851.754660951636);
  parameter_mix[BASE_A].push_back(-8144.431687170413);
  parameter_mix[BASE_A].push_back(53606.39082954489);
  parameter_mix[BASE_A].push_back(-38083.51243782253);
  parameter_mix[BASE_A].push_back(8293.47107993253);

  parameter_mix[BASE_C].push_back(11733.2828871599);
  parameter_mix[BASE_C].push_back(-38.76775400274115);
  parameter_mix[BASE_C].push_back(-19318.84666423464);
  parameter_mix[BASE_C].push_back(-4507.915522704176);
  parameter_mix[BASE_C].push_back(30576.57671286052);
  parameter_mix[BASE_C].push_back(-20132.46696910844);
  parameter_mix[BASE_C].push_back(3947.8727087996162);

  parameter_mix[BASE_G].push_back(19146.612417237808);
  parameter_mix[BASE_G].push_back(-102.37046638004914);
  parameter_mix[BASE_G].push_back(-38718.96478190546);
  parameter_mix[BASE_G].push_back(-13238.106081860074);
  parameter_mix[BASE_G].push_back(87309.07460288722);
  parameter_mix[BASE_G].push_back(-68364.82442984737);
  parameter_mix[BASE_G].push_back(16815.362401369);

  parameter_mix[BASE_T].push_back(17050.440260819163);
  parameter_mix[BASE_T].push_back(-76.33750600643376);
  parameter_mix[BASE_T].push_back(-31849.539096715005);
  parameter_mix[BASE_T].push_back(-9484.498992751434);
  parameter_mix[BASE_T].push_back(62881.895771680494);
  parameter_mix[BASE_T].push_back(-46531.948557759024);
  parameter_mix[BASE_T].push_back(10734.196329884822);

  parameter_mix[BASE_U].push_back(11904.095265219023);
  parameter_mix[BASE_U].push_back(-34.67511054915295);
  parameter_mix[BASE_U].push_back(-18842.275003104005);
  parameter_mix[BASE_U].push_back(-3993.1174764890684);
  parameter_mix[BASE_U].push_back(27663.625106762345);
  parameter_mix[BASE_U].push_back(-17577.387831701515);
  parameter_mix[BASE_U].push_back(3273.183903219142);

  parameter_vac[BB_PO3].push_back(1599.7962466063982);
  parameter_vac[BB_PO3].push_back(-0.2734304923675441);
  parameter_vac[BB_PO3].push_back(-1023.9448073303214);
  parameter_vac[BB_PO3].push_back(-28.78655166266909);
  parameter_vac[BB_PO3].push_back(426.4382937268249);
  parameter_vac[BB_PO3].push_back(-109.57771615730755);
  parameter_vac[BB_PO3].push_back(-10.244595559424265);

  parameter_vac[BB_PO2].push_back(960.8822037291127);
  parameter_vac[BB_PO2].push_back(-0.09312135742159174);
  parameter_vac[BB_PO2].push_back(-469.39643497461844);
  parameter_vac[BB_PO2].push_back(-9.779346709041985);
  parameter_vac[BB_PO2].push_back(162.1581550003227);
  parameter_vac[BB_PO2].push_back(-37.06686233305618);
  parameter_vac[BB_PO2].push_back(-4.695586672655664);

  parameter_vac[BB_DNA].push_back(3720.2522996838984);
  parameter_vac[BB_DNA].push_back(-5.4229642176938);
  parameter_vac[BB_DNA].push_back(-4800.897672711981);
  parameter_vac[BB_DNA].push_back(-597.2274673070993);
  parameter_vac[BB_DNA].push_back(4825.908840953665);
  parameter_vac[BB_DNA].push_back(-2451.397454446564);
  parameter_vac[BB_DNA].push_back(294.93071756645685);

  parameter_vac[BB_DNA_5].push_back(3843.234214262163);
  parameter_vac[BB_DNA_5].push_back(-6.423078416284452);
  parameter_vac[BB_DNA_5].push_back(-5112.1216386963115);
  parameter_vac[BB_DNA_5].push_back(-713.8373583426668);
  parameter_vac[BB_DNA_5].push_back(5547.545382516269);
  parameter_vac[BB_DNA_5].push_back(-2973.5659871174225);
  parameter_vac[BB_DNA_5].push_back(407.2789106630427);

  parameter_vac[BB_DNA_3].push_back(3843.234214262163);
  parameter_vac[BB_DNA_3].push_back(-6.268636561475645);
  parameter_vac[BB_DNA_3].push_back(-5103.192931218086);
  parameter_vac[BB_DNA_3].push_back(-693.8705734390547);
  parameter_vac[BB_DNA_3].push_back(5443.979645097035);
  parameter_vac[BB_DNA_3].push_back(-2871.4337477324893);
  parameter_vac[BB_DNA_3].push_back(377.3062915349831);

  parameter_vac[BB_RNA].push_back(4760.071443398374);
  parameter_vac[BB_RNA].push_back(-11.0990475402486);
  parameter_vac[BB_RNA].push_back(-6934.713566418421);
  parameter_vac[BB_RNA].push_back(-1256.5202524085096);
  parameter_vac[BB_RNA].push_back(9024.962587066922);
  parameter_vac[BB_RNA].push_back(-5386.842667932241);
  parameter_vac[BB_RNA].push_back(907.42947751372);

  parameter_vac[BB_RNA_5].push_back(4899.051406033406);
  parameter_vac[BB_RNA_5].push_back(-12.279240472628238);
  parameter_vac[BB_RNA_5].push_back(-7276.273570813667);
  parameter_vac[BB_RNA_5].push_back(-1400.9520839250868);
  parameter_vac[BB_RNA_5].push_back(9912.215823228895);
  parameter_vac[BB_RNA_5].push_back(-6079.2565270404075);
  parameter_vac[BB_RNA_5].push_back(1073.53428175472);

  parameter_vac[BB_RNA_3].push_back(4899.051406033406);
  parameter_vac[BB_RNA_3].push_back(-11.642804195148194);
  parameter_vac[BB_RNA_3].push_back(-7213.774619570996);
  parameter_vac[BB_RNA_3].push_back(-1317.4463949342964);
  parameter_vac[BB_RNA_3].push_back(9450.928929264686);
  parameter_vac[BB_RNA_3].push_back(-5643.856117200917);
  parameter_vac[BB_RNA_3].push_back(949.4698817407918);

  parameter_vac[BASE_A].push_back(4756.697028810353);
  parameter_vac[BASE_A].push_back(-12.158940746852812);
  parameter_vac[BASE_A].push_back(-7106.473423744205);
  parameter_vac[BASE_A].push_back(-1376.295184173137);
  parameter_vac[BASE_A].push_back(9747.132255557788);
  parameter_vac[BASE_A].push_back(-5900.026637038756);
  parameter_vac[BASE_A].push_back(1004.6226388342955);

  parameter_vac[BASE_C].push_back(3246.698975674651);
  parameter_vac[BASE_C].push_back(-6.125036521218687);
  parameter_vac[BASE_C].push_back(-4280.666521437201);
  parameter_vac[BASE_C].push_back(-684.0183580843932);
  parameter_vac[BASE_C].push_back(5077.062889522692);
  parameter_vac[BASE_C].push_back(-2870.3239317750963);
  parameter_vac[BASE_C].push_back(434.51551177863547);

  parameter_vac[BASE_G].push_back(5924.105658596052);
  parameter_vac[BASE_G].push_back(-23.097956587232552);
  parameter_vac[BASE_G].push_back(-10149.526301285418);
  parameter_vac[BASE_G].push_back(-2752.9166169522528);
  parameter_vac[BASE_G].push_back(18239.32985385683);
  parameter_vac[BASE_G].push_back(-12749.277858800957);
  parameter_vac[BASE_G].push_back(2715.354663787367);

  parameter_vac[BASE_T].push_back(4222.889713694404);
  parameter_vac[BASE_T].push_back(-12.15861456306705);
  parameter_vac[BASE_T].push_back(-6395.50289789404);
  parameter_vac[BASE_T].push_back(-1421.3942549301012);
  parameter_vac[BASE_T].push_back(9757.061008720135);
  parameter_vac[BASE_T].push_back(-6399.630933839126);
  parameter_vac[BASE_T].push_back(1258.9874225605438);

  parameter_vac[BASE_U].push_back(3247.251361465539);
  parameter_vac[BASE_U].push_back(-5.211020853261115);
  parameter_vac[BASE_U].push_back(-4125.165310360279);
  parameter_vac[BASE_U].push_back(-575.1860235473902);
  parameter_vac[BASE_U].push_back(4457.6562621371495);
  parameter_vac[BASE_U].push_back(-2368.7250146912766);
  parameter_vac[BASE_U].push_back(313.23997718445014);

  for(unsigned i=0; i<atoms.size(); ++i) {
    std::string Aname = pdb.getAtomName(atoms[i]);
    std::string Rname = pdb.getResidueName(atoms[i]);
    Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
    if(Rname=="ALA") {
      atoi[residue_atom[i]]=ALA;
    } else if(Rname=="ARG") {
      atoi[residue_atom[i]]=ARG;
    } else if(Rname=="ASN") {
      atoi[residue_atom[i]]=ASN;
    } else if(Rname=="ASP") {
      atoi[residue_atom[i]]=ASP;
    } else if(Rname=="CYS") {
      atoi[residue_atom[i]]=CYS;
    } else if(Rname=="CYX") {
      atoi[residue_atom[i]]=CYX;
    } else if(Rname=="GLN") {
      atoi[residue_atom[i]]=GLN;
    } else if(Rname=="GLU") {
      atoi[residue_atom[i]]=GLU;
    } else if(Rname=="GLY") {
      atoi[residue_atom[i]]=GLY;
    } else if(Rname=="HIS") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HID") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIP") {
      atoi[residue_atom[i]]=HIP;
      // CHARMM NAMING FOR PROTONATION STATES OF HISTIDINE
    } else if(Rname=="HSD") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSP") {
      atoi[residue_atom[i]]=HIP;
    } else if(Rname=="ILE") {
      atoi[residue_atom[i]]=ILE;
    } else if(Rname=="LEU") {
      atoi[residue_atom[i]]=LEU;
    } else if(Rname=="LYS") {
      atoi[residue_atom[i]]=LYS;
    } else if(Rname=="MET") {
      atoi[residue_atom[i]]=MET;
    } else if(Rname=="PHE") {
      atoi[residue_atom[i]]=PHE;
    } else if(Rname=="PRO") {
      atoi[residue_atom[i]]=PRO;
    } else if(Rname=="SER") {
      atoi[residue_atom[i]]=SER;
    } else if(Rname=="THR") {
      atoi[residue_atom[i]]=THR;
    } else if(Rname=="TRP") {
      atoi[residue_atom[i]]=TRP;
    } else if(Rname=="TYR") {
      atoi[residue_atom[i]]=TYR;
    } else if(Rname=="VAL") {
      atoi[residue_atom[i]]=VAL;
    }
    // NUCLEIC ACIDS
    // nucleobases are not automatically populated as an additional check on the health of the PDB.
    // RNA - G
    else if(Rname=="G") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G3
    } else if(Rname=="G3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G5
    } else if(Rname=="G5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - GT
    } else if(Rname=="GT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U
    } else if(Rname=="U") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U3
    } else if(Rname=="U3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U5
    } else if(Rname=="U5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - UT
    } else if(Rname=="UT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A
    } else if(Rname=="A") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A3
    } else if(Rname=="A3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A5
    } else if(Rname=="A5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - AT
    } else if(Rname=="AT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C
    } else if(Rname=="C") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C3
    } else if(Rname=="C3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C5
    } else if(Rname=="C5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - CT
    } else if(Rname=="CT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G
    } else if(Rname=="DG") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G3
    } else if(Rname=="DG3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G5
    } else if(Rname=="DG5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - GT
    } else if(Rname=="DGT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T
    } else if(Rname=="DT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T3
    } else if(Rname=="DT3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T5
    } else if(Rname=="DT5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - TT
    } else if(Rname=="DTT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A
    } else if(Rname=="DA") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A3
    } else if(Rname=="DA3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A5
    } else if(Rname=="DA5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - AT
    } else if(Rname=="DAT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C
    } else if(Rname=="DC") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C3
    } else if(Rname=="DC3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C5
    } else if(Rname=="DC5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - CT
    } else if(Rname=="DCT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
    } else {
      error("Residue not known: "+Rname);
    }
  }
}

void SAXS::getOnebeadparam_sansH(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_H, std::vector<std::vector<long double> > &parameter_mix_H, std::vector<std::vector<long double> > &parameter_solv_H) {
  parameter_solv_H[TRP].push_back(60737.60249988011);
  parameter_solv_H[TRP].push_back(-77.77344118516487);
  parameter_solv_H[TRP].push_back(-205962.80436572764);
  parameter_solv_H[TRP].push_back(-62014.18523271781);
  parameter_solv_H[TRP].push_back(680712.0512548896);
  parameter_solv_H[TRP].push_back(-681337.967312983);
  parameter_solv_H[TRP].push_back(211474.00338118695);

  parameter_solv_H[TYR].push_back(46250.803599880084);
  parameter_solv_H[TYR].push_back(-45.827646837514614);
  parameter_solv_H[TYR].push_back(-143872.94597686914);
  parameter_solv_H[TYR].push_back(-39049.51580628132);
  parameter_solv_H[TYR].push_back(441321.31246635393);
  parameter_solv_H[TYR].push_back(-434477.6826175059);
  parameter_solv_H[TYR].push_back(133179.21673104732);

  parameter_solv_H[PHE].push_back(42407.164900118914);
  parameter_solv_H[PHE].push_back(-159.20054549009086);
  parameter_solv_H[PHE].push_back(-123847.83591352346);
  parameter_solv_H[PHE].push_back(-41797.78884558073);
  parameter_solv_H[PHE].push_back(380283.87543872406);
  parameter_solv_H[PHE].push_back(-361432.81356389285);
  parameter_solv_H[PHE].push_back(107750.69385054041);

  parameter_solv_H[HIP].push_back(24473.473600119047);
  parameter_solv_H[HIP].push_back(-111.6412807124612);
  parameter_solv_H[HIP].push_back(-65826.17293437096);
  parameter_solv_H[HIP].push_back(-23305.902022201855);
  parameter_solv_H[HIP].push_back(194795.09953510275);
  parameter_solv_H[HIP].push_back(-180454.47859494278);
  parameter_solv_H[HIP].push_back(52699.36922660704);

  parameter_solv_H[ARG].push_back(34106.70239988039);
  parameter_solv_H[ARG].push_back(152.74468231268114);
  parameter_solv_H[ARG].push_back(-117086.46040369231);
  parameter_solv_H[ARG].push_back(-19664.37512726012);
  parameter_solv_H[ARG].push_back(364454.3721646173);
  parameter_solv_H[ARG].push_back(-382076.05102190404);
  parameter_solv_H[ARG].push_back(122775.83306003918);

  parameter_solv_H[LYS].push_back(32292.09000011877);
  parameter_solv_H[LYS].push_back(-111.97888350941923);
  parameter_solv_H[LYS].push_back(-91953.05212591762);
  parameter_solv_H[LYS].push_back(-30691.03615444628);
  parameter_solv_H[LYS].push_back(282092.82233263896);
  parameter_solv_H[LYS].push_back(-269503.6095978623);
  parameter_solv_H[LYS].push_back(80777.92760273012);

  parameter_solv_H[CYS].push_back(11352.902500119093);
  parameter_solv_H[CYS].push_back(-45.52255488723637);
  parameter_solv_H[CYS].push_back(-20925.086525675117);
  parameter_solv_H[CYS].push_back(-5662.681335644281);
  parameter_solv_H[CYS].push_back(38559.09602816144);
  parameter_solv_H[CYS].push_back(-27885.22747486708);
  parameter_solv_H[CYS].push_back(6280.148346561226);

  parameter_solv_H[CYX].push_back(10281.960000119348);
  parameter_solv_H[CYX].push_back(-42.315998754511);
  parameter_solv_H[CYX].push_back(-19328.174487327480);
  parameter_solv_H[CYX].push_back(-5523.191775626829);
  parameter_solv_H[CYX].push_back(38185.463172673335);
  parameter_solv_H[CYX].push_back(-28940.174693042034);
  parameter_solv_H[CYX].push_back(6925.390014187676);

  parameter_solv_H[ASP].push_back(13511.73760011933);
  parameter_solv_H[ASP].push_back(-59.92934247210642);
  parameter_solv_H[ASP].push_back(-25849.867077822244);
  parameter_solv_H[ASP].push_back(-7541.679510407563);
  parameter_solv_H[ASP].push_back(50760.93853987092);
  parameter_solv_H[ASP].push_back(-37677.89102528413);
  parameter_solv_H[ASP].push_back(8745.710458140105);

  parameter_solv_H[GLU].push_back(20443.280400119456);
  parameter_solv_H[GLU].push_back(-113.77513581661388);
  parameter_solv_H[GLU].push_back(-45587.79863958479);
  parameter_solv_H[GLU].push_back(-16187.534798976243);
  parameter_solv_H[GLU].push_back(112609.61802147207);
  parameter_solv_H[GLU].push_back(-93362.01894077536);
  parameter_solv_H[GLU].push_back(24519.546829431332);

  parameter_solv_H[ILE].push_back(27858.948100119596);
  parameter_solv_H[ILE].push_back(-159.27394962770595);
  parameter_solv_H[ILE].push_back(-61571.43025249802);
  parameter_solv_H[ILE].push_back(-21324.89659912433);
  parameter_solv_H[ILE].push_back(144070.7880009225);
  parameter_solv_H[ILE].push_back(-115021.84534734003);
  parameter_solv_H[ILE].push_back(28939.093300284097);

  parameter_solv_H[LEU].push_back(27858.948100119596);
  parameter_solv_H[LEU].push_back(-165.61872365361);
  parameter_solv_H[LEU].push_back(-62564.5706162518);
  parameter_solv_H[LEU].push_back(-22465.325666767214);
  parameter_solv_H[LEU].push_back(151616.7844050042);
  parameter_solv_H[LEU].push_back(-122905.60389771541);
  parameter_solv_H[LEU].push_back(31436.66201442061);

  parameter_solv_H[MET].push_back(25609.60090011981);
  parameter_solv_H[MET].push_back(-135.38816369794569);
  parameter_solv_H[MET].push_back(-67771.01548433342);
  parameter_solv_H[MET].push_back(-25228.91756660071);
  parameter_solv_H[MET].push_back(199649.92084565928);
  parameter_solv_H[MET].push_back(-182251.9246506795);
  parameter_solv_H[MET].push_back(52502.876819125115);

  parameter_solv_H[ASN].push_back(14376.010000119095);
  parameter_solv_H[ASN].push_back(-67.65587848183215);
  parameter_solv_H[ASN].push_back(-28302.877059425664);
  parameter_solv_H[ASN].push_back(-8577.444107282141);
  parameter_solv_H[ASN].push_back(57532.88704197217);
  parameter_solv_H[ASN].push_back(-43261.79974462857);
  parameter_solv_H[ASN].push_back(10186.450874679671);

  parameter_solv_H[PRO].push_back(16866.21690011944);
  parameter_solv_H[PRO].push_back(-70.84312112734995);
  parameter_solv_H[PRO].push_back(-31465.8423531932);
  parameter_solv_H[PRO].push_back(-8653.362744540535);
  parameter_solv_H[PRO].push_back(58032.27079924916);
  parameter_solv_H[PRO].push_back(-41521.001733021694);
  parameter_solv_H[PRO].push_back(9184.527523759205);

  parameter_solv_H[GLN].push_back(21503.289600119);
  parameter_solv_H[GLN].push_back(-121.3012777474678);
  parameter_solv_H[GLN].push_back(-50468.58503443957);
  parameter_solv_H[GLN].push_back(-18462.47495329696);
  parameter_solv_H[GLN].push_back(132718.42007501892);
  parameter_solv_H[GLN].push_back(-113787.20224345029);
  parameter_solv_H[GLN].push_back(30920.340813688686);

  parameter_solv_H[SER].push_back(9181.47240011935);
  parameter_solv_H[SER].push_back(-28.775273124392083);
  parameter_solv_H[SER].push_back(-15205.54229808512);
  parameter_solv_H[SER].push_back(-3377.785599913673);
  parameter_solv_H[SER].push_back(23345.562090489493);
  parameter_solv_H[SER].push_back(-15312.699787471944);
  parameter_solv_H[SER].push_back(3013.844610647712);

  parameter_solv_H[THR].push_back(15020.953600119403);
  parameter_solv_H[THR].push_back(-61.909951810375105);
  parameter_solv_H[THR].push_back(-27814.538986050964);
  parameter_solv_H[THR].push_back(-7532.222992514079);
  parameter_solv_H[THR].push_back(50586.29804970814);
  parameter_solv_H[THR].push_back(-35943.85986777198);
  parameter_solv_H[THR].push_back(7880.091610023207);

  parameter_solv_H[VAL].push_back(19647.628900119355);
  parameter_solv_H[VAL].push_back(-89.04968136833762);
  parameter_solv_H[VAL].push_back(-38050.10118919102);
  parameter_solv_H[VAL].push_back(-10921.421066774372);
  parameter_solv_H[VAL].push_back(72774.31277743122);
  parameter_solv_H[VAL].push_back(-52689.05168504517);
  parameter_solv_H[VAL].push_back(11806.48989635518);

  parameter_solv_H[ALA].push_back(7515.156100119273);
  parameter_solv_H[ALA].push_back(-20.226317591188526);
  parameter_solv_H[ALA].push_back(-11761.841775007797);
  parameter_solv_H[ALA].push_back(-2341.4903622033885);
  parameter_solv_H[ALA].push_back(16545.381259883452);
  parameter_solv_H[ALA].push_back(-10397.171546969075);
  parameter_solv_H[ALA].push_back(1921.5253045340198);

  parameter_solv_H[GLY].push_back(3594.002500119159);
  parameter_solv_H[GLY].push_back(-6.910832388009796);
  parameter_solv_H[GLY].push_back(-4937.3542895091905);
  parameter_solv_H[GLY].push_back(-785.4545979203357);
  parameter_solv_H[GLY].push_back(5852.853693316741);
  parameter_solv_H[GLY].push_back(-3391.2920205126734);
  parameter_solv_H[GLY].push_back(552.3278183161507);

  parameter_solv_H[HIS].push_back(22888.664100119073);
  parameter_solv_H[HIS].push_back(-133.86281863999585);
  parameter_solv_H[HIS].push_back(-57533.51412287858);
  parameter_solv_H[HIS].push_back(-21767.300111408193);
  parameter_solv_H[HIS].push_back(161255.15347073504);
  parameter_solv_H[HIS].push_back(-142176.65100718598);
  parameter_solv_H[HIS].push_back(39642.61507384587);

  parameter_mix_H[TRP].push_back(2974.6515001192306);
  parameter_mix_H[TRP].push_back(-18.361939022074825);
  parameter_mix_H[TRP].push_back(-7284.637435770752);
  parameter_mix_H[TRP].push_back(-2945.7969900381895);
  parameter_mix_H[TRP].push_back(21235.01878657283);
  parameter_mix_H[TRP].push_back(-18909.7406035548);
  parameter_mix_H[TRP].push_back(5324.324204245179);

  parameter_mix_H[TYR].push_back(2029.7362801192114);
  parameter_mix_H[TYR].push_back(-6.983186065527884);
  parameter_mix_H[TYR].push_back(-5041.996113037476);
  parameter_mix_H[TYR].push_back(-1744.5213085724158);
  parameter_mix_H[TYR].push_back(15329.420227814338);
  parameter_mix_H[TYR].push_back(-14648.322529889958);
  parameter_mix_H[TYR].push_back(4405.608657083287);

  parameter_mix_H[PHE].push_back(1704.6885401192117);
  parameter_mix_H[PHE].push_back(-10.077274979133408);
  parameter_mix_H[PHE].push_back(-3769.440088334303);
  parameter_mix_H[PHE].push_back(-1574.6255694551546);
  parameter_mix_H[PHE].push_back(10996.32497868798);
  parameter_mix_H[PHE].push_back(-9840.68281283696);
  parameter_mix_H[PHE].push_back(2792.098605716682);

  parameter_mix_H[HIP].push_back(1376.0462401192394);
  parameter_mix_H[HIP].push_back(-8.576320475413144);
  parameter_mix_H[HIP].push_back(-2796.8327726392167);
  parameter_mix_H[HIP].push_back(-1165.0473128576);
  parameter_mix_H[HIP].push_back(7495.063650365717);
  parameter_mix_H[HIP].push_back(-6331.20422098132);
  parameter_mix_H[HIP].push_back(1692.637366093312);

  parameter_mix_H[ARG].push_back(1280.940480119178);
  parameter_mix_H[ARG].push_back(-7.411214928783748);
  parameter_mix_H[ARG].push_back(-3747.6200569785033);
  parameter_mix_H[ARG].push_back(-1766.5282176004569);
  parameter_mix_H[ARG].push_back(14307.817638456267);
  parameter_mix_H[ARG].push_back(-14297.104122885643);
  parameter_mix_H[ARG].push_back(4450.526244207772);

  parameter_mix_H[LYS].push_back(570.7272001192143);
  parameter_mix_H[LYS].push_back(-5.371742288956095);
  parameter_mix_H[LYS].push_back(-1255.9868267793006);
  parameter_mix_H[LYS].push_back(-748.3071074443138);
  parameter_mix_H[LYS].push_back(4534.824932304509);
  parameter_mix_H[LYS].push_back(-4125.307867230812);
  parameter_mix_H[LYS].push_back(1178.781491068295);

  parameter_mix_H[CYS].push_back(410.21750011921864);
  parameter_mix_H[CYS].push_back(-0.7655651758449595);
  parameter_mix_H[CYS].push_back(-523.8897033718782);
  parameter_mix_H[CYS].push_back(-89.88478273744425);
  parameter_mix_H[CYS].push_back(655.3313542467919);
  parameter_mix_H[CYS].push_back(-407.87897719750896);
  parameter_mix_H[CYS].push_back(76.50541508448237);

  parameter_mix_H[CYX].push_back(466.237200119199);
  parameter_mix_H[CYX].push_back(-1.302082362010);
  parameter_mix_H[CYX].push_back(-667.565738985901);
  parameter_mix_H[CYX].push_back(-159.506437978088);
  parameter_mix_H[CYX].push_back(1085.648159448292);
  parameter_mix_H[CYX].push_back(-767.622943338598);
  parameter_mix_H[CYX].push_back(170.157274620163);

  parameter_mix_H[ASP].push_back(893.6531201192147);
  parameter_mix_H[ASP].push_back(-3.0756255172248874);
  parameter_mix_H[ASP].push_back(-1453.1760425275006);
  parameter_mix_H[ASP].push_back(-365.0424824614137);
  parameter_mix_H[ASP].push_back(2443.570600976796);
  parameter_mix_H[ASP].push_back(-1679.8996339740277);
  parameter_mix_H[ASP].push_back(352.33054461512455);

  parameter_mix_H[GLU].push_back(1075.4955601191884);
  parameter_mix_H[GLU].push_back(-6.917429973203965);
  parameter_mix_H[GLU].push_back(-2262.861870389347);
  parameter_mix_H[GLU].push_back(-909.8078514527992);
  parameter_mix_H[GLU].push_back(5841.923857549836);
  parameter_mix_H[GLU].push_back(-4784.620969556751);
  parameter_mix_H[GLU].push_back(1230.873134652953);

  parameter_mix_H[ILE].push_back(466.0127201192081);
  parameter_mix_H[ILE].push_back(-0.9323443258150218);
  parameter_mix_H[ILE].push_back(-576.7178005955719);
  parameter_mix_H[ILE].push_back(-103.03003361062478);
  parameter_mix_H[ILE].push_back(706.4269951176641);
  parameter_mix_H[ILE].push_back(-420.4412859632717);
  parameter_mix_H[ILE].push_back(71.53175726608731);

  parameter_mix_H[LEU].push_back(466.0127201192081);
  parameter_mix_H[LEU].push_back(-1.9793605752606065);
  parameter_mix_H[LEU].push_back(-718.3988478701591);
  parameter_mix_H[LEU].push_back(-227.36409339012113);
  parameter_mix_H[LEU].push_back(1389.2058254917304);
  parameter_mix_H[LEU].push_back(-990.0033118748643);
  parameter_mix_H[LEU].push_back(213.15736815883042);

  parameter_mix_H[MET].push_back(562.9855401192196);
  parameter_mix_H[MET].push_back(-3.7994094933771643);
  parameter_mix_H[MET].push_back(-1139.6331862451661);
  parameter_mix_H[MET].push_back(-516.6313269725724);
  parameter_mix_H[MET].push_back(3268.957245190869);
  parameter_mix_H[MET].push_back(-2809.178864807947);
  parameter_mix_H[MET].push_back(761.4832732100416);

  parameter_mix_H[ASN].push_back(828.7488001191887);
  parameter_mix_H[ASN].push_back(-2.1275493073799625);
  parameter_mix_H[ASN].push_back(-1222.248291388804);
  parameter_mix_H[ASN].push_back(-238.94210659613537);
  parameter_mix_H[ASN].push_back(1660.8322402171973);
  parameter_mix_H[ASN].push_back(-1008.7934996077323);
  parameter_mix_H[ASN].push_back(173.6082238625797);

  parameter_mix_H[PRO].push_back(578.4409801192146);
  parameter_mix_H[PRO].push_back(-0.5379505780909722);
  parameter_mix_H[PRO].push_back(-648.146493857212);
  parameter_mix_H[PRO].push_back(-56.67223895342921);
  parameter_mix_H[PRO].push_back(509.88860586987437);
  parameter_mix_H[PRO].push_back(-214.57871784725265);
  parameter_mix_H[PRO].push_back(11.99659463759968);

  parameter_mix_H[GLN].push_back(989.2334401191976);
  parameter_mix_H[GLN].push_back(-6.307760694331967);
  parameter_mix_H[GLN].push_back(-1971.7067150503622);
  parameter_mix_H[GLN].push_back(-791.333088386235);
  parameter_mix_H[GLN].push_back(4900.009768434847);
  parameter_mix_H[GLN].push_back(-3909.7740976374153);
  parameter_mix_H[GLN].push_back(975.4952613244343);

  parameter_mix_H[SER].push_back(426.39900011920196);
  parameter_mix_H[SER].push_back(-0.42304498358319664);
  parameter_mix_H[SER].push_back(-484.2066027682147);
  parameter_mix_H[SER].push_back(-45.38968988754228);
  parameter_mix_H[SER].push_back(401.3420977115044);
  parameter_mix_H[SER].push_back(-178.0861461692512);
  parameter_mix_H[SER].push_back(13.540349238730284);

  parameter_mix_H[THR].push_back(525.0470401191992);
  parameter_mix_H[THR].push_back(-0.7419102811534484);
  parameter_mix_H[THR].push_back(-652.7134808154495);
  parameter_mix_H[THR].push_back(-80.39481224407903);
  parameter_mix_H[THR].push_back(641.5487902728123);
  parameter_mix_H[THR].push_back(-320.4227667104819);
  parameter_mix_H[THR].push_back(36.03558531183942);

  parameter_mix_H[VAL].push_back(414.6228601192123);
  parameter_mix_H[VAL].push_back(-0.35889335250521337);
  parameter_mix_H[VAL].push_back(-453.11631644097474);
  parameter_mix_H[VAL].push_back(-36.402101097644284);
  parameter_mix_H[VAL].push_back(336.24049431626804);
  parameter_mix_H[VAL].push_back(-127.42235327515239);
  parameter_mix_H[VAL].push_back(0.8013280923923705);

  parameter_mix_H[ALA].push_back(285.21010011920816);
  parameter_mix_H[ALA].push_back(-0.1573012439142169);
  parameter_mix_H[ALA].push_back(-282.8945838800694);
  parameter_mix_H[ALA].push_back(-16.32030056827785);
  parameter_mix_H[ALA].push_back(178.065895049598);
  parameter_mix_H[ALA].push_back(-60.27423229179658);
  parameter_mix_H[ALA].push_back(-1.4695219304131588);

  parameter_mix_H[GLY].push_back(207.18720011920414);
  parameter_mix_H[GLY].push_back(-0.1036587134183235);
  parameter_mix_H[GLY].push_back(-185.70794948240638);
  parameter_mix_H[GLY].push_back(-11.008101039836257);
  parameter_mix_H[GLY].push_back(115.30600405624061);
  parameter_mix_H[GLY].push_back(-42.46629718037158);
  parameter_mix_H[GLY].push_back(0.9238928987070913);

  parameter_mix_H[HIS].push_back(1443.9117601192354);
  parameter_mix_H[HIS].push_back(-7.478618745973115);
  parameter_mix_H[HIS].push_back(-2715.0835155803215);
  parameter_mix_H[HIS].push_back(-918.5243015382779);
  parameter_mix_H[HIS].push_back(5821.6258431396);
  parameter_mix_H[HIS].push_back(-4415.32722209556);
  parameter_mix_H[HIS].push_back(1044.7044029209756);

  parameter_vac_H[TRP].push_back(36.42122511920832);
  parameter_vac_H[TRP].push_back(-0.36925500341767903);
  parameter_vac_H[TRP].push_back(-51.34228792835503);
  parameter_vac_H[TRP].push_back(-34.10021080004831);
  parameter_vac_H[TRP].push_back(132.647034983933);
  parameter_vac_H[TRP].push_back(-82.89152328934257);
  parameter_vac_H[TRP].push_back(13.029994092013231);

  parameter_vac_H[TYR].push_back(22.268961119209557);
  parameter_vac_H[TYR].push_back(-0.1995573892347673);
  parameter_vac_H[TYR].push_back(-36.54202179838511);
  parameter_vac_H[TYR].push_back(-23.820801043096694);
  parameter_vac_H[TYR].push_back(127.46799692275353);
  parameter_vac_H[TYR].push_back(-107.63783234245744);
  parameter_vac_H[TYR].push_back(28.180858902960775);

  parameter_vac_H[PHE].push_back(17.131321119209627);
  parameter_vac_H[PHE].push_back(-0.15766725674246446);
  parameter_vac_H[PHE].push_back(-19.19964432024534);
  parameter_vac_H[PHE].push_back(-12.34326882843138);
  parameter_vac_H[PHE].push_back(38.17216645824474);
  parameter_vac_H[PHE].push_back(-11.245288857407298);
  parameter_vac_H[PHE].push_back(-3.8114731300899343);

  parameter_vac_H[HIP].push_back(19.34240411920875);
  parameter_vac_H[HIP].push_back(-0.13533410292592593);
  parameter_vac_H[HIP].push_back(-25.924121027387276);
  parameter_vac_H[HIP].push_back(-12.36586927492752);
  parameter_vac_H[HIP].push_back(56.75268702111191);
  parameter_vac_H[HIP].push_back(-31.126240293638094);
  parameter_vac_H[HIP].push_back(2.749811579250848);

  parameter_vac_H[ARG].push_back(12.027024119209557);
  parameter_vac_H[ARG].push_back(-0.41927538341868287);
  parameter_vac_H[ARG].push_back(-22.137566939867042);
  parameter_vac_H[ARG].push_back(-43.22615008762667);
  parameter_vac_H[ARG].push_back(165.77466655520323);
  parameter_vac_H[ARG].push_back(-140.68664871425898);
  parameter_vac_H[ARG].push_back(36.67401195170306);

  parameter_vac_H[LYS].push_back(2.5217441192093717);
  parameter_vac_H[LYS].push_back(-0.0032825476242835413);
  parameter_vac_H[LYS].push_back(14.019071697737793);
  parameter_vac_H[LYS].push_back(7.8634074595069245);
  parameter_vac_H[LYS].push_back(-82.44639716451474);
  parameter_vac_H[LYS].push_back(94.32937851921197);
  parameter_vac_H[LYS].push_back(-32.324473823417);

  parameter_vac_H[CYS].push_back(3.705624880856525);
  parameter_vac_H[CYS].push_back(0.005214780840206113);
  parameter_vac_H[CYS].push_back(1.25680902661715);
  parameter_vac_H[CYS].push_back(0.5779209425501814);
  parameter_vac_H[CYS].push_back(-3.716408071089366);
  parameter_vac_H[CYS].push_back(2.3947518943233117);
  parameter_vac_H[CYS].push_back(-0.40204949737133333);

  parameter_vac_H[CYX].push_back(5.285401118868);
  parameter_vac_H[CYX].push_back(-0.006119528779);
  parameter_vac_H[CYX].push_back(-3.091212256902);
  parameter_vac_H[CYX].push_back(-0.679948780910);
  parameter_vac_H[CYX].push_back(4.495837313271);
  parameter_vac_H[CYX].push_back(-2.827133444940);
  parameter_vac_H[CYX].push_back(0.494583310914);

  parameter_vac_H[ASP].push_back(14.776336119209605);
  parameter_vac_H[ASP].push_back(-0.037351220316916435);
  parameter_vac_H[ASP].push_back(-18.556358387626286);
  parameter_vac_H[ASP].push_back(-4.1737354794552886);
  parameter_vac_H[ASP].push_back(28.424721213037405);
  parameter_vac_H[ASP].push_back(-17.51389895324883);
  parameter_vac_H[ASP].push_back(2.9729111724708597);

  parameter_vac_H[GLU].push_back(14.145121119208973);
  parameter_vac_H[GLU].push_back(-0.11468766098213011);
  parameter_vac_H[GLU].push_back(-26.272637652294613);
  parameter_vac_H[GLU].push_back(-13.769758826440151);
  parameter_vac_H[GLU].push_back(80.4575683578491);
  parameter_vac_H[GLU].push_back(-64.19346347075);
  parameter_vac_H[GLU].push_back(15.545440117656236);

  parameter_vac_H[ILE].push_back(1.9488158808808775);
  parameter_vac_H[ILE].push_back(0.05873132133874459);
  parameter_vac_H[ILE].push_back(12.032778845884135);
  parameter_vac_H[ILE].push_back(7.148416980612881);
  parameter_vac_H[ILE].push_back(-41.87377843832961);
  parameter_vac_H[ILE].push_back(33.96120749582283);
  parameter_vac_H[ILE].push_back(-8.362535852631256);

  parameter_vac_H[LEU].push_back(1.9488158808977816);
  parameter_vac_H[LEU].push_back(0.0778305500414777);
  parameter_vac_H[LEU].push_back(12.333370614594);
  parameter_vac_H[LEU].push_back(9.449427967560764);
  parameter_vac_H[LEU].push_back(-52.65457680603262);
  parameter_vac_H[LEU].push_back(44.681877289399615);
  parameter_vac_H[LEU].push_back(-11.460498338671227);

  parameter_vac_H[MET].push_back(3.0940808808117652);
  parameter_vac_H[MET].push_back(0.04903755678213222);
  parameter_vac_H[MET].push_back(8.981927022922049);
  parameter_vac_H[MET].push_back(8.654862771879014);
  parameter_vac_H[MET].push_back(-57.09889409156816);
  parameter_vac_H[MET].push_back(58.87704775164829);
  parameter_vac_H[MET].push_back(-18.60431990258862);

  parameter_vac_H[ASN].push_back(11.943936119209074);
  parameter_vac_H[ASN].push_back(-0.0005000836239497835);
  parameter_vac_H[ASN].push_back(-9.581236453763157);
  parameter_vac_H[ASN].push_back(0.16244025786232308);
  parameter_vac_H[ASN].push_back(2.9276580099749574);
  parameter_vac_H[ASN].push_back(2.133535783835143);
  parameter_vac_H[ASN].push_back(-1.5709968820975018);

  parameter_vac_H[PRO].push_back(4.9595288808229245);
  parameter_vac_H[PRO].push_back(0.017853932680793515);
  parameter_vac_H[PRO].push_back(4.5421559293101605);
  parameter_vac_H[PRO].push_back(2.008455612787203);
  parameter_vac_H[PRO].push_back(-12.444117841318494);
  parameter_vac_H[PRO].push_back(8.511723688836447);
  parameter_vac_H[PRO].push_back(-1.6337543903496765);

  parameter_vac_H[GLN].push_back(11.377129119208574);
  parameter_vac_H[GLN].push_back(-0.0674805307761209);
  parameter_vac_H[GLN].push_back(-16.56692720411458);
  parameter_vac_H[GLN].push_back(-6.477707440126834);
  parameter_vac_H[GLN].push_back(34.78232259512621);
  parameter_vac_H[GLN].push_back(-19.450886909938312);
  parameter_vac_H[GLN].push_back(1.944286925108988);

  parameter_vac_H[SER].push_back(4.95062488096605);
  parameter_vac_H[SER].push_back(0.004676435440506079);
  parameter_vac_H[SER].push_back(-0.1896653085608564);
  parameter_vac_H[SER].push_back(0.5142284931977218);
  parameter_vac_H[SER].push_back(-2.8946087252759893);
  parameter_vac_H[SER].push_back(2.1031239401634836);
  parameter_vac_H[SER].push_back(-0.38226443516361713);

  parameter_vac_H[THR].push_back(4.588163880808971);
  parameter_vac_H[THR].push_back(0.018587905993982613);
  parameter_vac_H[THR].push_back(3.5289861308270214);
  parameter_vac_H[THR].push_back(2.0780583604591567);
  parameter_vac_H[THR].push_back(-12.3802007068414);
  parameter_vac_H[THR].push_back(8.720986674116094);
  parameter_vac_H[THR].push_back(-1.683256475122275);

  parameter_vac_H[VAL].push_back(2.187440880853519);
  parameter_vac_H[VAL].push_back(0.028351524826584255);
  parameter_vac_H[VAL].push_back(8.36584512491955);
  parameter_vac_H[VAL].push_back(3.1686206615123926);
  parameter_vac_H[VAL].push_back(-19.81959917770108);
  parameter_vac_H[VAL].push_back(13.293003038570571);
  parameter_vac_H[VAL].push_back(-2.4595257726774125);

  parameter_vac_H[ALA].push_back(2.7060248808167935);
  parameter_vac_H[ALA].push_back(0.004618897267213416);
  parameter_vac_H[ALA].push_back(2.4990261487383947);
  parameter_vac_H[ALA].push_back(0.49579332664340864);
  parameter_vac_H[ALA].push_back(-3.850400071630347);
  parameter_vac_H[ALA].push_back(1.9501161562030942);
  parameter_vac_H[ALA].push_back(-0.18332582719788362);

  parameter_vac_H[GLY].push_back(2.985983880876256);
  parameter_vac_H[GLY].push_back(0.0005033131808079042);
  parameter_vac_H[GLY].push_back(-0.42250170279962684);
  parameter_vac_H[GLY].push_back(0.05620517453257455);
  parameter_vac_H[GLY].push_back(-0.16801962822020733);
  parameter_vac_H[GLY].push_back(0.23635459648780555);
  parameter_vac_H[GLY].push_back(-0.06585244715658795);

  parameter_vac_H[HIS].push_back(22.77198411920933);
  parameter_vac_H[HIS].push_back(-0.06607491006655417);
  parameter_vac_H[HIS].push_back(-27.277710268717247);
  parameter_vac_H[HIS].push_back(-5.674444390934355);
  parameter_vac_H[HIS].push_back(33.4059567406171);
  parameter_vac_H[HIS].push_back(-11.60826210271092);
  parameter_vac_H[HIS].push_back(-1.7359607560773076);

  // NUCLEIC ACIDS

  // BB_PO2-BB_PO3 H and D parameters are identical as there is no H or D in the bead.
  parameter_solv_H[BB_PO3].push_back(1464.5929001192262);
  parameter_solv_H[BB_PO3].push_back(-2.316908934494931);
  parameter_solv_H[BB_PO3].push_back(-1882.4844584696532);
  parameter_solv_H[BB_PO3].push_back(-258.8660305554736);
  parameter_solv_H[BB_PO3].push_back(2007.5216385943972);
  parameter_solv_H[BB_PO3].push_back(-1087.6423562424877);
  parameter_solv_H[BB_PO3].push_back(154.89236486681165);

  parameter_solv_H[BB_PO2].push_back(575.5201001192197);
  parameter_solv_H[BB_PO2].push_back(-0.6126595489733864);
  parameter_solv_H[BB_PO2].push_back(-623.3371092254899);
  parameter_solv_H[BB_PO2].push_back(-68.05795957022144);
  parameter_solv_H[BB_PO2].push_back(561.8052621243661);
  parameter_solv_H[BB_PO2].push_back(-283.39573309540276);
  parameter_solv_H[BB_PO2].push_back(35.550016980100295);

  parameter_solv_H[BB_DNA].push_back(21211.009600118316);
  parameter_solv_H[BB_DNA].push_back(-90.18805990529991);
  parameter_solv_H[BB_DNA].push_back(-39731.13373512149);
  parameter_solv_H[BB_DNA].push_back(-10920.373563712872);
  parameter_solv_H[BB_DNA].push_back(72882.21702424981);
  parameter_solv_H[BB_DNA].push_back(-51747.487078112776);
  parameter_solv_H[BB_DNA].push_back(11308.678429018755);

  parameter_solv_H[BB_DNA_5].push_back(22737.624100119025);
  parameter_solv_H[BB_DNA_5].push_back(-102.72714886664161);
  parameter_solv_H[BB_DNA_5].push_back(-43685.329677789734);
  parameter_solv_H[BB_DNA_5].push_back(-12564.25937409345);
  parameter_solv_H[BB_DNA_5].push_back(83454.87540484878);
  parameter_solv_H[BB_DNA_5].push_back(-60367.15652138887);
  parameter_solv_H[BB_DNA_5].push_back(13507.333729868991);

  parameter_solv_H[BB_DNA_3].push_back(22737.62410011902);
  parameter_solv_H[BB_DNA_3].push_back(-101.57816684452251);
  parameter_solv_H[BB_DNA_3].push_back(-43488.536705576145);
  parameter_solv_H[BB_DNA_3].push_back(-12345.056184958425);
  parameter_solv_H[BB_DNA_3].push_back(81963.52364114887);
  parameter_solv_H[BB_DNA_3].push_back(-58791.59443618196);
  parameter_solv_H[BB_DNA_3].push_back(13003.199362335583);

  parameter_solv_H[BB_RNA].push_back(23953.752900120977);
  parameter_solv_H[BB_RNA].push_back(-117.35779348824417);
  parameter_solv_H[BB_RNA].push_back(-47644.41735332843);
  parameter_solv_H[BB_RNA].push_back(-14641.556643789861);
  parameter_solv_H[BB_RNA].push_back(96893.48627050371);
  parameter_solv_H[BB_RNA].push_back(-72249.62534169303);
  parameter_solv_H[BB_RNA].push_back(16792.055521055358);

  parameter_solv_H[BB_RNA_5].push_back(25574.406400119024);
  parameter_solv_H[BB_RNA_5].push_back(-131.99642772933734);
  parameter_solv_H[BB_RNA_5].push_back(-52136.51404531251);
  parameter_solv_H[BB_RNA_5].push_back(-16682.14273917604);
  parameter_solv_H[BB_RNA_5].push_back(110278.01921639398);
  parameter_solv_H[BB_RNA_5].push_back(-83715.92027818544);
  parameter_solv_H[BB_RNA_5].push_back(19875.89133770605);

  parameter_solv_H[BB_RNA_3].push_back(25574.406400119027);
  parameter_solv_H[BB_RNA_3].push_back(-127.96875237036166);
  parameter_solv_H[BB_RNA_3].push_back(-51407.18391758439);
  parameter_solv_H[BB_RNA_3].push_back(-15922.900669975606);
  parameter_solv_H[BB_RNA_3].push_back(105078.5888910626);
  parameter_solv_H[BB_RNA_3].push_back(-78289.16276190645);
  parameter_solv_H[BB_RNA_3].push_back(18156.832143441192);

  parameter_solv_H[BASE_A].push_back(13282.562500119211);
  parameter_solv_H[BASE_A].push_back(-76.45124168404048);
  parameter_solv_H[BASE_A].push_back(-28376.06994108963);
  parameter_solv_H[BASE_A].push_back(-9972.910773722022);
  parameter_solv_H[BASE_A].push_back(65873.86341939073);
  parameter_solv_H[BASE_A].push_back(-52064.33492910885);
  parameter_solv_H[BASE_A].push_back(12931.608989412513);

  parameter_solv_H[BASE_C].push_back(10600.76160011891);
  parameter_solv_H[BASE_C].push_back(-49.1670871249108);
  parameter_solv_H[BASE_C].push_back(-20239.818742072875);
  parameter_solv_H[BASE_C].push_back(-6020.278780090207);
  parameter_solv_H[BASE_C].push_back(39632.13288981881);
  parameter_solv_H[BASE_C].push_back(-28954.779736165576);
  parameter_solv_H[BASE_C].push_back(6551.541109526305);

  parameter_solv_H[BASE_G].push_back(15470.384400119934);
  parameter_solv_H[BASE_G].push_back(-93.8013620200972);
  parameter_solv_H[BASE_G].push_back(-36188.29687013545);
  parameter_solv_H[BASE_G].push_back(-13717.685098209471);
  parameter_solv_H[BASE_G].push_back(95658.18473657136);
  parameter_solv_H[BASE_G].push_back(-81262.37811451119);
  parameter_solv_H[BASE_G].push_back(21841.903930943085);

  parameter_solv_H[BASE_T].push_back(17210.81610011936);
  parameter_solv_H[BASE_T].push_back(-93.10189802920208);
  parameter_solv_H[BASE_T].push_back(-36466.51927689957);
  parameter_solv_H[BASE_T].push_back(-12425.55615716932);
  parameter_solv_H[BASE_T].push_back(83847.427808925);
  parameter_solv_H[BASE_T].push_back(-66735.64997846584);
  parameter_solv_H[BASE_T].push_back(16757.3463987507);

  parameter_solv_H[BASE_U].push_back(10909.802500119395);
  parameter_solv_H[BASE_U].push_back(-46.17712672768298);
  parameter_solv_H[BASE_U].push_back(-20149.67695512526);
  parameter_solv_H[BASE_U].push_back(-5590.242961204435);
  parameter_solv_H[BASE_U].push_back(37169.2740983132);
  parameter_solv_H[BASE_U].push_back(-26475.631627167604);
  parameter_solv_H[BASE_U].push_back(5808.201015156168);

  parameter_mix_H[BB_PO3].push_back(143.5890401192106);
  parameter_mix_H[BB_PO3].push_back(-0.0679405156108208);
  parameter_mix_H[BB_PO3].push_back(-131.78648321068806);
  parameter_mix_H[BB_PO3].push_back(-7.222980065241985);
  parameter_mix_H[BB_PO3].push_back(79.67309464590994);
  parameter_mix_H[BB_PO3].push_back(-27.950095608460042);
  parameter_mix_H[BB_PO3].push_back(0.12790403369995257);

  parameter_mix_H[BB_PO2].push_back(80.12660011920252);
  parameter_mix_H[BB_PO2].push_back(-0.0278885551982023);
  parameter_mix_H[BB_PO2].push_back(-60.532194918222984);
  parameter_mix_H[BB_PO2].push_back(-2.976882903409687);
  parameter_mix_H[BB_PO2].push_back(33.30645116638125);
  parameter_mix_H[BB_PO2].push_back(-11.601573219761374);
  parameter_mix_H[BB_PO2].push_back(0.12551046492022422);

  parameter_mix_H[BB_DNA].push_back(712.7621601191935);
  parameter_mix_H[BB_DNA].push_back(-0.3228709821198571);
  parameter_mix_H[BB_DNA].push_back(-784.5118228559945);
  parameter_mix_H[BB_DNA].push_back(-27.196125702249613);
  parameter_mix_H[BB_DNA].push_back(410.0185035102729);
  parameter_mix_H[BB_DNA].push_back(-54.453513369320355);
  parameter_mix_H[BB_DNA].push_back(-44.85506789237683);

  parameter_mix_H[BB_DNA_5].push_back(625.175339965785);
  parameter_mix_H[BB_DNA_5].push_back(0.2691706617748245);
  parameter_mix_H[BB_DNA_5].push_back(-582.8721350420001);
  parameter_mix_H[BB_DNA_5].push_back(46.512408351374326);
  parameter_mix_H[BB_DNA_5].push_back(-58.93886949899108);
  parameter_mix_H[BB_DNA_5].push_back(307.29720336085046);
  parameter_mix_H[BB_DNA_5].push_back(-131.71996309259953);

  parameter_mix_H[BB_DNA_3].push_back(625.1753399401266);
  parameter_mix_H[BB_DNA_3].push_back(0.08763234414546289);
  parameter_mix_H[BB_DNA_3].push_back(-606.8067575087485);
  parameter_mix_H[BB_DNA_3].push_back(20.84427254872218);
  parameter_mix_H[BB_DNA_3].push_back(92.53523123608991);
  parameter_mix_H[BB_DNA_3].push_back(162.04688035654937);
  parameter_mix_H[BB_DNA_3].push_back(-89.13571774638052);

  parameter_mix_H[BB_RNA].push_back(936.9775801191857);
  parameter_mix_H[BB_RNA].push_back(-1.3233686929680253);
  parameter_mix_H[BB_RNA].push_back(-1212.1627155263773);
  parameter_mix_H[BB_RNA].push_back(-141.35324744384351);
  parameter_mix_H[BB_RNA].push_back(1155.8281658363026);
  parameter_mix_H[BB_RNA].push_back(-548.9340055857343);
  parameter_mix_H[BB_RNA].push_back(50.81734777881503);

  parameter_mix_H[BB_RNA_5].push_back(848.5355201165631);
  parameter_mix_H[BB_RNA_5].push_back(-0.49570338490120175);
  parameter_mix_H[BB_RNA_5].push_back(-976.1033073783973);
  parameter_mix_H[BB_RNA_5].push_back(-32.943532187986605);
  parameter_mix_H[BB_RNA_5].push_back(475.66177061923884);
  parameter_mix_H[BB_RNA_5].push_back(17.51955845824258);
  parameter_mix_H[BB_RNA_5].push_back(-96.74451972314769);

  parameter_mix_H[BB_RNA_3].push_back(848.5355201192122);
  parameter_mix_H[BB_RNA_3].push_back(-0.8301109354355396);
  parameter_mix_H[BB_RNA_3].push_back(-1019.9524389785406);
  parameter_mix_H[BB_RNA_3].push_back(-84.1388451424885);
  parameter_mix_H[BB_RNA_3].push_back(787.1277245040931);
  parameter_mix_H[BB_RNA_3].push_back(-294.67807432795627);
  parameter_mix_H[BB_RNA_3].push_back(-1.3214626461251089);

  parameter_mix_H[BASE_A].push_back(1504.9345001191857);
  parameter_mix_H[BASE_A].push_back(-3.5306888452552663);
  parameter_mix_H[BASE_A].push_back(-2234.3933572775572);
  parameter_mix_H[BASE_A].push_back(-380.0255208494757);
  parameter_mix_H[BASE_A].push_back(2726.27802432048);
  parameter_mix_H[BASE_A].push_back(-1490.8825754968443);
  parameter_mix_H[BASE_A].push_back(199.7501110740159);

  parameter_mix_H[BASE_C].push_back(939.8188801192172);
  parameter_mix_H[BASE_C].push_back(-1.489638186262854);
  parameter_mix_H[BASE_C].push_back(-1244.5515798554075);
  parameter_mix_H[BASE_C].push_back(-161.3972705672055);
  parameter_mix_H[BASE_C].push_back(1276.3568466722545);
  parameter_mix_H[BASE_C].push_back(-643.3057776225742);
  parameter_mix_H[BASE_C].push_back(72.75963113826273);

  parameter_mix_H[BASE_G].push_back(1768.434840119199);
  parameter_mix_H[BASE_G].push_back(-6.505347007077434);
  parameter_mix_H[BASE_G].push_back(-2919.3856777898427);
  parameter_mix_H[BASE_G].push_back(-701.2456464463938);
  parameter_mix_H[BASE_G].push_back(4464.594230284102);
  parameter_mix_H[BASE_G].push_back(-2733.138521295608);
  parameter_mix_H[BASE_G].push_back(458.1177706235891);

  parameter_mix_H[BASE_T].push_back(1179.3981001192033);
  parameter_mix_H[BASE_T].push_back(-3.2037849252756527);
  parameter_mix_H[BASE_T].push_back(-1821.255498763799);
  parameter_mix_H[BASE_T].push_back(-371.01993266441303);
  parameter_mix_H[BASE_T].push_back(2604.074226688971);
  parameter_mix_H[BASE_T].push_back(-1648.1965787713084);
  parameter_mix_H[BASE_T].push_back(307.2962186436368);

  parameter_mix_H[BASE_U].push_back(956.3442001192266);
  parameter_mix_H[BASE_U].push_back(-1.724458000760567);
  parameter_mix_H[BASE_U].push_back(-1287.9746970192687);
  parameter_mix_H[BASE_U].push_back(-192.74748379510373);
  parameter_mix_H[BASE_U].push_back(1459.0789258833893);
  parameter_mix_H[BASE_U].push_back(-810.0763075080915);
  parameter_mix_H[BASE_U].push_back(119.81810290248339);

  parameter_vac_H[BB_PO3].push_back(3.519375907888525);
  parameter_vac_H[BB_PO3].push_back(7.742660056553524e-05);
  parameter_vac_H[BB_PO3].push_back(-1.3856562746347367);
  parameter_vac_H[BB_PO3].push_back(0.00821183249969174);
  parameter_vac_H[BB_PO3].push_back(0.21213096729728745);
  parameter_vac_H[BB_PO3].push_back(0.032592855104331325);
  parameter_vac_H[BB_PO3].push_back(-0.02236149538438434);

  parameter_vac_H[BB_PO2].push_back(2.7889001116093275);
  parameter_vac_H[BB_PO2].push_back(-0.00011178884266113128);
  parameter_vac_H[BB_PO2].push_back(-1.1702605818380667);
  parameter_vac_H[BB_PO2].push_back(-0.011278044036819933);
  parameter_vac_H[BB_PO2].push_back(0.3214006584089025);
  parameter_vac_H[BB_PO2].push_back(-0.04097165983591666);
  parameter_vac_H[BB_PO2].push_back(-0.017525098100539722);

  parameter_vac_H[BB_DNA].push_back(5.987809026456476);
  parameter_vac_H[BB_DNA].push_back(9.945454528827912e-05);
  parameter_vac_H[BB_DNA].push_back(-1.1884708569330031);
  parameter_vac_H[BB_DNA].push_back(-0.007457733256362841);
  parameter_vac_H[BB_DNA].push_back(0.05666858781418339);
  parameter_vac_H[BB_DNA].push_back(-0.15158797629971757);
  parameter_vac_H[BB_DNA].push_back(0.11642340861329734);

  parameter_vac_H[BB_DNA_5].push_back(4.297328944539055);
  parameter_vac_H[BB_DNA_5].push_back(0.0014793971885106831);
  parameter_vac_H[BB_DNA_5].push_back(1.3961088365255605);
  parameter_vac_H[BB_DNA_5].push_back(0.08974639858979384);
  parameter_vac_H[BB_DNA_5].push_back(-1.5198099705167643);
  parameter_vac_H[BB_DNA_5].push_back(-0.12127122359433733);
  parameter_vac_H[BB_DNA_5].push_back(0.4134601046223601);

  parameter_vac_H[BB_DNA_3].push_back(4.297328886488132);
  parameter_vac_H[BB_DNA_3].push_back(0.0041802954281271905);
  parameter_vac_H[BB_DNA_3].push_back(1.6065462295705266);
  parameter_vac_H[BB_DNA_3].push_back(0.4399805535688805);
  parameter_vac_H[BB_DNA_3].push_back(-3.3806711791929804);
  parameter_vac_H[BB_DNA_3].push_back(1.6729551563628675);
  parameter_vac_H[BB_DNA_3].push_back(-0.10911063067909885);

  parameter_vac_H[BB_RNA].push_back(9.162728984394093);
  parameter_vac_H[BB_RNA].push_back(0.00019952321584579868);
  parameter_vac_H[BB_RNA].push_back(-4.744748946331966);
  parameter_vac_H[BB_RNA].push_back(0.025106563403946364);
  parameter_vac_H[BB_RNA].push_back(1.2302956694109803);
  parameter_vac_H[BB_RNA].push_back(0.12359062278096915);
  parameter_vac_H[BB_RNA].push_back(-0.1725633367685285);

  parameter_vac_H[BB_RNA_5].push_back(7.038408898671503);
  parameter_vac_H[BB_RNA_5].push_back(0.005106788424920148);
  parameter_vac_H[BB_RNA_5].push_back(-0.8981588221803118);
  parameter_vac_H[BB_RNA_5].push_back(0.4922588155214312);
  parameter_vac_H[BB_RNA_5].push_back(-2.6667853454023644);
  parameter_vac_H[BB_RNA_5].push_back(1.533316567240718);
  parameter_vac_H[BB_RNA_5].push_back(-0.07199604869737707);

  parameter_vac_H[BB_RNA_3].push_back(7.038408892621863);
  parameter_vac_H[BB_RNA_3].push_back(0.002993083907266898);
  parameter_vac_H[BB_RNA_3].push_back(-1.3626596831098492);
  parameter_vac_H[BB_RNA_3].push_back(0.3138856961130113);
  parameter_vac_H[BB_RNA_3].push_back(-1.684185014289391);
  parameter_vac_H[BB_RNA_3].push_back(1.1862168720864616);
  parameter_vac_H[BB_RNA_3].push_back(-0.1443894172417523);

  parameter_vac_H[BASE_A].push_back(42.62784088079008);
  parameter_vac_H[BASE_A].push_back(0.02302908536431516);
  parameter_vac_H[BASE_A].push_back(-33.22707177297222);
  parameter_vac_H[BASE_A].push_back(2.6853748424439834);
  parameter_vac_H[BASE_A].push_back(-1.6632902891624768);
  parameter_vac_H[BASE_A].push_back(11.905766349515268);
  parameter_vac_H[BASE_A].push_back(-4.547083454788805);

  parameter_vac_H[BASE_C].push_back(20.83009588079022);
  parameter_vac_H[BASE_C].push_back(0.017055822321768378);
  parameter_vac_H[BASE_C].push_back(-8.349634734370916);
  parameter_vac_H[BASE_C].push_back(1.9324634367723073);
  parameter_vac_H[BASE_C].push_back(-8.435199734060882);
  parameter_vac_H[BASE_C].push_back(8.272798368731268);
  parameter_vac_H[BASE_C].push_back(-1.986671440757263);

  parameter_vac_H[BASE_G].push_back(50.53788088079374);
  parameter_vac_H[BASE_G].push_back(0.024035597617780367);
  parameter_vac_H[BASE_G].push_back(-47.94916639302998);
  parameter_vac_H[BASE_G].push_back(3.143375731466498);
  parameter_vac_H[BASE_G].push_back(4.297009866708155);
  parameter_vac_H[BASE_G].push_back(15.855448505050578);
  parameter_vac_H[BASE_G].push_back(-7.827484135873966);

  parameter_vac_H[BASE_T].push_back(20.20502488079069);
  parameter_vac_H[BASE_T].push_back(0.033659966153300002);
  parameter_vac_H[BASE_T].push_back(-6.057999187718758);
  parameter_vac_H[BASE_T].push_back(4.146969282504351);
  parameter_vac_H[BASE_T].push_back(-20.664315319574357);
  parameter_vac_H[BASE_T].push_back(19.982178623201648);
  parameter_vac_H[BASE_T].push_back(-5.440921587349456);

  parameter_vac_H[BASE_U].push_back(20.958084119209754);
  parameter_vac_H[BASE_U].push_back(-0.005164660707148803);
  parameter_vac_H[BASE_U].push_back(-14.53831312442302);
  parameter_vac_H[BASE_U].push_back(-0.5276995756310442);
  parameter_vac_H[BASE_U].push_back(7.060900707522138);
  parameter_vac_H[BASE_U].push_back(-1.8988408480951036);
  parameter_vac_H[BASE_U].push_back(-0.215000567681094);

  for(unsigned i=0; i<atoms.size(); ++i) {
    std::string Aname = pdb.getAtomName(atoms[i]);
    std::string Rname = pdb.getResidueName(atoms[i]);
    Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
    if(Rname=="ALA") {
      atoi[residue_atom[i]]=ALA;
    } else if(Rname=="ARG") {
      atoi[residue_atom[i]]=ARG;
    } else if(Rname=="ASN") {
      atoi[residue_atom[i]]=ASN;
    } else if(Rname=="ASP") {
      atoi[residue_atom[i]]=ASP;
    } else if(Rname=="CYS") {
      atoi[residue_atom[i]]=CYS;
    } else if(Rname=="CYX") {
      atoi[residue_atom[i]]=CYX;
    } else if(Rname=="GLN") {
      atoi[residue_atom[i]]=GLN;
    } else if(Rname=="GLU") {
      atoi[residue_atom[i]]=GLU;
    } else if(Rname=="GLY") {
      atoi[residue_atom[i]]=GLY;
    } else if(Rname=="HIS") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HID") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIP") {
      atoi[residue_atom[i]]=HIP;
      // CHARMM NAMING FOR PROTONATION STATES OF HISTIDINE
    } else if(Rname=="HSD") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSP") {
      atoi[residue_atom[i]]=HIP;
    } else if(Rname=="ILE") {
      atoi[residue_atom[i]]=ILE;
    } else if(Rname=="LEU") {
      atoi[residue_atom[i]]=LEU;
    } else if(Rname=="LYS") {
      atoi[residue_atom[i]]=LYS;
    } else if(Rname=="MET") {
      atoi[residue_atom[i]]=MET;
    } else if(Rname=="PHE") {
      atoi[residue_atom[i]]=PHE;
    } else if(Rname=="PRO") {
      atoi[residue_atom[i]]=PRO;
    } else if(Rname=="SER") {
      atoi[residue_atom[i]]=SER;
    } else if(Rname=="THR") {
      atoi[residue_atom[i]]=THR;
    } else if(Rname=="TRP") {
      atoi[residue_atom[i]]=TRP;
    } else if(Rname=="TYR") {
      atoi[residue_atom[i]]=TYR;
    } else if(Rname=="VAL") {
      atoi[residue_atom[i]]=VAL;
    }
    // NUCLEIC ACIDS
    // nucleobases are not automatically populated as an additional check on the health of the PDB.
    // RNA - G
    else if(Rname=="G") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G3
    } else if(Rname=="G3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G5
    } else if(Rname=="G5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - GT
    } else if(Rname=="GT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U
    } else if(Rname=="U") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U3
    } else if(Rname=="U3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U5
    } else if(Rname=="U5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - UT
    } else if(Rname=="UT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A
    } else if(Rname=="A") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A3
    } else if(Rname=="A3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A5
    } else if(Rname=="A5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - AT
    } else if(Rname=="AT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C
    } else if(Rname=="C") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C3
    } else if(Rname=="C3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C5
    } else if(Rname=="C5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - CT
    } else if(Rname=="CT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G
    } else if(Rname=="DG") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G3
    } else if(Rname=="DG3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G5
    } else if(Rname=="DG5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - GT
    } else if(Rname=="DGT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T
    } else if(Rname=="DT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T3
    } else if(Rname=="DT3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T5
    } else if(Rname=="DT5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - TT
    } else if(Rname=="DTT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A
    } else if(Rname=="DA") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A3
    } else if(Rname=="DA3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A5
    } else if(Rname=="DA5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - AT
    } else if(Rname=="DAT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C
    } else if(Rname=="DC") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C3
    } else if(Rname=="DC3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C5
    } else if(Rname=="DC5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - CT
    } else if(Rname=="DCT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
    } else {
      error("Residue not known: "+Rname);
    }
  }
}

void SAXS::getOnebeadparam_sansD(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_D, std::vector<std::vector<long double> > &parameter_mix_D) {
  // parameter_solv is identical in SAXS/SANS_H/SANS_D since it depends exclusively on param_v. For that reason we kept param_solv only in SAXS and SANS_H.
  parameter_mix_D[TRP].push_back(8105.740500119327);
  parameter_mix_D[TRP].push_back(-41.785616935469804);
  parameter_mix_D[TRP].push_back(-25456.92790554363);
  parameter_mix_D[TRP].push_back(-10058.20599969184);
  parameter_mix_D[TRP].push_back(86171.76479108425);
  parameter_mix_D[TRP].push_back(-83227.63139882773);
  parameter_mix_D[TRP].push_back(25121.390436258724);

  parameter_mix_D[TYR].push_back(6059.530560118732);
  parameter_mix_D[TYR].push_back(-24.522695525705736);
  parameter_mix_D[TYR].push_back(-17180.858815360847);
  parameter_mix_D[TYR].push_back(-5990.1358528219325);
  parameter_mix_D[TYR].push_back(52936.46126637543);
  parameter_mix_D[TYR].push_back(-50150.0042622683);
  parameter_mix_D[TYR].push_back(14914.553672440441);

  parameter_mix_D[PHE].push_back(5563.404880119222);
  parameter_mix_D[PHE].push_back(-33.609784645922794);
  parameter_mix_D[PHE].push_back(-14576.935030777448);
  parameter_mix_D[PHE].push_back(-5759.170105553782);
  parameter_mix_D[PHE].push_back(43316.895956549866);
  parameter_mix_D[PHE].push_back(-39106.58694570862);
  parameter_mix_D[PHE].push_back(11143.375742877468);

  parameter_mix_D[HIP].push_back(3981.7108801192553);
  parameter_mix_D[HIP].push_back(-23.788371565946427);
  parameter_mix_D[HIP].push_back(-9471.73953776056);
  parameter_mix_D[HIP].push_back(-3690.3981577198365);
  parameter_mix_D[HIP].push_back(26365.958584217453);
  parameter_mix_D[HIP].push_back(-23067.58974902849);
  parameter_mix_D[HIP].push_back(6390.507451097114);

  parameter_mix_D[ARG].push_back(6279.489359881259);
  parameter_mix_D[ARG].push_back(1.2061878338083583);
  parameter_mix_D[ARG].push_back(-20305.413937989913);
  parameter_mix_D[ARG].push_back(-5621.666335222669);
  parameter_mix_D[ARG].push_back(67341.96785520067);
  parameter_mix_D[ARG].push_back(-68849.15464591733);
  parameter_mix_D[ARG].push_back(21773.0630363882);

  parameter_mix_D[LYS].push_back(5434.487400119193);
  parameter_mix_D[LYS].push_back(-29.32356328987909);
  parameter_mix_D[LYS].push_back(-14363.66155749977);
  parameter_mix_D[LYS].push_back(-5650.383128516514);
  parameter_mix_D[LYS].push_back(44573.73888236887);
  parameter_mix_D[LYS].push_back(-41515.980945300485);
  parameter_mix_D[LYS].push_back(12181.965046747513);

  parameter_mix_D[CYS].push_back(1519.4030001192032);
  parameter_mix_D[CYS].push_back(-3.564386334921097);
  parameter_mix_D[CYS].push_back(-2275.813167459516);
  parameter_mix_D[CYS].push_back(-409.54431591328125);
  parameter_mix_D[CYS].push_back(2969.5412742839258);
  parameter_mix_D[CYS].push_back(-1798.3157146799638);
  parameter_mix_D[CYS].push_back(314.568167888235);

  parameter_mix_D[CYX].push_back(1310.696400119220);
  parameter_mix_D[CYX].push_back(-2.919852579787);
  parameter_mix_D[CYX].push_back(-1902.283026713150);
  parameter_mix_D[CYX].push_back(-340.431267947190);
  parameter_mix_D[CYX].push_back(2480.025274590502);
  parameter_mix_D[CYX].push_back(-1529.188197179144);
  parameter_mix_D[CYX].push_back(278.926068515295);

  parameter_mix_D[ASP].push_back(1861.6998401191709);
  parameter_mix_D[ASP].push_back(-5.349780637260551);
  parameter_mix_D[ASP].push_back(-2960.36741510377);
  parameter_mix_D[ASP].push_back(-621.8270745040523);
  parameter_mix_D[ASP].push_back(4334.798300452934);
  parameter_mix_D[ASP].push_back(-2776.8560521554878);
  parameter_mix_D[ASP].push_back(527.9777182094936);

  parameter_mix_D[GLU].push_back(2861.6017201192253);
  parameter_mix_D[GLU].push_back(-13.146456903921809);
  parameter_mix_D[GLU].push_back(-5393.408563875243);
  parameter_mix_D[GLU].push_back(-1646.460570818364);
  parameter_mix_D[GLU].push_back(10884.544923253858);
  parameter_mix_D[GLU].push_back(-8159.923373048856);
  parameter_mix_D[GLU].push_back(1914.545660397314);

  parameter_mix_D[ILE].push_back(4288.585540119189);
  parameter_mix_D[ILE].push_back(-19.937215352880365);
  parameter_mix_D[ILE].push_back(-8324.540144463375);
  parameter_mix_D[ILE].push_back(-2431.835931316717);
  parameter_mix_D[ILE].push_back(16079.9912986194);
  parameter_mix_D[ILE].push_back(-11637.693060394462);
  parameter_mix_D[ILE].push_back(2600.8258068480495);

  parameter_mix_D[LEU].push_back(4288.585540119186);
  parameter_mix_D[LEU].push_back(-21.50343599461759);
  parameter_mix_D[LEU].push_back(-8479.703435720274);
  parameter_mix_D[LEU].push_back(-2647.8693829269596);
  parameter_mix_D[LEU].push_back(17297.18115838578);
  parameter_mix_D[LEU].push_back(-12826.972408323161);
  parameter_mix_D[LEU].push_back(2953.1262521615645);

  parameter_mix_D[MET].push_back(3561.6276801191552);
  parameter_mix_D[MET].push_back(-22.19323392975885);
  parameter_mix_D[MET].push_back(-8348.33907053846);
  parameter_mix_D[MET].push_back(-3323.053272414289);
  parameter_mix_D[MET].push_back(23153.238909304255);
  parameter_mix_D[MET].push_back(-20091.960440908682);
  parameter_mix_D[MET].push_back(5518.759669687693);

  parameter_mix_D[ASN].push_back(2326.5396001192003);
  parameter_mix_D[ASN].push_back(-8.634908921289112);
  parameter_mix_D[ASN].push_back(-4057.4552636749636);
  parameter_mix_D[ASN].push_back(-1032.743130124821);
  parameter_mix_D[ASN].push_back(6957.141592429445);
  parameter_mix_D[ASN].push_back(-4808.265318722317);
  parameter_mix_D[ASN].push_back(1016.3944815533755);

  parameter_mix_D[PRO].push_back(2471.1663601191985);
  parameter_mix_D[PRO].push_back(-6.360795284260088);
  parameter_mix_D[PRO].push_back(-3825.4533158429153);
  parameter_mix_D[PRO].push_back(-728.7164844824666);
  parameter_mix_D[PRO].push_back(5195.036303827973);
  parameter_mix_D[PRO].push_back(-3183.733716480742);
  parameter_mix_D[PRO].push_back(563.2376162754107);

  parameter_mix_D[GLN].push_back(3431.669280119236);
  parameter_mix_D[GLN].push_back(-19.412747205646166);
  parameter_mix_D[GLN].push_back(-7298.017973002134);
  parameter_mix_D[GLN].push_back(-2659.3014182337706);
  parameter_mix_D[GLN].push_back(17890.76595805173);
  parameter_mix_D[GLN].push_back(-14684.603067192957);
  parameter_mix_D[GLN].push_back(3814.338335151394);

  parameter_mix_D[SER].push_back(1423.885200119192);
  parameter_mix_D[SER].push_back(-2.586428606204385);
  parameter_mix_D[SER].push_back(-1966.7369507188134);
  parameter_mix_D[SER].push_back(-289.17277383434106);
  parameter_mix_D[SER].push_back(2209.478296043199);
  parameter_mix_D[SER].push_back(-1216.1521614944);
  parameter_mix_D[SER].push_back(177.0615931546754);

  parameter_mix_D[THR].push_back(2311.2364801191825);
  parameter_mix_D[THR].push_back(-6.258071321531929);
  parameter_mix_D[THR].push_back(-3656.295629081312);
  parameter_mix_D[THR].push_back(-716.4013890357804);
  parameter_mix_D[THR].push_back(5071.656317108832);
  parameter_mix_D[THR].push_back(-3125.8076789667816);
  parameter_mix_D[THR].push_back(555.9775741081131);

  parameter_mix_D[VAL].push_back(3041.128320119224);
  parameter_mix_D[VAL].push_back(-9.314034190716423);
  parameter_mix_D[VAL].push_back(-5075.684780220629);
  parameter_mix_D[VAL].push_back(-1070.7083380665008);
  parameter_mix_D[VAL].push_back(7455.654515006894);
  parameter_mix_D[VAL].push_back(-4701.19187164774);
  parameter_mix_D[VAL].push_back(863.4906179388547);

  parameter_mix_D[ALA].push_back(1187.65300011922);
  parameter_mix_D[ALA].push_back(-1.7011187932116822);
  parameter_mix_D[ALA].push_back(-1521.0113615359212);
  parameter_mix_D[ALA].push_back(-187.93745840575576);
  parameter_mix_D[ALA].push_back(1514.6745873304449);
  parameter_mix_D[ALA].push_back(-775.3890045113897);
  parameter_mix_D[ALA].push_back(96.41428177656567);

  parameter_mix_D[GLY].push_back(581.6349001192067);
  parameter_mix_D[GLY].push_back(-0.5877833598361395);
  parameter_mix_D[GLY].push_back(-640.0421286186524);
  parameter_mix_D[GLY].push_back(-64.58515074152534);
  parameter_mix_D[GLY].push_back(551.9509853583185);
  parameter_mix_D[GLY].push_back(-264.1522021146006);
  parameter_mix_D[GLY].push_back(28.36986478439301);

  parameter_mix_D[HIS].push_back(3648.812220119277);
  parameter_mix_D[HIS].push_back(-22.703075403555548);
  parameter_mix_D[HIS].push_back(-8260.235189881098);
  parameter_mix_D[HIS].push_back(-3190.3176569039265);
  parameter_mix_D[HIS].push_back(21589.074332364213);
  parameter_mix_D[HIS].push_back(-18108.640157613925);
  parameter_mix_D[HIS].push_back(4801.237639634437);

  parameter_vac_D[TRP].push_back(270.43802511921314);
  parameter_vac_D[TRP].push_back(-2.196022464340772);
  parameter_vac_D[TRP].push_back(-780.9546710244318);
  parameter_vac_D[TRP].push_back(-371.1573508312626);
  parameter_vac_D[TRP].push_back(2668.7678731652445);
  parameter_vac_D[TRP].push_back(-2478.2920954223678);
  parameter_vac_D[TRP].push_back(722.3731624901676);

  parameter_vac_D[TYR].push_back(198.471744119211);
  parameter_vac_D[TYR].push_back(-1.236792846228289);
  parameter_vac_D[TYR].push_back(-508.0448711054671);
  parameter_vac_D[TYR].push_back(-210.55908129481216);
  parameter_vac_D[TYR].push_back(1558.3884734212413);
  parameter_vac_D[TYR].push_back(-1418.36319255665);
  parameter_vac_D[TYR].push_back(407.21567613893296);

  parameter_vac_D[PHE].push_back(182.46606411921402);
  parameter_vac_D[PHE].push_back(-1.2708008333861447);
  parameter_vac_D[PHE].push_back(-424.50905926426054);
  parameter_vac_D[PHE].push_back(-177.97207825696387);
  parameter_vac_D[PHE].push_back(1180.839971941918);
  parameter_vac_D[PHE].push_back(-1004.004765231886);
  parameter_vac_D[PHE].push_back(269.34384064610344);

  parameter_vac_D[HIP].push_back(161.95107611920753);
  parameter_vac_D[HIP].push_back(-0.9661246983835707);
  parameter_vac_D[HIP].push_back(-332.04673226423995);
  parameter_vac_D[HIP].push_back(-125.41755194926544);
  parameter_vac_D[HIP].push_back(808.705672166199);
  parameter_vac_D[HIP].push_back(-648.8340711218191);
  parameter_vac_D[HIP].push_back(163.71251277400307);

  parameter_vac_D[ARG].push_back(289.0340011192071);
  parameter_vac_D[ARG].push_back(-1.4195753436279361);
  parameter_vac_D[ARG].push_back(-836.3864005546434);
  parameter_vac_D[ARG].push_back(-346.7081039129904);
  parameter_vac_D[ARG].push_back(2922.003491580559);
  parameter_vac_D[ARG].push_back(-2864.816533173085);
  parameter_vac_D[ARG].push_back(877.9525045072293);

  parameter_vac_D[LYS].push_back(228.64464111920753);
  parameter_vac_D[LYS].push_back(-1.686580749083617);
  parameter_vac_D[LYS].push_back(-544.8870548339771);
  parameter_vac_D[LYS].push_back(-252.11087773186324);
  parameter_vac_D[LYS].push_back(1693.784850493428);
  parameter_vac_D[LYS].push_back(-1514.2375008160348);
  parameter_vac_D[LYS].push_back(427.0713155512121);

  parameter_vac_D[CYS].push_back(50.836900116324315);
  parameter_vac_D[CYS].push_back(-0.040204572899665315);
  parameter_vac_D[CYS].push_back(-55.592868149339424);
  parameter_vac_D[CYS].push_back(-4.341359624977117);
  parameter_vac_D[CYS].push_back(41.55290573185214);
  parameter_vac_D[CYS].push_back(-17.248208429078456);
  parameter_vac_D[CYS].push_back(1.0736187172140528);

  parameter_vac_D[CYX].push_back(41.770369115535);
  parameter_vac_D[CYX].push_back(-0.019277246931);
  parameter_vac_D[CYX].push_back(-40.006821199463);
  parameter_vac_D[CYX].push_back(-2.056736901533);
  parameter_vac_D[CYX].push_back(23.707430747544);
  parameter_vac_D[CYX].push_back(-8.010813092204);
  parameter_vac_D[CYX].push_back(-0.023482540763);

  parameter_vac_D[ASP].push_back(64.12806411920792);
  parameter_vac_D[ASP].push_back(-0.08245818875074411);
  parameter_vac_D[ASP].push_back(-78.95500211069523);
  parameter_vac_D[ASP].push_back(-9.030157332821238);
  parameter_vac_D[ASP].push_back(74.72033164806712);
  parameter_vac_D[ASP].push_back(-36.71042192737952);
  parameter_vac_D[ASP].push_back(4.0989206257493676);

  parameter_vac_D[GLU].push_back(100.14004911920799);
  parameter_vac_D[GLU].push_back(-0.28685123265362006);
  parameter_vac_D[GLU].push_back(-152.44619103423773);
  parameter_vac_D[GLU].push_back(-32.99432901288321);
  parameter_vac_D[GLU].push_back(225.5853175183811);
  parameter_vac_D[GLU].push_back(-144.8489352831419);
  parameter_vac_D[GLU].push_back(27.49692658880534);

  parameter_vac_D[ILE].push_back(165.04540911921134);
  parameter_vac_D[ILE].push_back(-0.5061553029227089);
  parameter_vac_D[ILE].push_back(-275.1890586090823);
  parameter_vac_D[ILE].push_back(-57.288063177375356);
  parameter_vac_D[ILE].push_back(398.9780357099449);
  parameter_vac_D[ILE].push_back(-245.42678814428692);
  parameter_vac_D[ILE].push_back(42.72941025472001);

  parameter_vac_D[LEU].push_back(165.04540911921134);
  parameter_vac_D[LEU].push_back(-0.580034983510499);
  parameter_vac_D[LEU].push_back(-281.30910057877514);
  parameter_vac_D[LEU].push_back(-66.19427345166183);
  parameter_vac_D[LEU].push_back(445.19214155995115);
  parameter_vac_D[LEU].push_back(-287.0653610399624);
  parameter_vac_D[LEU].push_back(53.86626261066706);

  parameter_vac_D[MET].push_back(123.83238411920684);
  parameter_vac_D[MET].push_back(-0.7698672022751385);
  parameter_vac_D[MET].push_back(-251.2481622173618);
  parameter_vac_D[MET].push_back(-100.67742019193848);
  parameter_vac_D[MET].push_back(641.1563254731632);
  parameter_vac_D[MET].push_back(-524.8742634212379);
  parameter_vac_D[MET].push_back(135.36487813767542);

  parameter_vac_D[ASN].push_back(94.12880411921148);
  parameter_vac_D[ASN].push_back(-0.22986194121078912);
  parameter_vac_D[ASN].push_back(-138.78769705028003);
  parameter_vac_D[ASN].push_back(-25.896846049402594);
  parameter_vac_D[ASN].push_back(184.55609781654326);
  parameter_vac_D[ASN].push_back(-110.14043851975404);
  parameter_vac_D[ASN].push_back(18.388834098004153);

  parameter_vac_D[PRO].push_back(90.51619611920745);
  parameter_vac_D[PRO].push_back(-0.0977238494110807);
  parameter_vac_D[PRO].push_back(-109.43531311067846);
  parameter_vac_D[PRO].push_back(-10.592981104983805);
  parameter_vac_D[PRO].push_back(93.64863466237733);
  parameter_vac_D[PRO].push_back(-42.348197720920865);
  parameter_vac_D[PRO].push_back(3.5854078482704574);

  parameter_vac_D[GLN].push_back(136.91340111920806);
  parameter_vac_D[GLN].push_back(-0.7259026842220699);
  parameter_vac_D[GLN].push_back(-257.0347011897067);
  parameter_vac_D[GLN].push_back(-89.99600255417684);
  parameter_vac_D[GLN].push_back(570.3890595917421);
  parameter_vac_D[GLN].push_back(-438.8977029769549);
  parameter_vac_D[GLN].push_back(105.48846039376491);

  parameter_vac_D[SER].push_back(55.20490011583253);
  parameter_vac_D[SER].push_back(-0.038078030710377984);
  parameter_vac_D[SER].push_back(-58.79085960838952);
  parameter_vac_D[SER].push_back(-4.067364063406562);
  parameter_vac_D[SER].push_back(41.319899403658475);
  parameter_vac_D[SER].push_back(-15.865682241288962);
  parameter_vac_D[SER].push_back(0.5028409006168431);

  parameter_vac_D[THR].push_back(88.90604111920842);
  parameter_vac_D[THR].push_back(-0.11566717587697625);
  parameter_vac_D[THR].push_back(-114.4541243837681);
  parameter_vac_D[THR].push_back(-12.541537413808342);
  parameter_vac_D[THR].push_back(106.4974738790947);
  parameter_vac_D[THR].push_back(-50.15009912825225);
  parameter_vac_D[THR].push_back(4.719349514074467);

  parameter_vac_D[VAL].push_back(117.67910411920792);
  parameter_vac_D[VAL].push_back(-0.18187311248567883);
  parameter_vac_D[VAL].push_back(-162.8697844894754);
  parameter_vac_D[VAL].push_back(-19.769248288711825);
  parameter_vac_D[VAL].push_back(162.59270939168965);
  parameter_vac_D[VAL].push_back(-79.37261506441627);
  parameter_vac_D[VAL].push_back(8.230771959393175);

  parameter_vac_D[ALA].push_back(46.92250011448002);
  parameter_vac_D[ALA].push_back(-0.020339064649444412);
  parameter_vac_D[ALA].push_back(-44.41584945233503);
  parameter_vac_D[ALA].push_back(-2.1483754537886113);
  parameter_vac_D[ALA].push_back(25.713667829058785);
  parameter_vac_D[ALA].push_back(-8.222782061575268);
  parameter_vac_D[ALA].push_back(-0.2521732728817875);

  parameter_vac_D[GLY].push_back(23.532201119209795);
  parameter_vac_D[GLY].push_back(-0.00628609590047614);
  parameter_vac_D[GLY].push_back(-17.28421910139733);
  parameter_vac_D[GLY].push_back(-0.6641226821159686);
  parameter_vac_D[GLY].push_back(8.536119110048007);
  parameter_vac_D[GLY].push_back(-2.5438638688361466);
  parameter_vac_D[GLY].push_back(-0.11165675928832643);

  parameter_vac_D[HIS].push_back(145.41948111920982);
  parameter_vac_D[HIS].push_back(-0.8548328183368781);
  parameter_vac_D[HIS].push_back(-290.8653238004162);
  parameter_vac_D[HIS].push_back(-107.85375269366395);
  parameter_vac_D[HIS].push_back(685.7025818759361);
  parameter_vac_D[HIS].push_back(-538.2592043545858);
  parameter_vac_D[HIS].push_back(132.17357375729733);

  // NUCLEIC ACIDS
  parameter_mix_D[BB_PO3].push_back(223.2671801192072);
  parameter_mix_D[BB_PO3].push_back(-0.14452515213607267);
  parameter_mix_D[BB_PO3].push_back(-219.64134852678032);
  parameter_mix_D[BB_PO3].push_back(-15.527993497328728);
  parameter_mix_D[BB_PO3].push_back(153.27197635784856);
  parameter_mix_D[BB_PO3].push_back(-61.17793915482464);
  parameter_mix_D[BB_PO3].push_back(2.92608540200577);

  parameter_mix_D[BB_PO2].push_back(80.12660011920252);
  parameter_mix_D[BB_PO2].push_back(-0.02788855519820236);
  parameter_mix_D[BB_PO2].push_back(-60.53219491822279);
  parameter_mix_D[BB_PO2].push_back(-2.9768829034096806);
  parameter_mix_D[BB_PO2].push_back(33.30645116638123);
  parameter_mix_D[BB_PO2].push_back(-11.601573219761375);
  parameter_mix_D[BB_PO2].push_back(0.12551046492022438);

  parameter_mix_D[BB_DNA].push_back(2835.3195201193003);
  parameter_mix_D[BB_DNA].push_back(-7.954301723608173);
  parameter_mix_D[BB_DNA].push_back(-4509.325563460958);
  parameter_mix_D[BB_DNA].push_back(-909.1870692311344);
  parameter_mix_D[BB_DNA].push_back(6375.156903893768);
  parameter_mix_D[BB_DNA].push_back(-3956.4787847570715);
  parameter_mix_D[BB_DNA].push_back(708.9872879613656);

  parameter_mix_D[BB_DNA_5].push_back(3136.73358011921);
  parameter_mix_D[BB_DNA_5].push_back(-10.023435855160427);
  parameter_mix_D[BB_DNA_5].push_back(-5208.921666368173);
  parameter_mix_D[BB_DNA_5].push_back(-1160.4403539440214);
  parameter_mix_D[BB_DNA_5].push_back(7962.598421448727);
  parameter_mix_D[BB_DNA_5].push_back(-5149.059857691847);
  parameter_mix_D[BB_DNA_5].push_back(984.5217027570121);

  parameter_mix_D[BB_DNA_3].push_back(3136.73358011921);
  parameter_mix_D[BB_DNA_3].push_back(-9.618834865806274);
  parameter_mix_D[BB_DNA_3].push_back(-5164.249220443828);
  parameter_mix_D[BB_DNA_3].push_back(-1103.2721475326382);
  parameter_mix_D[BB_DNA_3].push_back(7633.46089052312);
  parameter_mix_D[BB_DNA_3].push_back(-4826.171688395644);
  parameter_mix_D[BB_DNA_3].push_back(888.1820863683546);

  parameter_mix_D[BB_RNA].push_back(3192.5955601188807);
  parameter_mix_D[BB_RNA].push_back(-11.475781582628308);
  parameter_mix_D[BB_RNA].push_back(-5486.264576931735);
  parameter_mix_D[BB_RNA].push_back(-1344.2878288415961);
  parameter_mix_D[BB_RNA].push_back(9035.26109892441);
  parameter_mix_D[BB_RNA].push_back(-6068.471909763036);
  parameter_mix_D[BB_RNA].push_back(1226.3696076463866);

  parameter_mix_D[BB_RNA_5].push_back(3512.1630401192215);
  parameter_mix_D[BB_RNA_5].push_back(-14.191020069433975);
  parameter_mix_D[BB_RNA_5].push_back(-6293.687102187508);
  parameter_mix_D[BB_RNA_5].push_back(-1689.3688494490984);
  parameter_mix_D[BB_RNA_5].push_back(11193.448566821942);
  parameter_mix_D[BB_RNA_5].push_back(-7806.9064399949375);
  parameter_mix_D[BB_RNA_5].push_back(1662.4594983069844);

  parameter_mix_D[BB_RNA_3].push_back(3512.1630401192215);
  parameter_mix_D[BB_RNA_3].push_back(-12.978118135595812);
  parameter_mix_D[BB_RNA_3].push_back(-6149.290195451877);
  parameter_mix_D[BB_RNA_3].push_back(-1515.8309761505627);
  parameter_mix_D[BB_RNA_3].push_back(10176.605450440278);
  parameter_mix_D[BB_RNA_3].push_back(-6813.250569884159);
  parameter_mix_D[BB_RNA_3].push_back(1366.823518955858);

  parameter_mix_D[BASE_A].push_back(2464.736500119229);
  parameter_mix_D[BASE_A].push_back(-12.127452038444783);
  parameter_mix_D[BASE_A].push_back(-4710.661256689607);
  parameter_mix_D[BASE_A].push_back(-1462.6964141954452);
  parameter_mix_D[BASE_A].push_back(9451.725575888277);
  parameter_mix_D[BASE_A].push_back(-6883.018479948857);
  parameter_mix_D[BASE_A].push_back(1540.1526599737797);

  parameter_mix_D[BASE_C].push_back(1797.2697601191685);
  parameter_mix_D[BASE_C].push_back(-5.963855532295215);
  parameter_mix_D[BASE_C].push_back(-2955.077717756034);
  parameter_mix_D[BASE_C].push_back(-689.4543508746372);
  parameter_mix_D[BASE_C].push_back(4665.914740532565);
  parameter_mix_D[BASE_C].push_back(-3051.4605913706982);
  parameter_mix_D[BASE_C].push_back(590.2201952719585);

  parameter_mix_D[BASE_G].push_back(2804.271480119049);
  parameter_mix_D[BASE_G].push_back(-16.928072935469974);
  parameter_mix_D[BASE_G].push_back(-5989.82519987899);
  parameter_mix_D[BASE_G].push_back(-2275.490326521775);
  parameter_mix_D[BASE_G].push_back(15007.832401865428);
  parameter_mix_D[BASE_G].push_back(-12287.520690325606);
  parameter_mix_D[BASE_G].push_back(3172.98306575258);

  parameter_mix_D[BASE_T].push_back(2545.0860001192113);
  parameter_mix_D[BASE_T].push_back(-10.975141620541738);
  parameter_mix_D[BASE_T].push_back(-4636.058358764447);
  parameter_mix_D[BASE_T].push_back(-1340.3746388296138);
  parameter_mix_D[BASE_T].push_back(8850.604320505428);
  parameter_mix_D[BASE_T].push_back(-6421.852532013674);
  parameter_mix_D[BASE_T].push_back(1443.371517335904);

  parameter_mix_D[BASE_U].push_back(1608.7389001192062);
  parameter_mix_D[BASE_U].push_back(-3.9816849036181434);
  parameter_mix_D[BASE_U].push_back(-2411.056432130769);
  parameter_mix_D[BASE_U].push_back(-451.8236361945487);
  parameter_mix_D[BASE_U].push_back(3220.4418252803644);
  parameter_mix_D[BASE_U].push_back(-1944.2515577994325);
  parameter_mix_D[BASE_U].push_back(332.9259542628691);

  parameter_vac_D[BB_PO3].push_back(8.508889119209273);
  parameter_vac_D[BB_PO3].push_back(-0.0010408625482164885);
  parameter_vac_D[BB_PO3].push_back(-5.656130990440752);
  parameter_vac_D[BB_PO3].push_back(-0.10748040057053611);
  parameter_vac_D[BB_PO3].push_back(2.1441246977168227);
  parameter_vac_D[BB_PO3].push_back(-0.3967083127147655);
  parameter_vac_D[BB_PO3].push_back(-0.10110003105909898);

  parameter_vac_D[BB_PO2].push_back(2.7889001116093284);
  parameter_vac_D[BB_PO2].push_back(-0.00011178884266113128);
  parameter_vac_D[BB_PO2].push_back(-1.1702605818380654);
  parameter_vac_D[BB_PO2].push_back(-0.011278044036819927);
  parameter_vac_D[BB_PO2].push_back(0.3214006584089024);
  parameter_vac_D[BB_PO2].push_back(-0.04097165983591666);
  parameter_vac_D[BB_PO2].push_back(-0.017525098100539684);

  parameter_vac_D[BB_DNA].push_back(94.75075611920529);
  parameter_vac_D[BB_DNA].push_back(-0.13973533952241124);
  parameter_vac_D[BB_DNA].push_back(-123.45402430039046);
  parameter_vac_D[BB_DNA].push_back(-15.19494522082691);
  parameter_vac_D[BB_DNA].push_back(123.34749914811465);
  parameter_vac_D[BB_DNA].push_back(-61.038507985345504);
  parameter_vac_D[BB_DNA].push_back(6.601587478585944);

  parameter_vac_D[BB_DNA_5].push_back(108.18080111920679);
  parameter_vac_D[BB_DNA_5].push_back(-0.2055953690887981);
  parameter_vac_D[BB_DNA_5].push_back(-150.7924892157235);
  parameter_vac_D[BB_DNA_5].push_back(-22.700459516383198);
  parameter_vac_D[BB_DNA_5].push_back(172.2599851655527);
  parameter_vac_D[BB_DNA_5].push_back(-93.4983124807692);
  parameter_vac_D[BB_DNA_5].push_back(12.867661230942868);

  parameter_vac_D[BB_DNA_3].push_back(108.18080111920537);
  parameter_vac_D[BB_DNA_3].push_back(-0.18263717534168372);
  parameter_vac_D[BB_DNA_3].push_back(-148.5918817744255);
  parameter_vac_D[BB_DNA_3].push_back(-19.90799847398835);
  parameter_vac_D[BB_DNA_3].push_back(157.55184203379557);
  parameter_vac_D[BB_DNA_3].push_back(-80.28471270058103);
  parameter_vac_D[BB_DNA_3].push_back(9.313712500298278);

  parameter_vac_D[BB_RNA].push_back(106.37859611922117);
  parameter_vac_D[BB_RNA].push_back(-0.2380766148121975);
  parameter_vac_D[BB_RNA].push_back(-153.74131338570024);
  parameter_vac_D[BB_RNA].push_back(-26.415436217574932);
  parameter_vac_D[BB_RNA].push_back(191.90585451112776);
  parameter_vac_D[BB_RNA].push_back(-109.61737794316868);
  parameter_vac_D[BB_RNA].push_back(16.663804191332204);

  parameter_vac_D[BB_RNA_5].push_back(120.58236111920618);
  parameter_vac_D[BB_RNA_5].push_back(-0.340258533619014);
  parameter_vac_D[BB_RNA_5].push_back(-186.08333929996334);
  parameter_vac_D[BB_RNA_5].push_back(-38.493337147644795);
  parameter_vac_D[BB_RNA_5].push_back(266.2262415641144);
  parameter_vac_D[BB_RNA_5].push_back(-164.73088478359585);
  parameter_vac_D[BB_RNA_5].push_back(29.07014157680879);

  parameter_vac_D[BB_RNA_3].push_back(120.5823611192099);
  parameter_vac_D[BB_RNA_3].push_back(-0.274146129206928);
  parameter_vac_D[BB_RNA_3].push_back(-179.24499182395388);
  parameter_vac_D[BB_RNA_3].push_back(-30.315729372259426);
  parameter_vac_D[BB_RNA_3].push_back(222.2645581367648);
  parameter_vac_D[BB_RNA_3].push_back(-125.13581171514033);
  parameter_vac_D[BB_RNA_3].push_back(18.350308154920107);

  parameter_vac_D[BASE_A].push_back(114.34024911921);
  parameter_vac_D[BASE_A].push_back(-0.4136665918383359);
  parameter_vac_D[BASE_A].push_back(-192.33138384655922);
  parameter_vac_D[BASE_A].push_back(-46.74428306691412);
  parameter_vac_D[BASE_A].push_back(312.9511030981905);
  parameter_vac_D[BASE_A].push_back(-199.6349962647333);
  parameter_vac_D[BASE_A].push_back(36.15938693202153);

  parameter_vac_D[BASE_C].push_back(76.17798411921166);
  parameter_vac_D[BASE_C].push_back(-0.1444475142707445);
  parameter_vac_D[BASE_C].push_back(-102.66873668949485);
  parameter_vac_D[BASE_C].push_back(-15.813768367725821);
  parameter_vac_D[BASE_C].push_back(119.63436338715553);
  parameter_vac_D[BASE_C].push_back(-64.22251971660583);
  parameter_vac_D[BASE_C].push_back(8.351952332828862);

  parameter_vac_D[BASE_G].push_back(127.08052911921965);
  parameter_vac_D[BASE_G].push_back(-0.7137457014712297);
  parameter_vac_D[BASE_G].push_back(-239.67686838772786);
  parameter_vac_D[BASE_G].push_back(-88.53661981200943);
  parameter_vac_D[BASE_G].push_back(556.7254485453866);
  parameter_vac_D[BASE_G].push_back(-432.0234649577737);
  parameter_vac_D[BASE_G].push_back(104.407200463848);

  parameter_vac_D[BASE_T].push_back(94.09000011920868);
  parameter_vac_D[BASE_T].push_back(-0.27147149980458524);
  parameter_vac_D[BASE_T].push_back(-143.65649702254174);
  parameter_vac_D[BASE_T].push_back(-30.861235738371892);
  parameter_vac_D[BASE_T].push_back(212.3643014774958);
  parameter_vac_D[BASE_T].push_back(-133.06675501066275);
  parameter_vac_D[BASE_T].push_back(23.951588200687073);

  parameter_vac_D[BASE_U].push_back(59.30540111665979);
  parameter_vac_D[BASE_U].push_back(-0.06146929846591808);
  parameter_vac_D[BASE_U].push_back(-67.43680950211682);
  parameter_vac_D[BASE_U].push_back(-6.625289749170134);
  parameter_vac_D[BASE_U].push_back(58.37012229348065);
  parameter_vac_D[BASE_U].push_back(-26.23044613101723);
  parameter_vac_D[BASE_U].push_back(2.061238351422343);

  for(unsigned i=0; i<atoms.size(); ++i) {
    std::string Aname = pdb.getAtomName(atoms[i]);
    std::string Rname = pdb.getResidueName(atoms[i]);
    Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
    if(Rname=="ALA") {
      atoi[residue_atom[i]]=ALA;
    } else if(Rname=="ARG") {
      atoi[residue_atom[i]]=ARG;
    } else if(Rname=="ASN") {
      atoi[residue_atom[i]]=ASN;
    } else if(Rname=="ASP") {
      atoi[residue_atom[i]]=ASP;
    } else if(Rname=="CYS") {
      atoi[residue_atom[i]]=CYS;
    } else if(Rname=="CYX") {
      atoi[residue_atom[i]]=CYX;
    } else if(Rname=="GLN") {
      atoi[residue_atom[i]]=GLN;
    } else if(Rname=="GLU") {
      atoi[residue_atom[i]]=GLU;
    } else if(Rname=="GLY") {
      atoi[residue_atom[i]]=GLY;
    } else if(Rname=="HIS") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HID") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HIP") {
      atoi[residue_atom[i]]=HIP;
      // CHARMM NAMING FOR PROTONATION STATES OF HISTIDINE
    } else if(Rname=="HSD") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSE") {
      atoi[residue_atom[i]]=HIS;
    } else if(Rname=="HSP") {
      atoi[residue_atom[i]]=HIP;
    } else if(Rname=="ILE") {
      atoi[residue_atom[i]]=ILE;
    } else if(Rname=="LEU") {
      atoi[residue_atom[i]]=LEU;
    } else if(Rname=="LYS") {
      atoi[residue_atom[i]]=LYS;
    } else if(Rname=="MET") {
      atoi[residue_atom[i]]=MET;
    } else if(Rname=="PHE") {
      atoi[residue_atom[i]]=PHE;
    } else if(Rname=="PRO") {
      atoi[residue_atom[i]]=PRO;
    } else if(Rname=="SER") {
      atoi[residue_atom[i]]=SER;
    } else if(Rname=="THR") {
      atoi[residue_atom[i]]=THR;
    } else if(Rname=="TRP") {
      atoi[residue_atom[i]]=TRP;
    } else if(Rname=="TYR") {
      atoi[residue_atom[i]]=TYR;
    } else if(Rname=="VAL") {
      atoi[residue_atom[i]]=VAL;
    }
    // NUCLEIC ACIDS
    // nucleobases are not automatically populated as an additional check on the health of the PDB.
    // RNA - G
    else if(Rname=="G") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G3
    } else if(Rname=="G3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G5
    } else if(Rname=="G5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - GT
    } else if(Rname=="GT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U
    } else if(Rname=="U") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U3
    } else if(Rname=="U3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U5
    } else if(Rname=="U5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - UT
    } else if(Rname=="UT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        atoi[residue_atom[i]]=BASE_U;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A
    } else if(Rname=="A") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A3
    } else if(Rname=="A3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A5
    } else if(Rname=="A5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - AT
    } else if(Rname=="AT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C
    } else if(Rname=="C") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C3
    } else if(Rname=="C3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        atoi[residue_atom[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C5
    } else if(Rname=="C5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - CT
    } else if(Rname=="CT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G
    } else if(Rname=="DG") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G3
    } else if(Rname=="DG3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G5
    } else if(Rname=="DG5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - GT
    } else if(Rname=="DGT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        atoi[residue_atom[i]]=BASE_G;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T
    } else if(Rname=="DT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T3
    } else if(Rname=="DT3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T5
    } else if(Rname=="DT5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - TT
    } else if(Rname=="DTT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        atoi[residue_atom[i]]=BASE_T;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A
    } else if(Rname=="DA") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A3
    } else if(Rname=="DA3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A5
    } else if(Rname=="DA5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - AT
    } else if(Rname=="DAT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        atoi[residue_atom[i]]=BASE_A;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C
    } else if(Rname=="DC") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C3
    } else if(Rname=="DC3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        atoi [residue_atom[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C5
    } else if(Rname=="DC5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        atoi[residue_atom[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - CT
    } else if(Rname=="DCT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        atoi [residue_atom[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        atoi[residue_atom[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        atoi[residue_atom[i]]=BASE_C;
      } else {
        error("Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
    } else {
      error("Residue not known: "+Rname);
    }
  }
}

double SAXS::calculateAFF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double rho) {
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

  param_a[H][0] = 0.493002;
  param_b[H][0] = 10.5109;
  param_c[H] = 0.003038;
  param_a[H][1] = 0.322912;
  param_b[H][1] = 26.1257;
  param_v[H] = 5.15;
  param_a[H][2] = 0.140191;
  param_b[H][2] = 3.14236;
  param_a[H][3] = 0.040810;
  param_b[H][3] = 57.7997;
  param_a[H][4] = 0.0;
  param_b[H][4] = 1.0;

  param_a[C][0] = 2.31000;
  param_b[C][0] = 20.8439;
  param_c[C] = 0.215600;
  param_a[C][1] = 1.02000;
  param_b[C][1] = 10.2075;
  param_v[C] = 16.44;
  param_a[C][2] = 1.58860;
  param_b[C][2] = 0.56870;
  param_a[C][3] = 0.86500;
  param_b[C][3] = 51.6512;
  param_a[C][4] = 0.0;
  param_b[C][4] = 1.0;

  param_a[N][0] = 12.2126;
  param_b[N][0] = 0.00570;
  param_c[N] = -11.529;
  param_a[N][1] = 3.13220;
  param_b[N][1] = 9.89330;
  param_v[N] = 2.49;
  param_a[N][2] = 2.01250;
  param_b[N][2] = 28.9975;
  param_a[N][3] = 1.16630;
  param_b[N][3] = 0.58260;
  param_a[N][4] = 0.0;
  param_b[N][4] = 1.0;

  param_a[O][0] = 3.04850;
  param_b[O][0] = 13.2771;
  param_c[O] = 0.250800 ;
  param_a[O][1] = 2.28680;
  param_b[O][1] = 5.70110;
  param_v[O] = 9.13;
  param_a[O][2] = 1.54630;
  param_b[O][2] = 0.32390;
  param_a[O][3] = 0.86700;
  param_b[O][3] = 32.9089;
  param_a[O][4] = 0.0;
  param_b[O][4] = 1.0;

  param_a[P][0] = 6.43450;
  param_b[P][0] = 1.90670;
  param_c[P] = 1.11490;
  param_a[P][1] = 4.17910;
  param_b[P][1] = 27.1570;
  param_v[P] = 5.73;
  param_a[P][2] = 1.78000;
  param_b[P][2] = 0.52600;
  param_a[P][3] = 1.49080;
  param_b[P][3] = 68.1645;
  param_a[P][4] = 0.0;
  param_b[P][4] = 1.0;

  param_a[S][0] = 6.90530;
  param_b[S][0] = 1.46790;
  param_c[S] = 0.866900;
  param_a[S][1] = 5.20340;
  param_b[S][1] = 22.2151;
  param_v[S] = 19.86;
  param_a[S][2] = 1.43790;
  param_b[S][2] = 0.25360;
  param_a[S][3] = 1.58630;
  param_b[S][3] = 56.1720;
  param_a[S][4] = 0.0;
  param_b[S][4] = 1.0;

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);

  double Iq0=0.;
  if( moldat ) {
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
        for(unsigned j=0; j<4; ++j) {
          Iq0 += param_a[index][j];
        }
        Iq0 = Iq0 -rho*param_v[index] + param_c[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  if(absolute) {
    Iq0 = 1;
  }

  return Iq0;
}

double SAXS::calculateAFFsans(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double deuter_conc) {
  std::map<std::string, unsigned> AA_map;
  AA_map["H"] = H;
  AA_map["C"] = C;
  AA_map["N"] = N;
  AA_map["O"] = O;
  AA_map["P"] = P;
  AA_map["S"] = S;

  std::vector<double> param_b;
  std::vector<double> param_v;

  param_b.resize(NTT);
  param_v.resize(NTT);

  param_b[H] = -0.374;
  param_v[H] = 5.15;
  // param_b[D] = 0.667;
  param_b[C] =  0.665;
  param_v[C] = 16.44;
  param_b[N] =  0.94;
  param_v[N] = 2.49;
  param_b[O] =  0.580;
  param_v[O] = 9.13;
  param_b[P] =  0.51;
  param_v[P] = 5.73;
  param_b[S] =  0.28;
  param_v[S] = 19.86;

  double solv_sc_length = 0.1*(param_b[O] + 2.*((1. - deuter_conc) * param_b[H] + deuter_conc * 0.667)); // per water electron (10 electrons)

  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);

  double Iq0=0.;
  if( moldat ) {
    // cycle over the atom types
    for(unsigned i=0; i<NTT; ++i) {
      double volr = std::pow(param_v[i], (2.0/3.0)) /(4. * M_PI);
      // cycle over q
      for(unsigned k=0; k<q_list.size(); ++k) {
        const double q = q_list[k];
        FF_tmp[k][i] = param_b[i];
        // subtract solvation: rho * v_i * EXP( (- v_i^(2/3) / (4pi)) * q^2  ) // since  D in Fraser 1978 is 2*s
        FF_tmp[k][i] -= solv_sc_length*rho*param_v[i]*std::exp(-volr*q*q);
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
        Iq0 += param_b[index]-solv_sc_length*rho*param_v[index];
      } else {
        error("Wrong atom type "+type_s+" from atom name "+name+"\n");
      }
    }
  } else {
    error("MOLINFO DATA not found\n");
  }
  if(absolute) {
    Iq0 = 1;
  }

  return Iq0;
}

std::map<std::string, std::vector<double> > SAXS::setupLCPOparam() {
  std::map<std::string, std::vector<double> > lcpomap;

  // We arbitrarily set OC1/OT1 as the charged oxygen.

  lcpomap["ALA_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["ALA_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["ALA_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["ALA_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_CB"] = { 1.7,  0.77887,  -0.28063,  -0.0012968,  0.00039328};
  lcpomap["ALA_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ALA_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["ALA_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["ALA_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["ASP_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["ASN_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["ARG_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["CYS_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

  lcpomap["CYX_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["CYX_CA"] = { 1.7,  0.23348,  -0.072627,  -0.00020079,  7.967e-05};
  lcpomap["CYX_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["CYX_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYX_CB"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["CYX_SG"] = { 1.9,  0.54581,  -0.19477,  -0.0012873,  0.00029247};
  lcpomap["CYX_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["CYX_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYX_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["CYX_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["CYX_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["GLU_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["GLN_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

  lcpomap["GLY_N"] = { 1.65,  0.41102,  -0.12254,  -7.5448e-05,  0.00011804};
  lcpomap["GLY_CA"] = { 1.7,  0.56482,  -0.19608,  -0.0010219,  0.0002658};
  lcpomap["GLY_C"] = { 1.7,  0.070344,  -0.019015,  -2.2009e-05,  1.6875e-05};
  lcpomap["GLY_O"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLY_OC1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLY_OC2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLY_OT1"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};
  lcpomap["GLY_OT2"] = { 1.6,  0.68563,  -0.1868,  -0.00135573,  0.00023743};
  lcpomap["GLY_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HIS_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HIE_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HSE_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HID_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HSD_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HIP_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["HSP_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["ILE_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["LEU_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["LYS_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["MET_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["PHE_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["PRO_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["SER_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["THR_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["TRP_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["TYR_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

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
  lcpomap["VAL_OXT"] = { 1.6,  0.88857,  -0.33421,  -0.0018683,  0.00049372};

  // nucleic acids - WARNING: ONLY AMBER (OL3-rna/OL15-dna) FORMAT

  lcpomap["A3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A3_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A3_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A3_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["A3_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A3_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A3_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A3_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["A3_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A3_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["A3_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["A5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A5_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A5_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A5_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["A5_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A5_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A5_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A5_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["A5_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A5_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["A5_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["AT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["AT_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["AT_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["AT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["AT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["AT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["AT_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["AT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["AT_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["AT_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["AT_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["AT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["AT_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["AT_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["AT_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["AT_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["AT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["AT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["AT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["AT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["AT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["AT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["AT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["AT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["AT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["AT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["A_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["A_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["A_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["A_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["A_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["A_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["A_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["A_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["A_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["A_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["A_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["C3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C3_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C3_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C3_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["C3_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C3_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["C3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["C3_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["C3_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["C3_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["C5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C5_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C5_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C5_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["C5_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C5_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["C5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["C5_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["C5_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["C5_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["CT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["CT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["CT_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["CT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["CT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["CT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["CT_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["CT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["CT_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["CT_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["CT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["CT_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["CT_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["CT_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["CT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["CT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["CT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["CT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["CT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["CT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["CT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["CT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["CT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["CT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["C_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["C_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["C_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["C_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["C_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["C_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["C_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["C_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["C_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["C_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["C_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["C_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DA3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA3_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA3_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA3_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA3_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA3_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA3_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA3_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DA3_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA3_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DA3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DA5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA5_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA5_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA5_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA5_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA5_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA5_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA5_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DA5_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA5_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DA5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA5_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DAT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DAT_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DAT_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DAT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DAT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DAT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DAT_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DAT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DAT_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DAT_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DAT_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DAT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DAT_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DAT_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DAT_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DAT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DAT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DAT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DAT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DAT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DAT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DAT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DAT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DAT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DAT_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DA_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA_C2"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DA_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DA_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DA_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DA_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA_N6"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DA_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DA_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DA_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DA_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DA_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DA_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DC3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC3_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC3_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC3_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC3_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC3_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DC3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DC3_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DC3_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DC3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DC5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC5_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC5_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC5_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC5_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC5_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DC5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DC5_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DC5_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DC5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DCT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DCT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05};
  lcpomap["DCT_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DCT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DCT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05};
  lcpomap["DCT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DCT_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DCT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DCT_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DCT_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05};
  lcpomap["DCT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DCT_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DCT_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DCT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DCT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DCT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DCT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DCT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DCT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DCT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DCT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DCT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DCT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DC_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DC_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DC_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DC_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DC_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DC_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DC_N4"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DC_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DC_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DC_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DC_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DC_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DG3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG3_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG3_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG3_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG3_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG3_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DG3_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG3_N2"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DG3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG3_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG3_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DG3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG3_O6"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DG3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG3_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DG5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG5_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG5_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG5_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG5_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG5_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DG5_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG5_N2"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DG5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG5_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG5_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DG5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG5_O6"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DG5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG5_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DGT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DGT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DGT_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DGT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DGT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DGT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DGT_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DGT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DGT_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DGT_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DGT_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DGT_N2"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DGT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DGT_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DGT_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DGT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DGT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DGT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DGT_O6"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DGT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DGT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DGT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DGT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DGT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DGT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DGT_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DG_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DG_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DG_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DG_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DG_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG_N2"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["DG_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DG_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DG_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DG_O6"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DG_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DG_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DG_P"]   = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DT3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT3_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT3_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT3_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT3_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DT3_C7"]  = { 1.7,  0.77887, -0.28063, -1.2968e-03, 3.9328e-04 };
  lcpomap["DT3_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DT3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DT3_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT3_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DT5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT5_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT5_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT5_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT5_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DT5_C7"]  = { 1.7,  0.77887, -0.28063, -1.2968e-03, 3.9328e-04 };
  lcpomap["DT5_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DT5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DT5_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT5_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DTT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DTT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DTT_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DTT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DTT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DTT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DTT_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DTT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DTT_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DTT_C7"]  = { 1.7,  0.77887, -0.28063, -1.2968e-03, 3.9328e-04 };
  lcpomap["DTT_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DTT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DTT_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DTT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DTT_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DTT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DTT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DTT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DTT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DTT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DTT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DTT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DTT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DTT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["DT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT_C2'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["DT_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["DT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["DT_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["DT_C7"]  = { 1.7,  0.77887, -0.28063, -1.2968e-03, 3.9328e-04 };
  lcpomap["DT_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["DT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["DT_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["DT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["DT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["DT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["DT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["G3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G3_C2"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G3_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G3_C4"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G3_C5"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["G3_C6"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G3_C8"] = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["G3_N1"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G3_N2"] = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["G3_N3"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G3_N7"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G3_N9"] = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["G3_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G3_O6"] = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["G3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["G5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G5_C2"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G5_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G5_C4"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G5_C5"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["G5_C6"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G5_C8"] = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["G5_N1"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G5_N2"] = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["G5_N3"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G5_N7"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G5_N9"] = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["G5_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G5_O6"] = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["G5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["GT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["GT_C2"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["GT_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["GT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["GT_C4"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["GT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["GT_C5"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["GT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["GT_C6"] = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["GT_C8"] = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["GT_N1"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["GT_N2"] = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["GT_N3"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["GT_N7"] = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["GT_N9"] = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["GT_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["GT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["GT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["GT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["GT_O6"] = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["GT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["GT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["GT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["GT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["GT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["GT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["GT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["G_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["G_C5"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["G_C6"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["G_C8"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["G_N1"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G_N2"]  = { 1.65,  0.73511, -0.22116, -8.9148e-04, 2.523e-04 };
  lcpomap["G_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G_N7"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["G_N9"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["G_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["G_O6"] = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["G_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["G_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["G_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["U3_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U3_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U3_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U3_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U3_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U3_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U3_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U3_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["U3_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U3_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["U3_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["U3_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U3_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U3_O3'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U3_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U3_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U3_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U3_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U3_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U3_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U3_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U3_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U3_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U3_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["U5_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U5_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U5_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U5_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U5_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U5_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U5_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U5_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["U5_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U5_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["U5_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["U5_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U5_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U5_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U5_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U5_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U5_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U5_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U5_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U5_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U5_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U5_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["UT_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["UT_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["UT_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["UT_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["UT_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["UT_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["UT_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["UT_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["UT_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["UT_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["UT_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["UT_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["UT_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["UT_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["UT_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["UT_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["UT_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["UT_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["UT_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["UT_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["UT_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["UT_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["UT_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["UT_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  lcpomap["U_C1'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U_C2"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U_C2'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U_C3'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U_C4"]  = { 1.7,  0.070344, -0.019015, -2.2009e-05, 1.6875e-05 };
  lcpomap["U_C4'"] = { 1.7,  0.23348, -0.072627, -2.0079e-04, 7.967e-05 };
  lcpomap["U_C5"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U_C5'"] = { 1.7,  0.56482, -0.19608, -1.0219e-03, 2.658e-04 };
  lcpomap["U_C6"]  = { 1.7,  0.51245, -0.15966, -1.9781e-04, 1.6392e-04 };
  lcpomap["U_N1"]  = { 1.65,  0.062577, -0.017874, -8.312e-05, 1.9849e-05 };
  lcpomap["U_N3"]  = { 1.65,  0.41102, -0.12254, -7.5448e-05, 1.1804e-04 };
  lcpomap["U_O2"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U_O2'"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U_O3'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U_O4"]  = { 1.6,  0.68563, -0.1868, -1.35573e-03, 2.3743e-04 };
  lcpomap["U_O4'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U_O5'"] = { 1.6,  0.49392, -0.16038, -1.5512e-04, 1.6453e-04 };
  lcpomap["U_OP1"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U_OP2"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U_OP3"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U_O1P"] = { 1.6,  0.77914, -0.25262, -1.6056e-03, 3.5071e-04 };
  lcpomap["U_O2P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U_O3P"] = { 1.6,  0.88857, -0.33421, -1.8683e-03, 4.9372e-04 };
  lcpomap["U_P"] = { 1.9,  0.03873,  -0.0089339, 8.3582e-06,  3.0381e-06};

  return lcpomap;
}

// assigns LCPO parameters to each atom reading from database
void SAXS::readLCPOparam(const std::vector<std::vector<std::string> > &AtomResidueName, unsigned natoms) {
  std::map<std::string, std::vector<double> > lcpomap = setupLCPOparam();

  for(unsigned i=0; i<natoms; ++i) {
    if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S' || (AtomResidueName[0][i][0]=='P'))) {
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
      if ((AtomResidueName[0][i][0]=='O') || (AtomResidueName[0][i][0]=='N') || (AtomResidueName[0][i][0]=='C') || (AtomResidueName[0][i][0]=='S') || (AtomResidueName[0][i][0]=='P')) {
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

void SAXS::resolution_function() {
  const unsigned numq = q_list.size();

  // only OpenMP because numq might be smaller than the number of ranks
  #pragma omp parallel for num_threads(OpenMP::getNumThreads())
  for (unsigned i=0; i<numq; i++) {
    double qi = q_list[i];
    double dq = 6*sigma_res[i]/(Nj-1);
    double sigma_sq = sigma_res[i]*sigma_res[i];
    double qstart = qi - 3*sigma_res[i];
    for (unsigned j=0; j<Nj; j++) {
      double qj = qstart + j*dq;
      double I0exp = i0e(qj*qi/sigma_sq);

      qj_list[i][j] = qj;
      Rij[i][j] = (qj/sigma_sq)*std::exp(-0.5*(qj - qi)*(qj - qi)/sigma_sq)*I0exp;
    }
  }
}

// i0e function from cephes
// compute I0(x) * exp (-x), with I0 being the modified Bessel function
// of first kind and zeroth order.
inline double SAXS::i0e(double x) {
  double y = 0.0;

  if (x < 0) {
    x = -x;
  }
  if (x <= 8.0) {
    y = (x/2.0) - 2.0;
    return chbevl(y, A);
  }

  return chbevl(32.0/x - 2.0, B) / sqrt(x);
}

double SAXS::chbevl(double x, const std::vector<double> &coeffs) {
  double b0, b1, b2;
  unsigned n = coeffs.size();

  b0 = coeffs[0];
  b1 = 0.0;
  b2 = 0.0;

  for (unsigned i = 1; i < n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + coeffs[i];
  }
  return 0.5 * (b0 - b2);
}

inline double SAXS::interpolation(std::vector<SplineCoeffs> &coeffs, double x) {
  unsigned s = 0;
  while ((x >= q_list[s+1]) && (s+1 < q_list.size()-1)) {
    s++;
  }

  double dx = x - coeffs[s].x;
  return coeffs[s].a + coeffs[s].b*dx + coeffs[s].c*dx*dx + coeffs[s].d*dx*dx*dx;
}

// natural bc cubic spline implementation from the Wikipedia algorithm
// modified from https://stackoverflow.com/a/19216702/3254658
std::vector<SAXS::SplineCoeffs> SAXS::spline_coeffs(std::vector<double> &x, std::vector<double> &y) {
  unsigned n = x.size()-1;
  std::vector<double> a;
  a.insert(a.begin(), y.begin(), y.end());
  std::vector<double> b(n);
  std::vector<double> d(n);
  std::vector<double> h;

  for(unsigned i=0; i<n; i++) {
    h.push_back(x[i+1]-x[i]);
  }

  std::vector<double> alpha;
  alpha.push_back(0);
  for(unsigned i=1; i<n; i++) {
    alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );
  }

  std::vector<double> c(n+1);
  std::vector<double> l(n+1);
  std::vector<double> mu(n+1);
  std::vector<double> z(n+1);
  l[0] = 1;
  mu[0] = 0;
  z[0] = 0;

  for(unsigned i=1; i<n; i++) {
    l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
    mu[i] = h[i]/l[i];
    z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
  }

  l[n] = 1;
  z[n] = 0;
  c[n] = 0;

  for(int j=n-1; j>=0; j--) {
    c[j] = z[j] - mu[j] * c[j+1];
    b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
    d[j] = (c[j+1]-c[j])/3/h[j];
  }

  std::vector<SplineCoeffs> output_set(n);
  for(unsigned i=0; i<n; i++) {
    output_set[i].a = a[i];
    output_set[i].b = b[i];
    output_set[i].c = c[i];
    output_set[i].d = d[i];
    output_set[i].x = x[i];
  }

  return output_set;
}


}//namespace isdb
}//namespace PLMD
