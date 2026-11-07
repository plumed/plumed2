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
#ifndef __PLUMED_isdb_SAXS_h
#define __PLUMED_isdb_SAXS_h

#include "MetainferenceBase.h"
#include "SAXSDebyeBackend.h"
#include "tools/PDB.h"
#include "tools/Random.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

namespace PLMD {
namespace isdb {

/// Implementation of SAXS and SANS.
/// The member functions are split over several files:
/// - SAXS.cpp: setup, ATOMISTIC form factors, Debye sum, metainference;
/// - SAXSMartini.cpp: MARTINI form factors and bead mapping;
/// - SAXSOnebead.cpp: ONEBEAD mapping, form factors (SAXS/SANS) and glycan/ion helpers;
/// - SAXSLCPO.cpp: LCPO parameters used for the SASA-based solvation correction.
class SAXS :
  public MetainferenceBase {
private:
  // ATOMISTIC X-ray/neutron form-factor slots. The first six are the original
  // biological elements; the AFF_* entries are monoatomic ions whose Cromer-Mann
  // and neutron parameters are supplied by getIonData().
  enum { H, C, N, O, P, S,
         AFF_NA, AFF_K, AFF_CL, AFF_CA, AFF_MG, AFF_ZN, AFF_FE2, AFF_FE3, AFF_MN,
         NTT
       };
  // Canonical ion identifiers, shared by the ATOMISTIC and ONEBEAD routes.
  // getIonData() and ionTypeFromResname() are keyed on these values. The AFF_*
  // (atomistic) slots above and the OB_* (onebead) slots below MUST be listed in
  // this same order so that affSlot()/obSlot() can map an ion id to its slot.
  enum { ION_NA, ION_K, ION_CL, ION_CA, ION_MG, ION_ZN, ION_FE2, ION_FE3, ION_MN, NIONS };
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
         DT_5TE, DT_TE3, DT_TE5,
         // Extra protein beads present only in the Martini 3 topology (MARTINI_FF=3)
         ALA_SC1, TRP_SC5, TYR_SC4,
         NMARTINI
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
         // Monosaccharide bead types (one bead per residue, exactly as for the
         // amino acids). GLYC_NAG/GLYC_NGA are the N-acetylated hexosamines
         // (GlcNAc/GalNAc); see getOnebeadparam() for the provenance of the
         // corresponding 21-parameter sets.
         GLYC_FUC, GLYC_MAN, GLYC_BMA, GLYC_GAL, GLYC_GLC, GLYC_NAG, GLYC_NGA, GLYC_SIA,
         // Monoatomic-ion bead types (one atom per bead). Their vac/solv/mix terms are
         // filled analytically by getOnebeadIonparam_*(); see ionTypeFromResname() for
         // residue-name detection. This list MUST stay in the same order as ION_*/AFF_*.
         OB_NA, OB_K, OB_CL, OB_CA, OB_MG, OB_ZN, OB_FE2, OB_FE3, OB_MN,
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
  // rows of the Debye sum computed by this rank: first_row_, first_row_+row_stride_, ...
  // (none if has_rows_ is false)
  unsigned first_row_;
  unsigned row_stride_;
  bool has_rows_;
  // whether the partial results of the ranks are summed
  bool reduce_;
  // the Debye sum is computed by backend_ (USEGPU)
  bool usegpu_;
  // accelerated Debye sum (USEGPU), provided by a plugin
  std::unique_ptr<SAXSDebyeBackend> backend_;
  // the form factors as a single array, ff_flat_[i*numq+k], for the backend
  std::vector<double> ff_flat_;
  bool onebead;
  bool resolution;
  bool chi2;
  bool chi2ensemble;
  int chi2ensemble_replicas;
  bool isFirstStep;
  unsigned nres;
  std::vector<unsigned> atoi;
  std::vector<unsigned> atoms_per_bead;
  std::vector<double>   atoms_masses;
  std::vector<double>   q_list;
  std::vector<double>   expint;
  std::vector<double>   experr;
  std::vector<double>   FF_rank;
  std::vector<std::vector<double> > FF_value_vacuum;
  std::vector<std::vector<double> > FF_value_solv;
  std::vector<std::vector<double> > FF_value_mixed;
  std::vector<std::vector<double> > FF_value;
  // SANS:
  std::vector<std::vector<double> > FF_value_vacuum_H;
  std::vector<std::vector<double> > FF_value_solv_H;
  std::vector<std::vector<double> > FF_value_mixed_H;
  std::vector<std::vector<double> > FF_value_vacuum_D;
  std::vector<std::vector<double> > FF_value_mixed_D;

  std::vector<std::vector<double> > LCPOparam;
  std::vector<unsigned> residue_atom;
  // ONEBEAD glycans. glyc_type[i] is the bead type (GLYC_*) of atom i, or -1 when
  // atom i does not belong to a monosaccharide. glyc_bridge[i] flags the glycan
  // oxygens that carry a glycosidic bond in the TEMPLATE, so that the LCPO SASA
  // treats them as ethers rather than as free hydroxyls. Both are filled by
  // getOnebeadMapping() and are indexed exactly like LCPOparam.
  std::vector<int>  glyc_type;
  std::vector<bool> glyc_bridge;

  double rho, rho_corr, sasa_cutoff;
  double deuter_conc;
  unsigned solv_stride;
  std::vector<double> Iq0_vac;
  std::vector<double> Iq0_solv;
  std::vector<double> Iq0_mix;
  double Iq0;
  // Sign of each bead's net scattering amplitude, +1 or -1.
  //
  // The single-bead approximation stores B(q) = F_vac + rho^2 F_solv - rho F_mix,
  // which is the orientational average of |sum_i A_i exp(i q r_i)|^2 over the
  // atoms of the bead with A_i = b_i - rho v_i. Taking sqrt(B) therefore gives
  // the MAGNITUDE of the bead amplitude and discards its sign. In the Debye sum
  // I(q) = sum_j sum_l F_j F_l sinc(q R_jl) the sign matters as soon as two
  // beads have different signs: a uniform flip cancels, a mixed one does not.
  //
  // This is reachable in SANS near a contrast-match point, and for a nucleus
  // with a negative coherent scattering length such as Mn-55. In SAXS every
  // bead has positive contrast against water, so the sign is always +1 and the
  // result is unchanged.
  //
  // Recomputed together with Iq0 because it depends on the solvent density, on
  // the solvation-layer correction, and in SANS on which residues are deuterated.
  std::vector<double> bead_sign;
  // Sign of A = sum_b - rho sum_v for one bead, from its q=0 terms.
  // At q=0: F_solv(0) = sum_v^2 and F_mix(0) = 2 sum_b sum_v, with sum_v > 0.
  // Hence 2 sum_v A = mix - 2 rho solv, and the two have the same sign.
  static double beadSign(double iq0_solv, double iq0_mix, double rho) {
    return (iq0_mix - 2.*rho*iq0_solv < 0.) ? -1. : 1.;
  }

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
  // the smearing with the resolution function is linear: smeared I(q_k) = sum_l smear_[k*numq+l] I(q_l)
  std::vector<double> smear_;
  // q values as a uniform grid q0+k*dq plus small deviations (see sincos_q)
  double qref0_;
  double qref_dq_;
  double qref_maxeps_;
  std::vector<double> qref_eps_;
  // random numbers for the deuteration of the solvent exposed residues (SANS)
  Random random_deut_;

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

  void calculate_cpu(std::vector<Vector> &pos, std::vector<Vector> &deriv, const bool reduce_deriv);
  // the partial Debye sum of the rows of this rank (without the self term), on the CPU
  void debye_cpu(const std::vector<Vector> &pos, std::vector<double> &sum, std::vector<Vector> &deriv);
  void sincos_q(const double r, double *s, double *c) const;
  void getMartiniFFparam(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter);
  void getMartiniFFparam_v3(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter);
  void getOnebeadparam(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac, std::vector<std::vector<long double> > &parameter_mix, std::vector<std::vector<long double> > &parameter_solv, const std::vector<unsigned> & residue_atoms);
  unsigned getOnebeadMapping(const PDB &pdb, const std::vector<AtomNumber> &atoms);
  double calculateAFF(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double rho);
  std::map<std::string, std::vector<double> > setupLCPOparam();
  void readLCPOparam(const std::vector<std::vector<std::string> > &AtomResidueName, unsigned natoms);
  // LCPO: atoms with LCPO parameters, those with an area, Verlet list of the
  // candidate neighbours (CSR) with the positions at which it was built, and the
  // neighbours (CSR) with the overlap term A(i,j) of each directed pair
  std::vector<unsigned> lcpo_atoms_;
  std::vector<unsigned> lcpo_area_atoms_;
  std::vector<Vector> lcpo_ref_pos_;
  std::vector<std::size_t> lcpo_cand_off_;
  std::vector<unsigned> lcpo_cand_idx_;
  std::vector<double> lcpo_cand_d_;
  std::vector<std::size_t> lcpo_off_;
  std::vector<unsigned> lcpo_idx_;
  std::vector<double> lcpo_term_;
  void lcpoBuildCandidates();
  void lcpoUpdateList();
  void sasa_calculate(std::vector<bool> &solv_res);
  // SANS:
  void getOnebeadparam_sansH(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_H, std::vector<std::vector<long double> > &parameter_mix_H, std::vector<std::vector<long double> > &parameter_solv_H);
  void getOnebeadparam_sansD(const PDB &pdb, const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &parameter_vac_D, std::vector<std::vector<long double> > &parameter_mix_D);
  double calculateAFFsans(const std::vector<AtomNumber> &atoms, std::vector<std::vector<long double> > &FF_tmp, const double deuter_conc);
  // Ions (shared by ATOMISTIC and ONEBEAD).
  // getIonData returns, for a canonical ION_* id: Cromer-Mann a[4]/b[4]/c (X-ray),
  // the displaced-solvent volume vol (A^3), the bound coherent neutron scattering
  // length neutron_b (10^-12 cm), and the atomic mass (g/mol).
  void getIonData(unsigned ion, std::vector<double> &cm_a, std::vector<double> &cm_b, double &cm_c, double &vol, double &neutron_b, double &mass) const;
  // Returns the canonical ION_* id for a residue name, or -1 if it is not a known ion.
  int ionTypeFromResname(const std::string &resname) const;
  // True for monosaccharides that are recognised but not parameterised.
  bool isUnsupportedGlycan(const std::string &resname) const;
  // Canonical name of a residue whose name encodes a glycosylation site.
  std::string onebeadResidueAlias(const std::string &resname) const;
  // True when the residue is represented by a single bead (amino acids, incl. CYX).
  bool isOnebeadProteinResidue(const PDB &pdb, const std::string &resname) const;
  // Returns the GLYC_* bead type for a monosaccharide residue name, or -1 if the
  // residue is not a supported glycan. hasN must be true when the residue carries
  // an amide nitrogen: that is what distinguishes GlcNAc from Glc and GalNAc from
  // Gal when the input file names an N-acetylhexosamine after its parent hexose.
  int glycanTypeFromResname(const std::string &resname, bool hasN) const;
  // Map a canonical ion id to its atomistic / onebead form-factor slot.
  inline unsigned affSlot(unsigned ion) const {
    return AFF_NA + ion;
  }
  inline unsigned obSlot(unsigned ion) const {
    return OB_NA + ion;
  }
  // Fill the onebead vac/solv/mix terms (and q=0 Iq0 values) for every ion bead.
  void getOnebeadIonparam_saxs(std::vector<std::vector<long double> > &FF_tmp_vac, std::vector<std::vector<long double> > &FF_tmp_mix, std::vector<std::vector<long double> > &FF_tmp_solv, std::vector<double> &iq0_vac, std::vector<double> &iq0_mix, std::vector<double> &iq0_solv) const;
  void getOnebeadIonparam_sans(std::vector<std::vector<long double> > &FF_tmp_vac_H, std::vector<std::vector<long double> > &FF_tmp_mix_H, std::vector<std::vector<long double> > &FF_tmp_solv_H, std::vector<std::vector<long double> > &FF_tmp_vac_D, std::vector<std::vector<long double> > &FF_tmp_mix_D, std::vector<double> &iq0_vac_H, std::vector<double> &iq0_mix_H, std::vector<double> &iq0_solv_H, std::vector<double> &iq0_vac_D, std::vector<double> &iq0_mix_D) const;
  void resolution_function();
  std::vector<SplineCoeffs> spline_coeffs(std::vector<double> &x, std::vector<double> &y);
  inline double interpolation(std::vector<SplineCoeffs> &coeffs, double x);
  inline double i0e(double x);
  double chbevl(double x, const std::vector<double> &coeffs);

public:
  static void registerKeywords( Keywords& keys );
  explicit SAXS(const ActionOptions&);
  void calculate() override;
};

}//namespace isdb
}//namespace PLMD

#endif
