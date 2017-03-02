/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017 The plumed team
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
#ifdef ENABLEC11 

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/SetupMolInfo.h"
#include "tools/OpenMP.h"

// These are only required when using the vectorized version
#include "emmintrin.h"
#include "immintrin.h"

#define INV_PI_SQRT_PI 0.179587122
#define KCAL_TO_KJ 4.184
#define ANG_TO_NM 0.1
#define ANG3_TO_NM3 0.001
#define INNER_LOOP_STRIDE 4

// WARNING: Here be dragons!
// Uncomment the next line for SIMD version, which is about 1.25 times faster
// clang will choke on this, gcc requires -mveclibabi=svml, icpc should work just fine
// #define AVXNLLOOP

using namespace std;

namespace PLMD {
    namespace colvar {

//+PLUMEDOC COLVAR IMPLICIT
/*

Calculate EEF1-SB solvation free energy for a group of atoms.

EEF1-SB is a solvent-accessible surface area based model, where the free energy of solvation is computed using a pairwise interaction term for non-hydrogen atoms:
\f[
    \Delta G^\mathrm{solv}_i = \Delta G^\mathrm{ref}_i - \sum_{j \neq i} f_i(r_{ij}) V_j
\f]
where \f$\Delta G^\mathrm{solv}_i\f$ is the free energy of solvation, \f$\Delta G^\mathrm{ref}_i\f$ is the reference solvation free energy, \f$V_j\f$ is the volume of atom \f$j\f$ and
\f[
    f_i(r) 4\pi r^2 = \frac{2}{\sqrt{\pi}} \frac{\Delta G^\mathrm{free}_i}{\lambda_i} \exp\left\{ - \frac{(r-R_i)^2}{\lambda^2_i}\right\}
\f]
where \f$\Delta G^\mathrm{free}_i\f$ is the solvation free energy of the isolated group, \f$\lambda_i\f$ is the correlation length equal to the width of the first solvation shell and \f$R_i\f$ is the van der Waals radius of atom \f$i\f$.

The output from this collective variable, the free energy of solvation, can be used with the \ref BIASVALUE keyword to provide implicit solvation to a system. All parameters are designed to be used with a modified CHARMM36 force field. It takes only non-hydrogen atoms as input, these can be conveniently specified using the \ref GROUP action with the NDX_GROUP parameter. To speed up the calculation, IMPLICIT internally uses a neighbourlist with a cutoff dependent on the type of atom (maximum of 1.95 nm). This cutoff can be extended further by using the NL_BUFFER keyword.

\par Examples
\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=peptide.pdb
WHOLEMOLECULES ENTITY0=1-111

# This allows us to select only non-hydrogen atoms
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H

# We extend the cutoff by 0.2 nm and update the neighbourlist every 10 steps
solv: IMPLICIT ATOMS=protein-h NL_STRIDE=10 NL_BUFFER=0.2

# Here we actually add our calculated energy back to the potential
bias: BIASVALUE ARG=solv

PRINT ARG=solv FILE=SOLV
\endverbatim
(see also \ref PRINT, \ref GROUP, \ref MOLINFO, \ref WHOLEMOLECULES)

*/
//+ENDPLUMEDOC

        class Implicit : public Colvar {
            private:
                bool pbc;
                bool tcorr;
                double buffer;
                unsigned stride;
                unsigned nl_update;
                vector<vector<unsigned> > nl;
                vector<vector<double> > parameter;
                map<string, map<string, string> > typemap;
                map<string, vector<double> > valuemap;
                void setupTypeMap();

            public:
                static void registerKeywords(Keywords& keys);
                explicit Implicit(const ActionOptions&);
                virtual void calculate();
                void update_neighb();
                void setupConstants(const vector<AtomNumber> &atoms, vector<vector<double> > &parameter);
        };

        PLUMED_REGISTER_ACTION(Implicit,"IMPLICIT")

        void Implicit::registerKeywords(Keywords& keys) {
            Colvar::registerKeywords(keys);
            componentsAreNotOptional(keys);
            useCustomisableComponents(keys);
            keys.add("atoms", "ATOMS", "The atoms to be included in the calculation, e.g. the whole protein.");
            keys.add("compulsory", "NL_BUFFER", "The buffer to the intrinsic cutoff used when calculating pairwise interactions.");
            keys.add("compulsory", "NL_STRIDE", "The frequency with which the neighbourlist is updated.");
            keys.addFlag("TEMP_CORRECTION", false, "Correct free energy of solvation constants for temperatures different from 298.15 K");
        }

        Implicit::Implicit(const ActionOptions&ao):
            PLUMED_COLVAR_INIT(ao),
            pbc(true),
            tcorr(false),
            buffer(0.1),
            stride(10),
            nl_update(0)
        {
            vector<AtomNumber> atoms;
            parseAtomList("ATOMS", atoms);
            const unsigned size = atoms.size();

            parseFlag("TEMP_CORRECTION", tcorr);
            parse("NL_BUFFER", buffer);
            parse("NL_STRIDE", stride);

            bool nopbc = !pbc;
            parseFlag("NOPBC", nopbc);
            pbc = !nopbc;

            checkRead();

            log << "  Bibliography " << plumed.cite("Bottaro S, Lindorff-Larsen K, Best R, J. Chem. Theory Comput. 9, 5641 (2013)")
                                     << plumed.cite("Lazaridis T, Karplus M, Proteins Struct. Funct. Genet. 35, 133 (1999)"); log << "\n";


            nl.resize(size);
            parameter.resize(size);
            setupConstants(atoms, parameter);

            addValueWithDerivatives();
            setNotPeriodic();
            requestAtoms(atoms);
        }

        void Implicit::update_neighb() {
            const double lower_c2 = 0.24 * 0.24; // this is the cut-off for bonded atoms
            const unsigned size = getNumberOfAtoms();
            const unsigned nt = OpenMP::getGoodNumThreads(nl);
            #pragma omp parallel num_threads(nt)
            {
                #pragma omp for
                for (unsigned i=0; i<size; ++i) {
                    const Vector posi = getPosition(i);

                    // Clear old values, reserve space for new indices
                    nl[i].clear();
                    nl[i].reserve(100);

                    // Loop through neighboring atoms, add the ones below cutoff
                    for (unsigned j=i+1; j<size; ++j) {
                        const double d2 = delta(posi, getPosition(j)).modulo2();
                        if (d2 < lower_c2 && j < i+14) {
                            // crude approximation for i-i+1/2 interactions,
                            // we want to exclude atoms separated by less than three bonds
                            continue;
                        }

                        double mlambda = parameter[i][5];
                        if (parameter[j][5] > mlambda) {
                            // We choose the maximum lambda value and use a more conservative cutoff
                            mlambda = parameter[j][5];
                        }

                        double c2 = (3. * mlambda + buffer) * (3. * mlambda + buffer);
                        if (d2 < c2 ) {
                            nl[i].push_back(j);
                        }
                    }
                }
            }
        }

#ifdef AVXNLLOOP
        void Implicit::calculate() {
            if (pbc) makeWhole();
            if (getExchangeStep()) nl_update = 0;
            if (nl_update == 0) update_neighb();

            const unsigned size = getNumberOfAtoms();
            vector<Vector> fedensity_deriv(size);
            double bias = 0.0;

            // Define constants (4 64bit doubles in one 256bit (AVX) register)
            const __m256d m_INV_PI_SQRT_PI = _mm256_set1_pd(INV_PI_SQRT_PI);
            const __m256d m_ones           = _mm256_set1_pd(1.0);
            const __m256d m_zero           = _mm256_set1_pd(0.0);

            // Main loop over all non-hydrogen atoms
            for (unsigned i = 0; i < size; ++i) {
                const Vector posi = getPosition(i);
                const double delta_g_ref = parameter[i][1];

                // Get the inner loop size from the neighbourlist
                const unsigned inner_size = nl[i].size();

                // size modulo 4 gives us the size of the remainder
                const unsigned inner_remainder = inner_size % INNER_LOOP_STRIDE;

                // Set free-energy density accumulator to 0
                __m256d m_fedensity = _mm256_set1_pd(0.0);

                // Main vectorized loop, we handle 4 values at once
                for (unsigned i_nl = 0; i_nl < (inner_size - inner_remainder); i_nl += INNER_LOOP_STRIDE) {
                    const unsigned j0 = nl[i][i_nl + 0];
                    const unsigned j1 = nl[i][i_nl + 1];
                    const unsigned j2 = nl[i][i_nl + 2];
                    const unsigned j3 = nl[i][i_nl + 3];

                    // Calculate distance vector between atom i and j
                    const Vector dist0 = delta(posi, getPosition(j0));
                    const Vector dist1 = delta(posi, getPosition(j1));
                    const Vector dist2 = delta(posi, getPosition(j2));
                    const Vector dist3 = delta(posi, getPosition(j3));

                    // Calculate distance between atom i and j
                    const __m256d m_rij               = _mm256_set_pd(dist3.modulo(),
                                                                      dist2.modulo(),
                                                                      dist1.modulo(),
                                                                      dist0.modulo());
                    // 1 / r_ij
                    const __m256d m_inv_rij           = _mm256_div_pd(m_ones, m_rij);
                    // 1 / r_ij^2
                    const __m256d m_inv_rij2          = _mm256_mul_pd(m_inv_rij, m_inv_rij);

                    // i-j interaction
                    {
                        // R_i
                        const __m256d m_vdwr          = _mm256_set1_pd(parameter[i][6]);
                        const __m256d m_lambda        = _mm256_set1_pd(parameter[i][5]);
                        // 1 / lambda
                        const __m256d m_inv_lambda    = _mm256_div_pd(m_ones, m_lambda);
                        // 1 / lambda^2
                        const __m256d m_inv_lambda2   = _mm256_mul_pd(m_inv_lambda, m_inv_lambda);
                        // Delta G^free
                        const __m256d m_delta_g_free  = _mm256_set1_pd(parameter[i][2]);
                        // V_j
                        const __m256d m_vdw_volume    = _mm256_set_pd(parameter[j3][0],
                                                                      parameter[j2][0],
                                                                      parameter[j1][0],
                                                                      parameter[j0][0]);

                        // r_ij - R_i
                        const __m256d m_rij_vdwr_diff = _mm256_sub_pd(m_rij, m_vdwr);

                        // exp(-(r_ij - R_i)^2 / lambda^2)
                              __m256d m_exponent      = _mm256_mul_pd(m_rij_vdwr_diff, m_rij_vdwr_diff);
                                      m_exponent      = _mm256_mul_pd(m_inv_lambda2, m_exponent);
                        // _mm256_exp_pd requires SVML
                        const __m256d m_expo          = _mm256_exp_pd(_mm256_sub_pd(m_zero, m_exponent));

                        // Delta G^free * V_j * expo / (lambda * r_ij^2 * pi * sqrt(pi))
                        const __m256d m_fact_a        = _mm256_mul_pd(m_delta_g_free, m_vdw_volume);
                        const __m256d m_fact_b        = _mm256_mul_pd(m_expo, m_INV_PI_SQRT_PI);
                        const __m256d m_fact          = _mm256_mul_pd(_mm256_mul_pd(m_fact_a, m_fact_b),
                                                                      _mm256_mul_pd(m_inv_rij2, m_inv_lambda));

                        // (1 / r_ij + (r_ij - R_i) / lambda^2) * 1 / r_ij * fact
                              __m256d m_deriv         = _mm256_add_pd(m_inv_rij, _mm256_mul_pd(m_rij_vdwr_diff, m_inv_lambda2));
                                      m_deriv         = _mm256_mul_pd(_mm256_mul_pd(m_inv_rij, m_fact), m_deriv);

                        // Sum up bias and derivatives
                        m_fedensity = _mm256_add_pd(m_fedensity, m_fact);
                        fedensity_deriv[i]  += m_deriv[0] * dist0;
                        fedensity_deriv[i]  += m_deriv[1] * dist1;
                        fedensity_deriv[i]  += m_deriv[2] * dist2;
                        fedensity_deriv[i]  += m_deriv[3] * dist3;
                        fedensity_deriv[j0] -= m_deriv[0] * dist0;
                        fedensity_deriv[j1] -= m_deriv[1] * dist1;
                        fedensity_deriv[j2] -= m_deriv[2] * dist2;
                        fedensity_deriv[j3] -= m_deriv[3] * dist3;
                    }

                    // j-i interaction
                    {
                        const __m256d m_vdwr              = _mm256_set_pd(parameter[j3][6],
                                                                          parameter[j2][6],
                                                                          parameter[j1][6],
                                                                          parameter[j0][6]);
                        const __m256d m_lambda            = _mm256_set_pd(parameter[j3][5],
                                                                          parameter[j2][5],
                                                                          parameter[j1][5],
                                                                          parameter[j0][5]);
                        const __m256d m_inv_lambda        = _mm256_div_pd(m_ones, m_lambda);
                        const __m256d m_inv_lambda2       = _mm256_mul_pd(m_inv_lambda, m_inv_lambda);
                        const __m256d m_delta_g_free      = _mm256_set_pd(parameter[j3][2],
                                                                          parameter[j2][2],
                                                                          parameter[j1][2],
                                                                          parameter[j0][2]);
                        const __m256d m_vdw_volume        = _mm256_set1_pd(parameter[i][0]);

                        const __m256d m_rij_vdwr_diff = _mm256_sub_pd(m_rij, m_vdwr);
                              __m256d m_exponent      = _mm256_mul_pd(m_rij_vdwr_diff, m_rij_vdwr_diff);
                                      m_exponent      = _mm256_mul_pd(m_inv_lambda2, m_exponent);
                        const __m256d m_expo          = _mm256_exp_pd(_mm256_sub_pd(m_zero, m_exponent));

                        const __m256d m_fact_a        = _mm256_mul_pd(m_delta_g_free, m_vdw_volume);
                        const __m256d m_fact_b        = _mm256_mul_pd(m_expo, m_INV_PI_SQRT_PI);
                        const __m256d m_fact          = _mm256_mul_pd(_mm256_mul_pd(m_fact_a, m_fact_b),
                                                                      _mm256_mul_pd(m_inv_rij2, m_inv_lambda));

                              __m256d m_deriv         = _mm256_add_pd(m_inv_rij, _mm256_mul_pd(m_rij_vdwr_diff, m_inv_lambda2));
                                      m_deriv         = _mm256_mul_pd(_mm256_mul_pd(m_inv_rij, m_fact), m_deriv);

                        m_fedensity = _mm256_add_pd(m_fedensity, m_fact);
                        fedensity_deriv[i]  += m_deriv[0] * dist0;
                        fedensity_deriv[i]  += m_deriv[1] * dist1;
                        fedensity_deriv[i]  += m_deriv[2] * dist2;
                        fedensity_deriv[i]  += m_deriv[3] * dist3;
                        fedensity_deriv[j0] -= m_deriv[0] * dist0;
                        fedensity_deriv[j1] -= m_deriv[1] * dist1;
                        fedensity_deriv[j2] -= m_deriv[2] * dist2;
                        fedensity_deriv[j3] -= m_deriv[3] * dist3;
                    }
                }

                // Sum up 4 elements of vector
                double fedensity = m_fedensity[0] + m_fedensity[1] + m_fedensity[2] + m_fedensity[3];

                // Remainder loop, like the generic version, but max 3 iterations
                for (unsigned i_nl = (inner_size - inner_remainder); i_nl < inner_size; ++i_nl) {
                    const unsigned j = nl[i][i_nl];
                    const Vector dist = delta(posi, getPosition(j));
                    const double rij = dist.modulo();
                    const double inv_rij = 1.0 / rij;
                    const double inv_rij2 = inv_rij * inv_rij;

                    // i-j interaction
                    if(rij < 3.*parameter[i][5])
                    {
                        const double delta_g_free = parameter[i][2];
                        const double lambda = parameter[i][5];
                        const double vdw_radius = parameter[i][6];
                        const double inv_lambda = 1.0 / lambda;
                        const double inv_lambda2 = inv_lambda * inv_lambda;
                        const double vdw_volume = parameter[j][0];

                        const double rij_vdwr_diff = rij - vdw_radius;
                        const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                        const double fact = delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                        const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                        fedensity += fact;
                        fedensity_deriv[i] += deriv * dist;
                        fedensity_deriv[j] -= deriv * dist;
                    }

                    // j-i interaction
                    if(rij < 3.*parameter[j][5])
                    {
                        const double delta_g_free = parameter[j][2];
                        const double lambda = parameter[j][5];
                        const double vdw_radius = parameter[j][6];
                        const double inv_lambda = 1.0 / lambda;
                        const double inv_lambda2 = inv_lambda * inv_lambda;
                        const double vdw_volume = parameter[i][0];

                        const double rij_vdwr_diff = rij - vdw_radius;
                        const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                        const double fact = delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                        const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                        fedensity += fact;
                        fedensity_deriv[i] += deriv * dist;
                        fedensity_deriv[j] -= deriv * dist;
                    }
                }
                bias += delta_g_ref - 0.5 * fedensity;
            }

            Tensor deriv_box;
            for (unsigned i=0; i<size; ++i) {
                setAtomsDerivatives(i, -fedensity_deriv[i]);
                deriv_box += Tensor(getPosition(i), -fedensity_deriv[i]);
            }
            setBoxDerivatives(-deriv_box);
            setValue(bias);

            // Keep track of the neighbourlist updates
            ++nl_update;
            if (nl_update == stride) {
                nl_update = 0;
            }
        }

#else
        void Implicit::calculate() {
            if(pbc) makeWhole();
            if(getExchangeStep()) nl_update = 0;
            if (nl_update == 0) {
                update_neighb();
            }

            const unsigned size=getNumberOfAtoms();
            vector<Vector> fedensity_deriv(size);
            double bias = 0.0;

            const unsigned nt = OpenMP::getGoodNumThreads(nl);
            #pragma omp parallel num_threads(nt)
            {
                #pragma omp for reduction(+:bias)
                for (unsigned i=0; i<size; ++i) {
                    const Vector posi = getPosition(i);
                    const double delta_g_ref = parameter[i][1];
                    double fedensity = 0.0;

                    // The pairwise interactions are unsymmetric, but we can get away with calculating the distance only once
                    for (unsigned i_nl=0; i_nl<nl[i].size(); ++i_nl) {
                        const unsigned j = nl[i][i_nl];
                        const Vector dist = delta(posi, getPosition(j));
                        const double rij = dist.modulo();
                        const double inv_rij = 1.0 / rij;
                        const double inv_rij2 = inv_rij * inv_rij;

                        // i-j interaction
                        if(rij < 3.*parameter[i][5])
                        {
                            const double delta_g_free = parameter[i][2];
                            const double lambda = parameter[i][5];
                            const double vdw_radius = parameter[i][6];
                            const double inv_lambda = 1.0 / lambda;
                            const double inv_lambda2 = inv_lambda * inv_lambda;
                            const double vdw_volume = parameter[j][0];

                            const double rij_vdwr_diff = rij - vdw_radius;
                            const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                            const double fact = delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                            const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                            // This is needed for correct box derivs
                            #pragma omp critical(deriv)
                            {
                                fedensity += fact;
                                fedensity_deriv[i] += deriv * dist;
                                fedensity_deriv[j] -= deriv * dist;
                            }
                        }

                        // j-i interaction
                        if(rij < 3.*parameter[j][5])
                        {
                            const double delta_g_free = parameter[j][2];
                            const double lambda = parameter[j][5];
                            const double vdw_radius = parameter[j][6];
                            const double inv_lambda = 1.0 / lambda;
                            const double inv_lambda2 = inv_lambda * inv_lambda;
                            const double vdw_volume = parameter[i][0];

                            const double rij_vdwr_diff = rij - vdw_radius;
                            const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                            const double fact = delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                            const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                            #pragma omp critical(deriv)
                            {
                                fedensity += fact;
                                fedensity_deriv[i] += deriv * dist;
                                fedensity_deriv[j] -= deriv * dist;
                            }
                        }
                    }
                    bias += delta_g_ref - 0.5 * fedensity;
                }
            }

            Tensor deriv_box;
            const unsigned ntd = OpenMP::getGoodNumThreads(fedensity_deriv);
            {
                for (unsigned i=0; i<size; ++i) {
                    setAtomsDerivatives(i, -fedensity_deriv[i]);
                    deriv_box += Tensor(getPosition(i), -fedensity_deriv[i]);
                }
            }

            setBoxDerivatives(-deriv_box);
            setValue(bias);

            // Keep track of the neighbourlist updates
            ++nl_update;
            if (nl_update == stride) {
                nl_update = 0;
            }
        }
#endif

        void Implicit::setupConstants(const vector<AtomNumber> &atoms, vector<vector<double> > &parameter) {
            setupTypeMap();
            vector<SetupMolInfo*> moldat = plumed.getActionSet().select<SetupMolInfo*>();
            if (moldat.size() == 1) {
                log << "  MOLINFO DATA found, using proper atom names\n";
                for(unsigned i=0; i<atoms.size(); ++i) {

                    // Get atom and residue names
                    string Aname = moldat[0]->getAtomName(atoms[i]);
                    string Rname = moldat[0]->getResidueName(atoms[i]);
                    string Atype = typemap[Rname][Aname];

                    // Check for terminal COOH or COO- (different atomtypes & parameters!)
                    if (moldat[0]->getAtomName(atoms[i]) == "OT1") {
                        // We create a temporary AtomNumber object to access future atoms
                        unsigned ai = atoms[i].index();
                        AtomNumber tmp_an;
                        tmp_an.setIndex(ai + 2);
                        if (moldat[0]->getAtomName(tmp_an) == "HT2") {
                            // COOH
                            Atype = "OB";
                        } else {
                            // COO-
                            Atype = "OC";
                        }
                    }
                    if (moldat[0]->getAtomName(atoms[i]) == "OT2") {
                        unsigned ai = atoms[i].index();
                        AtomNumber tmp_an;
                        tmp_an.setIndex(ai + 1);
                        if (moldat[0]->getAtomName(tmp_an) == "HT2") {
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
                    if (!isdigit(first)){
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
                        parameter[i] = valuemap.at(Atype);
                    } catch (exception &e) {
                        log << "Type: " << Atype << "  Name: " << Aname << "  Residue: " << Rname << "\n";
                        error("Invalid atom type!\n");
                    }

                    // Temperature correction
                    if (tcorr && parameter[i][1] > 0.0) {
                        const double t0 = 298.15;
                        const double delta_g_ref_t0 = parameter[i][1];
                        const double delta_h_ref_t0 = parameter[i][3];
                        const double delta_cp = parameter[i][4];
                        const double delta_s_ref_t0 = (delta_h_ref_t0 - delta_g_ref_t0) / t0;
                        const double t = plumed.getAtoms().getKbT() / plumed.getAtoms().getKBoltzmann();
                        parameter[i][1] -= delta_s_ref_t0 * (t - t0) - delta_cp * t * std::log(t / t0) + delta_cp * (t - t0);
                        parameter[i][2] *= parameter[i][1] / delta_g_ref_t0;
                    }
                }
            } else {
                error("MOLINFO DATA not found\n");
            }
        }

        void Implicit::setupTypeMap() {
            typemap = {
                {"ALA", {
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
                {"ARG", {
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
                {"ASN", {
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
                {"ASPP", {
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
                {"ASP", {
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
                {"CYS", {
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
                {"GLN", {
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
                {"GLUP", {
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
                {"GLU", {
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
                {"GLY", {
                            {"N",   "NH1"},
                            {"HN",  "H"  },
                            {"CA",  "CT2"},
                            {"HA1", "HB2"},
                            {"HA2", "HB2"},
                            {"C",   "C"  },
                            {"O",   "O"  }
                        }
                },
                {"HSD", {
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
                {"HIS", {
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
                {"HSE", {
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
                {"HSP", {
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
                {"ILE", {
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
                {"LEU", {
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
                {"LYS", {
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
                {"MET", {
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
                {"PHE", {
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
                {"PRO", {
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
                {"SER", {
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
                {"THR", {
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
                {"TRP", {
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
                {"TYR", {
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
                {"VAL", {
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

            // Volume ∆Gref ∆Gfree ∆H ∆Cp λ vdw_radius
            valuemap = {
                {"C", {
                          ANG3_TO_NM3 * 14.720,
                          KCAL_TO_KJ * 0.000,
                          KCAL_TO_KJ * 0.000,
                          KCAL_TO_KJ * 0.000,
                          KCAL_TO_KJ * 0.0,
                          ANG_TO_NM * 3.5,
                          0.20,
                      }
                },
                {"CD", {
                           ANG3_TO_NM3 * 14.720,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.0,
                           ANG_TO_NM * 3.5,
                           0.20,
                       }
                },
                {"CT1", {
                            ANG3_TO_NM3 * 11.507,
                            KCAL_TO_KJ * -0.187,
                            KCAL_TO_KJ * -0.187,
                            KCAL_TO_KJ * 0.876,
                            KCAL_TO_KJ * 0.0,
                            ANG_TO_NM * 3.5,
                            0.20,
                        }
                },
                {"CT2", {
                            ANG3_TO_NM3 * 18.850,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * -0.610,
                            KCAL_TO_KJ * 18.6,
                            ANG_TO_NM * 3.5,
                            0.20,
                        }
                },
                {"CT2A", {
                             ANG3_TO_NM3 * 18.666,
                             KCAL_TO_KJ * 0.372,
                             KCAL_TO_KJ * 0.372,
                             KCAL_TO_KJ * -0.610,
                             KCAL_TO_KJ * 18.6,
                             ANG_TO_NM * 3.5,
                             0.20,
                         }
                },
                {"CT3", {
                            ANG3_TO_NM3 * 27.941,
                            KCAL_TO_KJ * 1.089,
                            KCAL_TO_KJ * 1.089,
                            KCAL_TO_KJ * -1.779,
                            KCAL_TO_KJ * 35.6,
                            ANG_TO_NM * 3.5,
                            0.204,
                        }
                },
                {"CPH1", {
                             ANG3_TO_NM3 * 5.275,
                             KCAL_TO_KJ * 0.057,
                             KCAL_TO_KJ * 0.080,
                             KCAL_TO_KJ * -0.973,
                             KCAL_TO_KJ * 6.9,
                             ANG_TO_NM * 3.5,
                             0.18,
                         }
                },
                {"CPH2", {
                             ANG3_TO_NM3 * 11.796,
                             KCAL_TO_KJ * 0.057,
                             KCAL_TO_KJ * 0.080,
                             KCAL_TO_KJ * -0.973,
                             KCAL_TO_KJ * 6.9,
                             ANG_TO_NM * 3.5,
                             0.18,
                         }
                },
                {"CPT", {
                            ANG3_TO_NM3 * 4.669,
                            KCAL_TO_KJ * -0.890,
                            KCAL_TO_KJ * -0.890,
                            KCAL_TO_KJ * 2.220,
                            KCAL_TO_KJ * 6.9,
                            ANG_TO_NM * 3.5,
                            0.186,
                        }
                },
                {"CY", {
                           ANG3_TO_NM3 * 10.507,
                           KCAL_TO_KJ * -0.890,
                           KCAL_TO_KJ * -0.890,
                           KCAL_TO_KJ * 2.220,
                           KCAL_TO_KJ * 6.9,
                           ANG_TO_NM * 3.5,
                           0.199,
                       }
                },
                {"CP1", {
                            ANG3_TO_NM3 * 25.458,
                            KCAL_TO_KJ * -0.187,
                            KCAL_TO_KJ * -0.187,
                            KCAL_TO_KJ * 0.876,
                            KCAL_TO_KJ * 0.0,
                            ANG_TO_NM * 3.5,
                            0.227,
                        }
                },
                {"CP2", {
                            ANG3_TO_NM3 * 19.880,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * -0.610,
                            KCAL_TO_KJ * 18.6,
                            ANG_TO_NM * 3.5,
                            0.217,
                        }
                },
                {"CP3", {
                            ANG3_TO_NM3 * 26.731,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * 0.372,
                            KCAL_TO_KJ * -0.610,
                            KCAL_TO_KJ * 18.6,
                            ANG_TO_NM * 3.5,
                            0.217,
                        }
                },
                {"CC", {
                           ANG3_TO_NM3 * 16.539,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.000,
                           KCAL_TO_KJ * 0.0,
                           ANG_TO_NM * 3.5,
                           0.20,
                       }
                },
                {"CAI", {
                            ANG3_TO_NM3 * 18.249,
                            KCAL_TO_KJ * 0.057,
                            KCAL_TO_KJ * 0.057,
                            KCAL_TO_KJ * -0.973,
                            KCAL_TO_KJ * 6.9,
                            ANG_TO_NM * 3.5,
                            0.199,
                        }
                },
                {"CA", {
                           ANG3_TO_NM3 * 18.249,
                           KCAL_TO_KJ * 0.057,
                           KCAL_TO_KJ * 0.057,
                           KCAL_TO_KJ * -0.973,
                           KCAL_TO_KJ * 6.9,
                           ANG_TO_NM * 3.5,
                           0.199,
                       }
                },
                {"N", {
                          ANG3_TO_NM3 * 0.000,
                          KCAL_TO_KJ * -1.000,
                          KCAL_TO_KJ * -1.000,
                          KCAL_TO_KJ * -1.250,
                          KCAL_TO_KJ * 8.8,
                          ANG_TO_NM * 3.5,
                          0.185,
                      }
                },
                {"NR1", {
                            ANG3_TO_NM3 * 15.273,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -9.059,
                            KCAL_TO_KJ * -8.8,
                            ANG_TO_NM * 3.5,
                            0.185,
                        }
                },
                {"NR2", {
                            ANG3_TO_NM3 * 15.111,
                            KCAL_TO_KJ * -3.820,
                            KCAL_TO_KJ * -3.820,
                            KCAL_TO_KJ * -4.654,
                            KCAL_TO_KJ * -8.8,
                            ANG_TO_NM * 3.5,
                            0.185,
                        }
                },
                {"NR3", {
                            ANG3_TO_NM3 * 15.071,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -9.059,
                            KCAL_TO_KJ * -8.8,
                            ANG_TO_NM * 3.5,
                            0.185,
                        }
                },
                {"NH1", {
                            ANG3_TO_NM3 * 10.197,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -9.059,
                            KCAL_TO_KJ * -8.8,
                            ANG_TO_NM * 3.5,
                            0.185,
                        }
                },
                {"NH2", {
                            ANG3_TO_NM3 * 18.182,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -5.950,
                            KCAL_TO_KJ * -9.059,
                            KCAL_TO_KJ * -8.8,
                            ANG_TO_NM * 3.5,
                            0.185,
                        }
                },
                {"NH3", {
                            ANG3_TO_NM3 * 18.817,
                            KCAL_TO_KJ * -20.000,
                            KCAL_TO_KJ * -20.000,
                            KCAL_TO_KJ * -25.000,
                            KCAL_TO_KJ * -18.0,
                            ANG_TO_NM * 6.0,
                            0.185,
                        }
                },
                {"NC2", {
                            ANG3_TO_NM3 * 18.215,
                            KCAL_TO_KJ * -10.000,
                            KCAL_TO_KJ * -10.000,
                            KCAL_TO_KJ * -12.000,
                            KCAL_TO_KJ * -7.0,
                            ANG_TO_NM * 6.0,
                            0.185,
                        }
                },
                {"NY", {
                           ANG3_TO_NM3 * 12.001,
                           KCAL_TO_KJ * -5.950,
                           KCAL_TO_KJ * -5.950,
                           KCAL_TO_KJ * -9.059,
                           KCAL_TO_KJ * -8.8,
                           ANG_TO_NM * 3.5,
                           0.185,
                       }
                },
                {"NP", {
                           ANG3_TO_NM3 * 4.993,
                           KCAL_TO_KJ * -20.000,
                           KCAL_TO_KJ * -20.000,
                           KCAL_TO_KJ * -25.000,
                           KCAL_TO_KJ * -18.0,
                           ANG_TO_NM * 6.0,
                           0.185,
                       }
                },
                {"O", {
                          ANG3_TO_NM3 * 11.772,
                          KCAL_TO_KJ * -5.330,
                          KCAL_TO_KJ * -5.330,
                          KCAL_TO_KJ * -5.787,
                          KCAL_TO_KJ * -8.8,
                          ANG_TO_NM * 3.5,
                          0.170,
                      }
                },
                {"OB", {
                           ANG3_TO_NM3 * 11.694,
                           KCAL_TO_KJ * -5.330,
                           KCAL_TO_KJ * -5.330,
                           KCAL_TO_KJ * -5.787,
                           KCAL_TO_KJ * -8.8,
                           ANG_TO_NM * 3.5,
                           0.170,
                       }
                },
                {"OC", {
                           ANG3_TO_NM3 * 12.003,
                           KCAL_TO_KJ * -10.000,
                           KCAL_TO_KJ * -10.000,
                           KCAL_TO_KJ * -12.000,
                           KCAL_TO_KJ * -9.4,
                           ANG_TO_NM * 6.0,
                           0.170,
                       }
                },
                {"OH1", {
                            ANG3_TO_NM3 * 15.528,
                            KCAL_TO_KJ * -5.920,
                            KCAL_TO_KJ * -5.920,
                            KCAL_TO_KJ * -9.264,
                            KCAL_TO_KJ * -11.2,
                            ANG_TO_NM * 3.5,
                            0.177,
                        }
                },
                {"OS", {
                           ANG3_TO_NM3 * 6.774,
                           KCAL_TO_KJ * -2.900,
                           KCAL_TO_KJ * -2.900,
                           KCAL_TO_KJ * -3.150,
                           KCAL_TO_KJ * -4.8,
                           ANG_TO_NM * 3.5,
                           0.177,
                       }
                },
                {"S", {
                          ANG3_TO_NM3 * 20.703,
                          KCAL_TO_KJ * -3.240,
                          KCAL_TO_KJ * -3.240,
                          KCAL_TO_KJ * -4.475,
                          KCAL_TO_KJ * -39.9,
                          ANG_TO_NM * 3.5,
                          0.20,
                      }
                },
                {"SM", {
                           ANG3_TO_NM3 * 21.306,
                           KCAL_TO_KJ * -3.240,
                           KCAL_TO_KJ * -3.240,
                           KCAL_TO_KJ * -4.475,
                           KCAL_TO_KJ * -39.9,
                           ANG_TO_NM * 3.5,
                           0.197,
                       }
                }
            };
        }
    }
}
#endif
