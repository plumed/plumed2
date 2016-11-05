/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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

// TODO:
// delta_g_free temperature correction   [ ]
// try constant vdw_radius approximation [ ]
// dihedral correction                   [x]
// distance dependent dielectric         [x]
// neutralize ionic sidechains           [x]
// implement neighborlist                [x]
// cutoff distance                       [x]

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "core/SetupMolInfo.h"
#include "tools/OpenMP.h"

#define INV_PI_SQRT_PI 0.179587122

using namespace std;

namespace PLMD {
    namespace colvar {

//+PLUMEDOC COLVAR IMPLICIT
/*

Calculate EEF1-SB solvation free energy

*/
//+ENDPLUMEDOC

        class Implicit : public Colvar {
            private:
                bool pbc;
                double cutoff;
                unsigned stride;
                unsigned nl_update;
                vector<vector<unsigned> > nl;
                vector<vector<double> > parameter;
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
            keys.add("compulsory", "NL_CUTOFF", "The cutoff used when calculating pairwise interactions.");
            keys.add("compulsory", "NL_STRIDE", "The frequency with which the neighbourlist is updated.");
        }

        Implicit::Implicit(const ActionOptions&ao):
            PLUMED_COLVAR_INIT(ao),
            pbc(true),
            cutoff(1.0),
            stride(10),
            nl_update(0)
        {
            vector<AtomNumber> atoms;
            parseAtomList("ATOMS", atoms);
            const unsigned size = atoms.size();

            parse("NL_CUTOFF", cutoff);
            parse("NL_STRIDE", stride);

            bool nopbc = !pbc;
            parseFlag("NOPBC", nopbc);
            pbc = !nopbc;

            checkRead();

            nl.resize(size);

            parameter.resize(size);
            for (unsigned i=0; i<size; ++i) {
                parameter[i].resize(8, 0.0);
            }
            setupConstants(atoms, parameter);

            addValueWithDerivatives();
            setNotPeriodic();
            requestAtoms(atoms);
        }

        void Implicit::update_neighb() {
            const double c2 = cutoff*cutoff;
            const unsigned size=getNumberOfAtoms();
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
                        if (d2 < c2 && i != j) {
                            nl[i].push_back(j);
                        }
                    }
                }
            }
        }

        void Implicit::calculate() {
            if(pbc) makeWhole();

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
                        if (rij > cutoff) continue;
                        const double inv_rij = 1.0 / rij;
                        const double inv_rij2 = inv_rij * inv_rij;

                        // i-j interaction
                        {
                            const double delta_g_free = parameter[i][2];
                            const double lambda = parameter[i][5];
                            const double vdw_radius = parameter[i][6];
                            const double inv_lambda = 1.0 / lambda;
                            const double inv_lambda2 = inv_lambda * inv_lambda;
                            const double vdw_volume = parameter[j][0];

                            const double rij_vdwr_diff = rij - vdw_radius;
                            const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                            const double fact = -delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                            const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                            // This is needed for correct box derivs
                            #pragma omp critical(deriv)
                            {
                                fedensity += -fact;
                                fedensity_deriv[i] += deriv * dist;
                                fedensity_deriv[j] -= deriv * dist;
                            }
                        }

                        // j-i interaction
                        {
                            const double delta_g_free = parameter[j][2];
                            const double lambda = parameter[j][5];
                            const double vdw_radius = parameter[j][6];
                            const double inv_lambda = 1.0 / lambda;
                            const double inv_lambda2 = inv_lambda * inv_lambda;
                            const double vdw_volume = parameter[i][0];

                            const double rij_vdwr_diff = rij - vdw_radius;
                            const double expo = exp(-inv_lambda2 * rij_vdwr_diff * rij_vdwr_diff);
                            const double fact = -delta_g_free * vdw_volume * expo * INV_PI_SQRT_PI * inv_rij2 * inv_lambda;
                            const double deriv = inv_rij * fact * (inv_rij + rij_vdwr_diff * inv_lambda2);

                            #pragma omp critical(deriv)
                            {
                                fedensity += -fact;
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
            #pragma omp declare reduction(tensor_sum:Tensor: omp_out += omp_in) initializer(omp_priv = Tensor())
            #pragma omp parallel num_threads(ntd)
            {
                #pragma omp for reduction(tensor_sum:deriv_box)
                for (unsigned i=0; i<size; ++i) {
                    setAtomsDerivatives(i, fedensity_deriv[i]);
                    deriv_box += Tensor(getPosition(i), fedensity_deriv[i]);
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

        void Implicit::setupConstants(const vector<AtomNumber> &atoms, vector<vector<double> > &parameter) {
            vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
            if (moldat.size() == 1) {
                log << "  MOLINFO DATA found, using proper atom names\n";
                for(unsigned i=0; i<atoms.size(); ++i) {
                    string Aname = moldat[0]->getAtomName(atoms[i]);
                    string Rname = moldat[0]->getResidueName(atoms[i]);
                    // Volume ∆Gref ∆Gfree ∆H ∆Cp λ vdw_radius
                    if (Aname == "C") {
                        parameter[i][0] = 0.001 * 14.720;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CD") {
                        parameter[i][0] = 0.001 * 14.720;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CT1") {
                        parameter[i][0] = 0.001 * 11.507;
                        parameter[i][1] = 4.184 * -0.187;
                        parameter[i][2] = 4.184 * -0.187;
                        parameter[i][3] = 0.876;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CT2") {
                        parameter[i][0] = 0.001 * 18.850;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CT2A") {
                        parameter[i][0] = 0.001 * 18.666;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CT3") {
                        parameter[i][0] = 0.001 * 27.941;
                        parameter[i][1] = 4.184 * 1.089;
                        parameter[i][2] = 4.184 * 1.089;
                        parameter[i][3] = -1.779;
                        parameter[i][4] = 35.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CPH1") {
                        parameter[i][0] = 0.001 * 5.275;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.080;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CPH2") {
                        parameter[i][0] = 0.001 * 11.796;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.080;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CPT") {
                        parameter[i][0] = 0.001 * 4.669;
                        parameter[i][1] = 4.184 * -0.890;
                        parameter[i][2] = 4.184 * -0.890;
                        parameter[i][3] = 2.220;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CY") {
                        parameter[i][0] = 0.001 * 10.507;
                        parameter[i][1] = 4.184 * -0.890;
                        parameter[i][2] = 4.184 * -0.890;
                        parameter[i][3] = 2.220;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CP1") {
                        parameter[i][0] = 0.001 * 25.458;
                        parameter[i][1] = 4.184 * -0.187;
                        parameter[i][2] = 4.184 * -0.187;
                        parameter[i][3] = 0.876;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CP2") {
                        parameter[i][0] = 0.001 * 19.880;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CP3") {
                        parameter[i][0] = 0.001 * 26.731;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CC") {
                        parameter[i][0] = 0.001 * 16.539;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CAI") {
                        parameter[i][0] = 0.001 * 18.249;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.057;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "CA") {
                        parameter[i][0] = 0.001 * 18.249;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.057;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Aname == "N") {
                        parameter[i][0] = 0.001 * 0.000;
                        parameter[i][1] = 4.184 * -1.000;
                        parameter[i][2] = 4.184 * -1.000;
                        parameter[i][3] = -1.250;
                        parameter[i][4] = 8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NR1") {
                        parameter[i][0] = 0.001 * 15.273;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NR2") {
                        parameter[i][0] = 0.001 * 15.111;
                        parameter[i][1] = 4.184 * -3.820;
                        parameter[i][2] = 4.184 * -3.820;
                        parameter[i][3] = -4.654;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NR3") {
                        parameter[i][0] = 0.001 * 15.071;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NH1") {
                        parameter[i][0] = 0.001 * 10.197;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NH2") {
                        parameter[i][0] = 0.001 * 18.182;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NH3") {
                        parameter[i][0] = 0.001 * 18.817;
                        parameter[i][1] = 4.184 * -20.000;
                        parameter[i][2] = 4.184 * -20.000;
                        parameter[i][3] = -25.000;
                        parameter[i][4] = -18.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NC2") {
                        parameter[i][0] = 0.001 * 18.215;
                        parameter[i][1] = 4.184 * -10.000;
                        parameter[i][2] = 4.184 * -10.000;
                        parameter[i][3] = -12.000;
                        parameter[i][4] = -7.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NY") {
                        parameter[i][0] = 0.001 * 12.001;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "NP") {
                        parameter[i][0] = 0.001 * 4.993;
                        parameter[i][1] = 4.184 * -20.000;
                        parameter[i][2] = 4.184 * -20.000;
                        parameter[i][3] = -25.000;
                        parameter[i][4] = -18.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Aname == "O") {
                        parameter[i][0] = 0.001 * 11.772;
                        parameter[i][1] = 4.184 * -5.330;
                        parameter[i][2] = 4.184 * -5.330;
                        parameter[i][3] = -5.787;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Aname == "OB") {
                        parameter[i][0] = 0.001 * 11.694;
                        parameter[i][1] = 4.184 * -5.330;
                        parameter[i][2] = 4.184 * -5.330;
                        parameter[i][3] = -5.787;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Aname == "OC") {
                        parameter[i][0] = 0.001 * 12.003;
                        parameter[i][1] = 4.184 * -10.000;
                        parameter[i][2] = 4.184 * -10.000;
                        parameter[i][3] = -12.000;
                        parameter[i][4] = -9.4;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.152;
                    } else if (Aname == "OH1") {
                        parameter[i][0] = 0.001 * 15.528;
                        parameter[i][1] = 4.184 * -5.920;
                        parameter[i][2] = 4.184 * -5.920;
                        parameter[i][3] = -9.264;
                        parameter[i][4] = -11.2;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Aname == "OS") {
                        parameter[i][0] = 0.001 * 6.774;
                        parameter[i][1] = 4.184 * -2.900;
                        parameter[i][2] = 4.184 * -2.900;
                        parameter[i][3] = -3.150;
                        parameter[i][4] = -4.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Aname == "S") {
                        parameter[i][0] = 0.001 * 20.703;
                        parameter[i][1] = 4.184 * -3.240;
                        parameter[i][2] = 4.184 * -3.240;
                        parameter[i][3] = -4.475;
                        parameter[i][4] = -39.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.18;
                    } else if (Aname == "SM") {
                        parameter[i][0] = 0.001 * 21.306;
                        parameter[i][1] = 4.184 * -3.240;
                        parameter[i][2] = 4.184 * -3.240;
                        parameter[i][3] = -4.475;
                        parameter[i][4] = -39.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.18;
                    } else {
                        parameter[i][0] = 0.0;
                        parameter[i][1] = 0.0;
                        parameter[i][2] = 0.0;
                        parameter[i][3] = 0.0;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5; // We need to avoid division by zero
                        parameter[i][6] = 0.0;
                    }

                    // Charge corrections
                    if (Rname == "ARG") {
                        if (Aname == "CD") {          
                            parameter[i][7] = -0.30;
                        } else if (Aname == "HD1" || Aname == "HD2") {     
                            parameter[i][7] = 0.05;
                        } else if (Aname == "NE") {          
                            parameter[i][7] = -0.28;
                        } else if (Aname == "HE") {          
                            parameter[i][7] = 0.12;
                        } else if (Aname == "CZ") {          
                            parameter[i][7] = -0.20;
                        } else if (Aname == "NH1" || Aname == "NH2") {     
                            parameter[i][7] = -0.121;
                        } else if (Aname == "HH1" || Aname == "HH2") {     
                            parameter[i][7] = 0.2005;
                        }
                    } else if (Rname == "HSP") {
                        if (Aname == "CB") {          
                            parameter[i][7] = -0.10;
                        } else if (Aname == "HB1" || Aname == "HB2") {     
                            parameter[i][7] = 0.05;
                        } else if (Aname == "CD2") {         
                            parameter[i][7] = 0.05;
                        } else if (Aname == "HD2") {         
                            parameter[i][7] = 0.00;
                        } else if (Aname == "CG") {          
                            parameter[i][7] = 0.05;
                        } else if (Aname == "NE2" || Aname == "ND1") {     
                            parameter[i][7] = -0.55;
                        } else if (Aname == "HE2" || Aname == "HD1") {     
                            parameter[i][7] = 0.45;
                        } else if (Aname == "CE1") {         
                            parameter[i][7] = 0.10;
                        } else if (Aname == "HE1") {         
                            parameter[i][7] = 0.00;
                        }
                    } else if (Rname == "ASP") {
                        if (Aname == "CB") {          
                            parameter[i][7] = -0.28;
                        } else if (Aname == "HB1" || Aname == "HB2") {     
                            parameter[i][7] = 0.14;
                        } else if (Aname == "HB1" || Aname == "HB2") {     
                            parameter[i][7] = 0.14;
                        } else if (Aname == "HB1" || Aname == "HB2") {     
                            parameter[i][7] = 0.14;
                        }
                    } else if (Rname == "LYS") {
                        if (Aname == "CE") {          
                            parameter[i][7] = 0.00;
                        } else if (Aname == "HE1" || Aname == "HE2") {     
                            parameter[i][7] = 0.00;
                        } else if (Aname == "NZ") {          
                            parameter[i][7] = -0.90;
                        } else if (Aname == "HZ1" || Aname == "HZ2" || Aname == "HZ3") { 
                            parameter[i][7] = 0.30;
                        }
                    } else if (Rname == "GLU") {
                        if (Aname == "CG") {          
                            parameter[i][7] = -0.28;
                        } else if (Aname == "HG1" || Aname == "HG2") {     
                            parameter[i][7] = 0.14;
                        } else if (Aname == "CD") {          
                            parameter[i][7] = 1.00;
                        } else if (Aname == "OE1" || Aname == "OE2") {     
                            parameter[i][7] = -0.50;
                        }
                    } else if (Rname == "NTER") {
                        if (Aname == "N") {           
                            parameter[i][7] = -0.90;
                        } else if (Aname == "HT1" || Aname == "HT2" || Aname == "HT3") { 
                            parameter[i][7] = 0.20;
                        } else if (Aname == "HA") {          
                            parameter[i][7] = 0.10;
                        } else if (Aname == "CA") {          
                            parameter[i][7] = 0.20;
                        }
                    } else if (Rname == "GLP") {
                        if (Aname == "CG") {          
                            parameter[i][7] = -0.21;
                        } else if (Aname == "HG1" || Aname == "HG2") {     
                            parameter[i][7] = 0.09;
                        } else if (Aname == "CD") {          
                            parameter[i][7] = 0.75;
                        } else if (Aname == "OE1") {         
                            parameter[i][7] = -0.55;
                        } else if (Aname == "OE2") {         
                            parameter[i][7] = -0.71;
                        } else if (Aname == "HE2") {         
                            parameter[i][7] = 0.44;
                        }
                    } else if (Rname == "CTER") {
                        if (Aname == "C") {           
                            parameter[i][7] = 1.00;
                        } else if (Aname == "OT1" || Aname == "OT2") {     
                            parameter[i][7] = -0.50;
                        }
                    }
                }
            } else {
                error("MOLINFO DATA not found\n");
            }
        }
    }
}
