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
                map<string, map<string, string> > typemap;
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
                        double mlambda =  parameter[i][5];
                        if (parameter[j][5] > mlambda) mlambda = parameter[j][5];
                        if (rij > 3. * mlambda) continue;
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
                            {"CG",  "CPH"},
                            {"CE1", "CPH"},
                            {"HE1", "HR1"},
                            {"NE2", "NR2"},
                            {"CD2", "CPH"},
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
                            {"CG",  "CPH"},
                            {"CE1", "CPH"},
                            {"HE1", "HR1"},
                            {"NE2", "NR1"},
                            {"HE2", "H"  },
                            {"CD2", "CPH"},
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
                            {"CD2", "CPH"},
                            {"HD2", "HR1"},
                            {"CG",  "CPH"},
                            {"NE2", "NR3"},
                            {"HE2", "H"  },
                            {"ND1", "NR3"},
                            {"HD1", "H"  },
                            {"CE1", "CPH"},
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
                            {"O",    "O"  },
                        }
                }
            };
        }

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

                    // Volume ∆Gref ∆Gfree ∆H ∆Cp λ vdw_radius
                    if (Atype == "C") {
                        parameter[i][0] = 0.001 * 14.720;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CD") {
                        parameter[i][0] = 0.001 * 14.720;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CT1") {
                        parameter[i][0] = 0.001 * 11.507;
                        parameter[i][1] = 4.184 * -0.187;
                        parameter[i][2] = 4.184 * -0.187;
                        parameter[i][3] = 0.876;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CT2") {
                        parameter[i][0] = 0.001 * 18.850;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CT2A") {
                        parameter[i][0] = 0.001 * 18.666;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CT3") {
                        parameter[i][0] = 0.001 * 27.941;
                        parameter[i][1] = 4.184 * 1.089;
                        parameter[i][2] = 4.184 * 1.089;
                        parameter[i][3] = -1.779;
                        parameter[i][4] = 35.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CPH1") {
                        parameter[i][0] = 0.001 * 5.275;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.080;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CPH2") {
                        parameter[i][0] = 0.001 * 11.796;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.080;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CPT") {
                        parameter[i][0] = 0.001 * 4.669;
                        parameter[i][1] = 4.184 * -0.890;
                        parameter[i][2] = 4.184 * -0.890;
                        parameter[i][3] = 2.220;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CY") {
                        parameter[i][0] = 0.001 * 10.507;
                        parameter[i][1] = 4.184 * -0.890;
                        parameter[i][2] = 4.184 * -0.890;
                        parameter[i][3] = 2.220;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CP1") {
                        parameter[i][0] = 0.001 * 25.458;
                        parameter[i][1] = 4.184 * -0.187;
                        parameter[i][2] = 4.184 * -0.187;
                        parameter[i][3] = 0.876;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CP2") {
                        parameter[i][0] = 0.001 * 19.880;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CP3") {
                        parameter[i][0] = 0.001 * 26.731;
                        parameter[i][1] = 4.184 * 0.372;
                        parameter[i][2] = 4.184 * 0.372;
                        parameter[i][3] = -0.610;
                        parameter[i][4] = 18.6;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CC") {
                        parameter[i][0] = 0.001 * 16.539;
                        parameter[i][1] = 4.184 * 0.000;
                        parameter[i][2] = 4.184 * 0.000;
                        parameter[i][3] = 0.000;
                        parameter[i][4] = 0.0;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CAI") {
                        parameter[i][0] = 0.001 * 18.249;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.057;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "CA") {
                        parameter[i][0] = 0.001 * 18.249;
                        parameter[i][1] = 4.184 * 0.057;
                        parameter[i][2] = 4.184 * 0.057;
                        parameter[i][3] = -0.973;
                        parameter[i][4] = 6.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.17;
                    } else if (Atype == "N") {
                        parameter[i][0] = 0.001 * 0.000;
                        parameter[i][1] = 4.184 * -1.000;
                        parameter[i][2] = 4.184 * -1.000;
                        parameter[i][3] = -1.250;
                        parameter[i][4] = 8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NR1") {
                        parameter[i][0] = 0.001 * 15.273;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NR2") {
                        parameter[i][0] = 0.001 * 15.111;
                        parameter[i][1] = 4.184 * -3.820;
                        parameter[i][2] = 4.184 * -3.820;
                        parameter[i][3] = -4.654;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NR3") {
                        parameter[i][0] = 0.001 * 15.071;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NH1") {
                        parameter[i][0] = 0.001 * 10.197;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NH2") {
                        parameter[i][0] = 0.001 * 18.182;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NH3") {
                        parameter[i][0] = 0.001 * 18.817;
                        parameter[i][1] = 4.184 * -20.000;
                        parameter[i][2] = 4.184 * -20.000;
                        parameter[i][3] = -25.000;
                        parameter[i][4] = -18.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NC2") {
                        parameter[i][0] = 0.001 * 18.215;
                        parameter[i][1] = 4.184 * -10.000;
                        parameter[i][2] = 4.184 * -10.000;
                        parameter[i][3] = -12.000;
                        parameter[i][4] = -7.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NY") {
                        parameter[i][0] = 0.001 * 12.001;
                        parameter[i][1] = 4.184 * -5.950;
                        parameter[i][2] = 4.184 * -5.950;
                        parameter[i][3] = -9.059;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "NP") {
                        parameter[i][0] = 0.001 * 4.993;
                        parameter[i][1] = 4.184 * -20.000;
                        parameter[i][2] = 4.184 * -20.000;
                        parameter[i][3] = -25.000;
                        parameter[i][4] = -18.0;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.155;
                    } else if (Atype == "O") {
                        parameter[i][0] = 0.001 * 11.772;
                        parameter[i][1] = 4.184 * -5.330;
                        parameter[i][2] = 4.184 * -5.330;
                        parameter[i][3] = -5.787;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Atype == "OB") {
                        parameter[i][0] = 0.001 * 11.694;
                        parameter[i][1] = 4.184 * -5.330;
                        parameter[i][2] = 4.184 * -5.330;
                        parameter[i][3] = -5.787;
                        parameter[i][4] = -8.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Atype == "OC") {
                        parameter[i][0] = 0.001 * 12.003;
                        parameter[i][1] = 4.184 * -10.000;
                        parameter[i][2] = 4.184 * -10.000;
                        parameter[i][3] = -12.000;
                        parameter[i][4] = -9.4;
                        parameter[i][5] = 0.1 * 6.0;
                        parameter[i][6] = 0.152;
                    } else if (Atype == "OH1") {
                        parameter[i][0] = 0.001 * 15.528;
                        parameter[i][1] = 4.184 * -5.920;
                        parameter[i][2] = 4.184 * -5.920;
                        parameter[i][3] = -9.264;
                        parameter[i][4] = -11.2;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Atype == "OS") {
                        parameter[i][0] = 0.001 * 6.774;
                        parameter[i][1] = 4.184 * -2.900;
                        parameter[i][2] = 4.184 * -2.900;
                        parameter[i][3] = -3.150;
                        parameter[i][4] = -4.8;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.152;
                    } else if (Atype == "S") {
                        parameter[i][0] = 0.001 * 20.703;
                        parameter[i][1] = 4.184 * -3.240;
                        parameter[i][2] = 4.184 * -3.240;
                        parameter[i][3] = -4.475;
                        parameter[i][4] = -39.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.18;
                    } else if (Atype == "SM") {
                        parameter[i][0] = 0.001 * 21.306;
                        parameter[i][1] = 4.184 * -3.240;
                        parameter[i][2] = 4.184 * -3.240;
                        parameter[i][3] = -4.475;
                        parameter[i][4] = -39.9;
                        parameter[i][5] = 0.1 * 3.5;
                        parameter[i][6] = 0.18;
                    } else {
                        error("Invalid atom type!\n");
                    }
                }
            } else {
                error("MOLINFO DATA not found\n");
            }
        }
    }
}
