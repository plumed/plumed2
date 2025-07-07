/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) envsim 2023-2024 The code team
   (see the PEOPLE-envsim file at the root of the distribution for a list of names)

   This file is part of envsim code module.

   The envsim code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The envsim code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the envsim code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/PDB.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace envsim {

// N.B. In this equation $\tilde{k}_ {\chi_0}(\chi_0)=1$ space after underscore ensures correct rendering.
// I don't know why GAT

//+PLUMEDOC MCOLVAR ENVIRONMENTSIMILARITY
/*
Measure how similar the environment around atoms is to that found in some reference crystal structure.

The starting point for the definition of the CV is the local atomic density around an atom.
We consider an environment $\chi$ around this atom and we define the density by

$$
 \rho_{\chi}(\mathbf{r})=\sum\limits_{i\in\chi} \exp\left(- \frac{|r_i-r|^2} {2\sigma^2} \right),
$$

where $i$ runs over the neighbors in the environment $\chi$, $\sigma$ is a broadening parameter, and $r_i$ are the
coordinates of the neighbors relative to the central atom.
We now define a reference environment or template $\chi_0$ that contains $n$ reference positions $\{r^0_1,...,r^0_n\}$
that describe, for instance, the nearest neighbors in a given lattice.
$\sigma$ is set using the SIGMA keyword and $\chi_0$ is chosen with the CRYSTAL_STRUCTURE keyword.
If only the SPECIES keyword is given then the atoms defined there will be the central and neighboring atoms.
If instead the SPECIESA and SPECIESB keywords are given then SPECIESA determines the central atoms and SPECIESB the neighbors.

The environments $\chi$ and $\chi_0$ are compared using the kernel,

$$
 k_{\chi_0}(\chi)= \int d\mathbf{r} \rho_{\chi}(\mathbf{r}) \rho_{\chi_0}(\mathbf{r}) .
$$

Combining the two equations above and performing the integration analytically we obtain,

$$
 k_{\chi_0}(\chi)= \sum\limits_{i\in\chi} \sum\limits_{j\in\chi_0} \pi^{3/2} \sigma^3  \exp\left(- \frac{|\mathbf{r}_i-\mathbf{r}^0_j|^2} {4\sigma^2} \right).
$$

The kernel is finally normalized,

$$
 \tilde{k}_{\chi_0}(\chi)  = \frac{1}{n} \sum\limits_{i\in\chi} \sum\limits_{j\in\chi_0} \exp\left( - \frac{|\mathbf{r}_i-\mathbf{r}^0_j|^2} {4\sigma^2} \right),
$$

such that $\tilde{k}_ {\chi_0}(\chi_0)=1$.
The above kernel is computed for each atom in the SPECIES or SPECIESA keywords.
This quantity is a multicolvar so you can compute it for multiple atoms using a single PLUMED action and then compute
the average value for the atoms in your system, the number of atoms that have an $\tilde{k}_{\chi_0}$ value that is more that some target and
so on.

The kernel can be generalized to crystal structures described as a lattice with a basis of more than one atom.
In this case there is more than one type of environment.
We consider the case of $M$ environments $X = \chi_1,\chi_2,...,\chi_M$ and we define the kernel through a best match strategy:


$$
 \tilde{k}_X(\chi)= \frac{1}{\lambda} \log \left ( \sum\limits_{l=1}^{M}\exp \left (\lambda \: \tilde{k}_{\chi_l}(\chi) \right ) \right ).
$$

For a large enough $\lambda$ this expression will select the largest $\tilde{k}_{\chi_l}(\chi)$ with $\chi_l \in X$.
This approach can be used, for instance, to target the hexagonal closed packed (HCP keyword) or the diamond structure (DIAMOND keyword).

The CRYSTAL_STRUCTURE keyword can take the values SC (simple cubic), BCC (body centered cubic), FCC (face centered cubic),
HCP (hexagonal closed pack), DIAMOND (cubic diamond), and CUSTOM (user defined).
All options follow the same conventions as in the [lattice command](https://lammps.sandia.gov/doc/lattice.html) of [LAMMPS](https://lammps.sandia.gov/).
If a CRYSTAL_STRUCTURE other than CUSTOM is used, then the lattice constants have to be specified using the keyword LATTICE_CONSTANTS.
One value has to be specified for SC, BCC, FCC, and DIAMOND and two values have to be set for HCP (a and c lattice constants in that order).

If the CUSTOM option is used then the reference environments have to be specified by the user.
The reference environments are specified in pdb files containing the distance vectors from the central atom to the neighbors.
Make sure your PDB file is correctly formatted as explained in the documenation for [MOLINFO](MOLINFO.md)
If only one reference environment is specified then the filename should be given as argument of the keyword REFERENCE.
If instead several reference environments are given, then they have to be provided in separate pdb files and given as arguments for the
keywords REFERENCE_1, REFERENCE_2, etc.
If you have a reference crystal structure configuration you can use the [Environment Finder](https://github.com/PabloPiaggi/EnvironmentFinder) app to determine the reference environments that you should use.

If multiple chemical species are involved in the calculation, it is possible to provide the atom types (names) both for atoms in the reference environments and in the simulation box.
This information is provided in pdb files using the atom name field.
The comparison between environments is performed taking into account whether the atom names match.

### Examples

The following input calculates the ENVIRONMENTSIMILARITY kernel for 250 atoms in the system
using the BCC atomic environment as target, and then calculates and prints the average value
 for this quantity.

```plumed
es: ENVIRONMENTSIMILARITY SPECIES=1-250 SIGMA=0.05 LATTICE_CONSTANTS=0.423 CRYSTAL_STRUCTURE=BCC MEAN

PRINT ARG=es.mean FILE=COLVAR
```

If you want to use a different set of atoms in the environments to the atoms for which you are calculating
the ENVIRONMENTSIMILARITY kernel you use the SPECIESA and SPECIESB keywords as shown below:

```plumed
es: ENVIRONMENTSIMILARITY SPECIESA=1-100 SPECIESB=101-250 SIGMA=0.05 LATTICE_CONSTANTS=0.423 CRYSTAL_STRUCTURE=BCC MEAN

PRINT ARG=es.mean FILE=COLVAR
```

The next example compares the environments of the 96 selected atoms with a user specified reference
environment. The reference environment is contained in the env1.pdb file. Once the kernel is computed
 the average and the number of atoms with a kernel larger than 0.5 are computed.

```plumed
#SETTINGS INPUTFILES=regtest/envsim/rt-env-sim-atom-names-match/env1.pdb
es: ENVIRONMENTSIMILARITY ...
 SPECIES=1-288:3
 SIGMA=0.05
 CRYSTAL_STRUCTURE=CUSTOM
 REFERENCE=regtest/envsim/rt-env-sim-atom-names-match/env1.pdb
 MEAN
 MORE_THAN={RATIONAL R_0=0.5 NN=12 MM=24}
...

PRINT ARG=es.mean,es.morethan FILE=COLVAR
```

The next example is similar to the one above but in this case 4 reference environments are specified.
 Each reference environment is given in a separate pdb file.

```plumed
#SETTINGS INPUTFILES=regtest/envsim/rt-env-sim-atom-names-match/env1.pdb,regtest/envsim/rt-env-sim-atom-names-match/env2.pdb,regtest/envsim/rt-env-sim-atom-names-match/env3.pdb,regtest/envsim/rt-env-sim-atom-names-match/env4.pdb
es: ENVIRONMENTSIMILARITY ...
 SPECIES=1-288:3
 SIGMA=0.05
 CRYSTAL_STRUCTURE=CUSTOM
 REFERENCE_1=regtest/envsim/rt-env-sim-atom-names-match/env1.pdb
 REFERENCE_2=regtest/envsim/rt-env-sim-atom-names-match/env2.pdb
 REFERENCE_3=regtest/envsim/rt-env-sim-atom-names-match/env3.pdb
 REFERENCE_4=regtest/envsim/rt-env-sim-atom-names-match/env4.pdb
 MEAN
 MORE_THAN={RATIONAL R_0=0.5 NN=12 MM=24}
...

PRINT ARG=es.mean,es.morethan FILE=COLVAR
```

The following examples illustrates the use of pdb files to provide information about different chemical species:


```plumed
#SETTINGS INPUTFILES=regtest/envsim/rt-env-sim-custom-1env/env1.pdb,regtest/envsim/rt-env-sim-atom-names-match/IceIh-atom-names.pdb
es: ENVIRONMENTSIMILARITY ...
 SPECIES=1-384
 SIGMA=0.05
 CRYSTAL_STRUCTURE=CUSTOM
 REFERENCE=regtest/envsim/rt-env-sim-custom-1env/env1.pdb
 MEAN
 MORE_THAN={RATIONAL R_0=0.5 NN=12 MM=24}
 ATOM_NAMES_FILE=regtest/envsim/rt-env-sim-atom-names-match/IceIh-atom-names.pdb
...
```

In this case, all atoms are used as centers, but only neighbors of type O are taken into account.

*/
//+ENDPLUMEDOC

class EnvironmentSimilarity : public ActionShortcut {
private:
  std::vector<std::pair<unsigned,Vector> > getReferenceEnvironment( const PDB& pdb, const std::vector<std::string>& anames, double& maxdist );
public:
  static void registerKeywords( Keywords& keys );
  explicit EnvironmentSimilarity(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(EnvironmentSimilarity,"ENVIRONMENTSIMILARITY")

void EnvironmentSimilarity::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms-3","SPECIES","this keyword is used for colvars such as coordination number. In that context it specifies that plumed should calculate "
           "one coordination number for each of the atoms specified.  Each of these coordination numbers specifies how many of the "
           "other specified atoms are within a certain cutoff of the central atom.  You can specify the atoms here as another multicolvar "
           "action or using a MultiColvarFilter or ActionVolume action.  When you do so the quantity is calculated for those atoms specified "
           "in the previous multicolvar.  This is useful if you would like to calculate the Steinhardt parameter for those atoms that have a "
           "coordination number more than four for example");
  keys.add("atoms-4","SPECIESA","this keyword is used for colvars such as the coordination number.  In that context it species that plumed should calculate "
           "one coordination number for each of the atoms specified in SPECIESA.  Each of these cooordination numbers specifies how many "
           "of the atoms specifies using SPECIESB is within the specified cutoff.  As with the species keyword the input can also be specified "
           "using the label of another multicolvar");
  keys.add("atoms-4","SPECIESB","this keyword is used for colvars such as the coordination number.  It must appear with SPECIESA.  For a full explanation see "
           "the documentation for that keyword");
  keys.add("compulsory","CRYSTAL_STRUCTURE","FCC","Targeted crystal structure. Options are: "
           "SC: simple cubic, "
           "BCC: body center cubic, "
           "FCC: face centered cubic, "
           "HCP: hexagonal closed pack, "
           "DIAMOND: cubic diamond, "
           "CUSTOM: user defined "
           " ");
  keys.add("compulsory","LATTICE_CONSTANTS","Lattice constants. Two comma separated values for HCP, "
           "one value for all other CRYSTAL_STRUCTURES.");
  keys.add("compulsory","SIGMA","0.1","the width to use for the gaussian kernels");
  keys.add("compulsory","LCUTOFF","0.0001","any atoms separated by less than this tolerance should be ignored");
  keys.add("optional","REFERENCE","PDB files with relative distances from central atom.  Use this keyword if you are targeting a single reference environment.");
  keys.add("numbered","REFERENCE_","PDB files with relative distances from central atom. Each file corresponds to one template. Use these keywords if you are targeting more than one reference environment.");
  keys.add("compulsory","LAMBDA","100","Lambda parameter.  This is only used if you have more than one reference environment");
  keys.add("compulsory","CUTOFF","3","how many multiples of sigma would you like to consider beyond the maximum distance in the environment");
  keys.add("optional","ATOM_NAMES_FILE","PDB file with atom names for all atoms in SPECIES. Atoms in reference environments will be compared only if atom names match.");
  keys.setValueDescription("vector","the environmental similar parameter for each of the input atoms");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("GROUP");
  keys.needsAction("DISTANCE_MATRIX");
  keys.needsAction("ONES");
  keys.needsAction("CONSTANT");
  keys.needsAction("CUSTOM");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("COMBINE");
  keys.addDOI("10.1063/1.5102104");
}

EnvironmentSimilarity::EnvironmentSimilarity(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string atomNamesFile;
  parse("ATOM_NAMES_FILE",atomNamesFile);
  PDB atomnamepdb;
  if( !atomNamesFile.empty() && !atomnamepdb.read(atomNamesFile,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
    error("missing input file " + atomNamesFile);
  }

  double maxdist=0;
  std::vector<std::string> allspec(1);
  std::string crystal_structure;
  parse("CRYSTAL_STRUCTURE", crystal_structure);
  std::vector<std::vector<std::pair<unsigned,Vector> > > environments;
  if( crystal_structure=="CUSTOM" ) {
    if( !atomNamesFile.empty()  ) {
      allspec[0]=atomnamepdb.getAtomName(atomnamepdb.getAtomNumbers()[0]);
      unsigned natoms=atomnamepdb.getPositions().size();
      for(unsigned i=0; i<natoms; ++i) {
        bool found=false;
        for(unsigned j=0; j<allspec.size(); ++j) {
          if( allspec[j]==atomnamepdb.getAtomName(atomnamepdb.getAtomNumbers()[i] ) ) {
            found=true;
            break;
          }
        }
        if( !found ) {
          allspec.push_back( atomnamepdb.getAtomName(atomnamepdb.getAtomNumbers()[i]) );
        }
      }
    }
    std::string reffile;
    parse("REFERENCE",reffile);
    if( reffile.length()>0 ) {
      PDB pdb;
      pdb.read(reffile,plumed.usingNaturalUnits(),0.1/plumed.getUnits().getLength());
      environments.push_back( getReferenceEnvironment( pdb, allspec, maxdist ) );
      log.printf("  reading %d reference vectors from %s \n", environments[0].size(), reffile.c_str() );
    } else {
      for(unsigned int i=1;; i++) {
        PDB pdb;
        if( !parseNumbered("REFERENCE_",i,reffile) ) {
          break;
        }
        if( !pdb.read(reffile,usingNaturalUnits(),0.1/getUnits().getLength()) ) {
          error("missing input file " + reffile );
        }
        environments.push_back( getReferenceEnvironment( pdb, allspec, maxdist ) );
        log.printf("  Reference environment %d : reading %d reference vectors from %s \n", i, environments[i-1].size(), reffile.c_str() );
      }
    }
  } else {
    std::vector<double> lattice_constants;
    parseVector("LATTICE_CONSTANTS", lattice_constants);
    if (crystal_structure == "FCC") {
      if (lattice_constants.size() != 1) {
        error("Number of LATTICE_CONSTANTS arguments must be one for FCC");
      }
      environments.resize(1);
      environments[0].resize(12);
      environments[0][0]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+0.5,+0.0)*lattice_constants[0] );
      environments[0][1]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-0.5,+0.0)*lattice_constants[0] );
      environments[0][2]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-0.5,+0.0)*lattice_constants[0] );
      environments[0][3]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+0.5,+0.0)*lattice_constants[0] );
      environments[0][4]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+0.0,+0.5)*lattice_constants[0] );
      environments[0][5]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+0.0,-0.5)*lattice_constants[0] );
      environments[0][6]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+0.0,+0.5)*lattice_constants[0] );
      environments[0][7]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+0.0,-0.5)*lattice_constants[0] );
      environments[0][8]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.5,+0.5)*lattice_constants[0] );
      environments[0][9]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,-0.5,-0.5)*lattice_constants[0] );
      environments[0][10] = std::pair<unsigned,Vector>( 0, Vector(+0.0,-0.5,+0.5)*lattice_constants[0] );
      environments[0][11] = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.5,-0.5)*lattice_constants[0] );
      maxdist = std::sqrt(2)*lattice_constants[0]/2.;
    } else if (crystal_structure == "SC") {
      if (lattice_constants.size() != 1) {
        error("Number of LATTICE_CONSTANTS arguments must be one for SC");
      }
      environments.resize(1);
      environments[0].resize(6);
      environments[0][0]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+0.0,+0.0)*lattice_constants[0] );
      environments[0][1]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,+0.0,+0.0)*lattice_constants[0] );
      environments[0][2]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+1.0,+0.0)*lattice_constants[0] );
      environments[0][3]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,-1.0,+0.0)*lattice_constants[0] );
      environments[0][4]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.0,+1.0)*lattice_constants[0] );
      environments[0][5]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.0,-1.0)*lattice_constants[0] );
      maxdist = lattice_constants[0];
    } else if( crystal_structure == "BCC") {
      if (lattice_constants.size() != 1) {
        error("Number of LATTICE_CONSTANTS arguments must be one for BCC");
      }
      environments.resize(1);
      environments[0].resize(14);
      environments[0][0]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+0.5,+0.5)*lattice_constants[0] );
      environments[0][1]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-0.5,-0.5)*lattice_constants[0] );
      environments[0][2]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+0.5,+0.5)*lattice_constants[0] );
      environments[0][3]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-0.5,+0.5)*lattice_constants[0] );
      environments[0][4]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+0.5,-0.5)*lattice_constants[0] );
      environments[0][5]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-0.5,+0.5)*lattice_constants[0] );
      environments[0][6]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-0.5,-0.5)*lattice_constants[0] );
      environments[0][7]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+0.5,-0.5)*lattice_constants[0] );
      environments[0][8]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+0.0,+0.0)*lattice_constants[0] );
      environments[0][9]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+1.0,+0.0)*lattice_constants[0] );
      environments[0][10] = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.0,+1.0)*lattice_constants[0] );
      environments[0][11] = std::pair<unsigned,Vector>( 0, Vector(-1.0,+0.0,+0.0)*lattice_constants[0] );
      environments[0][12] = std::pair<unsigned,Vector>( 0, Vector(+0.0,-1.0,+0.0)*lattice_constants[0] );
      environments[0][13] = std::pair<unsigned,Vector>( 0, Vector(+0.0,+0.0,-1.0)*lattice_constants[0] );
      maxdist = lattice_constants[0];
    } else if (crystal_structure == "HCP") {
      if (lattice_constants.size() != 2) {
        error("Number of LATTICE_CONSTANTS arguments must be two for HCP");
      }
      environments.resize(2);
      environments[0].resize(12);
      environments[1].resize(12);
      double sqrt3=std::sqrt(3);
      environments[0][0]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[0][1]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[0][2]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[0][3]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[0][4]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+0.0,+0.0)      *lattice_constants[0] );
      environments[0][5]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,+0.0,+0.0)      *lattice_constants[0] );
      environments[0][6]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[0][7]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[0][8]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,-sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[0][9]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      environments[0][10] = std::pair<unsigned,Vector>( 0, Vector(-0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      environments[0][11] = std::pair<unsigned,Vector>( 0, Vector(+0.0,-sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      environments[1][0]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,+sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[1][1]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,+sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[1][2]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[1][3]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-sqrt3/2.0,+0.0)*lattice_constants[0] );
      environments[1][4]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+0.0,+0.0)      *lattice_constants[0] );
      environments[1][5]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,+0.0,+0.0)      *lattice_constants[0] );
      environments[1][6]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[1][7]  = std::pair<unsigned,Vector>( 0, Vector(-0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[1][8]  = std::pair<unsigned,Vector>( 0, Vector(+0.0,+sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1] );
      environments[1][9]  = std::pair<unsigned,Vector>( 0, Vector(+0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      environments[1][10] = std::pair<unsigned,Vector>( 0, Vector(-0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      environments[1][11] = std::pair<unsigned,Vector>( 0, Vector(+0.0,+sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1] );
      maxdist = lattice_constants[0];
    } else if (crystal_structure == "DIAMOND") {
      if (lattice_constants.size() != 1) {
        error("Number of LATTICE_CONSTANTS arguments must be one for DIAMOND");
      }
      environments.resize(2);
      environments[0].resize(4);
      environments[1].resize(4);
      environments[0][0]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+1.0,+1.0)*lattice_constants[0]/4.0 );
      environments[0][1]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,-1.0,+1.0)*lattice_constants[0]/4.0 );
      environments[0][2]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,-1.0,-1.0)*lattice_constants[0]/4.0 );
      environments[0][3]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,+1.0,-1.0)*lattice_constants[0]/4.0 );
      environments[1][0]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,-1.0,+1.0)*lattice_constants[0]/4.0 );
      environments[1][1]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,+1.0,+1.0)*lattice_constants[0]/4.0 );
      environments[1][2]  = std::pair<unsigned,Vector>( 0, Vector(+1.0,+1.0,-1.0)*lattice_constants[0]/4.0 );
      environments[1][3]  = std::pair<unsigned,Vector>( 0, Vector(-1.0,-1.0,-1.0)*lattice_constants[0]/4.0 );
      maxdist = std::sqrt(3)*lattice_constants[0]/4.0;
    } else {
      error( crystal_structure + " is not a valid input for keyword CRYSTAL_STRUCTURE");
    }
  }
  std::string matlab = getShortcutLabel() + "_cmat";
  double cutoff, sig;
  parse("SIGMA",sig);
  parse("CUTOFF",cutoff);
  std::string lcutoff;
  parse("LCUTOFF",lcutoff);
  std::string sig2;
  Tools::convert( sig*sig, sig2 );
  std::vector<std::vector<std::string> > funcstr(environments.size());
  std::string str_cutoff;
  Tools::convert( maxdist + cutoff*sig, str_cutoff );
  std::string str_natoms, xpos, ypos, zpos;
  Tools::convert( environments[0].size(), str_natoms );
  for(unsigned j=0; j<environments.size(); ++j) {
    funcstr[j].resize( allspec.size() );
    for(unsigned k=0; k<allspec.size(); ++k) {
      for(unsigned i=0; i<environments[j].size(); ++i) {
        if( environments[j][i].first!=k ) {
          continue ;
        }
        Tools::convert( environments[j][i].second[0], xpos );
        Tools::convert( environments[j][i].second[1], ypos );
        Tools::convert( environments[j][i].second[2], zpos );
        if( i==0 ) {
          funcstr[j][k] = "FUNC=(step(w-" + lcutoff + ")*step(" + str_cutoff + "-w)/" + str_natoms + ")*(exp(-((x-" + xpos + ")^2+(y-" + ypos + ")^2+(z-" + zpos + ")^2)/(4*" + sig2 + "))";
        } else {
          funcstr[j][k] += "+exp(-((x-" + xpos + ")^2+(y-" + ypos + ")^2+(z-" + zpos + ")^2)/(4*" + sig2 + "))";
        }
      }
      if( funcstr[j][k].length()>0 ) {
        funcstr[j][k] += ")";
      } else {
        funcstr[j][k] ="FUNC=0";
      }
    }
  }

  // Create the constact matrix
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  if( sp_str.length()>0 ) {
    readInputLine( matlab + ": DISTANCE_MATRIX COMPONENTS GROUP=" + sp_str + " CUTOFF=" + str_cutoff );
    readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + sp_str );
  } else {
    if( specA.length()==0 ) {
      error("no atoms were specified use SPECIES or SPECIESA+SPECIESB");
    }
    if( specB.length()==0 ) {
      error("no atoms were specified for SPECIESB");
    }
    readInputLine( matlab + ": DISTANCE_MATRIX COMPONENTS GROUPA=" + specA + " GROUPB=" + specB + " CUTOFF=" + str_cutoff );
    readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + specA );
  }

  // Make a vector containing all ones
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( matlab );
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  if( allspec.size()==1 ) {
    readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  } else {
    unsigned natoms=atomnamepdb.getPositions().size();
    unsigned firstneigh=0;
    if( sp_str.length()==0 ) {
      firstneigh = (av->copyOutput(0))->getShape()[0];
    }
    for(unsigned i=0; i<allspec.size(); ++i) {
      std::string onesstr="0";
      if( atomnamepdb.getAtomName(atomnamepdb.getAtomNumbers()[firstneigh])==allspec[i] ) {
        onesstr = "1";
      }
      for(unsigned j=firstneigh+1; j<natoms; ++j) {
        if( atomnamepdb.getAtomName(atomnamepdb.getAtomNumbers()[j])==allspec[i] ) {
          onesstr += ",1";
        } else {
          onesstr += ",0";
        }
      }
      readInputLine( getShortcutLabel() + "_ones_" + allspec[i] + ": CONSTANT VALUES=" + onesstr );
    }
  }

  std::string envargstr,varstr, maxfuncstr, lambda;
  if( funcstr.size()>1 ) {
    parse("LAMBDA",lambda);
  }
  // And now do the funcstr bit
  for(unsigned j=0; j<funcstr.size(); ++j) {
    std::string jnum;
    Tools::convert( j+1, jnum );
    if(j==0) {
      varstr = "v" + jnum;
      maxfuncstr = "(1/" + lambda + ")*log(exp(" + lambda + "*v1)";
      envargstr = getShortcutLabel() + "_env" + jnum;
    } else {
      varstr += ",v" + jnum;
      maxfuncstr += "+exp(" + lambda + "*v" + jnum + ")";
      envargstr += "," + getShortcutLabel() + "_env" + jnum;
    }
    // And coordination numbers
    if( allspec.size()>1 ) {
      std::string argnames;
      for(unsigned i=0; i<allspec.size(); ++i) {
        readInputLine( getShortcutLabel() + "_" + allspec[i] + "_matenv" + jnum + ": CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z," + matlab + ".w VAR=x,y,z,w PERIODIC=NO " + funcstr[j][i] );
        readInputLine( getShortcutLabel() + "_" + allspec[i] + "_env" + jnum + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_" + allspec[i] + "_matenv" + jnum + "," + getShortcutLabel() + "_ones_" + allspec[i] );
        if( i==0 ) {
          argnames = getShortcutLabel() + "_" + allspec[i] + "_env" + jnum;
        } else {
          argnames += "," + getShortcutLabel() + "_" + allspec[i] + "_env" + jnum;
        }
      }
      if( funcstr.size()==1) {
        readInputLine( getShortcutLabel() + ": COMBINE PERIODIC=NO ARG=" + argnames );
      } else {
        readInputLine( getShortcutLabel() + "_env" + jnum + ": COMBINE PERIODIC=NO ARG=" + argnames );
      }
    } else {
      readInputLine( getShortcutLabel() + "_matenv" + jnum + ": CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z," + matlab + ".w VAR=x,y,z,w PERIODIC=NO " + funcstr[j][0] );
      if( funcstr.size()==1) {
        readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_matenv" + jnum + "," + getShortcutLabel() + "_ones");
      } else {
        readInputLine( getShortcutLabel() + "_env" + jnum + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_matenv" + jnum + "," + getShortcutLabel() + "_ones");
      }
    }
  }
  // And get the maximum
  if( funcstr.size()>1 ) {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + envargstr + " PERIODIC=NO VAR=" + varstr + " FUNC=" + maxfuncstr + ")" );
  }
  // Read in all the shortcut stuff
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", keymap, this );
}

std::vector<std::pair<unsigned,Vector> > EnvironmentSimilarity::getReferenceEnvironment( const PDB& pdb, const std::vector<std::string>& anames,  double& maxdist ) {
  unsigned natoms = pdb.getPositions().size();
  std::vector<std::pair<unsigned,Vector> > env( natoms );
  for(unsigned i=0; i<natoms; ++i) {
    unsigned identity=0;
    for(unsigned j=1; j<anames.size(); ++j) {
      if( pdb.getAtomName(pdb.getAtomNumbers()[i])==anames[j] ) {
        identity=j;
        break;
      }
    }
    env[i] = std::pair<unsigned,Vector>( identity, pdb.getPositions()[i] );
    double dist = env[i].second.modulo();
    if( dist>maxdist ) {
      maxdist = dist;
    }
  }
  return env;
}

}
}
