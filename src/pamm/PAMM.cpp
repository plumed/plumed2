/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "tools/IFile.h"
#include "core/ActionSetup.h"

//+PLUMEDOC MCOLVARF PAMM
/*
Probabilistic analysis of molecular motifs.

Probabilistic analysis of molecular motifs (PAMM) was introduced in this paper \cite pamm.
The essence of this approach involves calculating some large set of collective variables
for a set of atoms in a short trajectory and fitting this data using a Gaussian Mixture Model.
The idea is that modes in these distributions can be used to identify features such as hydrogen bonds or
secondary structure types.

The assumption within this implementation is that the fitting of the Gaussian mixture model has been
done elsewhere by a separate code.  You thus provide an input file to this action which contains the
means, covariance matrices and weights for a set of Gaussian kernels, \f$\{ \phi \}\f$.  The values and
derivatives for the following set of quantities is then computed:

\f[
s_k = \frac{ \phi_k}{ \sum_i \phi_i }
\f]

Each of the \f$\phi_k\f$ is a Gaussian function that acts on a set of quantities calculated within
a \ref mcolv .  These might be \ref TORSIONS, \ref DISTANCES, \ref ANGLES or any one of the many
symmetry functions that are available within \ref mcolv actions.  These quantities are then inserted into
the set of \f$n\f$ kernels that are in the the input file.   This will be done for multiple sets of values
for the input quantities and a final quantity will be calculated by summing the above \f$s_k\f$ values or
some transformation of the above.  This sounds less complicated than it is and is best understood by
looking through the example given below.

\warning Mixing \ref mcolv actions that are periodic with variables that are not periodic has not been tested

\par Examples

In this example I will explain in detail what the following input is computing:

\plumedfile
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=M1d.pdb
psi: TORSIONS ATOMS1=@psi-2 ATOMS2=@psi-3 ATOMS3=@psi-4
phi: TORSIONS ATOMS1=@phi-2 ATOMS2=@phi-3 ATOMS3=@phi-4
p: PAMM DATA=phi,psi CLUSTERS=clusters.pamm MEAN1={COMPONENT=1} MEAN2={COMPONENT=2}
PRINT ARG=p.mean-1,p.mean-2 FILE=colvar
\endplumedfile

The best place to start our explanation is to look at the contents of the clusters.pamm file

\auxfile{clusters.pamm}
#! FIELDS height phi psi sigma_phi_phi sigma_phi_psi sigma_psi_phi sigma_psi_psi
#! SET multivariate von-misses
#! SET kerneltype gaussian
      2.97197455E-0001     -1.91983118E+0000      2.25029540E+0000      2.45960237E-0001     -1.30615381E-0001     -1.30615381E-0001      2.40239117E-0001
      2.29131448E-0002      1.39809354E+0000      9.54585380E-0002      9.61755708E-0002     -3.55657919E-0002     -3.55657919E-0002      1.06147253E-0001
      5.06676398E-0001     -1.09648066E+0000     -7.17867907E-0001      1.40523052E-0001     -1.05385552E-0001     -1.05385552E-0001      1.63290557E-0001
\endauxfile

This files contains the parameters of two two-dimensional Gaussian functions.  Each of these Gaussian kernels has a weight, \f$w_k\f$,
a vector that specifies the position of its center, \f$\mathbf{c}_k\f$, and a covariance matrix, \f$\Sigma_k\f$.  The \f$\phi_k\f$ functions that
we use to calculate our PAMM components are thus:

\f[
\phi_k = \frac{w_k}{N_k} \exp\left( -(\mathbf{s} - \mathbf{c}_k)^T \Sigma^{-1}_k (\mathbf{s} - \mathbf{c}_k) \right)
\f]

In the above \f$N_k\f$ is a normalization factor that is calculated based on \f$\Sigma\f$.  The vector \f$\mathbf{s}\f$ is a vector of quantities
that are calculated by the \ref TORSIONS actions.  This vector must be two dimensional and in this case each component is the value of a
torsion angle.  If we look at the two \ref TORSIONS actions in the above we are calculating the \f$\phi\f$ and \f$\psi\f$ backbone torsional
angles in a protein (Note the use of \ref MOLINFO to make specification of atoms straightforward).  We thus calculate the values of our
2 \f$ \{ \phi \} \f$  kernels 3 times.  The first time we use the \f$\phi\f$ and \f$\psi\f$ angles in the second residue of the protein,
the second time it is the \f$\phi\f$ and \f$\psi\f$ angles of the third residue of the protein and the third time it is the \f$\phi\f$ and \f$\psi\f$ angles
of the fourth residue in the protein.  The final two quantities that are output by the print command, p.mean-1 and p.mean-2, are the averages
over these three residues for the quantities:
\f[
s_1 = \frac{ \phi_1}{ \phi_1 + \phi_2 }
\f]
and
\f[
s_2 = \frac{ \phi_2}{ \phi_1 + \phi_2 }
\f]
There is a great deal of flexibility in this input.  We can work with, and examine, any number of components, we can use any set of collective variables
and compute these PAMM variables and we can transform the PAMM variables themselves in a large number of different ways when computing these sums.
*/
//+ENDPLUMEDOC

namespace PLMD {
namespace pamm {

class PAMM : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit PAMM(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PAMM,"PAMM")

void PAMM::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the vectors from which the pamm coordinates are calculated");
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the clusters");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.add("compulsory","KERNELS","all","which kernels are we computing the PAMM values for");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("KERNEL");
  keys.needsAction("COMBINE");
}

PAMM::PAMM(const ActionOptions& ao) :
  Action(ao),
  ActionShortcut(ao) {
  // Must get list of input value names
  std::vector<std::string> valnames;
  parseVector("ARG",valnames);
  // Create input values
  std::string argstr=" ARG=" + valnames[0];
  for(unsigned j=1; j<valnames.size(); ++j) {
    argstr += "," + valnames[j];
  }

  // Create actions to calculate all pamm kernels
  unsigned nkernels = 0;
  double h;
  std::string fname;
  parse("CLUSTERS",fname);
  IFile ifile;
  ifile.open(fname);
  ifile.allowIgnoredFields();
  for(unsigned k=0;; ++k) {
    if( !ifile.scanField("height",h) ) {
      break;
    }
    // Create a kernel for this cluster
    std::string num, wstr, ktype;
    Tools::convert( k+1, num );
    Tools::convert(h,wstr);
    ifile.scanField("kerneltype",ktype);
    readInputLine( getShortcutLabel() + "_kernel-" + num + ": KERNEL NORMALIZED" + argstr  + " NUMBER=" + num + " REFERENCE=" + fname + " WEIGHT=" + wstr + " TYPE=" + ktype );
    nkernels++;
    ifile.scanField();
  }
  ifile.close();

  // And add on the regularization
  std::string regparam;
  parse("REGULARISE",regparam);
  // Now combine all the PAMM objects with the regparam
  std::string paramstr, cinput = getShortcutLabel() + "_ksum: COMBINE PERIODIC=NO";
  for(unsigned k=0; k<nkernels; ++k) {
    std::string num;
    Tools::convert( k+1, num );
    if( k==0 ) {
      cinput += " ARG=";
      paramstr=" PARAMETERS=-" + regparam;
    } else {
      cinput += ",";
      paramstr += ",0";
    }
    cinput += getShortcutLabel() + "_kernel-" + num;
  }
  readInputLine( cinput + paramstr );

  // And now compute all the pamm kernels
  std::string kchoice;
  parse("KERNELS",kchoice);
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  if( kchoice=="all" ) {
    for(unsigned k=0; k<nkernels; ++k) {
      std::string num;
      Tools::convert( k+1, num );
      readInputLine( getShortcutLabel() + "-" + num + ": CUSTOM ARG=" + getShortcutLabel() + "_kernel-" + num + "," + getShortcutLabel() + "_ksum FUNC=x/y PERIODIC=NO");
      multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel() + "-" + num, getShortcutLabel() + "-" + num, "", keymap, this );
    }
  } else {
    std::vector<std::string> awords=Tools::getWords(kchoice,"\t\n ,");
    Tools::interpretRanges( awords );
    for(unsigned k=0; k<awords.size(); ++k) {
      readInputLine( getShortcutLabel() + "-" + awords[k] + ": CUSTOM ARG=" + getShortcutLabel() + "_kernel-" + awords[k] + "," + getShortcutLabel() + "_ksum FUNC=x/y PERIODIC=NO");
      multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel() + "-" + awords[k], getShortcutLabel() + "-" + awords[k], "", keymap, this );
    }
  }
}

}
}
