/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "tools/KernelFunctions.h"
#include "tools/IFile.h"
#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "PammObject.h"

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

class PAMM : public multicolvar::MultiColvarBase {
private:
  PammObject mypamm;
public:
  static void registerKeywords( Keywords& keys );
  explicit PAMM(const ActionOptions&);
/// We have to overwrite this here
  unsigned getNumberOfQuantities() const override;
/// Calculate the weight of this object ( average of input weights )
  using PLMD::multicolvar::MultiColvarBase::calculateWeight;
  void calculateWeight( multicolvar::AtomValuePack& myatoms );
/// Actually do the calculation
  double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const override;
/// This returns the position of the central atom
  Vector getCentralAtom();
/// Is the variable periodic
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(PAMM,"PAMM")

void PAMM::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","the multicolvars from which the pamm coordinates are calculated");
  keys.add("compulsory","CLUSTERS","the name of the file that contains the definitions of all the clusters");
  keys.add("compulsory","REGULARISE","0.001","don't allow the denominator to be smaller then this value");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("SUM"); keys.use("LESS_THAN"); keys.use("HISTOGRAM");
  keys.use("MIN"); keys.use("MAX"); keys.use("LOWEST"); keys.use("HIGHEST"); keys.use("ALT_MIN"); keys.use("BETWEEN"); keys.use("MOMENTS");
  keys.setComponentsIntroduction("When the label of this action is used as the input for a second you are not referring to a scalar quantity as you are in "
                                 "regular collective variables.  The label is used to reference the full set of quantities calculated by "
                                 "the action.  This is usual when using \\ref multicolvarfunction. Generally when doing this the set of PAMM variables "
                                 "will be referenced using the DATA keyword rather than ARG.\n\n"
                                 "This Action can be used to calculate the following scalar quantities directly from the underlying set of PAMM variables. "
                                 "These quantities are calculated by employing the keywords listed below and they can be referenced elsewhere in the input "
                                 "file by using this Action's label followed by a dot and the name of the quantity.  The particular PAMM variable that should "
                                 "be averaged in a MEAN command or transformed by a switching function in a LESS_THAN command is specified using the COMPONENT "
                                 "keyword. COMPONENT=1 refers to the PAMM variable in which the first kernel in your input file is on the numerator, COMPONENT=2 refers to "
                                 "PAMM variable in which the second kernel in the input file is on the numerator and so on.  The same quantity can be calculated "
                                 "multiple times for different PAMM components by a single PAMM action.  In this case the relevant keyword must appear multiple "
                                 "times on the input line followed by a numerical identifier i.e. MEAN1, MEAN2, ...  The quantities calculated when multiple "
                                 "MEAN commands appear on the input line can be reference elsewhere in the input file by using the name of the quantity followed "
                                 "followed by a numerical identifier e.g. <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 etc.  Alternatively, you can "
                                 "customize the labels of the quantities by using the LABEL keyword in the description of the keyword.");
  keys.remove("ALL_INPUT_SAME_TYPE");
}

PAMM::PAMM(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // This builds the lists
  buildSets();
  // Check for reasonable input
  for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
    if( getBaseMultiColvar(i)->getNumberOfQuantities()!=2 ) error("cannot use PAMM with " + getBaseMultiColvar(i)->getName() );
  }

  bool mixed=getBaseMultiColvar(0)->isPeriodic();
  std::vector<bool> pbc( getNumberOfBaseMultiColvars() );
  std::vector<std::string> valnames( getNumberOfBaseMultiColvars() );
  std::vector<std::string> min( getNumberOfBaseMultiColvars() ), max( getNumberOfBaseMultiColvars() );
  for(unsigned i=0; i<getNumberOfBaseMultiColvars(); ++i) {
    if( getBaseMultiColvar(i)->isPeriodic()!=mixed ) warning("mixing of periodic and aperiodic base variables in pamm is untested");
    pbc[i]=getBaseMultiColvar(i)->isPeriodic();
    if( pbc[i] ) getBaseMultiColvar(i)->retrieveDomain( min[i], max[i] );
    valnames[i]=getBaseMultiColvar(i)->getLabel();
  }

  double regulariser; parse("REGULARISE",regulariser);
  std::string errorstr, filename; parse("CLUSTERS",filename);
  mypamm.setup( filename, regulariser, valnames, pbc, min, max, errorstr );
  if( errorstr.length()>0 ) error( errorstr );
}

unsigned PAMM::getNumberOfQuantities() const {
  return 1 + mypamm.getNumberOfKernels();
}

void PAMM::calculateWeight( multicolvar::AtomValuePack& myatoms ) {
  unsigned nvars = getNumberOfBaseMultiColvars();
  // Weight of point is average of weights of input colvars?
  std::vector<double> tval(2); double ww=0;
  for(unsigned i=0; i<nvars; ++i) {
    getInputData( i, false, myatoms, tval ); ww+=tval[0];
  }
  myatoms.setValue( 0, ww / static_cast<double>( nvars ) );

  if(!doNotCalculateDerivatives() ) {
    double pref = 1.0 / static_cast<double>( nvars );
    for(unsigned ivar=0; ivar<nvars; ++ivar) {
      // Get the values of derivatives
      MultiValue& myder=getInputDerivatives( ivar, false, myatoms );
      for(unsigned j=0; j<myder.getNumberActive(); ++j) {
        unsigned jder=myder.getActiveIndex(j);
        myatoms.addDerivative( 0, jder, pref*myder.getDerivative( 0, jder ) );
      }
    }
  }
}

double PAMM::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  unsigned nvars = getNumberOfBaseMultiColvars();
  std::vector<std::vector<double> > tderiv( mypamm.getNumberOfKernels() );
  for(unsigned i=0; i<tderiv.size(); ++i) tderiv[i].resize( nvars );
  std::vector<double> tval(2), invals( nvars ), vals( mypamm.getNumberOfKernels() );

  for(unsigned i=0; i<nvars; ++i) {
    getInputData( i, false, myatoms, tval ); invals[i]=tval[1];
  }
  mypamm.evaluate( invals, vals, tderiv );

  // Now set all values other than the first one
  // This is because of some peverse choices in multicolvar
  for(unsigned i=1; i<vals.size(); ++i) myatoms.setValue( 1+i, vals[i] );

  if( !doNotCalculateDerivatives() ) {
    std::vector<double> mypref( 1 + vals.size() );
    for(unsigned ivar=0; ivar<nvars; ++ivar) {
      // Get the values of the derivatives
      MultiValue& myder = getInputDerivatives( ivar, false, myatoms );
      // And calculate the derivatives
      for(unsigned i=0; i<vals.size(); ++i) mypref[1+i] = tderiv[i][ivar];
      // This is basically doing the chain rule to get the final derivatives
      splitInputDerivatives( 1, 1, 1+vals.size(), ivar, mypref, myder, myatoms );
      // And clear the derivatives
      myder.clearAll();
    }
  }

  return vals[0];
}

Vector PAMM::getCentralAtom() {
  // Who knows how this should work
  plumed_error();
  return Vector(1.0,0.0,0.0);
}

}
}
