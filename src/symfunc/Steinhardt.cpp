/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "CoordinationNumbers.h"
#include "core/ActionShortcut.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR Q1
/*
Calculate 1st order Steinhardt parameters

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q3
/*
Calculate 3rd order Steinhardt parameters.

The 3rd order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, \f$i\f$ is complex vector whose components are
calculated using the following formula:

\f[
q_{3m}(i) = \frac{\sum_j \sigma( r_{ij} ) Y_{3m}(\mathbf{r}_{ij}) }{\sum_j \sigma( r_{ij} ) }
\f]

where \f$Y_{3m}\f$ is one of the 3rd order spherical harmonics so \f$m\f$ is a number that runs from \f$-3\f$ to
\f$+3\f$.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The Steinhardt parameters can be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

\f[
Q_3(i) = \sqrt{ \sum_{m=-3}^3 q_{3m}(i)^{*} q_{3m}(i) }
\f]

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. Furthermore, when
the keywords LESS_THAN, MIN, MAX, HISTOGRAM, MEAN and so on are used with this colvar it is the distribution of these normed quantities
that is investigated.

Other measures of order can be taken by averaging the components of the individual \f$q_3\f$ vectors individually or by taking dot products of
the \f$q_{3}\f$ vectors on adjacent atoms.  More information on these variables can be found in the documentation for \ref LOCAL_Q3,
\ref LOCAL_AVERAGE and \ref NLINKS.

\par Examples

The following command calculates the average Q3 parameter for the 64 atoms in a box of Lennard Jones and prints this
quantity to a file called colvar:

\plumedfile
Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 MEAN LABEL=q3
PRINT ARG=q3.mean FILE=colvar
\endplumedfile

The following command calculates the histogram of Q3 parameters for the 64 atoms in a box of Lennard Jones and prints these
quantities to a file called colvar:

\plumedfile
Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=q3
PRINT ARG=q3.* FILE=colvar
\endplumedfile

The following command could be used to measure the Q3 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na\f$^+\f$ ions followed by the 64 Cl\f$-\f$ ions.  Once again the average Q3 parameter is calculated and output to a
file called colvar

\plumedfile
Q3 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN LABEL=q3
PRINT ARG=q3.mean FILE=colvar
\endplumedfile

If you simply want to examine the values of the Q3 parameters for each of the atoms in your system you can do so by exploiting the
command \ref DUMPATOMS as shown in the example below.  The following output file will output a file in an extended xyz format
called q3.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q3 parameter, columns
6-12 will contain the real parts of the director of the \f$q_{3m}\f$ vector while columns 12-19 will contain the imaginary parts of this director.

\plumedfile
q3: Q3 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q3 ARG=q3_anorm FILE=q3.xyz
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q4
/*
Calculate fourth order Steinhardt parameters.

The fourth order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, \f$i\f$ is complex vector whose components are
calculated using the following formula:

\f[
q_{4m}(i) = \frac{\sum_j \sigma( r_{ij} ) Y_{4m}(\mathbf{r}_{ij}) }{\sum_j \sigma( r_{ij} ) }
\f]

where \f$Y_{4m}\f$ is one of the fourth order spherical harmonics so \f$m\f$ is a number that runs from \f$-4\f$ to
\f$+4\f$.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The Steinhardt parameters can be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

\f[
Q_4(i) = \sqrt{ \sum_{m=-4}^4 q_{4m}(i)^{*} q_{4m}(i) }
\f]

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. Furthermore, when
the keywords LESS_THAN, MIN, MAX, HISTOGRAM, MEAN and so on are used with this colvar it is the distribution of these normed quantities
that is investigated.

Other measures of order can be taken by averaging the components of the individual \f$q_4\f$ vectors individually or by taking dot products of
the \f$q_{4}\f$ vectors on adjacent atoms.  More information on these variables can be found in the documentation for \ref LOCAL_Q4,
\ref LOCAL_AVERAGE and \ref NLINKS.

\par Examples

The following command calculates the average Q4 parameter for the 64 atoms in a box of Lennard Jones and prints this
quantity to a file called colvar:

\plumedfile
Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 MEAN LABEL=q4
PRINT ARG=q4.mean FILE=colvar
\endplumedfile

The following command calculates the histogram of Q4 parameters for the 64 atoms in a box of Lennard Jones and prints these
quantities to a file called colvar:

\plumedfile
Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=q4
PRINT ARG=q4.* FILE=colvar
\endplumedfile

The following command could be used to measure the Q4 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na\f$^+\f$ ions followed by the 64 Cl\f$-\f$ ions.  Once again the average Q4 parameter is calculated and output to a
file called colvar

\plumedfile
Q4 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN LABEL=q4
PRINT ARG=q4.mean FILE=colvar
\endplumedfile

If you simply want to examine the values of the Q4 parameters for each of the atoms in your system you can do so by exploiting the
command \ref DUMPATOMS as shown in the example below.  The following output file will output a file in an extended xyz format
called q$.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q4 parameter, columns
6-15 will contain the real parts of the director of the \f$q_{6m}\f$ vector while columns 15-24 will contain the imaginary parts of this director.

\plumedfile
q4: Q4 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ATOMS=q4 ARG=q4_anorm FILE=q4.xyz
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR Q6
/*
Calculate sixth order Steinhardt parameters.

The sixth order Steinhardt parameters allow us to measure the degree to which the first coordination shell
around an atom is ordered.  The Steinhardt parameter for atom, \f$i\f$ is complex vector whose components are
calculated using the following formula:

\f[
q_{6m}(i) = \frac{\sum_j \sigma( r_{ij} ) Y_{6m}(\mathbf{r}_{ij}) }{\sum_j \sigma( r_{ij} ) }
\f]

where \f$Y_{6m}\f$ is one of the sixth order spherical harmonics so \f$m\f$ is a number that runs from \f$-6\f$ to
\f$+6\f$.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

The Steinhardt parameters can be used to measure the degree of order in the system in a variety of different ways.  The
simplest way of measuring whether or not the coordination sphere is ordered is to simply take the norm of the above vector i.e.

\f[
Q_6(i) = \sqrt{ \sum_{m=-6}^6 q_{6m}(i)^{*} q_{6m}(i) }
\f]

This norm is small when the coordination shell is disordered and larger when the coordination shell is ordered. Furthermore, when
the keywords LESS_THAN, MIN, MAX, HISTOGRAM, MEAN and so on are used with this colvar it is the distribution of these normed quantities
that is investigated.

Other measures of order can be taken by averaging the components of the individual \f$q_6\f$ vectors individually or by taking dot products of
the \f$q_{6}\f$ vectors on adjacent atoms.  More information on these variables can be found in the documentation for \ref LOCAL_Q6,
\ref LOCAL_AVERAGE and \ref NLINKS.

\par Examples

The following command calculates the average Q6 parameter for the 64 atoms in a box of Lennard Jones and prints this
quantity to a file called colvar:

\plumedfile
Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 MEAN LABEL=q6
PRINT ARG=q6.mean FILE=colvar
\endplumedfile

The following command calculates the histogram of Q6 parameters for the 64 atoms in a box of Lennard Jones and prints these
quantities to a file called colvar:

\plumedfile
Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=q6
PRINT ARG=q6.* FILE=colvar
\endplumedfile

The following command could be used to measure the Q6 parameters that describe the arrangement of chlorine ions around the
sodium atoms in sodium chloride.  The imagined system here is composed of 64 NaCl formula units and the atoms are arranged in the input
with the 64 Na\f$^+\f$ ions followed by the 64 Cl\f$-\f$ ions.  Once again the average Q6 parameter is calculated and output to a
file called colvar

\plumedfile
Q6 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN LABEL=q6
PRINT ARG=q6.mean FILE=colvar
\endplumedfile

If you simply want to examine the values of the Q6 parameters for each of the atoms in your system you can do so by exploiting the
command \ref DUMPATOMS as shown in the example below.  The following output file will output a file in an extended xyz format
called q6.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q6 parameter, columns
6-19 will contain the real parts of the director of the \f$q_{6m}\f$ vector while columns 20-33 will contain the imaginary parts of this director.

\plumedfile
q6: Q6 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPATOMS ARG=q6_anorm ATOMS=q6 FILE=q6.xyz
\endplumedfile

*/
//+ENDPLUMEDOC

class Steinhardt : public ActionShortcut {
private:
  std::string getSymbol( const int& m ) const ;
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const int& l, const std::string& sep, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  explicit Steinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Steinhardt,"Q1")
PLUMED_REGISTER_ACTION(Steinhardt,"Q3")
PLUMED_REGISTER_ACTION(Steinhardt,"Q4")
PLUMED_REGISTER_ACTION(Steinhardt,"Q6")

void Steinhardt::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.addFlag("LOWMEM",false,"this flag does nothing and is present only to ensure back-compatibility");
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","the norm of the mean vector");
  keys.needsAction("GROUP");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("SPHERICAL_HARMONIC");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("MEAN");
  keys.needsAction("SUM");
}

Steinhardt::Steinhardt( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  std::string sp_str, specA, specB;
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this );
  int l;
  std::string sph_input = getShortcutLabel() + "_sh: SPHERICAL_HARMONIC ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z," + getShortcutLabel() + "_mat.w";

  if( getName()=="Q1" ) {
    sph_input +=" L=1";
    l=1;
  } else if( getName()=="Q3" ) {
    sph_input += " L=3";
    l=3;
  } else if( getName()=="Q4" ) {
    sph_input += " L=4";
    l=4;
  } else if( getName()=="Q6" ) {
    sph_input += " L=6";
    l=6;
  } else {
    plumed_merror("invalid input");
  }
  readInputLine( sph_input );

  // Input for denominator (coord)
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_denom_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat.w," + getShortcutLabel() + "_denom_ones" );
  readInputLine( getShortcutLabel() + "_sp: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_sh.*," + getShortcutLabel() + "_denom_ones");

  // If we are doing VMEAN determine sum of vector components
  std::string snum;
  bool do_vmean;
  parseFlag("VMEAN",do_vmean);
  bool do_vsum;
  parseFlag("VSUM",do_vsum);
  if( do_vmean || do_vsum ) {
    // Divide all components by coordination numbers
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine( getShortcutLabel() + "_rmn-" + snum + ": CUSTOM ARG=" + getShortcutLabel() + "_sp.rm-" + snum + "," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      // Imaginary part
      readInputLine( getShortcutLabel() + "_imn-" + snum + ": CUSTOM ARG=" + getShortcutLabel() + "_sp.im-" + snum + "," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    }
  }

  if( do_vmean ) {
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine( getShortcutLabel() + "_rms-" + snum + ": MEAN ARG=" + getShortcutLabel() + "_rmn-" + snum + " PERIODIC=NO");
      // Imaginary part
      readInputLine( getShortcutLabel() + "_ims-" + snum + ": MEAN ARG=" + getShortcutLabel() + "_imn-" + snum + " PERIODIC=NO");
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", l, "_", "ms" );
  }
  if( do_vsum ) {
    for(int i=-l; i<=l; ++i) {
      snum = getSymbol( i );
      // Real part
      readInputLine( getShortcutLabel() + "_rmz-" + snum + ": SUM ARG=" + getShortcutLabel() + "_rmn-" + snum + " PERIODIC=NO");
      // Imaginary part
      readInputLine( getShortcutLabel() + "_imz-" + snum + ": SUM ARG=" + getShortcutLabel() + "_imn-" + snum + " PERIODIC=NO");
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", l, "_", "mz" );
  }

  // Now calculate the total length of the vector
  createVectorNormInput( getShortcutLabel() + "_sp", getShortcutLabel() + "_norm", l, ".", "m" );
  // And take average
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_norm," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

void Steinhardt::createVectorNormInput( const std::string& ilab, const std::string& olab, const int& l, const std::string& sep, const std::string& vlab ) {
  std::string arg_inp, norm_input = olab + "2: COMBINE PERIODIC=NO POWERS=2,2";
  std::string snum = getSymbol( -l );
  arg_inp = " ARG=" + ilab + sep + "r" + vlab + "-" + snum +"," + ilab + sep + "i" + vlab + "-" + snum;
  for(int i=-l+1; i<=l; ++i) {
    snum = getSymbol( i );
    arg_inp += "," + ilab + sep + "r" + vlab + "-" + snum + "," + ilab + sep + "i" + vlab + "-" + snum;
    norm_input += ",2,2";
  }
  readInputLine( norm_input + arg_inp );
  readInputLine( olab + ": CUSTOM ARG=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO");
}

std::string Steinhardt::getSymbol( const int& m ) const {
  if( m<0 ) {
    std::string num;
    Tools::convert( -1*m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num;
    Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

}
}

