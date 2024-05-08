/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVARF LOCAL_Q1
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_1 vector on the central atom and the q_3 vector on the atoms in the first coordination sphere.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q3
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_3 vector on the central atom and the q_3 vector on the atoms in the first coordination sphere.

The \ref Q3 command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average \ref Q3 parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

\ref LOCAL_AVERAGE and \ref NLINKS are variables that can be combined with the Steinhardt parameters allow to calculate variables that satisfy these requirements.
LOCAL_Q3 is another variable that can be used in these sorts of calculations. The LOCAL_Q3 parameter for a particular atom is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

\f[
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-3}^3 q_{3m}^{*}(i)q_{3m}(j) }{ \sum_j \sigma( r_{ij} ) }
\f]

where \f$q_{3m}(i)\f$ and \f$q_{3m}(j)\f$ are the 3rd order Steinhardt vectors calculated for atom \f$i\f$ and atom \f$j\f$ respectively and the asterisk denotes complex
conjugation.  The function
\f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set
so that it the function is equal to one when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms \f$i\f$ and \f$j\f$ and thus measures the degree to which the orientations of these
adjacent atoms is correlated.

\par Examples

The following command calculates the average value of the LOCAL_Q3 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

\plumedfile
Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q3
LOCAL_Q3 SPECIES=q3 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=lq3
PRINT ARG=lq3.mean FILE=colvar
\endplumedfile

The following input calculates the distribution of LOCAL_Q3 parameters at any given time and outputs this information to a file.

\plumedfile
Q3 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q3
LOCAL_Q3 SPECIES=q3 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=lq3
PRINT ARG=lq3.* FILE=colvar
\endplumedfile

The following calculates the LOCAL_Q3 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

\plumedfile
Q3 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q3a
Q3 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q3b

LOCAL_Q3 SPECIES=q3a,q3b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LOWMEM LABEL=w3
PRINT ARG=w3.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q4
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_4 vector on the central atom and the q_4 vector on the atoms in the first coordination sphere.

The \ref Q4 command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average \ref Q4 parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

\ref LOCAL_AVERAGE and \ref NLINKS are variables that can be combined with the Steinhardt parameters allow to calculate variables that satisfy these requirements.
LOCAL_Q4 is another variable that can be used in these sorts of calculations. The LOCAL_Q4 parameter for a particular atom is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

\f[
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-4}^4 q_{4m}^{*}(i)q_{4m}(j) }{ \sum_j \sigma( r_{ij} ) }
\f]

where \f$q_{4m}(i)\f$ and \f$q_{4m}(j)\f$ are the fourth order Steinhardt vectors calculated for atom \f$i\f$ and atom \f$j\f$ respectively and the asterisk denotes
complex conjugation.  The function
\f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set
so that it the function is equal to one when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms \f$i\f$ and \f$j\f$ and thus measures the degree to which the orientations of these
adjacent atoms is correlated.

\par Examples

The following command calculates the average value of the LOCAL_Q4 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

\plumedfile
Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q4
LOCAL_Q4 SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=lq4
PRINT ARG=lq4.mean FILE=colvar
\endplumedfile

The following input calculates the distribution of LOCAL_Q4 parameters at any given time and outputs this information to a file.

\plumedfile
Q4 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q4
LOCAL_Q4 SPECIES=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=lq4
PRINT ARG=lq4.* FILE=colvar
\endplumedfile

The following calculates the LOCAL_Q4 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

\plumedfile
Q4 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q4a
Q4 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q4b

LOCAL_Q4 SPECIES=q4a,q4b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LOWMEM LABEL=w4
PRINT ARG=w4.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q6
/*
Calculate the local degree of order around an atoms by taking the average dot product between the q_6 vector on the central atom and the q_6 vector on the atoms in the first coordination sphere.

The \ref Q6 command allows one to calculate one complex vectors for each of the atoms in your system that describe the degree of order in the coordination sphere
around a particular atom. The difficulty with these vectors comes when combining the order parameters from all of the individual atoms/molecules so as to get a
measure of the global degree of order for the system. The simplest way of doing this - calculating the average Steinhardt parameter - can be problematic. If one is
examining nucleation say only the order parameters for those atoms in the nucleus will change significantly when the nucleus forms. The order parameters for the
atoms in the surrounding liquid will remain pretty much the same. As such if one models a small nucleus embedded in a very large amount of solution/melt any
change in the average order parameter will be negligible. Substantial changes in the value of this average can be observed in simulations of nucleation but only
because the number of atoms is relatively small.

When the average \ref Q6 parameter is used to bias the dynamics a problems
can occur. These averaged coordinates cannot distinguish between the correct,
single-nucleus pathway and a concerted pathway in which all the atoms rearrange
themselves into their solid-like configuration simultaneously. This second type
of pathway would be impossible in reality because there is a large entropic
barrier that prevents concerted processes like this from happening.  However,
in the finite sized systems that are commonly simulated this barrier is reduced
substantially. As a result in simulations where average Steinhardt parameters
are biased there are often quite dramatic system size effects

If one wants to simulate nucleation using some form on
biased dynamics what is really required is an order parameter that measures:

- Whether or not the coordination spheres around atoms are ordered
- Whether or not the atoms that are ordered are clustered together in a crystalline nucleus

\ref LOCAL_AVERAGE and \ref NLINKS are variables that can be combined with the Steinhardt parameters allow to calculate variables that satisfy these requirements.
LOCAL_Q6 is another variable that can be used in these sorts of calculations. The LOCAL_Q6 parameter for a particular atom is a number that measures the extent to
which the orientation of the atoms in the first coordination sphere of an atom match the orientation of the central atom.  It does this by calculating the following
quantity for each of the atoms in the system:

\f[
 s_i = \frac{ \sum_j \sigma( r_{ij} ) \sum_{m=-6}^6 q_{6m}^{*}(i)q_{6m}(j) }{ \sum_j \sigma( r_{ij} ) }
\f]

where \f$q_{6m}(i)\f$ and \f$q_{6m}(j)\f$ are the sixth order Steinhardt vectors calculated for atom \f$i\f$ and atom \f$j\f$ respectively and the asterisk denotes
complex conjugation.  The function
\f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set
so that it the function is equal to one when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.  The sum in the numerator
of this expression is the dot product of the Steinhardt parameters for atoms \f$i\f$ and \f$j\f$ and thus measures the degree to which the orientations of these
adjacent atoms is correlated.

\par Examples

The following command calculates the average value of the LOCAL_Q6 parameter for the 64 Lennard Jones atoms in the system under study and prints this
quantity to a file called colvar.

\plumedfile
Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q6
LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=lq6
PRINT ARG=lq6.mean FILE=colvar
\endplumedfile

The following input calculates the distribution of LOCAL_Q6 parameters at any given time and outputs this information to a file.

\plumedfile
Q6 SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=q6
LOCAL_Q6 SPECIES=q6 SWITCH={RATIONAL D_0=1.3 R_0=0.2} HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=1.0 NBINS=20 SMEAR=0.1} LABEL=lq6
PRINT ARG=lq6.* FILE=colvar
\endplumedfile

The following calculates the LOCAL_Q6 parameters for atoms 1-5 only. For each of these atoms comparisons of the geometry of the coordination sphere
are done with those of all the other atoms in the system.  The final quantity is the average and is outputted to a file

\plumedfile
Q6 SPECIESA=1-5 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q6a
Q6 SPECIESA=6-64 SPECIESB=1-64 D_0=1.3 R_0=0.2 LABEL=q6b

LOCAL_Q6 SPECIES=q6a,q6b SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LOWMEM LABEL=w6
PRINT ARG=w6.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

class LocalSteinhardt : public ActionShortcut {
private:
  std::string getSymbol( const int& m ) const ;
  std::string getArgsForStack( const int& l, const std::string& lab ) const;
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalSteinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q1")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q3")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q4")
PLUMED_REGISTER_ACTION(LocalSteinhardt,"LOCAL_Q6")

void LocalSteinhardt::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","SPECIES","");
  keys.add("optional","SPECIESA","");
  keys.add("optional","SPECIESB","");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("LOWMEM",false,"this flag does nothing and is present only to ensure back-compatibility");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("CONTACT_MATRIX"); keys.needsAction("MATRIX_PRODUCT"); keys.needsAction("GROUP");
  keys.needsAction("ONES"); keys.needsAction("OUTER_PRODUCT"); keys.needsAction("VSTACK");
  keys.needsAction("CONCATENATE"); keys.needsAction("CUSTOM"); keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

std::string LocalSteinhardt::getSymbol( const int& m ) const {
  if( m<0 ) {
    std::string num; Tools::convert( -1*m, num );
    return "n" + num;
  } else if( m>0 ) {
    std::string num; Tools::convert( m, num );
    return "p" + num;
  }
  return "0";
}

std::string LocalSteinhardt::getArgsForStack( const int& l, const std::string& sp_lab ) const {
  std::string numstr; Tools::convert( l, numstr );
  std::string data_mat = " ARG=" + sp_lab + "_sp.rm-n" + numstr + "," + sp_lab + "_sp.im-n" + numstr;
  for(int i=-l+1; i<=l; ++i) {
    numstr = getSymbol( i );
    data_mat += "," + sp_lab + "_sp.rm-" + numstr + "," + sp_lab + "_sp.im-" + numstr;
  }
  return data_mat;
}

LocalSteinhardt::LocalSteinhardt(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  bool lowmem; parseFlag("LOWMEM",lowmem);
  if( lowmem ) warning("LOWMEM flag is deprecated and is no longer required for this action");
  // Get the Q value
  int l; Tools::convert( getName().substr(7), l);
  // Create a vector filled with ones
  std::string twolplusone; Tools::convert( 2*(2*l+1), twolplusone );
  readInputLine( getShortcutLabel() + "_uvec: ONES SIZE=" + twolplusone );
  // Read in species keyword
  std::string sp_str; parse("SPECIES",sp_str);
  std::string spa_str; parse("SPECIESA",spa_str);
  if( sp_str.length()>0 ) {
    // Create a group with these atoms
    readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + sp_str );
    std::vector<std::string> sp_lab = Tools::getWords(sp_str, "\t\n ,");
    // This creates the stash to hold all the vectors
    if( sp_lab.size()==1 ) {
      // The lengths of all the vectors in a vector
      readInputLine( getShortcutLabel() + "_nmat: OUTER_PRODUCT ARG=" + sp_lab[0] + "_norm," + getShortcutLabel() + "_uvec");
      // The unormalised vectors
      readInputLine( getShortcutLabel() + "_uvecs: VSTACK" + getArgsForStack( l, sp_lab[0] ) );
    } else {
      std::string len_vec = getShortcutLabel() + "_mags: CONCATENATE ARG=" + sp_lab[0] + "_norm";
      for(unsigned i=1; i<sp_lab.size(); ++i) len_vec += "," + sp_lab[i] + "_norm";
      // This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = getShortcutLabel() + "_uvecs: CONCATENATE";
      for(unsigned i=0; i<sp_lab.size(); ++i) {
        std::string snum; Tools::convert( i+1, snum );
        concat_str += " MATRIX" + snum + "1=" + getShortcutLabel() + "_uvecs" + snum;
        readInputLine( getShortcutLabel() + "_uvecs" + snum + ": VSTACK" + getArgsForStack( l, sp_lab[i] ) );
      }
      // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
      readInputLine( getShortcutLabel() + "_nmat: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_mags," + getShortcutLabel() + "_uvec");
      // The unormalised vectors
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecs: CUSTOM ARG=" + getShortcutLabel() + "_uvecs," + getShortcutLabel() + "_nmat FUNC=x/y PERIODIC=NO");
    // And transpose the matrix
    readInputLine( getShortcutLabel() + "_vecsT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecs" );
    std::string sw_str; parse("SWITCH",sw_str); readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_str + " SWITCH={" + sw_str + "}");
    // And the matrix of dot products
    readInputLine( getShortcutLabel() + "_dpmat: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_vecs," + getShortcutLabel() + "_vecsT" );
  } else if( spa_str.length()>0 ) {
    // Create a group with these atoms
    readInputLine( getShortcutLabel() + "_grp: GROUP ATOMS=" + spa_str );
    std::string spb_str; parse("SPECIESB",spb_str);
    if( spb_str.length()==0 ) plumed_merror("need both SPECIESA and SPECIESB in input");
    std::vector<std::string> sp_laba = Tools::getWords(spa_str, "\t\n ,");
    std::vector<std::string> sp_labb = Tools::getWords(spb_str, "\t\n ,");
    if( sp_laba.size()==1 ) {
      // The matrix that is used for normalising
      readInputLine( getShortcutLabel() + "_nmatA: OUTER_PRODUCT ARG=" +  sp_laba[0] + "_norm," + getShortcutLabel() + "_uvec");
      // The unormalised vectors
      readInputLine( getShortcutLabel() + "_uvecsA: VSTACK" + getArgsForStack( l, sp_laba[0] ) );
    } else {
      std::string len_vec = getShortcutLabel() + "_magsA: CONCATENATE ARG=" + sp_laba[0] + "_norm";
      for(unsigned i=1; i<sp_laba.size(); ++i) len_vec += "," + sp_laba[i] + "_norm";
      //  This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = getShortcutLabel() + "_uvecsA: CONCATENATE";
      for(unsigned i=0; i<sp_laba.size(); ++i) {
        std::string snum; Tools::convert( i+1, snum );
        concat_str += " MATRIX" + snum + "1=" + getShortcutLabel() + "_uvecsA" + snum;
        readInputLine( getShortcutLabel() + "_uvecsA" + snum + ": VSTACK" + getArgsForStack( l, sp_laba[i] ) );
      }
      // And the normalising matrix by taking the column vector of magnitudes and multiplying by the row vector of ones
      readInputLine( getShortcutLabel() + "_nmatA: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_magsA," + getShortcutLabel() + "_uvec");
      // The unormalised vector
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecsA: CUSTOM ARG=" + getShortcutLabel() + "_uvecsA," + getShortcutLabel() + "_nmatA FUNC=x/y PERIODIC=NO");
    // Now do second matrix
    if( sp_labb.size()==1 ) {
      readInputLine( getShortcutLabel() + "_nmatB: OUTER_PRODUCT ARG=" +  sp_laba[0] + "_norm," + getShortcutLabel() + "_uvec");
      readInputLine( getShortcutLabel() + "_uvecsBT: VSTACK" + getArgsForStack( l, sp_labb[0] ) );
      readInputLine( getShortcutLabel() + "_uvecsB: TRANSPOSE ARG=" + getShortcutLabel() + "_uvecsBT");
    } else {
      std::string len_vec = getShortcutLabel() + "_magsB: CONCATENATE ARG=" +  sp_labb[0] + "_norm";
      for(unsigned i=1; i<sp_labb.size(); ++i) len_vec += "," + sp_labb[i] + "_norm";
      //  This is the vector that contains all the magnitudes
      readInputLine( len_vec );
      std::string concat_str = getShortcutLabel() + "_uvecsB: CONCATENATE";
      for(unsigned i=0; i<sp_labb.size(); ++i) {
        std::string snum; Tools::convert( i+1, snum );
        concat_str += " MATRIX1" + snum + "=" + getShortcutLabel() + "_uvecsB" + snum;
        readInputLine( getShortcutLabel() + "_uvecsBT" + snum + ": VSTACK" + getArgsForStack( l, sp_labb[i] ) );
        readInputLine( getShortcutLabel() + "_uvecsB" + snum + ": TRANSPOSE ARG=" + getShortcutLabel() + "_uvecsBT" + snum );
      }
      // And the normalising matrix
      readInputLine( getShortcutLabel() + "_nmatB: OUTER_PRODUCT ARG=" + getShortcutLabel() + "_uvec," + getShortcutLabel() + "_magsB");
      // The unormalised vectors
      readInputLine( concat_str );
    }
    // Now normalise all the vectors by doing Hadammard "product" with normalising matrix
    readInputLine( getShortcutLabel() + "_vecsB: CUSTOM ARG=" + getShortcutLabel() + "_uvecsB," + getShortcutLabel() + "_nmatB FUNC=x/y PERIODIC=NO");
    std::string sw_str; parse("SWITCH",sw_str); readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + spa_str + " GROUPB=" + spb_str + " SWITCH={" + sw_str + "}");
    readInputLine( getShortcutLabel() + "_dpmat: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_vecsA," + getShortcutLabel() + "_vecsB" );
  }

  // Now create the product matrix
  readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() + "_dpmat FUNC=x*y PERIODIC=NO");
  // Now the sum of coordination numbers times the switching functions
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmap");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size; Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() +"_prod," + getShortcutLabel() +"_ones");
  // And just the sum of the coordination numbers
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() +"_ones");
  // And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_av: CUSTOM ARG=" + getShortcutLabel() + "," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_av", "", this );
}

}
}
