/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "Steinhardt.h"
#include "LocalSteinhardt.h"
#include "core/ActionRegister.h"

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
command \ref DUMPMULTICOLVAR as shown in the example below.  The following output file will output a file in an extended xyz format
called q6.xyz for each frame of the analyzed MD trajectory.  The first column in this file will contain a dummy name for each of the
atoms, columns 2-4 will then contain the x, y and z positions of the atoms, column 5 will contain the value of the Q6 parameter, columns
6-19 will contain the real parts of the director of the \f$q_{6m}\f$ vector while columns 20-33 will contain the imaginary parts of this director.

\plumedfile
q6: Q6 SPECIESA=1-64 SPECIESB=65-128 D_0=1.3 R_0=0.2 MEAN
DUMPMULTICOLVAR DATA=q6 FILE=q6.xyz
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

namespace PLMD {
namespace crystallization {

class Q6 : public Steinhardt {
public:
  static void registerKeywords( Keywords& keys );
  explicit Q6( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(Q6,"Q6")
typedef LocalSteinhardt<Q6> LOCAL_Q6;
PLUMED_REGISTER_ACTION(LOCAL_Q6,"LOCAL_Q6")

void Q6::registerKeywords( Keywords& keys ) {
  Steinhardt::registerKeywords( keys );
}

Q6::Q6(const ActionOptions& ao ):
  Action(ao),
  Steinhardt(ao)
{
  setAngularMomentum(6);

  normaliz.resize( 7 );
  normaliz[0] = sqrt( ( 13.0*720.0 ) / (4.0*pi*720.0) );
  normaliz[1] = -sqrt( ( 13.0*120.0 ) / (4.0*pi*5040) );
  normaliz[2] = sqrt( ( 13.0*24) / (4.0*pi*40320) );
  normaliz[3] = -sqrt( ( 13.0*6) / (4.0*pi*362880) );
  normaliz[4] = sqrt( (13.0*2) / (4.0*pi*3628800) );
  normaliz[5] = -sqrt( (13.0*1) / (4.0*pi*39916800) );
  normaliz[6] = sqrt( (13.0*1) / (4.0*pi*479001600) );

  coeff_poly.resize( 7 );
  coeff_poly[0]=-0.3125; coeff_poly[1]=0.0;
  coeff_poly[2]=6.5625; coeff_poly[3]=0.0;
  coeff_poly[4]=-19.6875; coeff_poly[5]=0.0;
  coeff_poly[6]=14.4375;
}

}
}

