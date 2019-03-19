/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q3
/*
Calculate the local degree of order around an atoms by taking the average dot product between the \f$q_3\f$ vector on the central atom and the \f$q_3\f$ vector on the atoms in the first coordination sphere.

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


namespace PLMD {
namespace crystallization {

class Q3 : public Steinhardt {
public:
  static void registerKeywords( Keywords& keys );
  explicit Q3( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(Q3,"Q3")
typedef LocalSteinhardt<Q3> LOCAL_Q3;
PLUMED_REGISTER_ACTION(LOCAL_Q3,"LOCAL_Q3")

void Q3::registerKeywords( Keywords& keys ) {
  Steinhardt::registerKeywords( keys );
}

Q3::Q3(const ActionOptions& ao ):
  Action(ao),
  Steinhardt(ao)
{
  setAngularMomentum(3);

// Spherical harmonics normalization:
// even =  sqrt ( ((2l+1)*(l-m)!) / (4*pi*(l+m)!) )
// odd  = -sqrt ( ((2l+1)*(l-m)!) / (4*pi*(l+m)!) )

  normaliz.resize( 4 );
  normaliz[0] = sqrt( ( 7.0*6.0 ) / (4.0*pi*6.0) );
  normaliz[1] = -sqrt( ( 7.0*2.0 ) / (4.0*pi*24.0) );
  normaliz[2] = sqrt( ( 7.0*1.0) / (4.0*pi*120.0) );
  normaliz[3] = -sqrt( ( 7.0*1.0) / (4.0*pi*720.0) );

// Legendre polynomial coefficients of order three

  coeff_poly.resize( 4 );
  coeff_poly[0]=0.0;
  coeff_poly[1]=-1.5;
  coeff_poly[2]=0.0;
  coeff_poly[3]=2.5;

}

}
}

