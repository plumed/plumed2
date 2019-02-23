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
#include "core/ActionRegister.h"
#include "tools/HistogramBead.h"
#include "MultiColvarFilter.h"

//+PLUMEDOC MTRANSFORMS MTRANSFORM_BETWEEN
/*
This action can be used to transform the colvar values calculated by a MultiColvar using a \ref histogrambead

In this action each colvar, \f$s_i\f$, calculated by MultiColvar is transformed by a \ref histogrambead function that
is equal to one if the colvar is within a certain range and which is equal to zero otherwise.  In other words, we
compute:

\f[
f_i = \int_a^b K\left( \frac{s-s_i}{w} \right)
\f]

where \f$a, b\f$ and \f$w\f$ are parameters.

It is important to understand the distinction between what is done here and what is done by \ref MFILTER_BETWEEN.
In \ref MFILTER_BETWEEN a weight, \f$w_i\f$ for the colvar is calculated using the \ref histogrambead.  If one calculates the
MEAN for \ref MFILTER_BETWEEN one is thus calculating:

\f[
\mu = \frac{ \sum_i f_i s_i }{ \sum_i f_i}
\f]

In this action by contrast the colvar is being transformed by the \ref histogrambead.  If one thus calculates a MEAN for
this action one computes:

\f[
\mu = \frac{ \sum_{i=1}^N f_i }{ N }
\f]

In other words, you are calculating the mean for the transformed colvar.

\par Examples

The following input gives an example of how a \ref MTRANSFORM_BETWEEN action can be used to duplicate
functionality that is elsewhere in PLUMED.

\plumedfile
DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 LABEL=d1
... DISTANCES
MTRANSFORM_BETWEEN DATA=d1 LOWER=1.0 UPPER=2.0 SMEAR=0.5
\endplumedfile

In this case you can achieve the same result by using:

\plumedfile
DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 BETWEEN={GAUSSIAN LOWER=1.0 UPPER=2.0}
... DISTANCES
\endplumedfile
(see \ref DISTANCES)

The advantage of \ref MTRANSFORM_BETWEEN comes, however, if you want to use transformed colvars as input
for \ref MULTICOLVARDENS

*/
//+ENDPLUMEDOC

//+PLUMEDOC MFILTERS MFILTER_BETWEEN
/*
This action can be used to filter the colvar values calculated by a \ref mcolv
so that one can compute the mean and so on for only those multicolvars within a certain range.

This action can be used to create a dynamic group of atom based on the value of a multicolvar.
In this action a multicolvar is within the dynamic group if its value lies in a particular range.
In actuality a weight, \f$w_i\f$  is ascribed to each colvar, \f$s_i\f$ calculated by a multicolvar
and this weight measures the degree to which a colvar is a member of the group.  This weight is
calculated using a \ref histogrambead so it is given by:

\f[
w_i = \int_a^b K\left( \frac{s - s_i}{w} \right)
\f]

where \f$a, b\f$ and \f$w\f$ are parameters.  If one calculates a function of the set of multicolvars
these weights are included in the calculation.  As such if one calculates the MEAN, \f$\mu\f$ of a filtered
multicolvar what is computed is the following:

\f[
\mu = \frac{ \sum_i w_i s_i }{ \sum_i w_i}
\f]

One is thus calculating the mean for those colvars that are within the range of interest.

\par Examples

The example shown below calculates the mean for those distances that are between 0 and 3 nm in length

\plumedfile
DISTANCES GROUPA=1 GROUPB=2-50 MEAN LABEL=d1
MFILTER_BETWEEN DATA=d1 LOWER=0 UPPER=3.0 SMEAR=0.0001 MEAN LABEL=d4
\endplumedfile

More complicated things can be done by using the label of a filter as input to a new multicolvar as shown
in the example below.  Here the coordination numbers of all atoms are computed.  The atoms with a coordination
number between 4 and 6 are then identified using the filter.  This reduced list of atoms is then used as input
to a second coordination number calculation.  This second coordination number thus measures the number of atoms
4-6 coordinated atoms each of the 4-6 coordination atoms is bound to.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-150 SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
cf: MFILTER_BETWEEN DATA=c1 LOWER=4 UPPER=6 SMEAR=0.5 LOWMEM
c2: COORDINATIONNUMBER SPECIES=cf SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0} MORE_THAN={RATIONAL D_0=2.0 R_0=0.1}
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class FilterBetween : public MultiColvarFilter {
private:
  HistogramBead hb;
public:
  static void registerKeywords( Keywords& keys );
  explicit FilterBetween(const ActionOptions& ao);
  double applyFilter( const double& val, double& df ) const ;
};

PLUMED_REGISTER_ACTION(FilterBetween,"MFILTER_BETWEEN")
PLUMED_REGISTER_ACTION(FilterBetween,"MTRANSFORM_BETWEEN")

void FilterBetween::registerKeywords( Keywords& keys ) {
  MultiColvarFilter::registerKeywords( keys );
  keys.add("compulsory","LOWER","the lower boundary for the range of interest");
  keys.add("compulsory","UPPER","the upper boundary for the range of interest");
  keys.add("compulsory","SMEAR","0.5","the amount by which to smear the value for kernel density estimation");
  keys.add("optional","BEAD","This keywords is used if you want to employ an alternative to the function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER and SMEAR keywords.");
}

FilterBetween::FilterBetween(const ActionOptions& ao):
  Action(ao),
  MultiColvarFilter(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("BEAD",sw);
  if( getPntrToMultiColvar()->isPeriodic() ) {
    std::string min, max; getPntrToMultiColvar()->retrieveDomain( min, max );
    double mlow, mhigh; Tools::convert( min,mlow ); Tools::convert( max, mhigh);
    hb.isPeriodic( mlow, mhigh );
  } else {
    hb.isNotPeriodic();
  }

  if(sw.length()>0) {
    hb.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading BEAD keyword : " + errors );
  } else {
    double l, u, s; std::string ll, uu, ss;
    parse("LOWER",l); parse("UPPER",u); parse("SMEAR",s);
    Tools::convert(l,ll); Tools::convert(u,uu); Tools::convert(s,ss);
    sw="GAUSSIAN LOWER=" + ll + " UPPER=" + uu + " SMEAR=" + ss;
    hb.set(sw,errors); plumed_massert(errors.length()==0,"problems with bead" + errors);
  }
  log.printf("  filtering colvar values and focussing only on those values in range %s\n",( hb.description() ).c_str() );

  checkRead();
}

double FilterBetween::applyFilter( const double& val, double& df ) const {
  double f = hb.calculate( val, df );
  return f;
}

}
}
