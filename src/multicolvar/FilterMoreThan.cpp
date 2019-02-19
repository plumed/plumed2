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
#include "tools/SwitchingFunction.h"
#include "MultiColvarFilter.h"

//+PLUMEDOC MTRANSFORMS MTRANSFORM_MORE
/*
This action can be used to transform the colvar values calculated by a \ref mcolv using one minus a \ref switchingfunction

In this action each colvar, \f$s_i\f$, calculated by \ref mcolv is transformed by a \ref switchingfunction function that
is equal to one if the colvar is greater than a certain target value and which is equal to zero otherwise.
It is important to understand the distinction between what is done here and what is done by \ref MFILTER_MORE.
In \ref MFILTER_MORE a weight, \f$w_i\f$ for the colvar is calculated using the \ref histogrambead.  If one calculates the
MEAN for \ref MFILTER_MORE one is thus calculating:

\f[
\mu = \frac{ \sum_i [1 - \sigma(s_i) ] s_i }{\sum_i [1 - \sigma(s_i)] }
\f]

where \f$\sigma\f$ is the \ref switchingfunction.  In this action by contrast the colvar is being transformed by the \ref switchingfunction.
If one thus calculates a MEAN for this action one computes:

\f[
\mu = \frac{ \sum_{i=1}^N 1 - \sigma(s_i) }{ N }
\f]

In other words, you are calculating the mean for the transformed colvar.

\par Examples

The following input gives an example of how a MTRANSFORM_MORE action can be used to duplicate
functionality that is elsewhere in PLUMED.

\plumedfile
DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 LABEL=d1
... DISTANCES
MTRANSFORM_MORE DATA=d1 SWITCH={GAUSSIAN D_0=1.5 R_0=0.00001}
\endplumedfile

In this case you can achieve the same result by using:

\plumedfile
DISTANCES ...
 GROUPA=1-10 GROUPB=11-20
 MORE_THAN={GAUSSIAN D_0=1.5 R_0=0.00001}
... DISTANCES
\endplumedfile
(see \ref DISTANCES)

The advantage of MTRANSFORM_MORE comes, however, if you want to use transformed colvars as input
for \ref MULTICOLVARDENS

*/
//+ENDPLUMEDOC

//+PLUMEDOC MFILTERS MFILTER_MORE
/*
This action can be used to filter the distribution of colvar values in a \ref mcolv
so that one can compute the mean and so on for only those multicolvars more than a tolerance.

This action can be used to create a dynamic group of atom based on the value of a multicolvar.
In this action a multicolvar is within the dynamic group if its value is greater than a target.
In actuality a weight, \f$w_i\f$  is ascribed to each colvar, \f$s_i\f$ calculated by a multicolvar
and this weight measures the degree to which a colvar is a member of the group.  This weight is
calculated using a \ref switchingfunction , \f$\sigma\f$ so it is given by:

\f[
w_i = 1 - \sigma(s_i)
\f]

If one calculates a function of the set of multicolvars
these weights are included in the calculation.  As such if one calculates the MEAN, \f$\mu\f$ of a filtered
multicolvar what is computed is the following:

\f[
\mu = \frac{ \sum_i w_i s_i }{ \sum_i w_i}
\f]

One is thus calculating the mean for those colvars that are greater than the target.

\par Examples

The example shown below calculates the mean for those distances that greater than 1.5 nm in length

\plumedfile
DISTANCES GROUPA=1 GROUPB=2-50 MEAN LABEL=d1
MFILTER_MORE DATA=d1 SWITCH={GAUSSIAN D_0=1.5 R_0=0.00001} MEAN LABEL=d4
\endplumedfile

More complicated things can be done by using the label of a filter as input to a new multicolvar as shown
in the example below.  Here the coordination numbers of all atoms are computed.  The atoms with a coordination
number greater than 2 are then identified using the filter.  This reduced list of atoms is then used as input
to a second coordination number calculation.  This second coordination number thus measures the number of
two-coordinated atoms that each of the two-coordinated atoms is bound to.

\plumedfile
1: COORDINATIONNUMBER SPECIES=1-150 SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
cf: MFILTER_MORE DATA=c1 SWITCH={RATIONAL D_0=2.0 R_0=0.1} LOWMEM
c2: COORDINATIONNUMBER SPECIES=cf SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0} MORE_THAN={RATIONAL D_0=2.0 R_0=0.1}
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class FilterMore : public MultiColvarFilter {
private:
  SwitchingFunction sf;
public:
  static void registerKeywords( Keywords& keys );
  explicit FilterMore(const ActionOptions& ao);
  double applyFilter( const double& val, double& df ) const ;
};

PLUMED_REGISTER_ACTION(FilterMore,"MFILTER_MORE")
PLUMED_REGISTER_ACTION(FilterMore,"MTRANSFORM_MORE")

void FilterMore::registerKeywords( Keywords& keys ) {
  MultiColvarFilter::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
}

FilterMore::FilterMore(const ActionOptions& ao):
  Action(ao),
  MultiColvarFilter(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    sf.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    sf.set(nn,mm,r_0,d_0);
  }
  log.printf("  filtering colvar values and focussing only on those more than %s\n",( sf.description() ).c_str() );

  checkRead();
}

double FilterMore::applyFilter( const double& val, double& df ) const {
  double f = 1.0 - sf.calculate( val, df ); df*=-val;
  return f;
}

}
}
