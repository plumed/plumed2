/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "vesselbase/ActionWithAveraging.h"
#include "core/ActionRegister.h"
#include "AverageVessel.h"

//+PLUMEDOC GRIDCALC AVERAGE
/*
Calculate the ensemble average of a collective variable

The ensemble average for a non-periodic, collective variable, \f$s\f$ is given by the following expression:

\f[
\langle s \rangle = \frac{ \sum_{t'=0}^t w(t') s(t') }{ \sum_{t'=0}^t w(t') }
\f]

Here the sum runs over a the trajectory and \f$s(t')\f$ is used to denote the value of the collective variable
at time \f$t'\f$.  The final quantity evaluated is a weighted
average as the weights, \f$w(t')\f$, allow us to negate the effect any bias might have on the region of phase space
sampled by the system.  This is discussed in the section of the manual on \ref Analysis.

When the variable is periodic (e.g. \ref TORSION) and has a value, \f$s\f$, in \f$a \le s \le b\f$ the ensemble average is evaluated using:

\f[
\langle s \rangle = a + \frac{b - a}{2\pi} \arctan \left[ \frac{ \sum_{t'=0}^t w(t') \sin\left( \frac{2\pi [s(t')-a]}{b - a} \right) }{ \sum_{t'=0}^t w(t') \cos\left( \frac{2\pi [s(t')-a]}{b - a} \right) } \right]
\f]

\par Examples

The following example calculates the ensemble average for the distance between atoms 1 and 2
and output this to a file called COLVAR.  In this example it is assumed that no bias is acting
on the system and that the weights, \f$w(t')\f$ in the formulas above can thus all be set equal
to one.

\plumedfile
d1: DISTANCE ATOMS=1,2
d1a: AVERAGE ARG=d1
PRINT ARG=d1a FILE=colvar STRIDE=100
\endplumedfile

The following example calculates the ensemble average for the torsional angle involving atoms 1, 2, 3 and 4.
At variance with the previous example this quantity is periodic so the second formula in the above introduction
is used to calculate the average.  Furthermore, by using the CLEAR keyword we have specified that block averages
are to be calculated.  Consequently, after 100 steps all the information acquired thus far in the simulation is
forgotten and the process of averaging is begun again.  The quantities output in the colvar file are thus the
block averages taken over the first 100 frames of the trajectory, the block average over the second 100 frames
of trajectory and so on.

\plumedfile
t1: TORSION ATOMS=1,2,3,4
t1a: AVERAGE ARG=t1 CLEAR=100
PRINT ARG=t1a FILE=colvar STRIDE=100
\endplumedfile

This third example incorporates a bias.  Notice that the effect the bias has on the ensemble average is removed by taking
advantage of the \ref REWEIGHT_BIAS method.  The final ensemble averages output to the file are thus block ensemble averages for the
unbiased canonical ensemble at a temperature of 300 K.

\plumedfile
t1: TORSION ATOMS=1,2,3,4
RESTRAINT ARG=t1 AT=pi KAPPA=100.
ww: REWEIGHT_BIAS TEMP=300
t1a: AVERAGE ARG=t1 LOGWEIGHTS=ww CLEAR=100
PRINT ARG=t1a FILE=colvar STRIDE=100
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class Average : public vesselbase::ActionWithAveraging {
private:
  AverageVessel* myaverage;
public:
  static void registerKeywords( Keywords& keys );
  explicit Average( const ActionOptions& );
  void calculate() override {}
  void apply() override {}
  void performOperations( const bool& from_update ) override;
  void finishAveraging() override;
  bool isPeriodic() override { return false; }
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
  void accumulateAverage( MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(Average,"AVERAGE")

void Average::registerKeywords( Keywords& keys ) {
  vesselbase::ActionWithAveraging::registerKeywords( keys ); keys.use("ARG");
  keys.remove("SERIAL"); keys.remove("LOWMEM");
}

Average::Average( const ActionOptions& ao ):
  Action(ao),
  ActionWithAveraging(ao)
{
  addValue(); // Create a value so that we can output the average
  if( getNumberOfArguments()!=1 ) error("only one quantity can be averaged at a time");
  std::string instring;
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max; getPntrToArgument(0)->getDomain(min,max);
    instring = "PERIODIC=" + min + "," + max; setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  // Create a vessel to hold the average
  vesselbase::VesselOptions da("myaverage","",-1,instring,this);
  Keywords keys; AverageVessel::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  auto average=Tools::make_unique<AverageVessel>(dar);
  myaverage = average.get();
  setAveragingAction( std::move(average), false );
}

void Average::performOperations( const bool& from_update ) {
  myaverage->accumulate( cweight, getArgument(0) );
}

void Average::accumulateAverage( MultiValue& myvals ) const {
  plumed_dbg_assert( myvals.getNumberOfValues()==3 );
  myaverage->accumulate( myvals.get(0), myvals.get(1) );
}

void Average::finishAveraging() {
  setValue( myaverage->getAverage() );
}

}
}
