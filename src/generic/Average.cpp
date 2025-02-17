/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"

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
namespace generic {

class Average : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Average( const ActionOptions& );
};

PLUMED_REGISTER_ACTION(Average,"AVERAGE")

void Average::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","scalar/grid","the quantity that is being averaged");
  keys.add("optional","LOGWEIGHTS","the logarithm of the quantity to use as the weights when calculating averages");
  keys.add("compulsory","STRIDE","1","the frequency with which to store data for averaging");
  keys.add("compulsory","CLEAR","0","the frequency with whihc to clear the data that is being averaged");
  keys.add("optional","NORMALIZATION","keyword for old version of the code that is there to maintain back compatibility only. Adding this keyword does nothing");
  keys.setValueDescription("scalar/grid","the value of the average");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("ACCUMULATE");
}

Average::Average( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao) {

  std::string lw;
  parse("LOGWEIGHTS",lw);
  std::string stride, clearstride;
  parse("STRIDE",stride);
  parse("CLEAR",clearstride);
  if( lw.length()>0 ) {
    readInputLine( getShortcutLabel() + "_wsum: COMBINE ARG=" + lw + " PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_weight: CUSTOM ARG=" + getShortcutLabel() + "_wsum FUNC=exp(x) PERIODIC=NO");
  } else {
    readInputLine( getShortcutLabel() + "_weight: ONES SIZE=1" );
  }

  std::vector<std::string> arg;
  parseVector("ARG",arg);
  if( arg.size()!=1 ) {
    error("should only be one argument to this action");
  }
  std::vector<Value*> vals;
  ActionWithArguments::interpretArgumentList( arg, plumed.getActionSet(), this, vals );

  readInputLine( getShortcutLabel() + "_denom: ACCUMULATE ARG=" + getShortcutLabel() + "_weight STRIDE=" + stride + " CLEAR=" + clearstride );
  if( vals[0]->isPeriodic() ) {
    std::string lbound, ubound, pfactor;
    vals[0]->getDomain( lbound, ubound );
    pfactor = "((" + ubound + "-" + lbound + ")/(pi+pi))";
    readInputLine( getShortcutLabel() + "_sin: CUSTOM ARG=" + arg[0] + "," + getShortcutLabel() + "_weight FUNC=y*sin((x-" + lbound + ")/" + pfactor + ") PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_cos: CUSTOM ARG=" + arg[0] + "," + getShortcutLabel() + "_weight FUNC=y*cos((x-" + lbound + ")/" + pfactor + ") PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sinsum: ACCUMULATE ARG=" + getShortcutLabel() + "_sin STRIDE=" + stride + " CLEAR=" + clearstride );
    readInputLine( getShortcutLabel() + "_cossum: ACCUMULATE ARG=" + getShortcutLabel() + "_cos STRIDE=" + stride + " CLEAR=" + clearstride );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_sinsum," + getShortcutLabel() + "_cossum," + getShortcutLabel() + "_denom FUNC=" + lbound + "+" + pfactor + "*atan2(x/z,y/z) PERIODIC=" + lbound +"," + ubound);
  } else {
    std::string normstr;
    parse("NORMALIZATION",normstr);
    if( normstr=="true" || normstr=="false" ) {
      warning("NORMALIZATION is deprecated. You are advised to take this out of input files in future and use the new syntax with ACCUMULATE for unormalized data rather than the shortcut AVERAGE");
    } else if( normstr.length()>0 ) {
      error("NORMALIZATION=" + normstr + " is not valid PLUMED input.  If you want an unormalised 'average' use ACCUMULATE");
    }
    readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + arg[0] + "," + getShortcutLabel() + "_weight FUNC=x*y PERIODIC=NO");
    if( normstr.length()==0 || normstr=="true" ) {
      readInputLine( getShortcutLabel() + "_numer: ACCUMULATE ARG=" + getShortcutLabel() + "_prod STRIDE=" + stride + " CLEAR=" + clearstride  );
      readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_numer," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    } else if( normstr=="false" ) {
      readInputLine( getShortcutLabel() + ": ACCUMULATE ARG=" + getShortcutLabel() + "_prod STRIDE=" + stride + " CLEAR=" + clearstride  );
    } else {
      plumed_error();
    }
  }
}

}
}
