/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Between.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BETWEEN
/*
Use a switching function to determine how many of the input variables are within a certain range.

If we have multiple instances of a variable we can estimate the probability density function
for that variable using a process called kernel density estimation:

\f[
P(s) = \sum_i K\left( \frac{s - s_i}{w} \right)
\f]

In this equation \f$K\f$ is a symmetric function that must integrate to one that is often
called a kernel function and \f$w\f$ is a smearing parameter.  From a probability density function calculated using
kernel density estimation we can calculate the number/fraction of values between an upper and lower
bound using:

\f[
w(s) = \int_a^b \sum_i K\left( \frac{s - s_i}{w} \right)
\f]

All the input to calculate a quantity like \f$w(s)\f$ is generally provided through a single
keyword that will have the following form:

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ SMEAR=\f$\frac{w}{b-a}\f$}

This will calculate the number of values between \f$a\f$ and \f$b\f$.  To calculate
the fraction of values you add the word NORM to the input specification.  If the
function keyword SMEAR is not present \f$w\f$ is set equal to \f$0.5(b-a)\f$. Finally,
type should specify one of the kernel types that is present in plumed. These are listed
in the table below:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td> TYPE </td> <td> FUNCTION </td>
</tr> <tr>
<td> GAUSSIAN </td> <td> \f$\frac{1}{\sqrt{2\pi}w} \exp\left( -\frac{(s-s_i)^2}{2w^2} \right)\f$ </td>
</tr> <tr>
<td> TRIANGULAR </td> <td> \f$ \frac{1}{2w} \left( 1. - \left| \frac{s-s_i}{w} \right| \right) \quad \frac{s-s_i}{w}<1 \f$ </td>
</tr>
</table>

Some keywords can also be used to calculate a discrete version of the histogram.  That
is to say the number of values between \f$a\f$ and \f$b\f$, the number of values between
\f$b\f$ and \f$c\f$ and so on.  A keyword that specifies this sort of calculation would look
something like

KEYWORD={TYPE UPPER=\f$a\f$ LOWER=\f$b\f$ NBINS=\f$n\f$ SMEAR=\f$\frac{w}{n(b-a)}\f$}

This specification would calculate the following vector of quantities:

\f[
w_j(s) = \int_{a + \frac{j-1}{n}(b-a)}^{a + \frac{j}{n}(b-a)} \sum_i K\left( \frac{s - s_i}{w} \right)
\f]


\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION BETWEEN_VECTOR
/*
Use a switching function to determine how many of the input components are within a certain range

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR BETWEEN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is oen when the input value is within a particular range

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Between> BetweenShortcut;
PLUMED_REGISTER_ACTION(BetweenShortcut,"BETWEEN")
typedef FunctionOfVector<Between> VectorBetween;
PLUMED_REGISTER_ACTION(VectorBetween,"BETWEEN_VECTOR")
typedef FunctionOfMatrix<Between> MatrixBetween;
PLUMED_REGISTER_ACTION(MatrixBetween,"BETWEEN_MATRIX")

void Between::registerKeywords(Keywords& keys) {
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER, SMEAR and KERNEL keywords.");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the input falls within a particular range and zero otherwise");
}

void Between::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    action->error("should only be one argument to between actions");
  }

  std::string str_min, str_max, tstr_min, tstr_max;
  bool isPeriodic = action->getPntrToArgument(0)->isPeriodic();
  if( isPeriodic ) {
    action->getPntrToArgument(0)->getDomain( str_min, str_max );
  }

  std::string hinput;
  action->parse("SWITCH",hinput);
  if(hinput.length()==0) {
    std::string low, up, sme;
    action->parse("LOWER",low);
    action->parse("UPPER",up);
    action->parse("SMEAR",sme);
    hinput = "GAUSSIAN LOWER=" + low + " UPPER=" + up + " SMEAR=" + sme;
  }
  std::string errors;
  hist.set( hinput, errors );
  if( errors.size()!=0 ) {
    action->error( errors );
  }
  action->log.printf("  %s \n", hist.description().c_str() );

  if( !isPeriodic ) {
    hist.isNotPeriodic();
  } else {
    double min;
    Tools::convert( str_min, min );
    double max;
    Tools::convert( str_max, max );
    hist.isPeriodic( min, max );
  }
}

void Between::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  vals[0] = hist.calculate( args[0], derivatives(0,0) );
}

}
}


