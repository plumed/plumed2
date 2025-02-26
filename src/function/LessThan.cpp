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
#include "LessThan.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION LESS_THAN
/*
Use a switching function to determine how many of the input variables are less than a certain cutoff.

Switching functions \f$s(r)\f$ take a minimum of one input parameter \f$r_0\f$.
For \f$r \le d_0 \quad s(r)=1.0\f$ while for \f$r > d_0\f$ the function decays smoothly to 0.
The various switching functions available in PLUMED differ in terms of how this decay is performed.

Where there is an accepted convention in the literature (e.g. \ref COORDINATION) on the form of the
switching function we use the convention as the default.  However, the flexibility to use different
switching functions is always present generally through a single keyword. This keyword generally
takes an input with the following form:

\verbatim
KEYWORD={TYPE <list of parameters>}
\endverbatim

The following table contains a list of the various switching functions that are available in PLUMED 2
together with an example input.

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td> TYPE </td> <td> FUNCTION </td> <td> EXAMPLE INPUT </td> <td> DEFAULT PARAMETERS </td>
</tr> <tr> <td>RATIONAL </td> <td>
\f$
s(r)=\frac{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{n} }{ 1 - \left(\frac{ r - d_0 }{ r_0 }\right)^{m} }
\f$
</td> <td>
{RATIONAL R_0=\f$r_0\f$ D_0=\f$d_0\f$ NN=\f$n\f$ MM=\f$m\f$}
</td> <td> \f$d_0=0.0\f$, \f$n=6\f$, \f$m=2n\f$ </td>
</tr> <tr>
<td> EXP </td> <td>
\f$
s(r)=\exp\left(-\frac{ r - d_0 }{ r_0 }\right)
\f$
</td> <td>
{EXP  R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr>
<td> GAUSSIAN </td> <td>
\f$
s(r)=\exp\left(-\frac{ (r - d_0)^2 }{ 2r_0^2 }\right)
\f$
</td> <td>
{GAUSSIAN R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr>
<td> SMAP </td> <td>
\f$
s(r) = \left[ 1 + ( 2^{a/b} -1 )\left( \frac{r-d_0}{r_0} \right)^a \right]^{-b/a}
\f$
</td> <td>
{SMAP R_0=\f$r_0\f$ D_0=\f$d_0\f$ A=\f$a\f$ B=\f$b\f$}
</td> <td> \f$d_0=0.0\f$ </td>
</tr> <tr>
<td> Q </td> <td>
\f$
s(r) = \frac{1}{1 + \exp(\beta(r_{ij} - \lambda r_{ij}^0))}
\f$
</td> <td>
{Q REF=\f$r_{ij}^0\f$ BETA=\f$\beta\f$ LAMBDA=\f$\lambda\f$ }
</td> <td> \f$\lambda=1.8\f$,  \f$\beta=50 nm^-1\f$ (all-atom)<br/>\f$\lambda=1.5\f$,  \f$\beta=50 nm^-1\f$ (coarse-grained)  </td>
</tr> <tr>
<td> CUBIC </td> <td>
\f$
s(r) = (y-1)^2(1+2y) \qquad \textrm{where} \quad y = \frac{r - r_1}{r_0-r_1}
\f$
</td> <td>
{CUBIC D_0=\f$r_1\f$ D_MAX=\f$r_0\f$}
</td> <td> </td>
</tr> <tr>
<td> TANH </td> <td>
\f$
s(r) = 1 - \tanh\left( \frac{ r - d_0 }{ r_0 } \right)
\f$
</td> <td>
{TANH R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> </td>
</tr> <tr>
<td> COSINUS </td> <td>
\f$s(r) =\left\{\begin{array}{ll}
   1                                                           & \mathrm{if } r \leq d_0 \\
   0.5 \left( \cos ( \frac{ r - d_0 }{ r_0 } \pi ) + 1 \right) & \mathrm{if } d_0 < r\leq d_0 + r_0 \\
   0                                                           & \mathrm{if } r > d_0 + r_0
  \end{array}\right.
\f$
</td> <td>
{COSINUS R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> </td>
</tr> <tr>
<td> CUSTOM </td> <td>
\f$
s(r) = FUNC
\f$
</td> <td>
{CUSTOM FUNC=1/(1+x^6) R_0=\f$r_0\f$ D_0=\f$d_0\f$}
</td> <td> </td>
</tr>
</table>

Notice that most commonly used rational functions are better optimized and might run faster.

Notice that for backward compatibility we allow using `MATHEVAL` instead of `CUSTOM`.
Also notice that if the a `CUSTOM` switching function only depends on even powers of `x` it can be
made faster by using `x2` as a variable. For instance
\verbatim
{CUSTOM FUNC=1/(1+x2^3) R_0=0.3}
\endverbatim
is equivalent to
\verbatim
{CUSTOM FUNC=1/(1+x^6) R_0=0.3}
\endverbatim
but runs faster. The reason is that there is an expensive square root calculation that can be optimized out.


\attention
With the default implementation CUSTOM is slower than other functions
(e.g., it is slower than an equivalent RATIONAL function by approximately a factor 2).
Checkout page \ref Lepton to see how to improve its performance.

For all the switching functions in the above table one can also specify a further (optional) parameter using the parameter
keyword D_MAX to assert that for \f$r>d_{\textrm{max}}\f$ the switching function can be assumed equal to zero.
In this case the function is brought smoothly to zero by stretching and shifting it.
\verbatim
KEYWORD={RATIONAL R_0=1 D_MAX=3}
\endverbatim
the resulting switching function will be
\f$
s(r) = \frac{s'(r)-s'(d_{max})}{s'(0)-s'(d_{max})}
\f$
where
\f$
s'(r)=\frac{1-r^6}{1-r^{12}}
\f$
Since PLUMED 2.2 this is the default. The old behavior (no stretching) can be obtained with the
NOSTRETCH flag. The NOSTRETCH keyword is only provided for backward compatibility and might be
removed in the future. Similarly, the STRETCH keyword is still allowed but has no effect.

Notice that switching functions defined with the simplified syntax are never stretched
for backward compatibility. This might change in the future.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION LESS_THAN_VECTOR
/*
Use a switching function to determine how many components of the vector are less than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR LESS_THAN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is one when the input value is smaller than a threshold

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<LessThan> LessThanShortcut;
PLUMED_REGISTER_ACTION(LessThanShortcut,"LESS_THAN")
typedef FunctionOfVector<LessThan> VectorLessThan;
PLUMED_REGISTER_ACTION(VectorLessThan,"LESS_THAN_VECTOR")
typedef FunctionOfMatrix<LessThan> MatrixLessThan;
PLUMED_REGISTER_ACTION(MatrixLessThan,"LESS_THAN_MATRIX")

void LessThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the input is less than a threshold");
}

void LessThan::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }


  std::string sw,errors;
  action->parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    action->parse("R_0",r0);
    if(r0<=0.0) {
      action->error("R_0 should be explicitly specified and positive");
    }
    action->parse("D_0",d0);
    action->parse("NN",nn);
    action->parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",squared);
  if( squared ) {
    action->log<<"  input quantity is square of quantity that switching function acts upon\n";
  }
}

void LessThan::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  if( squared ) {
    vals[0] = switchingFunction.calculateSqr( args[0], derivatives(0,0) );
  } else {
    vals[0] = switchingFunction.calculate( args[0], derivatives(0,0) );
  }
  derivatives(0,0) = args[0]*derivatives(0,0);
}

}
}


