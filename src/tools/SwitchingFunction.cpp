/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "SwitchingFunction.h"
#include "Tools.h"
#include "Keywords.h"
#include "OpenMP.h"
#include <vector>
#include <limits>
#include <algorithm>
#include <optional>

namespace PLMD {

//+PLUMEDOC INTERNAL switchingfunction
/*
Functions that measure whether values are less than a certain quantity.

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
   0                                                           & \mathrm{if } r < d_0 + r_0
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

*/
//+ENDPLUMEDOC

namespace switchContainers {

baseSwitch::baseSwitch(double D0,double DMAX, double R0, std::string_view name)
  : d0(D0),
    dmax(DMAX),
    dmax_2([](const double d) {
  if(d<std::sqrt(std::numeric_limits<double>::max())) {
    return  d*d;
  } else {
    return std::numeric_limits<double>::max();
  }
}(dmax)),
invr0(1.0/R0),
invr0_2(invr0*invr0),
mytype(name) {}

baseSwitch::~baseSwitch()=default;

double baseSwitch::calculate(const double distance, double& dfunc) const {
  double res = 0.0;//RVO!
  dfunc = 0.0;
  if(distance <= dmax) {
    res = 1.0;
    const double rdist = (distance-d0)*invr0;
    if(rdist > 0.0) {
      res = function(rdist,dfunc);
      //the following comments came from the original
      // this is for the chain rule (derivative of rdist):
      dfunc *= invr0;
      // for any future switching functions, be aware that multiplying invr0 is only
      // correct for functions of rdist = (r-d0)/r0.

      // this is because calculate() sets dfunc to the derivative divided times the
      // distance.
      // (I think this is misleading and I would like to modify it - GB)
      dfunc /= distance;
    }
    res=res*stretch+shift;
    dfunc*=stretch;
  }
  return res;
}

double baseSwitch::calculateSqr(double distance2,double&dfunc) const {
  double res= calculate(std::sqrt(distance2),dfunc);//RVO!
  return res;
}
double baseSwitch::get_d0() const {
  return d0;
}
double baseSwitch::get_r0() const {
  return 1.0/invr0;
}
double baseSwitch::get_dmax() const {
  return dmax;
}
double baseSwitch::get_dmax2() const {
  return dmax_2;
}
std::string baseSwitch::description() const {
  std::ostringstream ostr;
  ostr<<get_r0()
      <<".  Using "
      << mytype
      <<" switching function with parameters d0="<< d0
      << specificDescription();
  return ostr.str();
}
std::string baseSwitch::specificDescription() const {
  return "";
}
void baseSwitch::setupStretch() {
  if(dmax!=std::numeric_limits<double>::max()) {
    stretch=1.0;
    shift=0.0;
    double dummy;
    double s0=calculate(0.0,dummy);
    double sd=calculate(dmax,dummy);
    stretch=1.0/(s0-sd);
    shift=-sd*stretch;
  }
}
void baseSwitch::removeStretch() {
  stretch=1.0;
  shift=0.0;
}
template<int N, std::enable_if_t< (N >0), bool> = true, std::enable_if_t< (N %2 == 0), bool> = true>
    class fixedRational :public baseSwitch {
  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr << " nn=" << N << " mm=" <<N*2;
    return ostr.str();
  }
public:
  fixedRational(double D0,double DMAX, double R0)
    :baseSwitch(D0,DMAX,R0,"rational") {}

  template <int POW>
  static inline double doRational(const double rdist, double&dfunc, double result=0.0) {
    const double rNdist=Tools::fastpow<POW-1>(rdist);
    result=1.0/(1.0+rNdist*rdist);
    dfunc = -POW*rNdist*result*result;
    return result;
  }

  inline double function(double rdist,double&dfunc) const override {
    //preRes and preDfunc are passed already set
    dfunc=0.0;
    double result = doRational<N>(rdist,dfunc);
    return result;
  }

  double calculateSqr(double distance2,double&dfunc) const override {
    double result=0.0;
    dfunc=0.0;
    if(distance2 <= dmax_2) {
      const double rdist = distance2*invr0_2;
      result = doRational<N/2>(rdist,dfunc);
      dfunc*=2*invr0_2;
      // stretch:
      result=result*stretch+shift;
      dfunc*=stretch;
    }
    return result;

  }
};

//these enums are useful for clarifying the settings in the factory
//and the code is autodocumented ;)
enum class rationalPow:bool {standard, fast};
enum class rationalForm:bool {standard, simplified};

template<rationalPow isFast, rationalForm nis2m>
class rational : public baseSwitch {
protected:
  const int nn=6;
  const int mm=12;
  const double preRes;
  const double preDfunc;
  const double preSecDev;
  const int nnf;
  const int mmf;
  const double preDfuncF;
  const double preSecDevF;
  //I am using PLMD::epsilon to be certain to call the one defined in Tools.h
  static constexpr double moreThanOne=1.0+5.0e10*PLMD::epsilon;
  static constexpr double lessThanOne=1.0-5.0e10*PLMD::epsilon;

  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr << " nn=" << nn << " mm=" <<mm;
    return ostr.str();
  }
public:
  rational(double D0,double DMAX, double R0, int N, int M)
    :baseSwitch(D0,DMAX,R0,"rational"),
     nn(N),
     mm([](int m,int n) {
    if (m==0) {
      return n*2;
    } else {
      return m;
    }
  }(M,N)),
  preRes(static_cast<double>(nn)/mm),
  preDfunc(0.5*nn*(nn-mm)/static_cast<double>(mm)),
  //wolfram <3:lim_(x->1) d^2/(dx^2) (1 - x^N)/(1 - x^M) = (N (M^2 - 3 M (-1 + N) + N (-3 + 2 N)))/(6 M)
  preSecDev ((nn * (mm * mm - 3.0* mm * (-1 + nn ) + nn *(-3 + 2* nn )))/(6.0* mm )),
  nnf(nn/2),
  mmf(mm/2),
  preDfuncF(0.5*nnf*(nnf-mmf)/static_cast<double>(mmf)),
  preSecDevF((nnf* (mmf*mmf - 3.0* mmf* (-1 + nnf) + nnf*(-3 + 2* nnf)))/(6.0* mmf)) {}

  static inline double doRational(const double rdist, double&dfunc,double secDev, const int N,
                                  const int M,double result=0.0) {
    //the result and dfunc are assigned in the drivers for doRational
    //if(rdist>(1.0-100.0*epsilon) && rdist<(1.0+100.0*epsilon)) {
    //result=preRes;
    //dfunc=preDfunc;
    //} else {
    if constexpr (nis2m==rationalForm::simplified) {
      const double rNdist=Tools::fastpow(rdist,N-1);
      result=1.0/(1.0+rNdist*rdist);
      dfunc = -N*rNdist*result*result;
    } else {
      if(!((rdist > lessThanOne) && (rdist < moreThanOne))) {
        const double rNdist=Tools::fastpow(rdist,N-1);
        const double rMdist=Tools::fastpow(rdist,M-1);
        const double num = 1.0-rNdist*rdist;
        const double iden = 1.0/(1.0-rMdist*rdist);
        result = num*iden;
        dfunc = ((M*result*rMdist)-(N*rNdist))*iden;
      } else {
        //here I imply that the correct initialized are being passed to doRational
        const double x =(rdist-1.0);
        result = result+ x * ( dfunc + 0.5 * x * secDev);
        dfunc  = dfunc + x * secDev;
      }
    }
    return result;
  }
  inline double function(double rdist,double&dfunc) const override {
    //preRes and preDfunc are passed already set
    dfunc=preDfunc;
    double result = doRational(rdist,dfunc,preSecDev,nn,mm,preRes);
    return result;
  }

  double calculateSqr(double distance2,double&dfunc) const override {
    if constexpr (isFast==rationalPow::fast) {
      double result=0.0;
      dfunc=0.0;
      if(distance2 <= dmax_2) {
        const double rdist = distance2*invr0_2;
        dfunc=preDfuncF;
        result = doRational(rdist,dfunc,preSecDevF,nnf,mmf,preRes);
        dfunc*=2*invr0_2;
// stretch:
        result=result*stretch+shift;
        dfunc*=stretch;
      }
      return result;
    } else {
      double res= calculate(std::sqrt(distance2),dfunc);//RVO!
      return res;
    }
  }
};


template<int EXP,std::enable_if_t< (EXP %2 == 0), bool> = true>
std::optional<std::unique_ptr<baseSwitch>> fixedRationalFactory(double D0,double DMAX, double R0, int N) {
  if constexpr (EXP == 0) {
    return  std::nullopt;
  } else {
    if (N==EXP) {
      return PLMD::Tools::make_unique<switchContainers::fixedRational<EXP>>(D0,DMAX,R0);
    } else {
      return fixedRationalFactory<EXP-2>(D0,DMAX,R0,N);
    }
  }
}

std::unique_ptr<baseSwitch>
rationalFactory(double D0,double DMAX, double R0, int N, int M) {
  bool fast = N%2==0 && M%2==0 && D0==0.0;
  //if (M==0) M will automatically became 2*NN
  constexpr int highestPrecompiledPower=12;
  //precompiled rational
  if(((2*N)==M || M == 0) && fast && N<=highestPrecompiledPower) {
    auto tmp = fixedRationalFactory<highestPrecompiledPower>(D0,DMAX,R0,N);
    if(tmp) {
      return std::move(*tmp);
    }
    //else continue with the at runtime implementation
  }
  //template<bool isFast, bool n2m>
  //class rational : public baseSwitch
  if(2*N==M || M == 0) {
    if(fast) {
      //fast rational
      return PLMD::Tools::make_unique<switchContainers::rational<
             rationalPow::fast,rationalForm::simplified>>(D0,DMAX,R0,N,M);
    }
    return PLMD::Tools::make_unique<switchContainers::rational<
           rationalPow::standard,rationalForm::simplified>>(D0,DMAX,R0,N,M);
  }
  if(fast) {
    //fast rational
    return PLMD::Tools::make_unique<switchContainers::rational<
           rationalPow::fast,rationalForm::standard>>(D0,DMAX,R0,N,M);
  }
  return PLMD::Tools::make_unique<switchContainers::rational<
         rationalPow::standard,rationalForm::standard>>(D0,DMAX,R0,N,M);
}
//function =

class exponentialSwitch: public baseSwitch {
public:
  exponentialSwitch(double D0, double DMAX, double R0)
    :baseSwitch(D0,DMAX,R0,"exponential") {}
protected:
  inline double function(const double rdist,double&dfunc) const override {
    double result = std::exp(-rdist);
    dfunc=-result;
    return result;
  }
};

class gaussianSwitch: public baseSwitch {
public:
  gaussianSwitch(double D0, double DMAX, double R0)
    :baseSwitch(D0,DMAX,R0,"gaussian") {}
protected:
  inline double function(const double rdist,double&dfunc) const override {
    double result = std::exp(-0.5*rdist*rdist);
    dfunc=-rdist*result;
    return result;
  }
};

class fastGaussianSwitch: public baseSwitch {
public:
  fastGaussianSwitch(double /*D0*/, double DMAX, double /*R0*/)
    :baseSwitch(0.0,DMAX,1.0,"fastgaussian") {}
protected:
  inline double function(const double rdist,double&dfunc) const override {
    double result = std::exp(-0.5*rdist*rdist);
    dfunc=-rdist*result;
    return result;
  }
  inline double calculateSqr(double distance2,double&dfunc) const override {
    double result = 0.0;
    if(distance2>dmax_2) {
      dfunc=0.0;
    } else  {
      result = exp(-0.5*distance2);
      dfunc = -result;
      // stretch:
      result=result*stretch+shift;
      dfunc*=stretch;
    }
    return result;
  }
};

class smapSwitch: public baseSwitch {
  const int a=0;
  const int b=0;
  const double c=0.0;
  const double d=0.0;
protected:
  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr<<" a="<<a<<" b="<<b;
    return ostr.str();
  }
public:
  smapSwitch(double D0, double DMAX, double R0, int A, int B)
    :baseSwitch(D0,DMAX,R0,"smap"),
     a(A),
     b(B),
     c(std::pow(2., static_cast<double>(a)/static_cast<double>(b) ) - 1.0),
     d(-static_cast<double>(b) / static_cast<double>(a)) {}
protected:
  inline double function(const double rdist,double&dfunc) const override {

    const double sx=c*Tools::fastpow( rdist, a );
    double result=std::pow( 1.0 + sx, d );
    dfunc=-b*sx/rdist*result/(1.0+sx);
    return result;
  }
};

class cubicSwitch: public baseSwitch {
protected:
  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr<<" dmax="<<dmax;
    return ostr.str();
  }
public:
  cubicSwitch(double D0, double DMAX)
    :baseSwitch(D0,DMAX,DMAX-D0,"cubic") {
    //this operation should be already done!!
    // R0 = dmax - d0;
    // invr0 = 1/R0;
    // invr0_2 = invr0*invr0;
  }
  ~cubicSwitch()=default;
protected:
  inline double function(const double rdist,double&dfunc) const override {
    const double tmp1 = rdist - 1.0;
    const double tmp2 = 1.0+2.0*rdist;
    //double result = tmp1*tmp1*tmp2;
    dfunc = 2*tmp1*tmp2 + 2*tmp1*tmp1;
    return tmp1*tmp1*tmp2;
  }
};

class tanhSwitch: public baseSwitch {
public:
  tanhSwitch(double D0, double DMAX, double R0)
    :baseSwitch(D0,DMAX,R0,"tanh") {}
protected:
  inline double function(const double rdist,double&dfunc) const override {
    const double tmp1 = std::tanh(rdist);
    //was dfunc=-(1-tmp1*tmp1);
    dfunc = tmp1 * tmp1 - 1.0;
    //return result;
    return 1.0 - tmp1;
  }
};

class cosinusSwitch: public baseSwitch {
public:
  cosinusSwitch(double D0, double DMAX, double R0)
    :baseSwitch(D0,DMAX,R0,"cosinus") {}
protected:
  inline double function(const double rdist,double&dfunc) const override {
    double result = 0.0;
    dfunc=0.0;
    if(rdist<=1.0) {
// rdist = (r-r1)/(r2-r1) ; 0.0<=rdist<=1.0 if r1 <= r <=r2; (r2-r1)/(r2-r1)=1
      double rdistPI = rdist * PLMD::pi;
      result = 0.5 * (std::cos ( rdistPI ) + 1.0);
      dfunc = -0.5 * PLMD::pi * std::sin ( rdistPI ) * invr0;
    }
    return result;
  }
};

class nativeqSwitch: public baseSwitch {
  double beta = 50.0;  // nm-1
  double lambda = 1.8; // unitless
  double ref=0.0;
protected:
  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr<<" beta="<<beta<<" lambda="<<lambda<<" ref="<<ref;
    return ostr.str();
  }
  inline double function(const double rdist,double&dfunc) const override {
    return 0.0;
  }
public:
  nativeqSwitch(double D0, double DMAX, double R0, double BETA, double LAMBDA,double REF)
    :  baseSwitch(D0,DMAX,R0,"nativeq"),beta(BETA),lambda(LAMBDA),ref(REF) {}
  double calculate(const double distance, double& dfunc) const override {
    double res = 0.0;//RVO!
    dfunc = 0.0;
    if(distance<=dmax) {
      res = 1.0;
      if(distance > d0) {
        const double rdist = beta*(distance - lambda * ref);
        double exprdist=std::exp(rdist);
        res=1.0/(1.0+exprdist);
        /*2.9
        //need to see if this (5op+assign)
        //double exprmdist=1.0 + exprdist;
        //dfunc = - (beta *exprdist)/(exprmdist*exprmdist);
        //or this (5op but 2 divisions) is faster
        dfunc = - beta /(exprdist+ 2.0 +1.0/exprdist);
        //this cames from - beta * exprdist/(exprdist*exprdist+ 2.0 *exprdist +1.0)
        //dfunc *= invr0;
        dfunc /= distance;
        */
        //2.10
        dfunc = - beta /(exprdist+ 2.0 +1.0/exprdist) /distance;

        dfunc*=stretch;
      }
      res=res*stretch+shift;
    }
    return res;
  }
};

class leptonSwitch: public baseSwitch {
/// Lepton expression.
  class funcAndDeriv {
    lepton::CompiledExpression expression;
    lepton::CompiledExpression deriv;
    double* varRef=nullptr;
    double* varDevRef=nullptr;
  public:
    funcAndDeriv(const std::string &func) {
      lepton::ParsedExpression pe=lepton::Parser::parse(func).optimize(lepton::Constants());
      expression=pe.createCompiledExpression();
      std::string arg="x";

      {
        auto vars=expression.getVariables();
        bool found_x=std::find(vars.begin(),vars.end(),"x")!=vars.end();
        bool found_x2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();

        if(found_x2) {
          arg="x2";
        }
        if (vars.size()==0) {
// this is necessary since in some cases lepton thinks a variable is not present even though it is present
// e.g. func=0*x
          varRef=nullptr;
        } else if(vars.size()==1 && (found_x || found_x2)) {
          varRef=&expression.getVariableReference(arg);
        } else {
          plumed_error()
              <<"Please declare a function with only ONE argument that can only be x or x2. Your function is: "
              << func;
        }
      }

      lepton::ParsedExpression ped=lepton::Parser::parse(func).differentiate(arg).optimize(lepton::Constants());
      deriv=ped.createCompiledExpression();
      {
        auto vars=expression.getVariables();
        if (vars.size()==0) {
          varDevRef=nullptr;
        } else {
          varDevRef=&deriv.getVariableReference(arg);
        }
      }

    }
    funcAndDeriv (const funcAndDeriv& other):
      expression(other.expression),
      deriv(other.deriv) {
      std::string arg="x";

      {
        auto vars=expression.getVariables();
        bool found_x=std::find(vars.begin(),vars.end(),"x")!=vars.end();
        bool found_x2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();

        if(found_x2) {
          arg="x2";
        }
        if (vars.size()==0) {
          varRef=nullptr;
        } else if(vars.size()==1 && (found_x || found_x2)) {
          varRef=&expression.getVariableReference(arg);
        }// UB: I assume that the function is already correct
      }

      {
        auto vars=expression.getVariables();
        if (vars.size()==0) {
          varDevRef=nullptr;
        } else {
          varDevRef=&deriv.getVariableReference(arg);
        }
      }
    }

    funcAndDeriv& operator= (const funcAndDeriv& other) {
      if(this != &other) {
        expression = other.expression;
        deriv = other.deriv;
        std::string arg="x";

        {
          auto vars=expression.getVariables();
          bool found_x=std::find(vars.begin(),vars.end(),"x")!=vars.end();
          bool found_x2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();

          if(found_x2) {
            arg="x2";
          }
          if (vars.size()==0) {
            varRef=nullptr;
          } else if(vars.size()==1 && (found_x || found_x2)) {
            varRef=&expression.getVariableReference(arg);
          }// UB: I assume that the function is already correct
        }

        {
          auto vars=expression.getVariables();
          if (vars.size()==0) {
            varDevRef=nullptr;
          } else {
            varDevRef=&deriv.getVariableReference(arg);
          }
        }
      }
      return *this;
    }

    std::pair<double,double> operator()(double const x) const {
      //FAQ: why this works? this thing is const and you are modifying things!
      //Actually I am modifying something that is pointed at, not my pointers,
      //so I am not mutating the state of this!
      if(varRef) {
        *varRef=x;
      }
      if(varDevRef) {
        *varDevRef=x;
      }
      return std::make_pair(
               expression.evaluate(),
               deriv.evaluate());
    }

    auto& getVariables() const {
      return expression.getVariables();
    }
    auto& getVariables_derivative() const {
      return deriv.getVariables();
    }
  };
  /// Function for lepton
  std::string lepton_func;
  /// \warning Since lepton::CompiledExpression is mutable, a vector is necessary for multithreading!
  std::vector <funcAndDeriv> expressions{};
  /// Set to true if lepton only uses x2
  bool leptonx2=false;
protected:
  std::string specificDescription() const override {
    std::ostringstream ostr;
    ostr<<" func=" << lepton_func;
    return ostr.str();
  }
  inline double function(const double,double&) const override {
    return 0.0;
  }
public:
  leptonSwitch(double D0, double DMAX, double R0, const std::string & func)
    :baseSwitch(D0,DMAX,R0,"lepton"),
     lepton_func(func),
     expressions  (OpenMP::getNumThreads(), lepton_func) {
    //this is a bit odd, but it works
    auto vars=expressions[0].getVariables();
    leptonx2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();
  }

  double calculate(const double distance,double&dfunc) const override {
    double res = 0.0;//RVO!
    dfunc = 0.0;
    if(leptonx2) {
      res= calculateSqr(distance*distance,dfunc);
    } else {
      if(distance<=dmax) {
        res = 1.0;
        const double rdist = (distance-d0)*invr0;
        if(rdist > 0.0) {
          const unsigned t=OpenMP::getThreadNum();
          plumed_assert(t<expressions.size());
          std::tie(res,dfunc) = expressions[t](rdist);
          dfunc *= invr0;
          dfunc /= distance;
        }
        res=res*stretch+shift;
        dfunc*=stretch;
      }
    }
    return res;
  }

  double calculateSqr(const double distance2,double&dfunc) const override {
    double result =0.0;
    dfunc=0.0;
    if(leptonx2) {
      if(distance2<=dmax_2) {
        const unsigned t=OpenMP::getThreadNum();
        const double rdist_2 = distance2*invr0_2;
        plumed_assert(t<expressions.size());
        std::tie(result,dfunc) = expressions[t](rdist_2);
        // chain rule:
        dfunc*=2*invr0_2;
        // stretch:
        result=result*stretch+shift;
        dfunc*=stretch;
      }
    } else {
      result = calculate(std::sqrt(distance2),dfunc);
    }
    return result;
  }
};
} // namespace switchContainers

void SwitchingFunction::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","R_0","the value of R_0 in the switching function");
  keys.add("compulsory","D_0","0.0","the value of D_0 in the switching function");
  keys.add("optional","D_MAX","the value at which the switching function can be assumed equal to zero");
  keys.add("compulsory","NN","6","the value of n in the switching function (only needed for TYPE=RATIONAL)");
  keys.add("compulsory","MM","0","the value of m in the switching function (only needed for TYPE=RATIONAL); 0 implies 2*NN");
  keys.add("compulsory","A","the value of a in the switching function (only needed for TYPE=SMAP)");
  keys.add("compulsory","B","the value of b in the switching function (only needed for TYPE=SMAP)");
}

void SwitchingFunction::set(const std::string & definition,std::string& errormsg) {
  std::vector<std::string> data=Tools::getWords(definition);
#define CHECKandPARSE(datastring,keyword,variable,errormsg) \
  if(Tools::findKeyword(datastring,keyword) && !Tools::parse(datastring,keyword,variable))\
    errormsg="could not parse " keyword; //adiacent strings are automagically concatenated
#define REQUIREDPARSE(datastring,keyword,variable,errormsg) \
  if(!Tools::parse(datastring,keyword,variable))\
    errormsg=keyword " is required for " + name ; //adiacent strings are automagically concatenated

  if( data.size()<1 ) {
    errormsg="missing all input for switching function";
    return;
  }
  std::string name=data[0];
  data.erase(data.begin());
  double r0=0.0;
  double d0=0.0;
  double dmax=std::numeric_limits<double>::max();
  init=true;
  CHECKandPARSE(data,"D_0",d0,errormsg);
  CHECKandPARSE(data,"D_MAX",dmax,errormsg);

  bool dostretch=false;
  Tools::parseFlag(data,"STRETCH",dostretch); // this is ignored now
  dostretch=true;
  bool dontstretch=false;
  Tools::parseFlag(data,"NOSTRETCH",dontstretch); // this is ignored now
  if(dontstretch) {
    dostretch=false;
  }
  if(name=="CUBIC") {
    //cubic is the only switch type that only uses d0 and dmax
    function = PLMD::Tools::make_unique<switchContainers::cubicSwitch>(d0,dmax);
  } else {
    REQUIREDPARSE(data,"R_0",r0,errormsg);
    if(name=="RATIONAL") {
      int nn=6;
      int mm=0;
      CHECKandPARSE(data,"NN",nn,errormsg);
      CHECKandPARSE(data,"MM",mm,errormsg);
      function = switchContainers::rationalFactory(d0,dmax,r0,nn,mm);
    } else if(name=="SMAP") {
      int a=0;
      int b=0;
      //in the original a and b are "default=0",
      //but you divide by a and b during the initialization!
      //better an error message than an UB, so no default
      REQUIREDPARSE(data,"A",a,errormsg);
      REQUIREDPARSE(data,"B",b,errormsg);
      function = PLMD::Tools::make_unique<switchContainers::smapSwitch>(d0,dmax,r0,a,b);
    } else if(name=="Q") {
      double beta = 50.0;  // nm-1
      double lambda = 1.8; // unitless
      double ref;
      CHECKandPARSE(data,"BETA",beta,errormsg);
      CHECKandPARSE(data,"LAMBDA",lambda,errormsg);
      REQUIREDPARSE(data,"REF",ref,errormsg);
      //the original error message was not standard
      // if(!Tools::parse(data,"REF",ref))
      //   errormsg="REF (reference distaance) is required for native Q";
      function = PLMD::Tools::make_unique<switchContainers::nativeqSwitch>(d0,dmax,r0,beta,lambda,ref);
    } else if(name=="EXP") {
      function = PLMD::Tools::make_unique<switchContainers::exponentialSwitch>(d0,dmax,r0);
    } else if(name=="GAUSSIAN") {
      if ( r0==1.0 && d0==0.0 ) {
        function = PLMD::Tools::make_unique<switchContainers::fastGaussianSwitch>(d0,dmax,r0);
      } else {
        function = PLMD::Tools::make_unique<switchContainers::gaussianSwitch>(d0,dmax,r0);
      }
    } else if(name=="TANH") {
      function = PLMD::Tools::make_unique<switchContainers::tanhSwitch>(d0,dmax,r0);
    } else if(name=="COSINUS") {
      function = PLMD::Tools::make_unique<switchContainers::cosinusSwitch>(d0,dmax,r0);
    } else if((name=="MATHEVAL" || name=="CUSTOM")) {
      std::string func;
      Tools::parse(data,"FUNC",func);
      function = PLMD::Tools::make_unique<switchContainers::leptonSwitch>(d0,dmax,r0,func);
    } else {
      errormsg="cannot understand switching function type '"+name+"'";
    }
  }
#undef CHECKandPARSE
#undef REQUIREDPARSE

  if( !data.empty() ) {
    errormsg="found the following rogue keywords in switching function input : ";
    for(unsigned i=0; i<data.size(); ++i) {
      errormsg = errormsg + data[i] + " ";
    }
  }

  if(dostretch && dmax!=std::numeric_limits<double>::max()) {
    function->setupStretch();
  }
}

std::string SwitchingFunction::description() const {
  // if this error is necessary, something went wrong in the constructor
  //  plumed_merror("Unknown switching function type");
  return function->description();
}

double SwitchingFunction::calculateSqr(double distance2,double&dfunc)const {
  return function -> calculateSqr(distance2, dfunc);
}

double SwitchingFunction::calculate(double distance,double&dfunc)const {
  plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  double result=function->calculate(distance,dfunc);
  return result;
}

void SwitchingFunction::set(const int nn,int mm, const double r0, const double d0) {
  init=true;
  if(mm == 0) {
    mm = 2*nn;
  }
  double dmax=d0+r0*std::pow(0.00001,1./(nn-mm));
  function = switchContainers::rationalFactory(d0,dmax,r0,nn,mm);
  function->setupStretch();
}

double SwitchingFunction::get_r0() const {
  return function->get_r0();
}

double SwitchingFunction::get_d0() const {
  return function->get_d0();
}

double SwitchingFunction::get_dmax() const {
  return function->get_dmax();
}

double SwitchingFunction::get_dmax2() const {
  return function->get_dmax2();
}

}// Namespace PLMD
