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
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

#pragma GCC diagnostic error "-Wswitch"
/*
IMPORTANT NOTE FOR DEVELOPERS:

If you add a new type of switching function in this file please add documentation for your new switching function type in function/LessThan.cpp
*/

namespace PLMD {
namespace switchContainers {

std::string description(switchType type, const Data& data);

Switch::~Switch()=default;

template <typename SF>
class SwitchInterface :public Switch {
  switchType type;
  Data data;
public:
  explicit SwitchInterface(const std::pair <switchType,Data>& d):
    type(d.first),
    data(d.second) {}
  double calculate(double distance, double& dfunc) const override {
    auto [f,d] = SF::calculate(data,distance);
    dfunc=d;
    return f;
  }
  double calculateSqr(double distance2, double& dfunc) const override {
    auto [f,d] = SF::calculateSqr(data,distance2);
    dfunc=d;
    return f;
  }
  void setupStretch() override {
    if(data.dmax<std::numeric_limits<double>::max()) {
      data.stretch=1.0;
      data.shift=0.0;
      double s0=SF::calculate(data,0.0).first;
      double sd=SF::calculate(data,data.dmax).first;
      data.stretch=1.0/(s0-sd);
      data.shift=-sd*data.stretch;
    }
  }
  const Data& getData() const override {
    return data;
  }
  switchType getType() const override {
    return type;
  }
  std::string description() const override {
    return switchContainers::description(type,data);
  }
};

Data Data::init(const double D0,const double DMAX, const double R0) {
  Data data  ;
  data.d0=D0;
  data.dmax=DMAX;
  data.dmax_2=(data.dmax<std::sqrt(std::numeric_limits<double>::max()))
              ? data.dmax*data.dmax
              : std::numeric_limits<double>::max();
  data.invr0=1.0/R0;
  data.invr0_2=data.invr0*data.invr0;
  return data;
}

void Data::toACCDevice() const {
#pragma acc enter data copyin(this[0:1], d0, dmax, dmax_2, invr0, invr0_2, \
                              stretch, shift, nn, mm, preRes, preDfunc, preSecDev, \
                              nnf, mmf, preDfuncF, preSecDevF, a, b, c, d, beta, \
                              lambda, ref )
}
void Data::removeFromACCDevice() const {
#pragma acc exit data delete(ref, lambda, beta, d, c, b, a, preSecDevF, preDfuncF, \
                             mmf, nnf, preSecDev, preDfunc, preRes, mm, nn, shift, \
                             stretch, invr0_2, invr0, dmax_2, dmax, d0, this[0:1])
}

std::string typeToString(switchType type);
namespace switchContainersUtils {
template<class, class = void>
constexpr bool has_function_data = false;

//this verifies that T has a method gatherForces_custom that can be called with this signature
template<class T>
constexpr bool has_function_data <T, std::void_t<
decltype(T::function(
           std::declval<const Data& >(),
           std::declval<double >()
         ))> > = true;
} //namespace switchContainersUtils

using ValueDerivative=std::pair<double,double>;

inline ValueDerivative applystretch(const Data&data,double distance,ValueDerivative in) {
  in.first=in.first*data.stretch+data.shift;
  in.second*=data.stretch;
  in.second *= data.invr0;
  in.second /= distance;
  return in;
}

/// container for the actual switching function used by PLMD::SwitchingFunction
template<typename switching>
struct baseSwitch {
  static ValueDerivative calculate(const Data&data, const double distance) {
    const double rdist = (distance-data.d0)*data.invr0;
    if constexpr (switchContainersUtils::has_function_data<switching>) {
      return (distance > data.dmax) ? ValueDerivative{0.0,0.0}
             :
             (rdist > 0.0) ?
             applystretch(data,distance,switching::function(data,rdist))
             : ValueDerivative{data.stretch+data.shift,0.0};
    } else {
      return  (distance > data.dmax) ? ValueDerivative{0.0,0.0}
              :
              (rdist > 0.0) ?
              applystretch(data,distance,switching::function(rdist))
              : ValueDerivative{data.stretch+data.shift,0.0};
    }
    /* //not branchless (legacy) implementation
        double res = 0.0;
        double dfunc = 0.0;

        if(distance <= data.dmax) {
          res = 1.0;
          const double rdist = (distance-data.d0)*data.invr0;
          if(rdist > 0.0) {
            if constexpr (switchContainersUtils::has_function_data<switching>) {
              std::tie(res,dfunc) = switching::function(data,rdist);
            } else {
              std::tie(res,dfunc) = switching::function(rdist);
            }
            //the following comments came from the original
            // this is for the chain rule (derivative of rdist):
            dfunc *= data.invr0;
            // for any future switching functions, be aware that multiplying invr0 is only
            // correct for functions of rdist = (r-d0)/r0.

            // this is because calculate() sets dfunc to the derivative divided times the
            // distance.
            // (I think this is misleading and I would like to modify it - GB)
            dfunc /= distance;
          }
          res=res*data.stretch+data.shift;
          dfunc*=data.stretch;
        }
    return {res,dfunc};
    */
  }

  static ValueDerivative calculateSqr(const Data&data, double distance2) {
    return switching::calculate(data,std::sqrt(distance2));
  }
};
template<int N,
         std::enable_if_t< (N >0), bool> = true,
         std::enable_if_t< (N %2 == 0), bool> = true>
             struct fixedRational :public baseSwitch<fixedRational<N>> {

  template <int POW>
  static inline ValueDerivative doRational(const double rdist, double dfunc=0.0, double result=0.0) {
    const double rNdist=Tools::fastpow<POW-1>(rdist);
    result=1.0/(1.0+rNdist*rdist);
    dfunc = -POW*rNdist*result*result;
    return {result,dfunc};
  }

  static inline ValueDerivative function(const Data&, double rdist) {
    //preRes and preDfunc are passed already set
    return doRational<N>(rdist);
  }

  static ValueDerivative calculateSqr(const Data& data,double distance2) {
    double result=0.0;
    double dfunc=0.0;
    if(distance2 <= data.dmax_2) {
      const double rdist = distance2*data.invr0_2;
      std::tie(result,dfunc) = doRational<N/2>(rdist);
      dfunc*=2*data.invr0_2;
      // stretch:
      result=result*data.stretch+data.shift;
      dfunc*=data.stretch;
    }
    return {result,dfunc};
  }
};

//these enums are useful for clarifying the settings in the factory
//and the code is autodocumented ;)
enum class rationalPow:bool {standard, fast};
enum class rationalForm:bool {standard, simplified};
template<rationalPow isFast, rationalForm nis2m>
struct rational : public baseSwitch<rational<isFast,nis2m>> {
  //I am using PLMD::epsilon to be certain to call the one defined in Tools.h
  static constexpr double moreThanOne=1.0+5.0e10*PLMD::epsilon;
  static constexpr double lessThanOne=1.0-5.0e10*PLMD::epsilon;
  static std::pair <switchType,Data> init(double D0,double DMAX, double R0, int N, int M) {
    auto data= Data::init(D0,DMAX,R0);
    data.nn=N;
    data.mm = (M==0) ? N*2 : M;
    data.preRes=static_cast<double>(data.nn)/data.mm;
    data.preDfunc=0.5*data.nn*(data.nn-data.mm)/static_cast<double>(data.mm);
    //wolfram <3:lim_(x->1) d^2/(dx^2) (1 - x^N)/(1 - x^M) = (N (M^2 - 3 M (-1 + N) + N (-3 + 2 N)))/(6 M);
    data.preSecDev = (data.nn * (data.mm * data.mm - 3.0* data.mm * (-1 + data.nn ) + data.nn *(-3 + 2* data.nn )))/(6.0* data.mm );
    data.nnf = data.nn/2;
    data.mmf = data.mm/2;
    data.preDfuncF = 0.5*data.nnf*(data.nnf-data.mmf)/static_cast<double>(data.mmf);
    data.preSecDevF = (data.nnf* (data.mmf*data.mmf - 3.0* data.mmf* (-1 + data.nnf) + data.nnf*(-3 + 2* data.nnf)))/(6.0* data.mmf);

    const bool fast = N%2==0 && M%2==0 && D0==0.0;
    if(2*N==M || M == 0) {
      if(fast) {
        return {switchType::rationalSimpleFast,data};
      }
      return {switchType::rationalSimple,data};
    }
    if(fast) {
      return {switchType::rationalFast,data};
    }
    return {switchType::rational,data};
  }

  static inline ValueDerivative doRational(const double rdist,double secDev, const int N,
      const int M, double dfunc=0.0,double result=0.0) {
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
    return {result,dfunc};
  }
  static inline ValueDerivative function(const Data&data, double rdist) {
    //preRes and preDfunc are passed already set
    return doRational(rdist,data.preSecDev,data.nn,data.mm,data.preDfunc,data.preRes);
  }

  static ValueDerivative calculateSqr(const Data&data, double distance2) {
    if constexpr (isFast==rationalPow::fast) {
      double result=0.0;
      double dfunc=0.0;
      if(distance2 <= data.dmax_2) {
        const double rdist = distance2*data.invr0_2;
        std::tie(result,dfunc) =
          doRational(rdist,data.preSecDevF,data.nnf,data.mmf,data.preDfuncF,data.preRes);
        dfunc*=2*data.invr0_2;
// stretch:
        result=result*data.stretch+data.shift;
        dfunc*=data.stretch;
      }
      return {result,dfunc};
    } else {
      return baseSwitch<rational<isFast,nis2m>>::calculate(data,std::sqrt(distance2));
    }
  }
};

std::unique_ptr<Switch>
rationalFactory (double D0,double DMAX, double R0, int N, int M) {
  bool fast = N%2==0 && M%2==0 && D0==0.0;
  //if (M==0) M will automatically became 2*NN
  constexpr int highestPrecompiledPower=12;
  //precompiled rational
  if(((2*N)==M || M == 0) && fast && N<=highestPrecompiledPower) {
#define FIXEDRATIONALENUM(x) case x: \
return std::make_unique<SwitchInterface<fixedRational<x>>>( \
  std::pair <switchType,Data> {switchType::rationalfix##x,Data::init(D0,DMAX,R0)});
    switch (N) {
      FIXEDRATIONALENUM(12)
      FIXEDRATIONALENUM(10)
      FIXEDRATIONALENUM(8)
      FIXEDRATIONALENUM(6)
      FIXEDRATIONALENUM(4)
      FIXEDRATIONALENUM(2)
    default:
      break;
    }
#undef FIXEDRATIONALENUM
  }
  //continue with the 'at runtime implementation'
  auto data = rational<rationalPow::standard,rationalForm::standard>::init(D0,DMAX,R0,N,M);
  if(2*N==M || M == 0) {
    if(fast) {
      //fast rational
      return PLMD::Tools::make_unique<SwitchInterface<rational<
             rationalPow::fast,rationalForm::simplified>>>(data);
    }
    return PLMD::Tools::make_unique<SwitchInterface<rational<
           rationalPow::standard,rationalForm::simplified>>>(data);
  }
  if(fast) {
    //fast rational
    return PLMD::Tools::make_unique<SwitchInterface<rational<
           rationalPow::fast,rationalForm::standard>>>(data);
  }
  return PLMD::Tools::make_unique<SwitchInterface<rational<
         rationalPow::standard,rationalForm::standard>>>(data);
}

std::pair <switchType,Data>
initRational(double D0,double DMAX, double R0, int N, int M) {
  bool fast = N%2==0 && M%2==0 && D0==0.0;
  //if (M==0) M will automatically became 2*NN
  constexpr int highestPrecompiledPower=12;
  //precompiled rational
  if(((2*N)==M || M == 0) && fast && N<=highestPrecompiledPower) {
#define FIXEDRATIONALENUM(x) case x: \
return std::pair <switchType,Data> {switchType::rationalfix##x,Data::init(D0,DMAX,R0)};
    switch (N) {
      FIXEDRATIONALENUM(12)
      FIXEDRATIONALENUM(10)
      FIXEDRATIONALENUM(8)
      FIXEDRATIONALENUM(6)
      FIXEDRATIONALENUM(4)
      FIXEDRATIONALENUM(2)
    default:
      break;
    }
#undef FIXEDRATIONALENUM
  }
  //continue with the 'at runtime implementation'
  return rational<rationalPow::standard,rationalForm::standard>::init(D0,DMAX,R0,N,M);
}

struct exponentialSwitch: public baseSwitch<exponentialSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0) {
    return{switchType::exponential,
           Data::init(D0,DMAX,R0)};
  }
  static inline ValueDerivative function(const double rdist) {
    double result = std::exp(-rdist);
    return {result,-result};
  }
};

struct gaussianSwitch: public baseSwitch<gaussianSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0) {
    return{switchType::gaussian,
           Data::init(D0,DMAX,R0)};
  }
  static inline ValueDerivative function(const double rdist) {
    double result = std::exp(-0.5*rdist*rdist);
    return {result,-rdist*result};
  }
};

struct fastgaussianSwitch: public baseSwitch<fastgaussianSwitch> {
  static std::pair <switchType,Data> init(
    const double DMAX) {
    return{switchType::fastgaussian,
           Data::init(0.0,DMAX,1.0)};
  }
  static inline ValueDerivative function(const double rdist) {
    double result = std::exp(-0.5*rdist*rdist);
    return {result,-rdist*result};
  }
  static inline ValueDerivative calculateSqr(const Data& data,const double distance2) {
    double result = 0.0;
    double dfunc = 0.0;
    if(distance2<data.dmax_2) {
      result=1.0;
      if(distance2 >0.0) {
        result = exp(-0.5*distance2);
        dfunc = -result;
        //I have to multiply and then divide dfunc by sqrt(distance2)
        //by omitting that I get a funny wrong value of the derivative in 0
        //by not omitting that I get a funny nan becasue I divide by 0
        //hence the extra if
        // stretch:
        result=result*data.stretch+data.shift;
        dfunc*=data.stretch;
      }
    }
    return {result,dfunc};
  }
};

struct smapSwitch: public baseSwitch<smapSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0,
    const int A,
    const int B) {
    auto data =Data::init(D0,DMAX,R0);
    data.a= A;
    data.b= B;
    data.c= std::pow(2., static_cast<double>(data.a)/static_cast<double>(data.b) ) - 1.0;
    data.d= -static_cast<double>(data.b) / static_cast<double>(data.a);
    return{switchType::smap,data};
  }

  static inline ValueDerivative function(const Data& data,const double rdist) {
    const double sx=data.c*Tools::fastpow( rdist, data.a );
    double result=std::pow( 1.0 + sx, data.d );
    double dfunc=-data.b*sx/rdist*result/(1.0+sx);
    return {result,dfunc};
  }
};

struct cubicSwitch: public baseSwitch<cubicSwitch> {
  static std::pair <switchType,Data> init(const double D0, const double DMAX) {
    return{switchType::cubic,
           Data::init(D0,DMAX,DMAX-D0)};
  }
  static inline ValueDerivative function(const double rdist) {
    const double tmp1 = rdist - 1.0;
    const double tmp2 = 1.0+2.0*rdist;
    //double result = tmp1*tmp1*tmp2;
    double dfunc = 2*tmp1*tmp2 + 2*tmp1*tmp1;
    return {tmp1*tmp1*tmp2,dfunc};
  }
};

struct tanhSwitch: public baseSwitch<tanhSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0) {
    return{switchType::tanh,
           Data::init(D0,DMAX,R0)};
  }
  static inline ValueDerivative function(const double rdist) {
    const double tmp1 = std::tanh(rdist);
    //was dfunc=-(1-tmp1*tmp1);
    double dfunc = tmp1 * tmp1 - 1.0;
    //return result;
    return {1.0 - tmp1,dfunc };
  }
};

struct cosinusSwitch: public baseSwitch<cosinusSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0) {
    return{switchType::cosinus,
           Data::init(D0,DMAX,R0)};
  }
  static inline ValueDerivative function(const Data& data,const double rdist) {
    double result = 0.0;
    double dfunc=0.0;
    if(rdist<=1.0) {
// rdist = (r-r1)/(r2-r1) ; 0.0<=rdist<=1.0 if r1 <= r <=r2; (r2-r1)/(r2-r1)=1
      double rdistPI = rdist * PLMD::pi;
      result = 0.5 * (std::cos ( rdistPI ) + 1.0);
      dfunc = -0.5 * PLMD::pi * std::sin ( rdistPI );
    }
    return {result,dfunc};
  }
};

struct nativeqSwitch: public baseSwitch<nativeqSwitch> {
  static std::pair <switchType,Data> init(
    const double D0,
    const double DMAX,
    const double R0,
    const double BETA, // nm-1
    const double LAMBDA,// unitless
    const double REF) {
    auto data=Data::init(D0,DMAX,R0);
    data.beta=BETA;
    data.lambda=LAMBDA;
    data.ref=REF;
    return{switchType::nativeq,data};
  }
  static inline ValueDerivative calculate(const Data& data,const double distance) {
    double res = 0.0;
    double dfunc = 0.0;
    if(distance<=data.dmax) {
      res = 1.0;
      if(distance > data.d0) {
        const double rdist = data.beta*(distance - data.lambda * data.ref);
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
        dfunc = - data.beta /(exprdist+ 2.0 +1.0/exprdist) /distance;

        dfunc*=data.stretch;
      }
      res=res*data.stretch+data.shift;
    }
    return {res,dfunc};
  }
};

class leptonSwitch {
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
  Data data;
public:
  leptonSwitch(double D0, double DMAX, double R0, const std::string & func)
    :lepton_func(func),
     expressions(OpenMP::getNumThreads(), lepton_func),
     data(Data::init(D0,DMAX,R0)) {
    //this is a bit odd, but it works
    auto vars=expressions[0].getVariables();
    leptonx2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();
  }

  leptonSwitch(const leptonSwitch& other)
    :lepton_func(other.lepton_func),
     expressions(OpenMP::getNumThreads(), lepton_func),
     data(other.data) {
    //this is a bit odd, but it works
    auto vars=expressions[0].getVariables();
    leptonx2=std::find(vars.begin(),vars.end(),"x2")!=vars.end();
  }

  inline double calculate(const double distance,double&dfunc) const  {
    double res = 0.0;
    dfunc = 0.0;
    if(leptonx2) {
      res= calculateSqr(distance*distance,dfunc);
    } else {
      if(distance<=data.dmax) {
        res = 1.0;
        const double rdist = (distance-data.d0)*data.invr0;
        if(rdist > 0.0) {
          const unsigned t=OpenMP::getThreadNum();
          plumed_assert(t<expressions.size());
          std::tie(res,dfunc) = expressions[t](rdist);
          dfunc *= data.invr0;
          dfunc /= distance;
        }
        res=res*data.stretch+data.shift;
        dfunc*=data.stretch;
      }
    }
    return res;
  }

  double calculateSqr(const double distance2,double&dfunc) const  {
    double result =0.0;
    dfunc=0.0;
    if(leptonx2) {
      if(distance2<=data.dmax_2) {
        const unsigned t=OpenMP::getThreadNum();
        const double rdist_2 = distance2*data.invr0_2;
        plumed_assert(t<expressions.size());
        std::tie(result,dfunc) = expressions[t](rdist_2);
        // chain rule:
        dfunc*=2*data.invr0_2;
        // stretch:
        result=result*data.stretch+data.shift;
        dfunc*=data.stretch;
      }
    } else {
      result = calculate(std::sqrt(distance2),dfunc);
    }
    return result;
  }
  const Data & getData() const {
    return data;
  }
  void setupStretch() {
    if(data.dmax!=std::numeric_limits<double>::max()) {
      data.stretch=1.0;
      data.shift=0.0;
      double dummy;
      double s0=calculate(0.0,dummy);
      double sd=calculate(data.dmax,dummy);
      data.stretch=1.0/(s0-sd);
      data.shift=-sd*data.stretch;
    }
  }
  void removeStretch() {
    data.stretch=1.0;
    data.shift=0.0;
  }
  std::string description() const {
    using namespace switchContainers;
    std::ostringstream ostr;
    ostr<<1.0/data.invr0
        <<".  Using "
        << typeToString(switchType::lepton)
        <<" switching function with parameters d0="<< data.d0
        <<" func=" << lepton_func;
    return ostr.str();
  }
};

//call to calculate with no inheritance
ValueDerivative calculate(const switchType type,
                          const Data& data,
                          const double rdist) {
#define SWITCHCALL(x) case switchType::x: return x##Switch::calculate(data, rdist);
#define RATCALL(x) case switchType::rationalfix##x:return fixedRational<x>::calculate(data, rdist);
  switch (type) {
    RATCALL(12)
    RATCALL(10)
    RATCALL(8)
    RATCALL(6)
    RATCALL(4)
    RATCALL(2)
  case switchType::rational:
    return rational<rationalPow::standard,rationalForm::standard>::calculate(data, rdist);
  case switchType::rationalFast:
    return rational<rationalPow::fast,rationalForm::standard>::calculate(data, rdist);
  case switchType::rationalSimple:
    return rational<rationalPow::standard,rationalForm::simplified>::calculate(data, rdist);
  case switchType::rationalSimpleFast:
    return rational<rationalPow::fast,rationalForm::simplified>::calculate(data, rdist);
    SWITCHCALL(exponential)
    SWITCHCALL(gaussian)
    SWITCHCALL(fastgaussian)
    SWITCHCALL(smap)
    SWITCHCALL(cubic)
    SWITCHCALL(tanh)
    SWITCHCALL(cosinus)
    SWITCHCALL(nativeq)
  default:
    break;
  }
#undef SWITCHCALL
#undef RATCALL
  return {0.0,0.0};
}

//call to calculateSqr with no inheritance
ValueDerivative calculateSqr(const switchType type,
                             const Data& data,
                             const double rdist2) {
#define SWITCHCALL(x) case switchType::x: return x##Switch::calculateSqr(data, rdist2);
#define RATCALL(x) case switchType::rationalfix##x:return fixedRational<x>::calculateSqr(data, rdist2);
  switch (type) {
    RATCALL(12)
    RATCALL(10)
    RATCALL(8)
    RATCALL(6)
    RATCALL(4)
    RATCALL(2)
  case switchType::rational:
    return rational<rationalPow::standard,rationalForm::standard>::calculateSqr(data, rdist2);
  case switchType::rationalFast:
    return rational<rationalPow::fast,rationalForm::standard>::calculateSqr(data, rdist2);
  case switchType::rationalSimple:
    return rational<rationalPow::standard,rationalForm::simplified>::calculateSqr(data, rdist2);
  case switchType::rationalSimpleFast:
    return rational<rationalPow::fast,rationalForm::simplified>::calculateSqr(data, rdist2);
    SWITCHCALL(exponential)
    SWITCHCALL(gaussian)
    SWITCHCALL(fastgaussian)
    SWITCHCALL(smap)
    SWITCHCALL(cubic)
    SWITCHCALL(tanh)
    SWITCHCALL(cosinus)
    SWITCHCALL(nativeq)
  default:
    break;
  }
#undef SWITCHCALL
#undef RATCALL
  return {0.0,0.0};
}

//call to setupStretch with no inheritance
void setupStretch(switchType type, Data& data) {
  if(data.dmax!=std::numeric_limits<double>::max()) {
    data.stretch=1.0;
    data.shift=0.0;
    double s0=calculate(type,data,0.0).first;
    double sd=calculate(type,data,data.dmax).first;
    data.stretch=1.0/(s0-sd);
    data.shift=-sd*data.stretch;
  }
}

void removeStretch(Data& data) {
  data.stretch=1.0;
  data.shift=0.0;
}

std::string typeToString(switchType type) {
#define DEFAULTPRINT(x) case switchType::x: return #x;
  switch (type) {
  case switchType::rationalfix12:
  case switchType::rationalfix10:
  case switchType::rationalfix8:
  case switchType::rationalfix6:
  case switchType::rationalfix4:
  case switchType::rationalfix2:
  case switchType::rationalFast:
  case switchType::rationalSimple:
  case switchType::rationalSimpleFast:
    DEFAULTPRINT(rational)
    DEFAULTPRINT(exponential)
    DEFAULTPRINT(gaussian)
    DEFAULTPRINT(fastgaussian)
    DEFAULTPRINT(smap)
    DEFAULTPRINT(cubic)
    DEFAULTPRINT(tanh)
    DEFAULTPRINT(cosinus)
    DEFAULTPRINT(nativeq)
    DEFAULTPRINT(lepton)
  case switchType::not_initialized:
    return "not initialized!";
  }
#undef DEFAULTPRINT
  return "";
}

std::string description(switchType type, const Data& data) {
  std::ostringstream ostr;
  ostr<<1.0/data.invr0
      <<".  Using "
      << typeToString(type)
      <<" switching function with parameters d0="<< data.d0;
  switch (type) {
#define RATFIX(N) case switchType::rationalfix##N: \
    ostr<< " " << " nn=" << N << " mm=" <<N*2; \
    break;
    RATFIX(12)
    RATFIX(10)
    RATFIX(8)
    RATFIX(6)
    RATFIX(4)
    RATFIX(2)
#undef RATFIX
  case switchType::rational:
    ostr << " nn=" << data.nn << " mm=" << data.mm;
    break;
  case switchType::smap:
    ostr<<" a="<<data.a<<" b="<<data.b;
    break;
  case switchType::cubic:
    ostr<<" dmax="<<data.dmax;
    break;
  case switchType::nativeq:
    ostr<<" beta="<<data.beta<<" lambda="<<data.lambda<<" ref="<<data.ref;
    break;
  default:

    break;
  }
  return ostr.str();
}

class SwitchInterface_lepton: public Switch {
  leptonSwitch function;
public:
  SwitchInterface_lepton(double D0, double DMAX, double R0, const std::string & func):
    function(D0,DMAX,R0,func) {}
  SwitchInterface_lepton(const SwitchInterface_lepton& other):
    function(other.function) {}
  double calculate(double distance, double& dfunc) const override {
    return function.calculate(distance,dfunc);
  }
  double calculateSqr(double distance2, double& dfunc) const override {
    return function.calculateSqr(distance2,dfunc);
  }
  void setupStretch() override {
    function.setupStretch();
  }
  const Data& getData() const override {
    return function.getData();
  }
  switchType getType() const override {
    return switchType::lepton;
  }
  std::string description() const override {
    return function.description();
  }
};

} // namespace switchContainers

SwitchingFunction::SwitchingFunction()=default;
SwitchingFunction::~SwitchingFunction()=default;

SwitchingFunction::SwitchingFunction(const SwitchingFunction& other) {
  copyFunction(other);
}

SwitchingFunction::SwitchingFunction(SwitchingFunction&& other):
  function(std::move(other.function)) {}

SwitchingFunction& SwitchingFunction::operator=(const SwitchingFunction& other) {
  if (this != &other) {
    copyFunction(other);
  }
  return *this;
}

SwitchingFunction& SwitchingFunction::operator=(SwitchingFunction&& other) {
  function.reset();
  function=std::move(other.function);
  return *this;
}

void SwitchingFunction::copyFunction(const SwitchingFunction& other) {
  if (this == &other) {
    return; // nothing to do
  }
  function.reset();
  if (!other.function) {
    //now both have the function uninitialized
    return; // nothing to copy
  }
  using namespace switchContainers;
  const auto settings = std::make_pair(other.function->getType(),other.function->getData());
#define SWITCHCALL(x) case switchType::x: \
  function = std::make_unique<SwitchInterface<x##Switch>>(settings); \
  break;
#define RATCALL(x) case switchType::rationalfix##x:\
  function = std::make_unique<SwitchInterface<fixedRational<x>>>(settings); \
  break;
  switch (settings.first) {
    RATCALL(12)
    RATCALL(10)
    RATCALL(8)
    RATCALL(6)
    RATCALL(4)
    RATCALL(2)
  case switchType::rational:
    function = std::make_unique<SwitchInterface<rational<rationalPow::standard,rationalForm::standard>>>(settings);
    break;
  case switchType::rationalFast:
    function = std::make_unique<SwitchInterface<rational<rationalPow::fast,rationalForm::standard>>>(settings);
    break;
  case switchType::rationalSimple:
    function = std::make_unique<SwitchInterface<rational<rationalPow::standard,rationalForm::simplified>>>(settings);
    break;
  case switchType::rationalSimpleFast:
    function = std::make_unique<SwitchInterface<rational<rationalPow::fast,rationalForm::simplified>>>(settings);
    break;
    SWITCHCALL(exponential)
    SWITCHCALL(gaussian)
    SWITCHCALL(fastgaussian)
    SWITCHCALL(smap)
    SWITCHCALL(cubic)
    SWITCHCALL(tanh)
    SWITCHCALL(cosinus)
    SWITCHCALL(nativeq)
  case switchType::lepton:
    function = std::make_unique<SwitchInterface_lepton>(*dynamic_cast<SwitchInterface_lepton*>(other.function.get()));
  case switchType::not_initialized:
    break;
  }
#undef SWITCHCALL
#undef RATCALL
}


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
  function.reset();
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
  double d0=0.0;
  double dmax=std::numeric_limits<double>::max();
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
  using namespace switchContainers;
  if(name=="CUBIC") {
    //cubic is the only switch type that only uses d0 and dmax
    function= std::make_unique<SwitchInterface<cubicSwitch>>(
                cubicSwitch::init(d0,dmax));
  } else {
    double r0=0.0;
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
      function= std::make_unique<SwitchInterface<smapSwitch>>(
                  smapSwitch::init(d0,dmax,r0,a,b));
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
      function= std::make_unique<SwitchInterface<nativeqSwitch>>(
                  nativeqSwitch::init(d0,dmax,r0,beta,lambda,ref));
    } else if(name=="EXP") {
      function= std::make_unique<SwitchInterface<exponentialSwitch>>(
                  exponentialSwitch::init(d0,dmax,r0));
    } else if(name=="GAUSSIAN") {
      if ( r0==1.0 && d0==0.0 ) {
        function= std::make_unique<SwitchInterface<fastgaussianSwitch>>(
                    fastgaussianSwitch::init(dmax));
      } else {
        function= std::make_unique<SwitchInterface<gaussianSwitch>>(
                    gaussianSwitch::init(d0,dmax,r0));
      }
    } else if(name=="TANH") {
      function= std::make_unique<SwitchInterface<tanhSwitch>>(
                  tanhSwitch::init(d0,dmax,r0));
    } else if(name=="COSINUS") {
      function= std::make_unique<SwitchInterface<cosinusSwitch>>(
                  cosinusSwitch::init(d0,dmax,r0));
    } else if((name=="MATHEVAL" || name=="CUSTOM")) {
      std::string func;
      Tools::parse(data,"FUNC",func);
      function= std::make_unique<SwitchInterface_lepton>(d0,dmax,r0,func);
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
  return function->calculateSqr( distance2, dfunc);
}

double SwitchingFunction::calculate(double distance,double&dfunc)const {
  plumed_massert(function,"you are trying to use an unset SwitchingFunction");
  return function->calculate( distance, dfunc);
}

void SwitchingFunction::set(const int nn,int mm, const double r0, const double d0) {
  function.reset();
  if(mm == 0) {
    mm = 2*nn;
  }
  double dmax=d0+r0*std::pow(0.00001,1./(nn-mm));
  function = switchContainers::rationalFactory(d0,dmax,r0,nn,mm);
  function->setupStretch();
}

double SwitchingFunction::get_r0() const {
  return 1.0/function->getData().invr0;
}

double SwitchingFunction::get_d0() const {
  return function->getData().d0;
}

double SwitchingFunction::get_dmax() const {
  return function->getData().dmax;
}

double SwitchingFunction::get_dmax2() const {
  return function->getData().dmax_2;
}

void SwitchingFunctionAccelerable::registerKeywords( Keywords& keys ) {
  SwitchingFunction::registerKeywords(keys);
}

void SwitchingFunctionAccelerable::set(const std::string & definition,std::string& errormsg) {
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
  double d0=0.0;
  double dmax=std::numeric_limits<double>::max();
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
  using namespace switchContainers;
  if(name=="CUBIC") {
    //cubic is the only switch type that only uses d0 and dmax
    std::tie(type,switchData) = cubicSwitch::init(d0,dmax);
  } else {
    double r0=0.0;
    REQUIREDPARSE(data,"R_0",r0,errormsg);
    if(name=="RATIONAL") {
      int nn=6;
      int mm=0;
      CHECKandPARSE(data,"NN",nn,errormsg);
      CHECKandPARSE(data,"MM",mm,errormsg);
      std::tie(type,switchData) = switchContainers::initRational(d0,dmax,r0,nn,mm);

    } else if(name=="SMAP") {
      int a=0;
      int b=0;
      //in the original a and b are "default=0",
      //but you divide by a and b during the initialization!
      //better an error message than an UB, so no default
      REQUIREDPARSE(data,"A",a,errormsg);
      REQUIREDPARSE(data,"B",b,errormsg);
      std::tie(type,switchData) = smapSwitch::init(d0,dmax,r0,a,b);
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
      std::tie(type,switchData) = nativeqSwitch::init(d0,dmax,r0,beta,lambda,ref);
    } else if(name=="EXP") {
      std::tie(type,switchData) = exponentialSwitch::init(d0,dmax,r0);
    } else if(name=="GAUSSIAN") {
      if ( r0==1.0 && d0==0.0 ) {
        std::tie(type,switchData) = fastgaussianSwitch::init(dmax);
      } else {
        std::tie(type,switchData) = gaussianSwitch::init(d0,dmax,r0);
      }
    } else if(name=="TANH") {
      std::tie(type,switchData) = tanhSwitch::init(d0,dmax,r0);
    } else if(name=="COSINUS") {
      std::tie(type,switchData) = cosinusSwitch::init(d0,dmax,r0);
    } else if((name=="MATHEVAL" || name=="CUSTOM")) {
      std::string func;
      Tools::parse(data,"FUNC",func);
      errormsg="the custom switching function is not comatible with acceleration";
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
    setupStretch(type,switchData);
  }
}

std::string SwitchingFunctionAccelerable::description() const {
  // if this error is necessary, something went wrong in the constructor
  //  plumed_merror("Unknown switching function type");
  return switchContainers::description(type,switchData);
}

double SwitchingFunctionAccelerable::calculateSqr(double distance2,double&dfunc)const {
  double result;
  std::tie(result, dfunc) = switchContainers::calculateSqr(type,switchData,distance2);
  return result;
}

double SwitchingFunctionAccelerable::calculate(double distance,double&dfunc)const {
  //massert do not go with openacc
  // plumed_massert(init,"you are trying to use an unset SwitchingFunction");
  double result;
  std::tie(result, dfunc) = switchContainers::calculate(type,switchData,distance);
  return result;
}

void SwitchingFunctionAccelerable::set(const int nn,int mm, const double r0, const double d0) {
  if(mm == 0) {
    mm = 2*nn;
  }
  double dmax=d0+r0*std::pow(0.00001,1./(nn-mm));
  std::tie(type,switchData) = switchContainers::initRational(d0,dmax,r0,nn,mm);
  setupStretch(type,switchData);
}

double SwitchingFunctionAccelerable::get_r0() const {
  return 1.0/switchData.invr0;
}

double SwitchingFunctionAccelerable::get_d0() const {
  return switchData.d0;
}

double SwitchingFunctionAccelerable::get_dmax() const {
  return switchData.dmax;
}

double SwitchingFunctionAccelerable::get_dmax2() const {
  return switchData.dmax_2;
}

void SwitchingFunctionAccelerable::toACCDevice() const {
#pragma acc enter data copyin(this[0:1],type)
  switchData.toACCDevice();
}

void SwitchingFunctionAccelerable::removeFromACCDevice() const {
  switchData.removeFromACCDevice();
#pragma acc exit data delete(type,this[0:1])
}
}// Namespace PLMD
