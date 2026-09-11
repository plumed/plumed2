/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_function_Combine_h
#define __PLUMED_function_Combine_h
#include "FunctionSetup.h"
#include "tools/Keywords.h"
#include "tools/View.h"
#include "core/ActionWithArguments.h"

#include <numeric>
namespace PLMD {
namespace function {

class Combine  {
  struct component {
    double coefficient{0.0};
    double parameter{0.0};
    double power{0.0};
    double max_minus_min{0.0};
    double inv_max_minus_min{0.0};
    bool periodic{false};
  };
  std::vector<component> components{};
public:
  std::size_t ncomponents{0};
  component * cmps=nullptr;
  void update() {
    cmps = components.data();
    ncomponents = components.size();
  }
  Combine() = default;
  ~Combine() = default;
  Combine(const Combine&x):
    components(x.components) {
    update();
  }
  Combine(Combine&&x):
    components(std::move(x.components)) {
    update();
  }
  Combine &operator=(const Combine&x) {
    if (this!=&x) {
      components=x.components;
      update();
    }
    return *this;
  }
  Combine &operator=(Combine&&x) {
    if (this!=&x) {
      components=std::move(x.components);
      update();
    }
    return *this;
  }
  static void registerKeywords(Keywords& keys);
  static void read( Combine& func,
                    ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const Combine& func,
                    bool noderiv,
                    const View<const double>& args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], ncomponents, cmps[0:ncomponents])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(cmps[0:ncomponents], ncomponents, this[0:1])
  }
#endif // __PLUMED_HAS_OPENACC
};

void Combine::registerKeywords(Keywords& keys) {
  keys.use("PERIODIC");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","PARAMETERS","0.0","the parameters of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficients so that in total they are equal to one");
  keys.setValueDescription("scalar/vector/matrix","a linear combination");
}

void Combine::read( Combine& func, ActionWithArguments* action,
                    FunctionOptions& options ) {
  unsigned nargs = action->getNumberOfArguments();
  ActionWithVector* av=dynamic_cast<ActionWithVector*>(action);
  if(av && av->getNumberOfMasks()>0) {
    nargs = nargs - av->getNumberOfMasks();
  }
  std::vector<double> coefficients(nargs);
  action->parseVector("COEFFICIENTS",coefficients);

  if(coefficients.size()!=nargs) {
    action->error("Size of COEFFICIENTS array should be the same as number for arguments");
  }
  std::vector<double> parameters(nargs);
  action->parseVector("PARAMETERS",parameters);
  if(parameters.size()!=nargs) {
    action->error("Size of PARAMETERS array should be the same as number for arguments");
  }
  std::vector<double> powers(nargs);
  action->parseVector("POWERS",powers);
  if(powers.size()!=nargs) {
    action->error("Size of POWERS array should be the same as number for arguments");
  }

  bool normalize;
  action->parseFlag("NORMALIZE",normalize);
  if(normalize) {
    double n=std::accumulate(coefficients.begin(),coefficients.end(),0.0);
    std::transform(coefficients.begin(),coefficients.end(),coefficients.begin(),
    [n](double x) {
      return x/n;
    });
  }

  action->log.printf("  with coefficients:");
  func.components.resize( nargs );
  func.update();
  for(unsigned i=0; i<func.components.size(); i++) {
    func.components[i].coefficient = coefficients[i];
    func.components[i].parameter = parameters[i];
    func.components[i].power = powers[i];
    action->log.printf(" %f",func.components[i].coefficient);
  }
  action->log.printf("\n  with parameters:");
  for(unsigned i=0; i<func.components.size(); i++) {
    action->log.printf(" %f",func.components[i].parameter);
  }
  action->log.printf("\n  and powers:");
  for(unsigned i=0; i<func.components.size(); i++) {
    action->log.printf(" %f",func.components[i].power);
  }
  action->log.printf("\n");
  // Store periodicity stuff
  for(unsigned i=0; i<nargs; ++i) {
    if( (action->getPntrToArgument(i))->isPeriodic() ) {
      func.components[i].periodic = true;
      std::string min, max;
      (action->getPntrToArgument(i))->getDomain( min, max );
      double dmin, dmax;
      Tools::convert( min, dmin );
      Tools::convert( max, dmax );
      func.components[i].max_minus_min = dmax - dmin;
      func.components[i].inv_max_minus_min = 1.0 /  func.components[i].max_minus_min ;
    }
  }

}

void Combine::calc( const Combine& func,
                    bool noderiv,
                    const View<const double>& args,
                    FunctionOutput& funcout ) {
  funcout.values[0]=0.0;
  for(unsigned i=0; i<func.ncomponents; ++i) {
    const auto& cmp = func.cmps[i];
    double cv = args[i] - cmp.parameter;
    if( cmp.periodic ) {
      cv = cmp.max_minus_min*Tools::pbc( cv*cmp.inv_max_minus_min );
    }
    funcout.values[0] += cmp.coefficient*pow( cv, cmp.power );
    if( !noderiv ) {
      funcout.derivs[0][i] = cmp.coefficient*cmp.power*pow(cv,cmp.power-1.0);
    }
  }
}

} // namespace function
} // namespace PLMD
#endif //__PLUMED_function_Combine_h
