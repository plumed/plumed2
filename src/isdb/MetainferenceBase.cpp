/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2026 The plumed team
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
#include "MetainferenceBase.h"

namespace PLMD {
namespace isdb {

void MetainferenceBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("optional","ARG","scalar","the labels of the values from which the function is calculated");
  keys.addFlag("DOSCORE",false,"activate metainference");
  MetainferenceCore::registerKeywords(keys);
  keys.addOutputComponent("score",        "default",      "scalar","the Metainference score");
}

MetainferenceBase::MetainferenceBase(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  doscore_(false),
  narg(0),
  valueScore(nullptr),
  iselect(0) {
  parseFlag("DOSCORE", doscore_);
  core_=std::make_unique<MetainferenceCore>(*this,doscore_);
  // reweight implies a different number of arguments (the latest one must always be the bias)
  if(core_->doReweight()&&getNumberOfArguments()!=1) {
    error("To REWEIGHT one must provide one single bias as an argument");
  }
}

MetainferenceBase::~MetainferenceBase() {
}

void MetainferenceBase::Initialise(const unsigned input) {
  setNarg(input);
  if(narg!=parameters.size()) {
    std::string num1;
    Tools::convert(parameters.size(),num1);
    std::string num2;
    Tools::convert(narg,num2);
    std::string msg = "The number of experimental values " + num1 +" must be the same of the calculated values " + num2;
    error(msg);
  }

  addComponentWithDerivatives("score");
  componentIsNotPeriodic("score");
  valueScore=getPntrToComponent("score");

  core_->setup(parameters, *this, [this](const std::string& name) {
    return getPntrToComponent(name);
  });
}

void MetainferenceBase::Selector() {
  iselect = core_->getSelector();
}

double MetainferenceBase::getScore() {
  const double ene = core_->calculate(iselect, core_->doReweight() ? getArgument(0) : 0.);
  if(core_->doReweight()) {
    setArgDerivatives(valueScore, core_->getBiasDer());
  }
  return ene;
}

void MetainferenceBase::writeStatus() {
  core_->writeStatus();
}

void MetainferenceBase::update() {
  // write the status file
  core_->update();
}

}
}
