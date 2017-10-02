/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "BasisFunctions.h"
#include "TargetDistribution.h"
#include "VesBias.h"
#include "VesTools.h"
#include "GridIntegrationWeights.h"

#include "tools/Grid.h"
#include "tools/Tools.h"


namespace PLMD {
namespace ves {

void BasisFunctions::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  keys.add("compulsory","ORDER","The order of the basis function expansion.");
  keys.add("compulsory","MINIMUM","The minimum of the interval on which the basis functions are defined.");
  keys.add("compulsory","MAXIMUM","The maximum of the interval on which the basis functions are defined.");
  keys.add("hidden","NGRID_POINTS","The number of grid points used for numerical integrals");
  keys.addFlag("DEBUG_INFO",false,"Print out more detailed information about the basis set. Useful for debugging.");
  keys.addFlag("NUMERICAL_INTEGRALS",false,"Calculate basis function integral for the uniform distribution numerically. Useful for debugging.");
}


BasisFunctions::BasisFunctions(const ActionOptions&ao):
  Action(ao),
  print_debug_info_(false),
  has_been_set(false),
  description_("Undefined"),
  type_("Undefined"),
  norder_(0),
  nbasis_(1),
  bf_label_prefix_("f"),
  bf_labels_(nbasis_,"f0"),
  periodic_(false),
  interval_bounded_(true),
  interval_intrinsic_min_str_("1.0"),
  interval_intrinsic_max_str_("-1.0"),
  interval_intrinsic_min_(1.0),
  interval_intrinsic_max_(-1.0),
  interval_intrinsic_range_(0.0),
  interval_intrinsic_mean_(0.0),
  interval_min_str_(""),
  interval_max_str_(""),
  interval_min_(0.0),
  interval_max_(0.0),
  interval_range_(0.0),
  interval_mean_(0.0),
  argT_derivf_(1.0),
  numerical_uniform_integrals_(false),
  nbins_(1001),
  uniform_integrals_(nbasis_,0.0),
  vesbias_pntr_(NULL),
  action_pntr_(NULL)
{
  bf_keywords_.push_back(getName());
  if(keywords.exists("ORDER")) {
    parse("ORDER",norder_); addKeywordToList("ORDER",norder_);
  }
  nbasis_=norder_+1;
  //
  std::string str_imin; std::string str_imax;
  if(keywords.exists("MINIMUM") && keywords.exists("MAXIMUM")) {
    parse("MINIMUM",str_imin); addKeywordToList("MINIMUM",str_imin);
    parse("MAXIMUM",str_imax); addKeywordToList("MAXIMUM",str_imax);
  }
  else {
    str_imin = "-1.0";
    str_imax = "1.0";
  }
  interval_min_str_ = str_imin;
  interval_max_str_ = str_imax;
  if(!Tools::convert(str_imin,interval_min_)) {
    plumed_merror(getName()+": cannot convert the value given in MINIMUM to a double");
  }
  if(!Tools::convert(str_imax,interval_max_)) {
    plumed_merror(getName()+": cannot convert the value given in MAXIMUM to a double");
  }
  if(interval_min_>interval_max_) {plumed_merror(getName()+": MINIMUM and MAXIMUM are not correctly defined");}
  //
  parseFlag("DEBUG_INFO",print_debug_info_);
  if(keywords.exists("NUMERICAL_INTEGRALS")) {
    parseFlag("NUMERICAL_INTEGRALS",numerical_uniform_integrals_);
  }
  if(keywords.exists("NGRID_POINTS")) {
    parse("NGRID_POINTS",nbins_);
  }
  // log.printf(" %s \n",getKeywordString().c_str());

}


void BasisFunctions::setIntrinsicInterval(const double interval_intrinsic_min_in, const double interval_intrinsic_max_in) {
  interval_intrinsic_min_ = interval_intrinsic_min_in;
  interval_intrinsic_max_ = interval_intrinsic_max_in;
  VesTools::convertDbl2Str(interval_intrinsic_min_,interval_intrinsic_min_str_);
  VesTools::convertDbl2Str(interval_intrinsic_max_,interval_intrinsic_max_str_);
  plumed_massert(interval_intrinsic_min_<interval_intrinsic_max_,"setIntrinsicInterval: intrinsic intervals are not defined correctly");
}


void BasisFunctions::setIntrinsicInterval(const std::string& interval_intrinsic_min_str_in, const std::string& interval_intrinsic_max_str_in) {
  interval_intrinsic_min_str_ = interval_intrinsic_min_str_in;
  interval_intrinsic_max_str_ = interval_intrinsic_max_str_in;
  if(!Tools::convert(interval_intrinsic_min_str_,interval_intrinsic_min_)) {
    plumed_merror("setIntrinsicInterval: cannot convert string value given for the minimum of the intrinsic interval to a double");
  }
  if(!Tools::convert(interval_intrinsic_max_str_,interval_intrinsic_max_)) {
    plumed_merror("setIntrinsicInterval: cannot convert string value given for the maximum of the intrinsic interval to a double");
  }
  plumed_massert(interval_intrinsic_min_<interval_intrinsic_max_,"setIntrinsicInterval: intrinsic intervals are not defined correctly");
}


void BasisFunctions::setInterval(const double interval_min_in, const double interval_max_in) {
  interval_min_ = interval_min_in;
  interval_max_ = interval_max_in;
  VesTools::convertDbl2Str(interval_min_,interval_min_str_);
  VesTools::convertDbl2Str(interval_max_,interval_max_str_);
  plumed_massert(interval_min_<interval_max_,"setInterval: intervals are not defined correctly");
}


void BasisFunctions::setInterval(const std::string& interval_min_str_in, const std::string& interval_max_str_in) {
  interval_min_str_ = interval_min_str_in;
  interval_max_str_ = interval_max_str_in;
  if(!Tools::convert(interval_min_str_,interval_min_)) {
    plumed_merror("setInterval: cannot convert string value given for the minimum of the interval to a double");
  }
  if(!Tools::convert(interval_max_str_,interval_max_)) {
    plumed_merror("setInterval: cannot convert string value given for the maximum of the interval to a double");
  }
  plumed_massert(interval_min_<interval_max_,"setInterval: intervals are not defined correctly");
}


void BasisFunctions::setupInterval() {
  // if(!intervalBounded()){plumed_merror("setupInterval() only works for bounded interval");}
  interval_intrinsic_range_ = interval_intrinsic_max_-interval_intrinsic_min_;
  interval_intrinsic_mean_  = 0.5*(interval_intrinsic_max_+interval_intrinsic_min_);
  interval_range_ = interval_max_-interval_min_;
  interval_mean_  = 0.5*(interval_max_+interval_min_);
  argT_derivf_ = interval_intrinsic_range_/interval_range_;
}


void BasisFunctions::setupLabels() {
  for(unsigned int i=0; i < nbasis_; i++) {
    std::string is; Tools::convert(i,is);
    bf_labels_[i]=bf_label_prefix_+is+"(s)";
  }
}


void BasisFunctions::setupUniformIntegrals() {
  numerical_uniform_integrals_=true;
  numericalUniformIntegrals();
}


void BasisFunctions::setupBF() {
  if(interval_intrinsic_min_>interval_intrinsic_max_) {plumed_merror("setupBF: default intervals are not correctly set");}
  setupInterval();
  setupLabels();
  if(bf_labels_.size()==1) {plumed_merror("setupBF: the labels of the basis functions are not correct.");}
  if(!numerical_uniform_integrals_) {setupUniformIntegrals();}
  else {numericalUniformIntegrals();}
  if(uniform_integrals_.size()==1) {plumed_merror("setupBF: the integrals of the basis functions is not correct.");}
  if(type_=="Undefined") {plumed_merror("setupBF: the type of the basis function is not defined.");}
  if(description_=="Undefined") {plumed_merror("setupBF: the description of the basis function is not defined.");}
  has_been_set=true;
  printInfo();
}


void BasisFunctions::printInfo() const {
  if(!has_been_set) {plumed_merror("the basis set has not be setup correctly");}
  log.printf("  One-dimensional basis set\n");
  log.printf("   Description: %s\n",description_.c_str());
  log.printf("   Type: %s\n",type_.c_str());
  if(periodic_) {log.printf("   The basis functions are periodic\n");}
  log.printf("   Order of basis set: %u\n",norder_);
  log.printf("   Number of basis functions: %u\n",nbasis_);
  // log.printf("   Interval of basis set: %f to %f\n",interval_min_,interval_max_);
  log.printf("   Interval of basis set: %s to %s\n",interval_min_str_.c_str(),interval_max_str_.c_str());
  log.printf("   Description of basis functions:\n");
  for(unsigned int i=0; i < nbasis_; i++) {log.printf("    %2u       %10s\n",i,bf_labels_[i].c_str());}
  //
  if(print_debug_info_) {
    log.printf("  Debug information:\n");
    // log.printf("   Default interval of basis set: [%f,%f]\n",interval_intrinsic_min_,interval_intrinsic_max_);
    log.printf("   Intrinsic interval of basis set: [%s,%s]\n",interval_intrinsic_min_str_.c_str(),interval_intrinsic_max_str_.c_str());
    log.printf("   Intrinsic interval of basis set: range=%f,  mean=%f\n",interval_intrinsic_range_,interval_intrinsic_mean_);
    // log.printf("   Defined interval of basis set: [%f,%f]\n",interval_min_,interval_max_);
    log.printf("   Defined interval of basis set: [%s,%s]\n",interval_min_str_.c_str(),interval_max_str_.c_str());
    log.printf("   Defined interval of basis set: range=%f,  mean=%f\n",interval_range_,interval_mean_);
    log.printf("   Derivative factor due to interval translation: %f\n",argT_derivf_);
    log.printf("   Integral of basis functions over the interval:\n");
    if(numerical_uniform_integrals_) {log.printf("   Note: calculated numerically\n");}
    for(unsigned int i=0; i < nbasis_; i++) {log.printf("    %2u       %16.10f\n",i,uniform_integrals_[i]);}
    log.printf("   --------------------------\n");
  }
}


void BasisFunctions::linkVesBias(VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);
}


void BasisFunctions::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


void BasisFunctions::numericalUniformIntegrals() {
  std::vector<std::string> grid_min(1); grid_min[0]=intervalMinStr();
  std::vector<std::string> grid_max(1); grid_max[0]=intervalMaxStr();
  std::vector<unsigned int> grid_bins(1); grid_bins[0]=nbins_;
  std::vector<Value*> arguments(1); arguments[0]= new Value(NULL,"arg",false);
  if(arePeriodic()) {arguments[0]->setDomain(intervalMinStr(),intervalMaxStr());}
  else {arguments[0]->setNotPeriodic();}
  Grid* uniform_grid = new Grid("uniform",arguments,grid_min,grid_max,grid_bins,false,false);
  //
  double inverse_normalization = 1.0/(intervalMax()-intervalMin());
  for(Grid::index_t l=0; l<uniform_grid->getSize(); l++) {
    uniform_grid->setValue(l,inverse_normalization);
  }
  uniform_integrals_ = numericalTargetDistributionIntegralsFromGrid(uniform_grid);
  delete arguments[0]; arguments.clear();
  delete uniform_grid;
}


std::vector<double> BasisFunctions::numericalTargetDistributionIntegralsFromGrid(const Grid* grid_pntr) const {
  plumed_massert(grid_pntr!=NULL,"the grid is not defined");
  plumed_massert(grid_pntr->getDimension()==1,"the target distribution grid should be one dimensional");
  //
  std::vector<double> targetdist_integrals(nbasis_,0.0);
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(grid_pntr);

  for(Grid::index_t k=0; k < grid_pntr->getSize(); k++) {
    double arg = grid_pntr->getPoint(k)[0];
    std::vector<double> bf_values(nbasis_);
    std::vector<double> bf_derivs(nbasis_);
    bool inside=true;
    double argT=0.0;
    getAllValues(arg,argT,inside,bf_values,bf_derivs);
    for(unsigned int i=0; i < nbasis_; i++) {
      targetdist_integrals[i] += (integration_weights[k] * grid_pntr->getValue(k)) * bf_values[i];
    }
  }
  // assume that the first function is the constant
  bool inside=true;
  double argT=0.0;
  targetdist_integrals[0] = getValue(0.0,0,argT,inside);
  return targetdist_integrals;
}


std::vector<double> BasisFunctions::getTargetDistributionIntegrals(const TargetDistribution* targetdist_pntr) const {
  if(targetdist_pntr==NULL) {
    return getUniformIntegrals();
  }
  else {
    Grid* targetdist_grid = targetdist_pntr->getTargetDistGridPntr();
    return numericalTargetDistributionIntegralsFromGrid(targetdist_grid);
  }
}


std::string BasisFunctions::getKeywordString() const {
  std::string str_keywords=bf_keywords_[0];
  for(unsigned int i=1; i<bf_keywords_.size(); i++) {str_keywords+=" "+bf_keywords_[i];}
  return str_keywords;
}


double BasisFunctions::getValue(const double arg, const unsigned int n, double& argT, bool& inside_range) const {
  plumed_massert(n<numberOfBasisFunctions(),"getValue: n is outside range of the defined order of the basis set");
  inside_range=true;
  std::vector<double> tmp_values(numberOfBasisFunctions());
  std::vector<double> tmp_derivs(numberOfBasisFunctions());
  getAllValues(arg, argT, inside_range, tmp_values, tmp_derivs);
  return tmp_values[n];
}


void BasisFunctions::getAllValuesNumericalDerivs(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // use forward difference, unless very close to the boundary
  double delta = sqrt(epsilon);
  if((arg+delta)>intervalMax()) {
    delta *= -1.0;
  }
  inside_range=true;
  std::vector<double> values_delta(numberOfBasisFunctions());
  std::vector<double> derivs_dummy(numberOfBasisFunctions());
  getAllValues(arg+delta, argT, inside_range, values_delta, derivs_dummy);
  getAllValues(arg, argT, inside_range, values, derivs_dummy);
  for(unsigned int i=0; i<numberOfBasisFunctions(); i++) {
    derivs[i] = (values_delta[i]-values[i])/delta;
  }
}


void BasisFunctions::getMultipleValue(const std::vector<double>& args, std::vector<double>& argsT, std::vector<std::vector<double> >& values, std::vector<std::vector<double> >& derivs, const bool numerical_deriv) const {
  argsT.resize(args.size());
  values.clear();
  derivs.clear();
  for(unsigned int i=0; i<args.size(); i++) {
    std::vector<double> tmp_values(getNumberOfBasisFunctions());
    std::vector<double> tmp_derivs(getNumberOfBasisFunctions());
    bool inside_interval=true;
    if(!numerical_deriv) {
      getAllValues(args[i],argsT[i],inside_interval,tmp_values,tmp_derivs);
    } else {
      getAllValuesNumericalDerivs(args[i],argsT[i],inside_interval,tmp_values,tmp_derivs);
    }
    values.push_back(tmp_values);
    derivs.push_back(tmp_derivs);
  }
}


void BasisFunctions::writeBasisFunctionsToFile(OFile& ofile_values, OFile& ofile_derivs, const std::string& min_in, const std::string& max_in, unsigned int nbins_in, const bool ignore_periodicity, const std::string& output_fmt, const bool numerical_deriv) const {
  std::vector<std::string> min(1); min[0]=min_in;
  std::vector<std::string> max(1); max[0]=max_in;
  std::vector<unsigned int> nbins(1); nbins[0]=nbins_in;
  std::vector<Value*> value_pntr(1);
  value_pntr[0]= new Value(NULL,"arg",false);
  if(arePeriodic() && !ignore_periodicity) {value_pntr[0]->setDomain(intervalMinStr(),intervalMaxStr());}
  else {value_pntr[0]->setNotPeriodic();}
  Grid args_grid = Grid("grid",value_pntr,min,max,nbins,false,false);

  std::vector<double> args(args_grid.getSize(),0.0);
  for(unsigned int i=0; i<args.size(); i++) {
    args[i] = args_grid.getPoint(i)[0];
  }
  std::vector<double> argsT;
  std::vector<std::vector<double> > values;
  std::vector<std::vector<double> > derivs;

  ofile_values.addConstantField("bf_keywords").printField("bf_keywords","{"+getKeywordString()+"}");
  ofile_derivs.addConstantField("bf_keywords").printField("bf_keywords","{"+getKeywordString()+"}");

  ofile_values.addConstantField("min").printField("min",intervalMinStr());
  ofile_values.addConstantField("max").printField("max",intervalMaxStr());

  ofile_derivs.addConstantField("min").printField("min",intervalMinStr());
  ofile_derivs.addConstantField("max").printField("max",intervalMaxStr());

  ofile_values.addConstantField("nbins").printField("nbins",static_cast<int>(args_grid.getNbin()[0]));
  ofile_derivs.addConstantField("nbins").printField("nbins",static_cast<int>(args_grid.getNbin()[0]));

  if(arePeriodic()) {
    ofile_values.addConstantField("periodic").printField("periodic","true");
    ofile_derivs.addConstantField("periodic").printField("periodic","true");
  }
  else {
    ofile_values.addConstantField("periodic").printField("periodic","false");
    ofile_derivs.addConstantField("periodic").printField("periodic","false");
  }

  getMultipleValue(args,argsT,values,derivs,numerical_deriv);
  ofile_values.fmtField(output_fmt);
  ofile_derivs.fmtField(output_fmt);
  for(unsigned int i=0; i<args.size(); i++) {
    ofile_values.printField("arg",args[i]);
    ofile_derivs.printField("arg",args[i]);
    for(unsigned int k=0; k<getNumberOfBasisFunctions(); k++) {
      ofile_values.printField(getBasisFunctionLabel(k),values[i][k]);
      ofile_derivs.printField("d_"+getBasisFunctionLabel(k),derivs[i][k]);
    }
    ofile_values.printField();
    ofile_derivs.printField();
  }
  ofile_values.fmtField();
  ofile_derivs.fmtField();

  delete value_pntr[0]; value_pntr.clear();

}


}
}
