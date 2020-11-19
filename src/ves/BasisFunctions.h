/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
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
#ifndef __PLUMED_ves_BasisFunctions_h
#define __PLUMED_ves_BasisFunctions_h

#include "core/Action.h"

#include <vector>
#include <string>
#include <cmath>


#define PLUMED_VES_BASISFUNCTIONS_INIT(ao) BasisFunctions(ao)

namespace PLMD {

/**
\ingroup INHERIT
Abstract base class for implenting new 1D basis sets.
*/

class Action;
class Grid;

namespace ves {

class VesBias;
class TargetDistribution;

class BasisFunctions :
  public Action
{
private:
  // print extra info about the basis set
  bool print_debug_info_;
  // to check if the basis set has been defined
  bool has_been_set;
  // description of the basis set
  std::string description_;
  // the type of the basis set
  std::string type_;
  // the maximum order of the basis functions
  unsigned int norder_;
  // the total number of basis functions
  unsigned int nbasis_;
  // the keywords used to invoke the basis set
  std::vector<std::string> bf_keywords_;
  // prefix for the basis function labels
  std::string bf_label_prefix_;
  // label of each basis function
  std::vector<std::string> bf_labels_;
  // if the basis functions are periodic or not
  bool periodic_;
  // if the basis functions are defined on a bounded interval or not
  bool interval_bounded_;
  // the intrinsic interval of the basis functions
  std::string interval_intrinsic_min_str_;
  std::string interval_intrinsic_max_str_;
  double interval_intrinsic_min_;
  double interval_intrinsic_max_;
  double interval_intrinsic_range_;
  double interval_intrinsic_mean_;
  // the defined (translated) interval of the basis functions
  std::string interval_min_str_;
  std::string interval_max_str_;
  double interval_min_;
  double interval_max_;
  double interval_range_;
  double interval_mean_;
  // the derivative term in the chain rule coming from the translation of the interval
  double argT_derivf_;
  // calculate numerically the integrals of the basis functions over the intervals
  bool numerical_uniform_integrals_;
  unsigned int nbins_;
  // the integrals of the basis functions over the interval on which they are defined
  std::vector <double> uniform_integrals_;
  //
  VesBias* vesbias_pntr_;
  Action* action_pntr_;
  //
  void getAllValuesNumericalDerivs(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;

protected:
  // setup various stuff
  void setupBF();
  void setupInterval();
  void setNumericalIntegrationBins(const unsigned int nbins) {nbins_=nbins;}
  void numericalUniformIntegrals();
  std::vector<double> numericalTargetDistributionIntegralsFromGrid(const Grid*) const ;
  virtual void setupLabels();
  virtual void setupUniformIntegrals();
  template<typename T>
  void addKeywordToList(const std::string&, const T);
  template<typename T>
  void addKeywordToList(const std::string&, const std::vector<T>&);
  void addKeywordToList(const std::string&, const bool);
  //
  void setPeriodic() {periodic_=true;}
  void setNonPeriodic() {periodic_=false;}
  void setIntervalBounded() {interval_bounded_=true;}
  void setIntervalNonBounded() {interval_bounded_=false;}
  void setType(const std::string& type_in) {type_=type_in;}
  void setDescription(const std::string& description_in) {description_=description_in;}
  //
  void setNumberOfBasisFunctions(const unsigned int);
  void setOrder(const unsigned int norder_in) {norder_=norder_in;}
  void setIntrinsicInterval(const double, const double);
  void setIntrinsicInterval(const std::string&, const std::string&);
  void setInterval(const double, const double);
  void setInterval(const std::string&, const std::string&);
  //
  double intrinsicIntervalMin() const {return interval_intrinsic_min_;}
  double intrinsicIntervalMax() const {return interval_intrinsic_max_;}
  std::string intrinsicIntervalMinStr() const {return interval_intrinsic_min_str_;}
  std::string intrinsicIntervalMaxStr() const {return interval_intrinsic_max_str_;}
  //
  void setUniformIntegral(const unsigned int, const double);
  void setUniformIntegrals(const std::vector<double>&);
  void setAllUniformIntegralsToZero();
  //
  void setLabelPrefix(const std::string&);
  void setLabel(const unsigned int, const std::string&);
  void setLabels(const std::vector<std::string>&);

public:
  static void registerKeywords(Keywords&);
  explicit BasisFunctions(const ActionOptions&ao);
  bool hasBeenSet() const {return has_been_set;}
  std::string getType() const {return type_;}
  std::string getDescription() const {return description_;}
  unsigned int getOrder() const {return norder_;}
  unsigned int getNumberOfBasisFunctions() const {return nbasis_;}
  unsigned int numberOfBasisFunctions() const {return nbasis_;}
  unsigned int getSize() const {return nbasis_;}
  bool arePeriodic() const {return periodic_;}
  bool intervalBounded() const {return interval_bounded_;}
  double intervalMin() const {return interval_min_;}
  double intervalMax() const {return interval_max_;}
  double intervalRange() const {return interval_range_;}
  double intervalMean() const {return interval_mean_;}
  double intervalDerivf() const {return argT_derivf_;}
  std::string intervalMinStr() const {return interval_min_str_;}
  std::string intervalMaxStr() const {return interval_max_str_;}
  std::vector<double> getUniformIntegrals() const {return uniform_integrals_;}
  std::vector<double> getTargetDistributionIntegrals(const TargetDistribution*) const;
  //
  std::vector<std::string> getKeywordList() const {return bf_keywords_;}
  std::string getKeywordString() const;
  //
  std::string getBasisFunctionLabel(const unsigned int index) const {return bf_labels_[index];}
  std::vector<std::string> getBasisFunctionLabels() const {return bf_labels_;}
  //
  void linkVesBias(VesBias*);
  void linkAction(Action*);
  VesBias* getPntrToVesBias() const;
  Action* getPntrToAction() const;
  //
  double translateArgument(const double, bool&) const;
  double checkIfArgumentInsideInterval(const double, bool&) const;
  //
  void apply() override {};
  void calculate() override {};
  // calculate the value for the n-th basis function
  double getValue(const double, const unsigned int, double&, bool&) const;
  // calculate the values for all basis functions
  virtual void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const = 0;
  //virtual void get2ndDerivatives(const double, std::vector<double>&)=0;
  void printInfo() const;
  //
  void getMultipleValue(const std::vector<double>&, std::vector<double>&, std::vector<std::vector<double> >&, std::vector<std::vector<double> >&, const bool numerical_deriv=false) const;
  void writeBasisFunctionsToFile(OFile&, OFile&, const std::string& min_in, const std::string& max_in, unsigned int nbins=1000, const bool ignore_periodicity=false, const std::string& output_fmt="%15.8f", const bool numerical_deriv=false) const;
};


inline
void BasisFunctions::setNumberOfBasisFunctions(const unsigned int nbasis_in) {
  nbasis_=nbasis_in;
  bf_labels_.assign(nbasis_,"");
  uniform_integrals_.assign(nbasis_,0.0);
}


inline
VesBias* BasisFunctions::getPntrToVesBias() const {
  plumed_massert(vesbias_pntr_!=NULL,"the VES bias has not been linked");
  return vesbias_pntr_;
}


inline
Action* BasisFunctions::getPntrToAction() const {
  plumed_massert(action_pntr_!=NULL,"the action has not been linked");
  return action_pntr_;
}


inline
void BasisFunctions::setUniformIntegral(const unsigned index, const double value) {
  uniform_integrals_[index] = value;
}


inline
void BasisFunctions::setUniformIntegrals(const std::vector<double>& uniform_integrals_in) {
  plumed_assert(uniform_integrals_in.size()==nbasis_);
  uniform_integrals_ = uniform_integrals_in;
}


inline
void BasisFunctions::setAllUniformIntegralsToZero() {
  uniform_integrals_.assign(nbasis_,0.0);
}

inline
void BasisFunctions::setLabelPrefix(const std::string& bf_label_prefix_in) {
  bf_label_prefix_ = bf_label_prefix_in;
}


inline
void BasisFunctions::setLabel(const unsigned int index, const std::string& label) {
  bf_labels_[index] = label;
}


inline
void BasisFunctions::setLabels(const std::vector<std::string>& bf_labels_in) {
  bf_labels_ = bf_labels_in;
}


inline
double BasisFunctions::translateArgument(const double arg, bool& inside_interval) const {
  // NOTE: only works for symmetric intrinsic intervals
  inside_interval=true;
  double argT = (arg-interval_mean_)*argT_derivf_;
  if(argT < interval_intrinsic_min_) {
    inside_interval=false;
    argT=interval_intrinsic_min_;
  }
  else if(argT > interval_intrinsic_max_) {
    inside_interval=false;
    argT=interval_intrinsic_max_;
  }
  return argT;
}


inline
double BasisFunctions::checkIfArgumentInsideInterval(const double arg, bool& inside_interval) const {
  inside_interval=true;
  double argT = arg;
  if(arg < interval_min_) {
    inside_interval=false;
    argT=interval_min_;
  }
  else if(arg > interval_max_) {
    inside_interval=false;
    argT=interval_max_;
  }
  return argT;
}



template<typename T>
void BasisFunctions::addKeywordToList(const std::string& keyword, const T value) {
  std::string str_value;
  Tools::convert(value,str_value);
  bf_keywords_.push_back(keyword+"="+str_value);
}


template<typename T>
void BasisFunctions::addKeywordToList(const std::string& keyword, const std::vector<T>& values) {
  std::string str_value;
  std::string str_keywordvalues;
  Tools::convert(values[0],str_value);
  str_keywordvalues = keyword + "=" + str_value;
  for(unsigned int i=1; i<values.size(); i++) {
    Tools::convert(values[i],str_value);
    str_keywordvalues += "," + str_value;
  }
  bf_keywords_.push_back(str_keywordvalues);
}


inline
void BasisFunctions::addKeywordToList(const std::string& keyword, const bool value) {
  if(value) {bf_keywords_.push_back(keyword);}
}


}
}

#endif
