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

#include "BasisFunctions.h"
#include "VesTools.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_COMBINED
/*
Combining other basis functions types

\par Examples

Here we define both Fourier cosine and sine expansions of order 10,
each with 11 basis functions, which are combined. This results
in a total number of 21 basis functions as only the constant from
is bf_cos is used.
\plumedfile
bf_cos: BF_COSINE MINIMUM=-pi MAXIMUM=+pi ORDER=10
bf_sin: BF_SINE   MINIMUM=-pi MAXIMUM=+pi ORDER=10
bf_comb: BF_COMBINED BASIS_FUNCTIONS=bf_cos,bf_sin
\endplumedfile
In principle this is the same as using BF_FOURIER with
ORDER=10 but with different ordering of the basis functions.
Note that the order used in BASIS_FUNCTIONS matters for the ordering
of the basis functions, using BASIS_FUNCTIONS=bf_sin,bf_cos would
results in a different order of the basis functions.
This should be kept in mind when restarting from previous
coefficients.






*/
//+ENDPLUMEDOC

class BF_Combined : public BasisFunctions {
  std::vector<BasisFunctions*> basisf_pntrs_;
  void setupLabels() override;
  void setupUniformIntegrals() override;
  // void getBFandValueIndices(const unsigned int, unsigned int&, unsigned int&) const;
public:
  static void registerKeywords(Keywords&);
  explicit BF_Combined(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Combined,"BF_COMBINED")


void BF_Combined::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.remove("ORDER");
  keys.remove("MAXIMUM");
  keys.remove("MINIMUM");
  keys.remove("NUMERICAL_INTEGRALS");
  keys.remove("NGRID_POINTS");
  keys.add("compulsory","BASIS_FUNCTIONS","Labels of the basis functions that should be combined. Note that the order used matters for the ordering of the basis functions. This needs to be kept in mind when restarting from previous coefficients.");
}


BF_Combined::BF_Combined(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  basisf_pntrs_(0)
{
  std::vector<std::string> basisf_labels;
  parseVector("BASIS_FUNCTIONS",basisf_labels); addKeywordToList("BASIS_FUNCTIONS",basisf_labels);
  if(basisf_labels.size()==1) {
    plumed_merror("using only one basis function in BF_COMBINED does not make sense");
  }
  std::string error_msg = "";
  basisf_pntrs_ = VesTools::getPointersFromLabels<BasisFunctions*>(basisf_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword BASIS_FUNCTIONS of "+getName()+": "+error_msg);}

  unsigned int nbasisf_total_ = 1;
  bool periodic = true;
  for(unsigned int i=0; i<basisf_pntrs_.size(); i++) {
    nbasisf_total_ += basisf_pntrs_[i]->getNumberOfBasisFunctions() - 1;
    if(basisf_pntrs_[i]->intervalMinStr()!=basisf_pntrs_[0]->intervalMinStr() || basisf_pntrs_[i]->intervalMaxStr()!=basisf_pntrs_[0]->intervalMaxStr()) {
      plumed_merror("all the basis functions to be combined should have same MINIMUM and MAXIMUM");
    }
    if(!basisf_pntrs_[i]->arePeriodic()) {periodic=false;}
  }
  setOrder(nbasisf_total_-1);
  setNumberOfBasisFunctions(nbasisf_total_);
  setInterval(basisf_pntrs_[0]->intervalMinStr(),basisf_pntrs_[0]->intervalMaxStr());
  setIntrinsicInterval("-1.0","+1.0");
  if(periodic) {setPeriodic();}
  else {setNonPeriodic();}
  setIntervalBounded();
  setType("combined");
  setDescription("Combined");
  setupBF();
  checkRead();
}


// void BF_Combined::getBFandValueIndices(const unsigned int n, unsigned int& bf_index, unsigned int& value_index) const {
//   bf_index = 0; value_index = 0;
//   if(n==0){
//     bf_index = 0;
//     value_index = 0;
//     return;
//   }
//   else{
//     unsigned int r=1;
//     for(unsigned int i=0; i<basisf_pntrs_.size(); i++){
//       for(unsigned int l=1; l<basisf_pntrs_[i]->numberOfBasisFunctions(); l++){
//         if(r==n){
//           bf_index = i;
//           value_index = l;
//           return;
//         }
//         r++;
//       }
//     }
//   }
// }


void BF_Combined::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // first BF, constant, argT, and inside_range taken from here
  unsigned int r=0;
  std::vector<double> values_tmp(basisf_pntrs_[0]->numberOfBasisFunctions(),0.0);
  std::vector<double> derivs_tmp(basisf_pntrs_[0]->numberOfBasisFunctions(),0.0);
  basisf_pntrs_[0]->getAllValues(arg,argT,inside_range,values_tmp,derivs_tmp);
  for(unsigned int l=0; l<basisf_pntrs_[0]->numberOfBasisFunctions(); l++) {
    values[r] = values_tmp[l];
    derivs[r] = derivs_tmp[l];
    r++;
  }
  // other BF
  for(unsigned int i=1; i<basisf_pntrs_.size(); i++) {
    values_tmp.assign(basisf_pntrs_[i]->numberOfBasisFunctions(),0.0);
    derivs_tmp.assign(basisf_pntrs_[i]->numberOfBasisFunctions(),0.0);
    double dummy_dbl; bool dummy_bool=true;
    basisf_pntrs_[i]->getAllValues(arg,dummy_dbl,dummy_bool,values_tmp,derivs_tmp);
    for(unsigned int l=1; l<basisf_pntrs_[i]->numberOfBasisFunctions(); l++) {
      values[r] = values_tmp[l];
      derivs[r] = derivs_tmp[l];
      r++;
    }
  }
}


void BF_Combined::setupLabels() {
  setLabel(0,basisf_pntrs_[0]->getBasisFunctionLabel(0));
  unsigned int r=1;
  for(unsigned int i=0; i<basisf_pntrs_.size(); i++) {
    for(unsigned int l=1; l<basisf_pntrs_[i]->numberOfBasisFunctions(); l++) {
      setLabel(r,basisf_pntrs_[i]->getBasisFunctionLabel(l));
      r++;
    }
  }
}


void BF_Combined::setupUniformIntegrals() {
  setUniformIntegral(0,basisf_pntrs_[0]->getUniformIntegrals()[0]);
  unsigned int r=1;
  for(unsigned int i=0; i<basisf_pntrs_.size(); i++) {
    std::vector<double> uniform_tmp = basisf_pntrs_[i]->getUniformIntegrals();
    for(unsigned int l=1; l<basisf_pntrs_[i]->numberOfBasisFunctions(); l++) {
      setUniformIntegral(r,uniform_tmp[l]);
      r++;
    }
  }
}


}
}
