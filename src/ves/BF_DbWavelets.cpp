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

#include "core/ActionRegister.h"
#include "tools/Grid.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_DB_WAVELETS
/*
Daubechies Wavelets as basis functions

order: number of vanishing moments

Support (of the scaling function) is then from 0 to 2*order - 1

Number of basis functions is therefore also 2*order - 1

\par Examples

\par Test

*/
//+ENDPLUMEDOC

class BF_DbWavelets : public BasisFunctions {
  std::unique_ptr<Grid> Wavelet_Grid_;
  virtual void setupLabels();
public:
  static void registerKeywords( Keywords&);
  explicit BF_DbWavelets(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
  void setup_Wavelet_Grid();
};


PLUMED_REGISTER_ACTION(BF_DbWavelets,"BF_DB_WAVELETS")


void BF_DbWavelets::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  // why is this removed?
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_DbWavelets::BF_DbWavelets(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions((getOrder()*2)-1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  setNonPeriodic();
  setIntervalBounded();
  setup_Wavelet_Grid();
  setType("daubechies_wavelets");
  setDescription("Daubechies Wavelets (minimum phase type)");
  setLabelPrefix("k");
  setupBF();
  checkRead();
}


void BF_DbWavelets::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    values[i] = 0.0;
    derivs[i] = 0.0;
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


// label according to positions?
void BF_DbWavelets::setupLabels() {
  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert(i,is);
    setLabel(i,"i"+is);
  }
}

// Creates and fills the Grid with the Wavelet values
void BF_DbWavelets::setup_Wavelet_Grid() {
  for (int i=0; i<10; ++i) {
    std::vector<double> derivtest = {0.3};
    double gridval = i*0.5;
    Wavelet_Grid_->addValueAndDerivatives(i, gridval, derivtest);
  }
  Ofile wv_gridfile;

  wv_gridfile.open("wv_griddump");
  Wavelet_Grid_->writeToFile(wv_gridfile)
}


}
}
