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
#include "DbWaveletGrid.h"
#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {


//+PLUMEDOC VES_BASISF BF_DB_WAVELETS
/*
Daubechies Wavelets as basis functions

Note: at the moment only the scaling function and not the wavelet function is used.
It should nevertheless form an orthogonal basis set and will be needed for multiscale.
The Wavelet function can be easily implemented by an additional matrix multiplication and a translation of the position axis.

order: number of vanishing moments

Support (of the scaling function) is then from 0 to 2*order - 1

Number of basis functions is therefore also 2*order - 1 + 1 constant


Method of construction: Strang, Nguyen - Vector cascade algorithm

\par Examples

\par Test

*/
//+ENDPLUMEDOC


class BF_DbWavelets : public BasisFunctions {
  // Grid that holds the Wavelet values and its derivative
  std::unique_ptr<DbWaveletGrid> Wavelet_Grid_;
  virtual void setupLabels();

public:
  static void registerKeywords( Keywords&);
  explicit BF_DbWavelets(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_DbWavelets,"BF_DB_WAVELETS")


void BF_DbWavelets::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","GRID_SIZE","The number of grid bins of the Wavelet function. Because of the used construction algorithm this value will be used as guiding value only, while the true number will be \"(ORDER*2 - 1) * 2**n\" with the smallest n such that the grid is at least as large as the specified number."); // Change the commentary a bit?
  keys.addFlag("DUMP_WAVELET_GRID", false, "If this flag is set the grid with the wavelet values will be written to a file called \"wavelet_grid.data\". Default is false.");
  // why is this removed?
  keys.remove("NUMERICAL_INTEGRALS");
}


BF_DbWavelets::BF_DbWavelets(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions((getOrder()*2));
  setIntrinsicInterval("0",std::to_string(getNumberOfBasisFunctions()-1));
  setNonPeriodic();
  setIntervalBounded();
  unsigned gridsize = 1000;
  parse("GRID_SIZE", gridsize);
  if(gridsize!=1000) {addKeywordToList("GRID_SIZE",gridsize);}
  Wavelet_Grid_.reset(new DbWaveletGrid(log, getOrder(), gridsize));
  bool dump_wavelet_grid=false;
  parseFlag("DUMP_WAVELET_GRID", dump_wavelet_grid);
  if (dump_wavelet_grid) {
    OFile wavelet_gridfile;
    wavelet_gridfile.open("wavelet_grid.data");
    Wavelet_Grid_->writeToFile(wavelet_gridfile);
  }
  setType("daubechies_wavelets");
  setDescription("Daubechies Wavelets (maximum phase type)");
  setLabelPrefix("k");
  setupBF();
  checkRead();
}


void BF_DbWavelets::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double x = arg - ((i-1)/intervalDerivf()); // shift argument by scaled i
    if (x < intervalMin() || x > intervalMax()) { // Wavelets are 0 outside the defined range
      values[i] = 0.0; derivs[i] = 0.0;
    }
    else {
      // declaring vectors and calling first a function to get the index is a bit cumbersome and might be slow
      std::vector<double> temp_deriv;
      std::vector<double> x_vec {x};
      values[i] = Wavelet_Grid_->getValueAndDerivatives(x_vec, temp_deriv);
      derivs[i] = temp_deriv[0];
    }
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
  return;
}


// label according to positions?
void BF_DbWavelets::setupLabels() {
  setLabel(0,"const");
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert((i-1)/intervalDerivf(),is);
    setLabel(i,"i = "+is);
  }
}


}
}
