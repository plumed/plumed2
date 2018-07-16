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
#include "../lapack/lapack.h"
#include "tools/Grid.h"
#include "tools/Matrix.h"
#include <unordered_map>

namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_DB_WAVELETS
/*
Daubechies Wavelets as basis functions

order: number of vanishing moments

Support (of the scaling function) is then from 0 to 2*order - 1

Number of basis functions is therefore also 2*order - 1


Method of construction: Strang, Nguyen - Vector cascade algorithm

\par Examples

\par Test

*/
//+ENDPLUMEDOC

class BF_DbWavelets : public BasisFunctions {
  // Grid that holds the Wavelet values and its derivative
  std::unique_ptr<Grid> Wavelet_Grid_;
  virtual void setupLabels();
  void setup_Wavelet_Grid(const unsigned recursion_number);
  static std::vector<double> get_filter_coefficients(const unsigned order);
  static void setup_Matrices(Matrix<double> &M0, Matrix<double> &M1, const std::vector<double> h);
  std::vector<double> get_eigenvector(const Matrix<double> &A, const double eigenvalue);
  std::vector<double> calc_integer_values(const Matrix<double> &M, const int deriv);
  std::unordered_map<std::string, std::vector<double>> cascade(std::vector<Matrix<double>>& Matvec, const std::vector<double>& values_at_integers, unsigned recursion_number, unsigned bins_per_int, unsigned derivnum);

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
  setNumberOfBasisFunctions((getOrder()*2)-1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  setNonPeriodic();
  setIntervalBounded();
  // gridbins_ defines the number of recursion steps at the Wavelet construction. Maybe change the name? Set by Keyword, how to combine both intuitively for later developers/users?
  unsigned gridsize = 1000;
  parse("GRID_SIZE", gridsize);
  if(gridsize!=1000) {addKeywordToList("GRID_SIZE",gridsize);}
  setup_Wavelet_Grid(gridsize);
  bool dump_wavelet_grid=false;
  parseFlag("DUMP_WAVELET_GRID", dump_wavelet_grid);
  if (dump_wavelet_grid) {
    OFile wavelet_gridfile;
    wavelet_gridfile.open("wavelet_grid.data");
    Wavelet_Grid_->writeToFile(wavelet_gridfile);
  }
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
  if(!inside_range) {for(size_t i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


// label according to positions?
void BF_DbWavelets::setupLabels() {
  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert(i,is);
    setLabel(i,"i"+is);
  }
}


// Creates and fills the Grid with the Wavelet values
void BF_DbWavelets::setup_Wavelet_Grid(const unsigned gridsize) {
  // NumberOfBasisFunctions is equal to the maximum support of the Wavelet scaling function
  std::string maxsupport; Tools::convert(getNumberOfBasisFunctions(), maxsupport);
  // determine needed recursion depth for specified size
  unsigned recursion_number = 0;
  while (getNumberOfBasisFunctions()*(1<<recursion_number) < gridsize) recursion_number++;
  unsigned bins_per_int = 1<<recursion_number;
  // define number of grid bins
  const std::vector<unsigned> gridbins {getNumberOfBasisFunctions() * (bins_per_int)};
  // set up Grid
  Wavelet_Grid_.reset(new Grid("db_wavelet", {"position"}, {"0"}, {maxsupport}, gridbins, false, true, true, {false}, {"0."}, {"0."}));
  //log.printf("\nProperties of Grid:\n Size: %d\nHasderivs: %d\nPeriodic: %d\n\n",Wavelet_Grid_->getSize(), Wavelet_Grid_->hasDerivatives(), Wavelet_Grid_->getIsPeriodic().at(0));


  // Filter coefficients
  std::vector<double> h_coeffs = get_filter_coefficients(getOrder());
  // Matrices for the cascade
  Matrix<double> M0, M1;
  // maybe change the function to return the vector directly --> less memory
  // or alternatively use a vector of pointers
  setup_Matrices(M0, M1, h_coeffs);
  std::vector<Matrix<double>> Matvec {M0, M1};

  // get the values at integers
  std::vector<double> values_at_integers = calc_integer_values(Matvec[0], 0);
  std::vector<double> derivs_at_integers = calc_integer_values(Matvec[0], 1);


  // do the cascade algorithm
  std::unordered_map<std::string, std::vector<double>> valuesmap = cascade(Matvec, values_at_integers, recursion_number, bins_per_int, 0);
  std::unordered_map<std::string, std::vector<double>> derivsmap = cascade(Matvec, derivs_at_integers, recursion_number, bins_per_int, 1);

  // Fill the Grid with the values of the unordered maps
  // this is somewhat complicated… not sure if the unordered_map way is the best way for c++
  for (const auto& value_element: valuesmap) {
    // get decimal of binary key
    int decimal = std::stoi(value_element.first, nullptr, 2);
    // corresponding iterator of deriv 
    auto deriv_iter = derivsmap.find(value_element.first);
    // calculate first grid element (this looks too complicated)
    unsigned first_grid_element = decimal * 1<<(recursion_number - value_element.first.length());
    if (value_element.first == "0") {
      log << "Integer values: ";
      for (auto ele : value_element.second) log << ele << " ";
    }
    for (unsigned j=0; j<value_element.second.size(); ++j) {
      // derivative(s) has to be in vector
      std::vector<double> deriv {deriv_iter->second.at(j)};
      Wavelet_Grid_->setValueAndDerivatives(first_grid_element + bins_per_int*j, value_element.second[j], deriv);
    }
  }

}


// returns the filter coefficients, at the moment simply a lookup table (copied from Mathematica)
std::vector<double> BF_DbWavelets::get_filter_coefficients(const unsigned order) {
  std::vector<double> h;
  switch(order) {
  case 4:
    h = { 0.16290171402564917413726008653520,
          0.50547285754591443144136667385662,
          0.44610006912337981158583475265059,
          -0.019787513117822321547681334086982,
          -0.13225358368451986802584241445340,
          0.021808150237088626328869997057481,
          0.023251800535490882302747575267605,
          -0.0074934946651807362225553368271239
        };
    break;
  case 6:
    h = { 0.078871216001450708360703821762941,
          0.34975190703761783105607104498368,
          0.53113187994086898454751440735467,
          0.22291566146501775627367242887952,
          -0.15999329944606139494194938783162,
          -0.091759032030147576133204962992087,
          0.068944046487372298805285738485360,
          0.019461604854164664143361603520751,
          -0.022331874165094534628441049888253,
          0.00039162557614857788770574331926167,
          0.0033780311814639378568864701169004,
          -0.00076176690280125322760585771112944
        };
    break;
  default:
    plumed_merror("Specified order currently not implemented");
  }
  return h;
}

// Fills the coefficient matrices needed for the cascade algorithm
void BF_DbWavelets::setup_Matrices(Matrix<double> &M0, Matrix<double> &M1, const std::vector<double> h_coeffs) {
  const int N = h_coeffs.size() -1;
  M0.resize(N,N); M1.resize(N,N);
  for (int i = 0; i < N; ++i) { // not very elegant, maybe change this later
    for (int j = 0; j < N; ++j) {
      int shift = 2*i -j;
      if (0 <= shift && shift <= N) {
        M0[i][j] = 2 * h_coeffs[2*i -j];
      }
      if (-1 <= shift && shift <= N -1) {
        M1[i][j] = 2 * h_coeffs[2*i -j + 1];
      }
    }
  }
}


// calculates the values of the Wavelet or its derivatives at the Integer points
std::vector<double> BF_DbWavelets::calc_integer_values(const Matrix<double> &M, const int deriv) {
  // corresponding eigenvalue of the matrix
  double eigenvalue = pow(0.5, deriv);
  std::vector<double> values = get_eigenvector(M, eigenvalue);

  // normalization of the eigenvector
  double normfactor = 0.;
  // i=0 is always 0; for deriv > 1 there is an additional factorial term missing
  for (unsigned i=1; i<values.size(); ++i) {
    normfactor += values[i] * pow(-i, deriv);
  }
  normfactor = 1/normfactor;
  for (auto& value : values) {
    value *= normfactor;
  }
  log << "Values at determination:\n";
  for (auto& value : values) {
    log << value << " ";
  }
  log << "\n";

  return values;
}


// maybe move this to the tools/matrix.h file?
// get eigenvector of square matrix A corresponding to some eigenvalue via SVD decomposition; this works reliably only for singular eigenvalues
std::vector<double> BF_DbWavelets::get_eigenvector(const Matrix<double> &A, const double eigenvalue) {
  // mostly copied from tools/matrix.h
  int info, N = A.ncols(); // ncols == nrows
  std::vector<double> da(N*N);
  std::vector<double> S(N);
  std::vector<double> U(N*N);
  std::vector<double> VT(N*N);
  std::vector<int> iwork(8*N);

  // Transfer the matrix to the local array and substract eigenvalue
  for (int i=0; i<N; ++i) for (int j=0; j<N; ++j) {
      da[i*N+j]=static_cast<double>( A(j,i) );
      if (i==j) da[i*N+j] -= eigenvalue;
    }

  // This optimizes the size of the work array used in lapack singular value decomposition
  int lwork=-1;
  std::vector<double> work(1);
  plumed_lapack_dgesdd( "A", &N, &N, da.data(), &N, S.data(), U.data(), &N, VT.data(), &N, work.data(), &lwork, iwork.data(), &info );

  // Retrieve correct sizes for work and reallocate
  lwork=(int) work[0]; work.resize(lwork);

  // This does the singular value decomposition
  plumed_lapack_dgesdd( "A", &N, &N, da.data(), &N, S.data(), U.data(), &N, VT.data(), &N, work.data(), &lwork, iwork.data(), &info );

  // fill eigenvector with last column of VT
  std::vector<double> eigenvector;
  for (int i=0; i<N; ++i) eigenvector.push_back(VT[N-1 + i*N]);

  return eigenvector;
}


// Calculate the values of the Wavelet or its derivative at 
std::unordered_map<std::string, std::vector<double>> BF_DbWavelets::cascade(std::vector<Matrix<double>>& Matvec, const std::vector<double>& values_at_integers, unsigned recursion_number, unsigned bins_per_int, unsigned derivnum) {
  // map for the cascade scheme with binary representation of the decimal part as keys
  std::unordered_map<std::string, std::vector<double>> binarymap;
  binarymap.reserve(bins_per_int);
  // vector to store the binary representation of all the decimal parts
  std::vector<std::string> binaryvec;
  // vector used as result of the matrix multiplications
  std::vector<double> new_values; // better name?!

  // multiply matrices by 2 if derivatives are calculated
  if (derivnum != 0) for (auto& M : Matvec) M *= 2;

  // fill the first two datasets by hand
  binarymap["0"] = values_at_integers;
  mult(Matvec[1], values_at_integers, new_values);
  binarymap["1"] = new_values;
  binaryvec.push_back("1");

  for (unsigned i=1; i<recursion_number; ++i) {
    std::vector<std::string> new_binaryvec;
    for (auto binary : binaryvec) {
      for (int k=0; k<2; ++k) {
        // prepend the new bit
        std::string new_binary = std::to_string(k) + binary;
        //log << "Current binary: " + new_binary + "\n";
        mult(Matvec[k], binarymap[binary], new_values);
        binarymap.insert(std::pair<std::string, std::vector<double>>(new_binary, new_values));
        new_binaryvec.push_back(new_binary);
      }
    }
    binaryvec = new_binaryvec;
  }
  return binarymap;
}

}
}
