/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2021 The VES code team
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

#include "WaveletGrid.h"
#include "tools/Exception.h"
#include "tools/Grid.h"
#include "tools/Matrix.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace PLMD {
namespace ves {


// construction of Wavelet grid according to the Daubechies-Lagarias method
// see Strang, Nguyen "Wavelets and Filter Banks" chapter 6.3
std::unique_ptr<Grid> WaveletGrid::setupGrid(const unsigned order, unsigned gridsize, const bool use_mother_wavelet, const Type type) {
  plumed_assert(order>=1) << "Wavelet order has to be a positive integer";
  // calculate the grid properties of the scaling grid
  // the range of the grid is from 0 to maxsupport
  unsigned maxsupport = order*2 -1;
  // determine needed recursion depth for specified size
  unsigned recursion_number = 0;
  while (maxsupport*(1U<<recursion_number) < gridsize) { recursion_number++; }
  // set "true" gridsize
  unsigned bins_per_int = 1U<<recursion_number;
  gridsize = maxsupport*bins_per_int;

  // Filter coefficients
  std::vector<double> h_coeffs = getFilterCoefficients(order, true, type);
  // Vector with the Matrices M0 and M1 for the cascade
  std::vector<Matrix<double>> h_Matvec = setupMatrices(h_coeffs);
  std::vector<Matrix<double>> g_Matvec; // only filled if needed for wavelet

  // get the values at integers
  std::vector<double> values_at_integers = calcIntegerValues(h_Matvec[0], 0);
  std::vector<double> derivs_at_integers = calcIntegerValues(h_Matvec[0], 1);

  std::string gridvarname; // stores the name of the grid variable
  if (use_mother_wavelet) { // get the highpass filter coefficients as well
    std::vector<double> g_coeffs = getFilterCoefficients(order, false, type);
    g_Matvec = setupMatrices(g_coeffs);
    gridvarname = typeToString(type,true)+std::to_string(order)+"_psi";
  }
  else {
    gridvarname = typeToString(type,true)+std::to_string(order)+"_phi";
  }

  // Set up the grid with correct properties
  auto grid = std::unique_ptr<Grid>(new Grid(gridvarname, {"position"}, {"0"},
  {std::to_string(maxsupport)}, {gridsize}, false, true, {false}, {"0."}, {"0."}));

  BinaryMap values = cascade(h_Matvec, g_Matvec, values_at_integers, recursion_number, bins_per_int, 0, use_mother_wavelet);
  BinaryMap derivs = cascade(h_Matvec, g_Matvec, derivs_at_integers, recursion_number, bins_per_int, 1, use_mother_wavelet);

  fillGridFromMaps(grid, values, derivs);

  return grid;
}



std::vector<Matrix<double>> WaveletGrid::setupMatrices(const std::vector<double>& coeffs) {
  Matrix<double> M0, M1;
  const int N = coeffs.size() -1;
  M0.resize(N,N); M1.resize(N,N);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      int shift = 2*i -j;
      if (0 <= shift && shift <= N) {
        M0[i][j] = coeffs[2*i -j]; // normalization is already included in coeffs
      }
      if (-1 <= shift && shift <= N -1) {
        M1[i][j] = coeffs[2*i -j + 1];
      }
    }
  }
  return std::vector<Matrix<double>> {M0, M1};
}


std::vector<double> WaveletGrid::calcIntegerValues(const Matrix<double> &M, const int deriv) {
  // corresponding eigenvalue of the matrix
  double eigenvalue = std::pow(0.5, deriv);
  std::vector<double> values = getEigenvector(M, eigenvalue);

  // normalization of the eigenvector
  double normfactor = 0.;
  // i=0 is always 0; for derivs > 1 an additional factorial would have to be added
  for (unsigned i=1; i<values.size(); ++i) {
    normfactor += values[i] * std::pow(-i, deriv);
  }
  normfactor = 1/normfactor;
  for (auto& value : values) {
    value *= normfactor;
  }

  return values;
}


// maybe move this to the tools/matrix.h file?
// this works reliably only for singular eigenvalues
//
std::vector<double> WaveletGrid::getEigenvector(const Matrix<double> &A, const double eigenvalue) {
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
  lwork=static_cast<int>(work[0]); work.resize(lwork);

  // This does the singular value decomposition
  plumed_lapack_dgesdd( "A", &N, &N, da.data(), &N, S.data(), U.data(), &N, VT.data(), &N, work.data(), &lwork, iwork.data(), &info );

  // fill eigenvector with last column of VT
  std::vector<double> eigenvector(N);
  for (int i=0; i<N; ++i) eigenvector[i] = VT[N-1 + i*N];

  return eigenvector;
}


WaveletGrid::BinaryMap WaveletGrid::cascade(std::vector<Matrix<double>>& h_Matvec, std::vector<Matrix<double>>& g_Matvec, const std::vector<double>& values_at_integers, const unsigned recursion_number, const unsigned bins_per_int, const unsigned derivnum, const bool use_mother_wavelet) {
  BinaryMap scaling_map, wavelet_map;
  scaling_map.reserve(bins_per_int);
  // vector to store the binary representation of all the decimal parts
  std::vector<std::string> binary_vec;
  // vector used as result of the matrix multiplications
  std::vector<double> temp_values;

  // multiply matrices by 2 if derivatives are calculated (assumes ascending order)
  if (derivnum != 0) for (auto& M : h_Matvec) M *= 2;

  if (use_mother_wavelet) {
    wavelet_map.reserve(bins_per_int);
    if (derivnum != 0) for (auto& M : g_Matvec) M *= 2;
  }

  // fill the first two datasets by hand
  scaling_map["0"] = values_at_integers;
  mult(h_Matvec[1], values_at_integers, temp_values);
  scaling_map["1"] = temp_values;

  if (use_mother_wavelet) {
    mult(g_Matvec[0], values_at_integers, temp_values);
    wavelet_map["0"] = temp_values;
    mult(g_Matvec[1], values_at_integers, temp_values);
    wavelet_map["1"] = temp_values;
  }

  // now do the cascade
  binary_vec.emplace_back("1");
  for (unsigned i=1; i<recursion_number; ++i) {
    std::vector<std::string> new_binary_vec;
    for (const auto& binary : binary_vec) {
      for (unsigned k=0; k<2; ++k) {
        // prepend the new bit
        std::string new_binary = std::to_string(k) + binary;
        mult(h_Matvec[k], scaling_map[binary], temp_values);
        scaling_map.insert(std::pair<std::string, std::vector<double>>(new_binary, temp_values));
        if (use_mother_wavelet) {
          mult(g_Matvec[k], scaling_map[binary], temp_values);
          wavelet_map.insert(std::pair<std::string, std::vector<double>>(new_binary, temp_values));
        }
        new_binary_vec.push_back(new_binary);
      }
    }
    binary_vec = new_binary_vec;
  }

  return use_mother_wavelet ? wavelet_map : scaling_map;
}


// Fill the Grid with the values of the unordered maps
void WaveletGrid::fillGridFromMaps(std::unique_ptr<Grid>& grid, const BinaryMap& values_map, const BinaryMap& derivs_map) {
  unsigned bins_per_int = values_map.size();
  // this is somewhat complicated… not sure if the unordered_map way is the best way for c++
  for (const auto& value_iter : values_map) {
    // get decimal of binary key
    unsigned decimal = std::stoi(value_iter.first, nullptr, 2);
    // corresponding iterator of deriv
    auto deriv_iter = derivs_map.find(value_iter.first);
    // calculate first grid element
    unsigned first_grid_element = decimal * (bins_per_int >> value_iter.first.length());
    for (unsigned i=0; i<value_iter.second.size(); ++i) {
      // derivative has to be passed as vector
      std::vector<double> deriv {deriv_iter->second[i]};
      grid->setValueAndDerivatives(first_grid_element + bins_per_int*i, value_iter.second[i], deriv);
    }
  }
}


WaveletGrid::Type WaveletGrid::stringToType(std::string& type_str) {
  std::unordered_map<std::string, Type> typemap = { {"DAUBECHIES", Type::db}, {"SYMLETS", Type::sym} };
  try { return typemap.at(type_str); }
  catch(const std::out_of_range& e) {plumed_merror("The specified wavelet type "+type_str+" is not implemented."); }
}


std::string WaveletGrid::typeToString(Type type, bool abbrev) {
  std::string type_str;
  switch(type) {
  case Type::db:
    type_str = abbrev ? "Db" : "DAUBECHIES";
    break;
  case Type::sym:
    type_str = abbrev ? "Sym" : "SYMLETS";
    break;
  }
  return type_str;
}


}
}
