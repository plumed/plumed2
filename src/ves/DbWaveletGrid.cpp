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

#include "DbWaveletGrid.h"

namespace PLMD {
namespace ves {


void DbWaveletGrid::setup_Grid(const unsigned order, const unsigned gridsize) {
  // helper variable, find recursion number
  unsigned bins_per_int = gridsize/(order*2-1);
  // check if it is a power of two
  if ((bins_per_int & (bins_per_int -1)) != 0) {
    plumed_merror("Trying to set up wavelet grid with bad size");
  }
  unsigned recursion_number = log2(bins_per_int);

  // Filter coefficients
  std::vector<double> h_coeffs = get_filter_coefficients(order);
  // Vector with the Matrices M0 and M1 for the cascade
  std::vector<Matrix<double>> Matvec = setup_Matrices(h_coeffs);
  //log << "Got coeffs and matrices\n";

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
    for (unsigned j=0; j<value_element.second.size(); ++j) {
      // derivative has to be passed as vector
      std::vector<double> deriv {deriv_iter->second.at(j)};
      this->setValueAndDerivatives(first_grid_element + bins_per_int*j, value_element.second[j], deriv);
    }
  }
}


std::vector<Matrix<double>> DbWaveletGrid::setup_Matrices(const std::vector<double>& h_coeffs) {
  Matrix<double> M0, M1;
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
  return std::vector<Matrix<double>> {M0, M1};
}


std::vector<double> DbWaveletGrid::calc_integer_values(const Matrix<double> &M, const int deriv) {
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

  return values;
}


// maybe move this to the tools/matrix.h file?
// this works reliably only for singular eigenvalues
//
std::vector<double> DbWaveletGrid::get_eigenvector(const Matrix<double> &A, const double eigenvalue) {
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


std::unordered_map<std::string, std::vector<double>> DbWaveletGrid::cascade(std::vector<Matrix<double>>& Matvec, const std::vector<double>& values_at_integers, unsigned recursion_number, unsigned bins_per_int, unsigned derivnum) {
  // map of all values with binary representation of the decimal part as keys
  std::unordered_map<std::string, std::vector<double>> binarymap;
  binarymap.reserve(bins_per_int);
  // vector to store the binary representation of all the decimal parts
  std::vector<std::string> binaryvec;
  // vector used as result of the matrix multiplications
  std::vector<double> new_values; // better name?!

  // multiply matrices by 2 if derivatives are calculated (assumes ascending order)
  if (derivnum != 0) for (auto& M : Matvec) M *= 2;

  // fill the first two datasets by hand
  binarymap["0"] = values_at_integers;
  mult(Matvec[1], values_at_integers, new_values);
  binarymap["1"] = new_values;

  // now do the cascade
  binaryvec.push_back("1");
  for (unsigned i=1; i<recursion_number; ++i) {
    std::vector<std::string> new_binaryvec;
    for (auto binary : binaryvec) {
      for (int k=0; k<2; ++k) {
        // prepend the new bit
        std::string new_binary = std::to_string(k) + binary;
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
