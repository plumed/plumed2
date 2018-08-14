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
#ifndef __PLUMED_ves_DbWaveletsGrid_h
#define __PLUMED_ves_DbWaveletsGrid_h

#include "../lapack/lapack.h"
#include <memory>
#include <unordered_map>
#include <vector>

namespace PLMD {

template <typename T>
class Matrix;

namespace ves {

class Grid;


// factory class to set up a Grid with DbWavelets
class DbWaveletGrid {
private:
  // unordered map with binary representation as string and corresponding values
  using BinaryMap = std::unordered_map<std::string, std::vector<double>>;


  // lookup function for the filter coefficients
  static std::vector<double> get_filter_coefficients(const unsigned order, bool lowpass);
  // Fills the coefficient matrices needed for the cascade algorithm
  static std::vector<Matrix<double>> setup_Matrices(const std::vector<double>& h);
  // calculates the values of the Wavelet or its derivatives at the Integer points
  static std::vector<double> calc_integer_values(const Matrix<double>& M, const int deriv);
  // get eigenvector of square matrix A corresponding to some eigenvalue via SVD decomposition
  static std::vector<double> get_eigenvector(const Matrix<double>& A, const double eigenvalue);
  // calculate the values of the Wavelet or its derivative via the vector cascade algorithm
  static BinaryMap cascade(std::vector<Matrix<double>>& h_Matvec, std::vector<Matrix<double>>& g_Matvec,
      const std::vector<double>& values_at_integers, unsigned recursion_number, unsigned bins_per_int,
      unsigned derivnum, bool dowavelet);
  static void fill_grid_from_map(std::unique_ptr<Grid>& grid, const BinaryMap& valuesmap,
      const BinaryMap& derivsmap);
public:
  // construct either a grid with the scaling function or the wavelet function, and its first derivative
 static std::unique_ptr<Grid> setup_Grid(const unsigned order, unsigned gridsize, bool dowavelet);
};


}
}

#endif
