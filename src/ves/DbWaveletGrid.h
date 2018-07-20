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
#include "tools/Grid.h"
#include "tools/Log.h"
#include "tools/Matrix.h"
#include <unordered_map>

namespace PLMD {
namespace ves {


// Grid that holds the values and 1st derivative of Daubechies wavelets
// mustn't be an inherited class at the moment, but that way new interpolation options can be added easily
class DbWaveletGrid : public Grid {
  void setup_Grid(const unsigned order, const unsigned gridsize);
  // helper functions to set up Grid
  // lookup function for the filter coefficients
  static std::vector<double> get_filter_coefficients(const unsigned order);
  // Fills the coefficient matrices needed for the cascade algorithm
  static std::vector<Matrix<double>> setup_Matrices(const std::vector<double>& h);
  // calculates the values of the Wavelet or its derivatives at the Integer points
  static std::vector<double> calc_integer_values(const Matrix<double>& M, const int deriv);
  // get eigenvector of square matrix A corresponding to some eigenvalue via SVD decomposition
  static std::vector<double> get_eigenvector(const Matrix<double>& A, const double eigenvalue);
  // calculate the values of the Wavelet or its derivative via the vector cascade algorithm
  static std::unordered_map<std::string, std::vector<double>> cascade(std::vector<Matrix<double>>& Matvec, const std::vector<double>& values_at_integers, unsigned recursion_number, unsigned bins_per_int, unsigned derivnum);

public:
  // constructor that directly creates Grid with right properties
  DbWaveletGrid(const unsigned order, const unsigned gridsize) : Grid("db_wavelet", {"position"}, {"0"}, {std::to_string(order*2-1)}, {gridsize}, false, true, true, {false}, {"0."}, {"0."}) {
    this->setup_Grid(order, gridsize);
  };

};


}
}

#endif
