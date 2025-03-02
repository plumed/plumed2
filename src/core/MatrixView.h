/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_core_MatrixView_h
#define __PLUMED_core_MatrixView_h

#include "Value.h"

namespace PLMD {

class MatrixView {
public:
  std::size_t start;
  std::vector<std::size_t> shape;
  std::size_t ncols;
  std::vector<std::size_t> bookeeping;
/// Pass the data from the value to the matrix
  void setup( std::size_t s, const Value* myval );
/// Get an element of the matrix
  static double getElement( std::size_t irow, std::size_t jcol, const MatrixView& mat, std::vector<double>& data );
/// Determine if there is an element
  static bool hasElement( std::size_t irow, std::size_t jcol, const MatrixView& mat, std::size_t& ind );
};

} // namespace PLMD
#endif // __PLUMED_tools_MatrixView_h
