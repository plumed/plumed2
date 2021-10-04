/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"
#include "function/Sum.h"
#include "function/LessThan.h"
#include "function/MoreThan.h"
#include "function/Between.h"
#include "function/Custom.h"
#include "function/Combine.h"

namespace PLMD {
namespace matrixtools {

typedef FunctionOfMatrix<function::Sum> MatrixSum;
PLUMED_REGISTER_ACTION(MatrixSum,"SUM_MATRIX")
typedef FunctionOfMatrix<function::LessThan> MatrixLessThan;
PLUMED_REGISTER_ACTION(MatrixLessThan,"LESS_THAN_MATRIX")
typedef FunctionOfMatrix<function::MoreThan> MatrixMoreThan;
PLUMED_REGISTER_ACTION(MatrixMoreThan,"MORE_THAN_MATRIX")
typedef FunctionOfMatrix<function::Between> MatrixBetween;
PLUMED_REGISTER_ACTION(MatrixBetween,"BETWEEN_MATRIX")
typedef FunctionOfMatrix<function::Custom> MatrixCustom;
PLUMED_REGISTER_ACTION(MatrixCustom,"CUSTOM_MATRIX")
PLUMED_REGISTER_ACTION(MatrixCustom,"MATHEVAL_MATRIX")
typedef FunctionOfMatrix<function::Combine> MatrixCombine;
PLUMED_REGISTER_ACTION(MatrixCombine,"COMBINE_MATRIX")

}
}
