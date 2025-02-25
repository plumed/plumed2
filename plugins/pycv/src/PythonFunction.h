/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 Daniele Rapetti

The pycv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The pycv module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_pycv_PythonFunction_h
#define __PLUMED_pycv_PythonFunction_h
#include "ActionWithPython.h"
#include "plumed/function/Function.h"
#include <pybind11/numpy.h>
namespace PLMD {

class NeighborList;

namespace pycv {
class PythonFunction :
  public function::Function,
  public ActionWithPython {
  static constexpr auto PYCV_DEFAULTINIT="plumedInit";
  static constexpr auto PYCV_DEFAULTCALCULATE="plumedCalculate";
  ::pybind11::module_ pyModule {};
  ::pybind11::object pyCalculate{};
  void calculateMultiComponent(pybind11::object &);
  void readReturn(const pybind11::object &, Value* );
public:
  explicit PythonFunction(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};
} // namespace pycv
} // namespace PLMD
#endif //__PLUMED_pycv_PythonFunction_h
