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
#ifndef __PLUMED_pycv_PythonCVInterface_h
#define __PLUMED_pycv_PythonCVInterface_h
#include "ActionWithPython.h"

#include "plumed/colvar/Colvar.h"

namespace PLMD {

class NeighborList;

namespace pycv {

///TODO: manual "you have to specify ATOMS=something for default atoms"
///TODO: add interface to pbc
///TODO: the topology can be assumed fixed and done on the go at each run by loading the pdb in the python code
class PythonCVInterface : public Colvar, public ActionWithPython {
  static constexpr auto PYCV_NOTIMPLEMENTED="PYCV_NOTIMPLEMENTED";
  static constexpr auto PYCV_DEFAULTINIT="plumedInit";
  static constexpr auto PYCV_DEFAULTCALCULATE="plumedCalculate";

  std::unique_ptr<NeighborList> nl{nullptr};

  ::pybind11::module_ pyModule {};
  ::pybind11::object pyCalculate{};
  ::pybind11::object pyPrepare;
  ::pybind11::object pyUpdate;

  bool pbc=false;
  bool hasPrepare = false;
  bool hasUpdate = false;
  bool invalidateList = true;
  bool firsttime = true;
  void calculateMultiComponent(pybind11::object &);
  void readReturn(const pybind11::object &, Value* );
public:
  ::pybind11::dict dataContainer {};
  explicit PythonCVInterface(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
// active methods:
  void calculate() override;
  void prepare() override;
  void update() override;

  NeighborList& getNL();
  ///redefinition of parseAtomList to avoid confict between plumed.dat and python options
  void pyParseAtomList(const char* key, const ::pybind11::dict &initDict, std::vector<AtomNumber> &);
};

} // namespace pycv
} // namespace PLMD
#endif //__PLUMED_pycv_PythonCVInterface_h
