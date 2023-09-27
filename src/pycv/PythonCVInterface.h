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
#include "PythonPlumedBase.h"

#include "colvar/Colvar.h"

namespace PLMD {
namespace pycv {

///TODO: manual "you have to specify ATOMS=something for default atoms"
///TODO: add interface to pbc
///TODO: the topology can be assumed fixed and done on the go at each run by loading the pdb in the python code
class PythonCVInterface : public Colvar,
  public PythonPlumedBase {
  static constexpr auto PYCV_NOTIMPLEMENTED="PYCV_NOTIMPLEMENTED";
  std::string import;
  std::string calculate_function;
  std::string prepare_function{PYCV_NOTIMPLEMENTED};

  std::vector<std::string> components;
  int ncomponents;
  int natoms;
  bool pbc=false;
  bool has_prepare{false};

  void check_dim(py::array_t<pycv_t>);
  void calculateSingleComponent(py::object &);
  void calculateMultiComponent(py::object &);

public:
  explicit PythonCVInterface(const ActionOptions&);
// active methods:
  void calculate() override;
  void prepare() override;
  static void registerKeywords( Keywords& keys );
};

} // namespace pycv
} // namespace PLMD
