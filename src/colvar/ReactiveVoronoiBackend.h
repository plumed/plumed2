/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 Reactive Voronoi contributors
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
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_colvar_ReactiveVoronoiBackend_h
#define __PLUMED_colvar_ReactiveVoronoiBackend_h

#include "core/RegisterBase.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PLMD {
class Pbc;
namespace colvar {

/// Optional backend for the soft assignment and its contracted derivatives.
class ReactiveVoronoiBackend {
public:
  /// Indices refer to the center and assigned selections of this context.
  struct InvalidPair : std::runtime_error {
    unsigned center;
    unsigned assigned;
    InvalidPair(unsigned c, unsigned a, const char* message):
      std::runtime_error(message), center(c), assigned(a) {}
  };
  virtual ~ReactiveVoronoiBackend() = default;
  virtual int device() const = 0;
  /// Queues the soft assignment (e.g. on a GPU) without waiting for it, see Action::launch().
  /// The positions can be reused as soon as this returns. finishForward() gives the occupancy.
  virtual void forward(const double* positions, const Pbc& box, bool pbc, double kappa) = 0;
  /// The same with a candidate list. The backend can reuse the uploaded candidate list for the same version.
  virtual void forwardNeighbors(const double* positions, const Pbc& box, bool pbc, double kappa,
                                double cutoff, const std::vector<int>& offsets,
                                const std::vector<int>& centers, unsigned long version) = 0;
  /// Waits for the last forward() or forwardNeighbors() and writes the occupancy of each center.
  /// Invalid geometry is reported here (InvalidPair) or by the forward call.
  virtual void finishForward(double* occupancy) = 0;
  virtual void backward(const double* sensitivity, double* derivatives, double* virial) = 0;
};

using ReactiveVoronoiBackendFactory = std::unique_ptr<ReactiveVoronoiBackend>(*)(int centers,
                                      int assigned, int device);
using ReactiveVoronoiBackendRegister = RegisterBase<ReactiveVoronoiBackendFactory>;

ReactiveVoronoiBackendRegister& reactiveVoronoiBackendRegister();

template<class Backend>
class ReactiveVoronoiBackendRegistration {
  static std::unique_ptr<ReactiveVoronoiBackend> create(int centers, int assigned, int device) {
    return std::make_unique<Backend>(centers,assigned,device);
  }
  ReactiveVoronoiBackendRegister::ID id;
public:
  explicit ReactiveVoronoiBackendRegistration(const std::string& name):
    id(reactiveVoronoiBackendRegister().add(name,create,reinterpret_cast<const void*>(&create))) {}
  ~ReactiveVoronoiBackendRegistration() {
    reactiveVoronoiBackendRegister().remove(id);
  }
  ReactiveVoronoiBackendRegistration(const ReactiveVoronoiBackendRegistration&) = delete;
  ReactiveVoronoiBackendRegistration& operator=(const ReactiveVoronoiBackendRegistration&) = delete;
};

}
}

#endif
