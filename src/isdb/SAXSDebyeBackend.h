/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#ifndef __PLUMED_isdb_SAXSDebyeBackend_h
#define __PLUMED_isdb_SAXSDebyeBackend_h

#include "core/RegisterBase.h"
#include "tools/Vector.h"

#include <memory>
#include <string>
#include <vector>

namespace PLMD {
namespace isdb {

/// Interface of an accelerated (e.g. GPU) implementation of the Debye sum of SAXS/SANS.
///
/// The backends live in plugins (see plugins/openaccSAXS) and are selected
/// with the USEGPU keyword of SAXS/SANS. All the rest of the calculation
/// (form factors, self term, resolution function, metainference) stays in SAXS.
///
/// For n particles at positions r_i, with form factors f_i(q_k), a backend computes
/// \f[
///   I_k = \sum_{i} \sum_{j\neq i} f_i(q_k) f_j(q_k) \frac{\sin(q_k r_{ij})}{q_k r_{ij}}
/// \f]
/// (the self term \f$\sum_i f_i^2\f$ is added by SAXS) and its derivatives with
/// respect to the positions. The work is split over the ranks: each rank computes
/// the contributions of the rows i = first, first+stride, ... and SAXS sums the
/// results of all the ranks.
class SAXSDebyeBackend {
public:
  virtual ~SAXSDebyeBackend() = default;
  /// Computes the contributions of the rows i = first, first+stride, ... < pos.size().
  /// ff[i*q.size()+k] is the form factor of particle i at q[k].
  /// The partial intensities are added to sum (q.size() elements).
  /// If deriv is not null, the derivatives of I_k with respect to r_i are added to
  /// (*deriv)[k*pos.size()+i] (q.size()*pos.size() elements).
  /// The derivatives are kept for the next call to contract().
  virtual void compute(const std::vector<Vector>& pos,
                       const std::vector<double>& ff,
                       const std::vector<double>& q,
                       unsigned first, unsigned stride,
                       std::vector<double>& sum,
                       std::vector<Vector>* deriv) = 0;
  /// Adds to forces[i] (pos.size() elements) the derivatives of the last compute()
  /// contracted with the weights w (q.size() elements): sum_k w[k] dI_k/dr_i.
  virtual void contract(const std::vector<double>& w,
                        std::vector<Vector>& forces) = 0;
};

/// Content of the register of the backends
struct SAXSDebyeBackendPointers {
  /// creates a backend for a given device (-1: default device)
  typedef std::unique_ptr<SAXSDebyeBackend>(*creator_pointer)(int deviceid);
  creator_pointer create;
};

/// Register of the backends of SAXS/SANS.
class SAXSDebyeBackendRegister :
  public RegisterBase<SAXSDebyeBackendPointers> {
public:
  /// Names of the available backends
  std::vector<std::string> getBackendNames() const;
  /// Creates a backend. If name is empty and only one backend is available, that one is used.
  std::unique_ptr<SAXSDebyeBackend> create(const std::string& name, int deviceid) const;
};

/// The register of the backends of SAXS/SANS
SAXSDebyeBackendRegister& saxsDebyeBackendRegister();

/// Registers a backend while it is in scope: to be used in a plugin as
/// `static SAXSDebyeBackendRegistration<MyBackend> reg("name");`
/// where MyBackend has a constructor taking the device id.
template<typename Backend>
class SAXSDebyeBackendRegistration {
  SAXSDebyeBackendRegister::ID id;
  static std::unique_ptr<SAXSDebyeBackend> create(int deviceid) {
    return std::make_unique<Backend>(deviceid);
  }
public:
  explicit SAXSDebyeBackendRegistration(const std::string& name):
    id(saxsDebyeBackendRegister().add(name,SAXSDebyeBackendPointers{create})) {}
  ~SAXSDebyeBackendRegistration() {
    saxsDebyeBackendRegister().remove(id);
  }
  SAXSDebyeBackendRegistration(const SAXSDebyeBackendRegistration&) = delete;
  SAXSDebyeBackendRegistration& operator=(const SAXSDebyeBackendRegistration&) = delete;
};

}
}

#endif
