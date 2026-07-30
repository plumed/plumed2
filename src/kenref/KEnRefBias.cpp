/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023-2026 of Amr Alhossary and the KEnRef Authors.

The kenref module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The kenref module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------

This module is only compiled when PLUMED is configured with --enable-kenref,
in which case it links against the external KEnRef core library
(https://github.com/Smith-Group/KEnRef), which is distributed separately under
the BSD 3-Clause License. That licence is permissive and imposes no additional
restriction on this module or on PLUMED; see LICENSE_CORE.txt in the KEnRef
distribution for its text.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*
 * KEnRefBias.cpp — PLUMED kenref module (frozen frame)
 *
 * KEnRefBias is a PLUMED Bias/ActionAtomistic that restrains kinetic-ensemble observables (interproton
 * relaxation / order parameters) computed by the KEnRef library. It implements kenref::EngineAdapter
 * and delegates the whole per-step refinement pipeline (no-jump -> Kabsch fit -> replica gather ->
 * energy/derivative -> scatter -> inverse-fit -> force) to kenref::KEnRefDriver, implementing here only
 * the handful of PLUMED-specific callbacks (positions, box, force application, inter-replica MPI).
 *
 * This file holds the STABLE parts of the action (declaration is in KEnRefBias.h via the KEnRef repo):
 * registerKeywords, the PLUMED<->Eigen glue, calculate(), the lock/unlock boilerplate, and the action
 * registration. The one-time constructor (model + sub-indexing + driver setup) evolves with the KEnRef
 * model abstraction, so it is hosted in the KEnRef repository and compiled in via KEnRefBias_setup.cpp.
 */

// The whole module is optional: PLUMED can be configured with --enable-modules=all on a machine that
// has no KEnRef installed, in which case configure leaves __PLUMED_HAS_KENREF undefined and this
// translation unit compiles to nothing. Same pattern as src/pytorch and the ttsketch/ITensor module.
#ifdef __PLUMED_HAS_KENREF

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"

#include "core/ModelRegistry.h"   // registerKeywords iterates the model registry's schemas

// Unprefixed: kenref_plumed.pc puts the directory holding this header directly on the include path.
// A "plumedinterface/..." form would look to plumedcheck like a cross-module include, and
// plumedinterface is a KEnRef directory, not a PLUMED module.
#include "KEnRefBias.h"

#include <chrono>
#include <set>

// Written as two nested namespaces rather than `namespace PLMD::kenref` because plumedcheck matches
// each of them line-by-line ("^ *namespace *PLMD" and "^ *namespace *kenref"); the C++17 nested form
// satisfies only the first. This is also the form every other module uses.
namespace PLMD {
namespace kenref {

//+PLUMEDOC BIAS KENREF
/*
Restrains an ensemble of replicas against NMR observables (kinetic ensemble refinement).

KENREF computes an NMR observable as an average over the whole ensemble of replicas, compares it with
experimental data, and applies the resulting restraint force to the atoms. Unlike most biases it acts on
coordinates rather than on collective variables, so it takes **no** `ARG`.

The energy model is selected by name with `MODEL`:

| `MODEL`    | restrains                                          |
|:-----------|:---------------------------------------------------|
| `SIGMA`    | cross-relaxation rates (interproton NOE build-up)  |
| `PLATEAUS` | NOE plateau values                                 |
| `RELAX`    | longitudinal / transverse relaxation               |

Because the observable is an ensemble average, `RELAX` needs at least two replicas (run with `--multi`);
`SIGMA` and `PLATEAUS` also work with a single replica.

`EXP_DATA_FOLDER` holds the experimental data and the atom-pair definitions. `ATOMNAME_MAPPING` is a PDB
that maps atom names to indices, and is also used as the reference structure when `REF` is omitted. When
`FIT_TO_REFERENCE` is set the coordinates are Kabsch-fitted onto the reference, using `GUIDE_ATOMS`,
before the observable is computed, and the derivatives are rotated back afterwards. `SATURATE_FORCES`
clamps the resulting force to `MAX_FORCE`.

The action outputs its restraint energy as `.energy` (equal to `.bias`) and, as a diagnostic, the RMSD
from the reference as `.rmsd`.

This module needs the external KEnRef library and is not compiled by default; see the module
documentation for how to build it.

## Examples

Restraining a fragment of GB3 against cross-relaxation rates:

```plumed
#SETTINGS MOLFILE=regtest/kenref/rt-kenref-sigma/gb3_frag.pdb
kenref: KENREF ...
  MODEL=SIGMA
  K=1.0
  N=0.25
  PROTON_MHZ=700.0
  EXP_DATA_FOLDER=regtest/kenref/rt-kenref-sigma/
  REF=regtest/kenref/rt-kenref-sigma/gb3_frag.pdb
  ATOMNAME_MAPPING=regtest/kenref/rt-kenref-sigma/gb3_frag.pdb
  GUIDE_ATOMS=1,5,18,20,22,35,37,39,56,58,60,78,80,82,97
  MAX_FORCE=999
  FIT_TO_REFERENCE
  SATURATE_FORCES
...

PRINT ARG=kenref.bias,kenref.energy,kenref.rmsd FILE=colvar
```

*/
//+ENDPLUMEDOC

// cppcheck-suppress unknownMacro
// cppcheck cannot expand PLUMED_REGISTER_ACTION when checking this module, so it cannot parse the
// line below. The reason is specific to a default-off module: `make links` only creates a module's
// core/ and tools/ header symlinks when the module is ENABLED, and in the codecheck job kenref is
// off (its optional dependency is absent), so "core/ActionRegister.h" does not resolve -- cppcheck
// reports that missing include too, as (information). Since the macro is undefined, the
// registration is an unparseable construct and cppcheck flags it.
//
// A forward declaration of the class does NOT help, and was tested: the unknown entity is the
// MACRO, not the class name. Verified with cppcheck 2.17.1 by supplying a stub header that defines
// PLUMED_REGISTER_ACTION, whereupon the diagnostic disappears with no other change.
PLUMED_REGISTER_ACTION(KEnRefBias, "KENREF")

// ============================================================
//  registerKeywords
// ============================================================
void KEnRefBias::registerKeywords(Keywords &keys) {
  Bias::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.setDisplayName("KENREF");

  // KENREF is atom-based: it biases on coordinates, never on a CV, so it takes no ARG. PLUMED
  // 2.11-dev made ARG compulsory for every Bias, which forced input files to carry an empty
  // `ARG=` on master while 2.10 required it to be ABSENT -- the same input could not serve both.
  // Dropping the keyword makes one input file work everywhere. Guarded by exists() so this stays
  // correct on branches whose Bias never added ARG in the first place.
  if (keys.exists("ARG")) {
    keys.remove("ARG");
  }

  // ---- General (framework-owned) keywords ----
  keys.add("compulsory", "MODEL", "SIGMA", "The energy model to use (e.g. SIGMA or PLATEAUS)");
  keys.add("compulsory", "K", "1.0", "Force constant");
  keys.add("compulsory", "N", "0.25", "Power scaling factor");
  keys.add("optional", "MAX_FORCE", "Maximum force magnitude (default 9999)");
  keys.add("atoms", "GUIDE_ATOMS", "Atoms used for alignment to reference");
  keys.add("optional", "REF", "Reference structure PDB for alignment");
  keys.add("compulsory", "ATOMNAME_MAPPING", "undefined",
           "PDB file with atom-name -> atom-index mapping (also used as reference if REF is omitted)");
  keys.addFlag("FIT_TO_REFERENCE", false, "Fit coordinates to reference before calculating energy");
  keys.addFlag("SATURATE_FORCES", false, "Clamp forces to MAX_FORCE");

  // ---- Model-specific keywords (the union across every registered KEnRef model, optional; the
  //      selected model validates required-ness in buildCache). Deduplicated by key. ----
  ::kenref::bootstrapModels();
  std::set<std::string> seen;
  for (const auto &[name, schema]: ::kenref::ModelRegistry<KEnRef_Real_t>::allSchemas()) {
    for (const auto &spec: schema.specs()) {
      if (seen.insert(spec.key).second) {
        keys.add("optional", spec.key, spec.doc);
      }
    }
  }

  // ---- Output components ----
  keys.addOutputComponent("energy", "default", "scalar", "Total KEnRef restraint energy");
  keys.addOutputComponent("rmsd", "default", "scalar", "RMSD from reference (after fitting)");

  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

// ============================================================
//  Helpers: current PLUMED positions (nm) -> Angstrom Eigen
// ============================================================
CoordsMatrixType<KEnRef_Real_t> KEnRefBias::getGuideAtomsX() const {
  const std::vector<Vector> &pos = getPositions();
  const int nGuide = static_cast<int>(guideAtoms_.size());
  CoordsMatrixType<KEnRef_Real_t> gx(nGuide, 3);
  for (int i = 0; i < nGuide; ++i) {
    const int li = serial_to_localIdx_.at(guideAtoms_[i].serial());
    gx(i, 0) = to<KEnRef_Real_t>(pos[li][0]) * 10;
    gx(i, 1) = to<KEnRef_Real_t>(pos[li][1]) * 10;
    gx(i, 2) = to<KEnRef_Real_t>(pos[li][2]) * 10;
  }
  return gx;
}

void KEnRefBias::fillSubAtomsX(CoordsMatrixType<KEnRef_Real_t> &out) const {
  const std::vector<Vector> &pos = getPositions();
  for (int i = 0; i < out.rows(); ++i) {
    const int serial = sub0Id_to_global1Serial_[i];
    const int li = serial_to_localIdx_.at(serial);
    out(i, 0) = to<KEnRef_Real_t>(pos[li][0]) * 10; // nm -> Angstrom
    out(i, 1) = to<KEnRef_Real_t>(pos[li][1]) * 10;
    out(i, 2) = to<KEnRef_Real_t>(pos[li][2]) * 10;
  }
}

// ============================================================
//  kenref::EngineAdapter implementation (PLUMED-specific I/O + replica MPI)
// ============================================================
std::optional<std::string> KEnRefBias::getRawParam(const std::string &key) const {
  const auto it = modelParams_.find(key);
  if (it == modelParams_.end() || it->second.empty()) {
    return std::nullopt;
  }
  return it->second;
}

int KEnRefBias::numOmpThreads() const {
  return static_cast<int>(OpenMP::getNumThreads());
}

void KEnRefBias::getLocalModelX(int /*localModel*/, CoordsMatrixType<KEnRef_Real_t> &guideX,
                                CoordsMatrixType<KEnRef_Real_t> &subX,
                                Eigen::Matrix<KEnRef_Real_t, 3, 3> &box) const {
  guideX = getGuideAtomsX();
  subX.resize(numSubAtoms_, 3);
  fillSubAtomsX(subX);
  // PLUMED box Tensor (nm) -> Eigen 3x3 (raw; the driver scales it to Angstrom via toAngstrom).
  const Tensor &b = getPbc().getBox();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      box(i, j) = to<KEnRef_Real_t>(b[i][j]);
    }
}

void KEnRefBias::addLocalModelDerivatives(int /*localModel*/, const CoordsMatrixType<KEnRef_Real_t> &derivs) {
  // derivs are already inverse-fitted, unit-scaled and saturated by the driver. KEnRef returns
  // dE/dx; the applied force is -dE/dx, so we negate.
  for (int i = 0; i < derivs.rows(); ++i) {
    const int serial = sub0Id_to_global1Serial_[i];
    const int li = serial_to_localIdx_.at(serial);
    const Vector d(static_cast<double>(derivs(i, 0)),
                   static_cast<double>(derivs(i, 1)),
                   static_cast<double>(derivs(i, 2)));
    addForce(getValueIndices(getAbsoluteIndex(li)), -d);
  }
}

void KEnRefBias::gatherFittedSubAtomsX(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &localFitted,
                                       std::vector<CoordsMatrixType<KEnRef_Real_t>> &all) const {
  if (!isMultiSim_) {
    all = localFitted; // single replica owns the whole ensemble
    return;
  }
  const auto &fitted = localFitted[0];
  const int subSize = static_cast<int>(fitted.size()); // rows*3
  multi_sim_comm.Gather(const_cast<KEnRef_Real_t *>(fitted.data()), subSize,
                        allSimulationsSubAtomsX_buffer_.data(), subSize, 0);
  if (simulationIndex_ == 0) {
    const Eigen::Index nSub = fitted.rows();
    all.resize(numSimulations_);
    for (int m = 0; m < numSimulations_; ++m) { // zero-copy view materialised by the assignment
      all[m] = CoordsMapTypeConst<KEnRef_Real_t>(&allSimulationsSubAtomsX_buffer_[m * nSub * 3], nSub, 3);
    }
  }
}

void KEnRefBias::scatterModelDerivatives(const std::vector<CoordsMatrixType<KEnRef_Real_t>> &allPerModel,
    std::vector<CoordsMatrixType<KEnRef_Real_t>> &localPerModel) const {
  if (!isMultiSim_) {
    localPerModel = allPerModel; // single replica
    return;
  }
  const int subSize = numSubAtoms_ * 3;
  if (simulationIndex_ == 0)
    for (int m = 0; m < numSimulations_; ++m) {
      std::copy_n(allPerModel[m].data(), subSize, &allDerivatives_buffer_[m * subSize]);
    }
  multi_sim_comm.Scatter(allDerivatives_buffer_.data(), subSize, derivatives_buffer_.data(), subSize, 0);
  localPerModel.resize(1);
  localPerModel[0] = CoordsMapType<KEnRef_Real_t>(derivatives_buffer_.data(), numSubAtoms_, 3);
}

// ============================================================
//  calculate()
// ============================================================
void KEnRefBias::calculate() {
  const auto begin = std::chrono::high_resolution_clock::now();
  const long step = getStep();

  // Refresh per-step replica state (read by the EngineAdapter callbacks the driver invokes).
  isMultiSim_ = Communicator::plumedHasMPI() && Communicator::initialized() && multi_sim_comm.Get_size() > 1;
  numSimulations_ = isMultiSim_ ? multi_sim_comm.Get_size() : 1;
  simulationIndex_ = isMultiSim_ ? multi_sim_comm.Get_rank() : 0;

  // Run the whole engine-agnostic pipeline (no-jump -> fit -> gather -> compute -> scatter ->
  // inverse-fit -> unit-scale -> saturate -> addForce) through the adapter callbacks above.
  const KEnRef_Real_t energy = driver_->step(*this, step % 10 == 0);

  // ---- RMSD (crude diagnostic; current guide atoms vs reference, Angstrom -> nm) ----
  double rmsd = 0.0;
  if (guideAtomsReferenceCoords_.rows() > 0) {
    const CoordsMatrixType<KEnRef_Real_t> diff = getGuideAtomsX() - guideAtomsReferenceCoords_;
    rmsd = std::sqrt(diff.rowwise().squaredNorm().mean()) * 0.1;
  }

  // ---- output components ----
  const double bias_value = energy;
  getPntrToComponent("bias")->set(bias_value);
  getPntrToComponent("energy")->set(bias_value);
  getPntrToComponent("rmsd")->set(rmsd);
  setBias(bias_value);

  const auto end = std::chrono::high_resolution_clock::now();
  calculate_time_ += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
  if (simulationIndex_ == 0 && step % 10 == 0) {
    std::cout << "  Step " << step << "  KEnRef energy = " << energy << std::endl;
  }
}

// ============================================================
//  Ambiguity resolution (lock both Atoms and Arguments/CVs)
// ============================================================
void KEnRefBias::lockRequests() {
  ActionAtomistic::lockRequests();
  Bias::lockRequests();
}

void KEnRefBias::unlockRequests() {
  ActionAtomistic::unlockRequests();
  Bias::unlockRequests();
}

void KEnRefBias::calculateNumericalDerivatives(ActionWithValue *a) {
  ActionAtomistic::calculateNumericalDerivatives(a);
}

} // namespace kenref
} // namespace PLMD

#endif // __PLUMED_HAS_KENREF
