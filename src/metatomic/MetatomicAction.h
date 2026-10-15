/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2024-2026 The METATOMIC-PLUMED team
(see the PEOPLE-METATOMIC file at the root of this folder for a list of names)

See https://docs.metatensor.org/metatomic/ for more information about the
metatomic package that this module allows you to call from PLUMED.

This file is part of METATOMIC-PLUMED module.

The METATOMIC-PLUMED module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The METATOMIC-PLUMED module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the METATOMIC-PLUMED module. If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_metatomic_MetatomicAction_h
#define __PLUMED_metatomic_MetatomicAction_h

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"

/*INDENT-OFF*/
namespace PLMD {
namespace metatomic {

/// The METATOMIC action. It owns the PLUMED-side state (atoms, output value)
/// and forwards the model evaluation to one of two backends: TorchScript models
/// through metatomic-torch, or plugin models through metatomic-core.
class MetatomicPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatomicPlumedAction(const ActionOptions&);
    ~MetatomicPlumedAction() override;

    void calculate() override;
    void apply() override;
    unsigned getNumberOfDerivatives() override;

    class Backend {
    public:
        virtual ~Backend() = default;
        virtual void calculate() = 0;
        virtual void apply() = 0;
    };

private:
    friend class TorchBackend;
    friend class CoreBackend;

    /// Parse SPECIES/SPECIES_TO_TYPES, request all atoms in order, and return
    /// the metatomic type of each atom.
    std::vector<int32_t> parseAtomicTypes();
    /// Parse SELECTED_ATOMS, checking them against the number of atoms
    std::vector<AtomNumber> parseSelectedAtoms(size_t n_atoms);
    /// Current PLUMED box, with near-zero cell vectors set to zero and the
    /// corresponding directions marked as non-periodic.
    std::pair<Tensor, std::array<bool, 3>> cellAndPeriodicity();
    /// Register the PLUMED value holding `n_samples x n_properties` outputs
    void addOutputValue(unsigned n_samples, unsigned n_properties);

    /// Store row-major `n_samples x n_properties` values, with output sample
    /// `i` going to row `rows[i]` of the PLUMED value.
    void setOutputValues(const std::vector<unsigned>& rows, const std::vector<double>& values);
    /// Forces acting on the PLUMED value, in the same layout as `setOutputValues`
    std::vector<double> outputForces(const std::vector<unsigned>& rows) const;

    unsigned n_samples_ = 0;
    unsigned n_properties_ = 0;

    std::unique_ptr<Backend> backend_;
};

/// Create the TorchScript backend, or error if PLUMED was built without it
std::unique_ptr<MetatomicPlumedAction::Backend> makeTorchBackend(MetatomicPlumedAction& action);
/// Create the metatomic-core backend, or error if PLUMED was built without it
std::unique_ptr<MetatomicPlumedAction::Backend> makeCoreBackend(
    MetatomicPlumedAction& action,
    const std::vector<std::string>& plugin_libraries,
    const std::string& plugin
);

} // namespace metatomic
} // namespace PLMD

#endif
