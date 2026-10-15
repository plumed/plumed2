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

#include "MetatomicAction.h"

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

//+PLUMEDOC COLVAR METATOMIC
/*
Use arbitrary machine learning models as collective variables.

!!! note

    This action requires `libmetatomic` (from metatomic-core) to be enabled and
    found. Check the instructions in [the module page](module_metatomic.md) for
    more information on how to enable this module.

This action enables the use of fully custom machine learning models — based on
the [metatomic] models interface — as collective variables in PLUMED. Such
machine learning model are typically written and customized using Python code,
and then exported to run within PLUMED as [TorchScript], which is a subset of
Python that can be executed by the C++ torch library.

Metatomic offers a way to define such models and pass data from PLUMED (or any
other simulation engine) to the model and back. For more information on how to
define such model, have a look at the [corresponding tutorials][mta_tutorials],
or at the code in `regtest/metatomic/`. Each of the Python scripts in this
directory defines a custom machine learning CV that can be used with PLUMED.

## Examples

The following input shows how you can call metatomic and evaluate the model that
is described in the file `custom_cv.pt` from PLUMED.

```plumed
metatomic_cv: METATOMIC ...
    MODEL=custom_cv.pt

    SPECIES1=1-26
    SPECIES2=27-62
    SPECIES3=63-76
    SPECIES_TO_TYPES=6,1,8
...
```

The numbered `SPECIES` labels are used to indicate the list of atoms that belong
to each atomic species in the system. The `SPECIES_TO_TYPE` keyword then
provides information on the atom type for each species. The first number here is
the atomic type of the atoms that have been specified using the `SPECIES1` flag,
the second number is the atomic number of the atoms that have been specified
using the `SPECIES2` flag and so on.

TorchScript models are executed through metatomic-core: PLUMED loads the Python
plugin of metatomic-core and wraps the model in `metatomic_torch.legacy.LegacyModel`,
so the Python environment found by the plugin needs `metatomic-torch` (and any
package providing TorchScript extensions used by the model). Gradients of
`system` features are computed with torch's automatic differentiation, while
`atom` features can be used as CVs, but can not be biased.

`METATOMIC` action also accepts the following options:

- `EXTENSIONS_DIRECTORY` should be the path to a directory containing
  TorchScript extensions (as shared libraries) that are required to load and
  execute the model. This matches the `collect_extensions` argument to
  `AtomisticModel.export` in Python.
- `DEVICE` can only be `cpu`, other devices are not supported yet;
- `CHECK_CONSISTENCY` can be used to enable internal consistency checks;
- `SELECTED_ATOMS` can be used to signal the metatomic models that it should
  only run its calculation for the selected subset of atoms. The model still
  need to know about all the atoms in the system (through the `SPECIES`
  keyword); but this can be used to reduce the calculation cost. Note that the
  indices of the selected atoms should start at 1 in the PLUMED input file, but
  they will be translated to start at 0 when given to the model (i.e. in
  Python/TorchScript, the `forward` method will receive a `selected_atoms` which
  starts at 0)

Here is another example with all the possible keywords:

```plumed
soap: METATOMIC ...
    MODEL=path/to/model.pt
    EXTENSIONS_DIRECTORY=path/to/extensions/
    DEVICE=cpu
    CHECK_CONSISTENCY

    SPECIES1=1-10
    SPECIES2=11-20
    SPECIES_TO_TYPES=8,13

    # only run the calculation for the Aluminium (type 13) atoms,
    # still including all other atoms as potential neighbors.
    SELECTED_ATOMS=11-20
...
```

## Plugin models

Models can also be loaded through any other plugin of metatomic-core. Unless
`MODEL` is a TorchScript file (ending in `.pt`) and `PLUGIN` is not given,
`MODEL` is forwarded to the plugins, and does not need to be a file:

- `PLUGIN_LIBRARY` is the path to a shared library containing a metatomic
  plugin, which is registered before loading the model. Use `PLUGIN_LIBRARY1`,
  `PLUGIN_LIBRARY2`, ... to register more than one plugin. If it is not given,
  the Python plugin installed next to the metatomic library used by PLUMED is
  registered;
- `PLUGIN` is the name of the plugin that should load `MODEL`. If it is not
  given, metatomic tries all registered plugins.

For PLUMED to apply forces, the model must provide `positions` and `strain`
gradients of its `feature` output explicitly. The pinned `metatomic.TorchModel`
does not compute gradients for the standard `feature` output, for either
per-system or per-atom samples. These gradients follow
the metatensor gradient layout and store derivatives of the feature itself.
`EXTENSIONS_DIRECTORY` can only be used with TorchScript models, and
metatomic-core's consistency checks do not accept explicit gradients on
`feature`, so `CHECK_CONSISTENCY` only works for models that are used without
biasing forces.

The simplest way to write such a model is a Python script that subclasses
`metatomic.NumpyModel` and defines a module-level `model` variable. Such scripts
are loaded with the Python plugin shipped in the `metatomic-core` Python package
(its path is given by `metatomic.utils.python_plugin_path`), which uses the
active Python environment. Have a look at `regtest/metatomic/rt-core/distance.py`
for a complete example.

```plumed
cv: METATOMIC ...
    PLUGIN=python
    MODEL=my_model.py

    SPECIES1=1-10
    SPECIES_TO_TYPES=6
...
```

Pair lists are given to such Python models as metatensor C++ arrays, and the
model script should register them with
`metatensor.register_external_data_wrapper("metatensor::SimpleDataArray",
metatensor.ExternalCpuArray)` before reading `system.pairs(...)`.

[TorchScript]: https://pytorch.org/docs/stable/jit.html [metatomic]:
https://docs.metatensor.org/metatomic/ [mta_tutorials]:
https://docs.metatensor.org/metatomic/latest/examples/ [features_output]:
https://docs.metatensor.org/metatomic/latest/outputs/features.html

*/
//+ENDPLUMEDOC

/*INDENT-OFF*/
namespace PLMD {
namespace metatomic {

MetatomicPlumedAction::MetatomicPlumedAction(const ActionOptions& options):
    Action(options),
    ActionAtomistic(options),
    ActionWithValue(options)
{
#ifndef __PLUMED_HAS_LIBMETATOMIC
    this->error(
        "Can not use metatomic action without the corresponding libraries. \n"
        "Make sure to configure with `--enable-libmetatomic` and that the "
        "metatomic and metatensor libraries are found"
    );
#else
    std::vector<std::string> plugin_libraries;
    std::string library;
    this->parse("PLUGIN_LIBRARY", library);
    if (!library.empty()) {
        plugin_libraries.push_back(library);
    }
    for (int i = 1; this->parseNumbered("PLUGIN_LIBRARY", i, library); i++) {
        plugin_libraries.push_back(library);
    }

    std::string plugin;
    this->parse("PLUGIN", plugin);

    backend_ = makeCoreBackend(*this, plugin_libraries, plugin);
#endif
}

MetatomicPlumedAction::~MetatomicPlumedAction() = default;

void MetatomicPlumedAction::calculate() {
    backend_->calculate();
}

void MetatomicPlumedAction::apply() {
    backend_->apply();
}

unsigned MetatomicPlumedAction::getNumberOfDerivatives() {
    // gradients w.r.t. positions (3 x N values) + gradients w.r.t. strain (9 values)
    return 3 * this->getNumberOfAtoms() + 9;
}

std::vector<int32_t> MetatomicPlumedAction::parseAtomicTypes() {
    std::vector<int32_t> atomic_types;
    std::vector<int32_t> species_to_types;
    this->parseVector("SPECIES_TO_TYPES", species_to_types);
    bool has_custom_types = !species_to_types.empty();

    std::vector<AtomNumber> all_atoms;
    this->parseAtomList("SPECIES", all_atoms);

    size_t n_species = 0;
    if (all_atoms.empty()) {
        // first parse each of the 'SPECIES' entry
        std::vector<std::vector<AtomNumber>> atoms_per_species;
        int i = 0;
        while (true) {
            i += 1;
            auto atoms = std::vector<AtomNumber>();
            this->parseAtomList("SPECIES", i, atoms);

            if (atoms.empty()) {
                break;
            }

            int32_t type = i;
            if (has_custom_types) {
                if (species_to_types.size() < static_cast<size_t>(i)) {
                    this->error(
                        "SPECIES_TO_TYPES is too small, it should have one entry "
                        "for each species (we have at least " + std::to_string(i) +
                        " species and " + std::to_string(species_to_types.size()) +
                        "entries in SPECIES_TO_TYPES)"
                    );
                }

                type = species_to_types[static_cast<size_t>(i - 1)];
            }

            log.printf("  atoms with type %d are: ", type);
            for(unsigned j=0; j<atoms.size(); j++) {
                log.printf("%u ", atoms[j].serial());
            }
            log.printf("\n");

            n_species += 1;
            atoms_per_species.emplace_back(std::move(atoms));
        }

        size_t n_atoms = 0;
        for (const auto& atoms: atoms_per_species) {
            n_atoms += atoms.size();
        }

        // then fill the atomic_types as required
        atomic_types.resize(n_atoms, 0);
        i = 0;
        for (const auto& atoms: atoms_per_species) {
            i += 1;

            int32_t type = i;
            if (has_custom_types) {
                type = species_to_types[static_cast<size_t>(i - 1)];
            }

            for (const auto& atom: atoms) {
                atomic_types[atom.index()] = type;
            }
        }
    } else {
        n_species = 1;

        int32_t type = 1;
        if (has_custom_types) {
            type = species_to_types[0];
        }
        atomic_types.resize(all_atoms.size(), type);
    }

    if (has_custom_types && species_to_types.size() != n_species) {
        this->warning(
            "SPECIES_TO_TYPES contains more entries (" +
            std::to_string(species_to_types.size()) +
            ") than there where species (" + std::to_string(n_species) + ")"
        );
    }

    // request atoms in order
    all_atoms.clear();
    for (size_t i=0; i<atomic_types.size(); i++) {
        all_atoms.push_back(AtomNumber::index(i));
    }
    this->requestAtoms(all_atoms);

    return atomic_types;
}

std::vector<AtomNumber> MetatomicPlumedAction::parseSelectedAtoms(size_t n_atoms) {
    std::vector<AtomNumber> selected_atoms;
    this->parseAtomList("SELECTED_ATOMS", selected_atoms);
    for (const auto& atom: selected_atoms) {
        if (atom.index() >= n_atoms) {
            this->error(
                "Values in metatomic's SELECTED_ATOMS should be between 1 "
                "and the number of atoms (" + std::to_string(n_atoms) + "), "
                "got " + std::to_string(atom.serial()));
        }
    }
    return selected_atoms;
}

std::pair<Tensor, std::array<bool, 3>> MetatomicPlumedAction::cellAndPeriodicity() {
    auto cell = this->getPbc().getBox();
    auto periodic = std::array<bool, 3>{true, true, true};

    // make sure the cell and pbc argument agree with each other
    for (unsigned i=0; i<3; i++) {
        auto norm = cell.getRow(i).modulo();
        if (norm < 1e-9) {
            if (norm > 1e-30) {
                this->warning(
                    "the cell vector " + std::string(1, "ABC"[i]) + " has a very small norm (" +
                    std::to_string(norm) + "), this direction will be treated as non periodic. "
                    "If this is intentional, ensure the vector is exactly zero to silence this warning"
                );
            }
            cell.setRow(i, Vector(0, 0, 0));
            periodic[i] = false;
        }
    }

    return {cell, periodic};
}

void MetatomicPlumedAction::addOutputValue(unsigned n_samples, unsigned n_properties) {
    n_samples_ = n_samples;
    n_properties_ = n_properties;

    if (n_samples_ == 1 && n_properties_ == 1) {
        log.printf("  the output of this model is a scalar\n");

        this->addValue();
    } else if (n_samples_ == 1) {
        log.printf("  the output of this model is 1x%d vector\n", n_properties_);

        this->addValue({n_properties_});
    } else if (n_properties_ == 1) {
        log.printf("  the output of this model is %dx1 vector\n", n_samples_);

        this->addValue({n_samples_});
    } else {
        log.printf("  the output of this model is a %dx%d matrix\n", n_samples_, n_properties_);

        this->addValue({n_samples_, n_properties_});
        this->getPntrToComponent(0)->reshapeMatrixStore(n_properties_);
    }

    this->setNotPeriodic();
}

void MetatomicPlumedAction::setOutputValues(const std::vector<unsigned>& rows, const std::vector<double>& values) {
    Value* value = this->getPntrToComponent(0);
    if (n_samples_ == 1 && n_properties_ == 1) {
        value->set(values[0]);
        return;
    }

    for (unsigned i=0; i<n_samples_; i++) {
        for (unsigned j=0; j<n_properties_; j++) {
            value->set(rows[i] * n_properties_ + j, values[i * n_properties_ + j]);
        }
    }
}

std::vector<double> MetatomicPlumedAction::outputForces(const std::vector<unsigned>& rows) const {
    const Value* value = this->getConstPntrToComponent(0);
    if (n_samples_ == 1 && n_properties_ == 1) {
        return {value->getForce()};
    }

    auto forces = std::vector<double>(
        static_cast<size_t>(n_samples_) * static_cast<size_t>(n_properties_)
    );
    for (unsigned i=0; i<n_samples_; i++) {
        for (unsigned j=0; j<n_properties_; j++) {
            forces[i * n_properties_ + j] = value->getForce(rows[i] * n_properties_ + j);
        }
    }
    return forces;
}

void MetatomicPlumedAction::registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
    ActionAtomistic::registerKeywords(keys);
    ActionWithValue::registerKeywords(keys);

    keys.add("compulsory", "MODEL", "path to the exported TorchScript model, or model passed to the metatomic plugin");
    keys.add("optional", "EXTENSIONS_DIRECTORY", "path to the directory containing TorchScript extensions to load");
    keys.add("optional", "DEVICE", "device to use for the calculations, only 'cpu' is supported");

    keys.add("numbered", "PLUGIN_LIBRARY", "shared library containing a metatomic-core plugin to register before loading the model");
    keys.add("optional", "PLUGIN", "name of the metatomic-core plugin used to load the model");

    keys.addFlag("CHECK_CONSISTENCY", false, "should we enable internal consistency checks when executing the model");

    keys.add("numbered", "SPECIES", "the indices of atoms in each PLUMED species");
    keys.reset_style("SPECIES", "atoms");

    keys.add("optional", "SELECTED_ATOMS", "subset of atoms that should be used for the calculation");
    keys.reset_style("SELECTED_ATOMS", "atoms");

    keys.add("optional", "SPECIES_TO_TYPES", "mapping from PLUMED SPECIES to metatomic's atom types");

    keys.add("optional", "VARIANT", "which variant of the 'features' output to pick");

    keys.addOutputComponent("outputs", "default", "scalar/vector/matrix", "collective variable created by the metatomic model");
    keys.setValueDescription("scalar/vector/matrix", "collective variable created by the metatomic model");
}

PLUMED_REGISTER_ACTION(MetatomicPlumedAction, "METATOMIC")

} // namespace metatomic
} // namespace PLMD
