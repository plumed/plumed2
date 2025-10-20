/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2024 The METATOMIC-PLUMED team
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

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

//+PLUMEDOC COLVAR METATOMIC
/*
Use arbitrary machine learning models as collective variables.

!!! note

    This action requires `libmetatomic` to be enabled and found. Check the
    instructions in [the module page](module_metatomic.md) for more information
    on how to enable this module.

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

`METATOMIC` action also accepts the following options:

- `EXTENSIONS_DIRECTORY` should be the path to a directory containing
  TorchScript extensions (as shared libraries) that are required to load and
  execute the model. This matches the `collect_extensions` argument to
  `AtomisticModel.export` in Python.
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
    DEVICE=cuda
    CHECK_CONSISTENCY

    SPECIES1=1-10
    SPECIES2=11-20
    SPECIES_TO_TYPES=8,13

    # only run the calculation for the Aluminium (type 13) atoms,
    # still including all other atoms as potential neighbors.
    SELECTED_ATOMS=11-20
...
```

[TorchScript]: https://pytorch.org/docs/stable/jit.html [metatomic]:
https://docs.metatensor.org/metatomic/ [mta_tutorials]:
https://docs.metatensor.org/metatomic/latest/examples/ [features_output]:
https://docs.metatensor.org/metatomic/latest/outputs/features.html

*/
//+ENDPLUMEDOC

/*INDENT-OFF*/
#if !defined(__PLUMED_HAS_LIBMETATOMIC) || !defined(__PLUMED_HAS_LIBTORCH)

namespace PLMD { namespace metatomic {
class MetatomicPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatomicPlumedAction(const ActionOptions& options):
        Action(options),
        ActionAtomistic(options),
        ActionWithValue(options)
    {
        throw std::runtime_error(
            "Can not use metatomic action without the corresponding libraries. \n"
            "Make sure to configure with `--enable-libmetatomic --enable-libtorch` "
            "and that the corresponding libraries are found"
        );
    }

    void calculate() override {}
    void apply() override {}
    unsigned getNumberOfDerivatives() override {return 0;}
};

}} // namespace PLMD::metatomic

#else

#include <type_traits>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wimplicit-float-conversion"
#pragma GCC diagnostic ignored "-Wimplicit-int-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <torch/script.h>
#include <torch/version.h>
#include <torch/cuda.h>
#if TORCH_VERSION_MAJOR >= 2
#include <torch/mps.h>
#endif

#pragma GCC diagnostic pop

#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>

#include "vesin.h"


namespace PLMD {
namespace metatomic {

// We will cast Vector/Tensor to pointers to arrays and doubles, so let's make
// sure this is legal to do
static_assert(std::is_standard_layout<PLMD::Vector>::value);
static_assert(sizeof(PLMD::Vector) == sizeof(std::array<double, 3>));
static_assert(alignof(PLMD::Vector) == alignof(std::array<double, 3>));

static_assert(std::is_standard_layout<PLMD::Tensor>::value);
static_assert(sizeof(PLMD::Tensor) == sizeof(std::array<std::array<double, 3>, 3>));
static_assert(alignof(PLMD::Tensor) == alignof(std::array<std::array<double, 3>, 3>));

class MetatomicPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatomicPlumedAction(const ActionOptions&);

    void calculate() override;
    void apply() override;
    unsigned getNumberOfDerivatives() override;

private:
    // fill this->system_ according to the current PLUMED data
    void createSystem();
    // compute a neighbor list following metatomic format, using data from PLUMED
    metatensor_torch::TensorBlock computeNeighbors(
        metatomic_torch::NeighborListOptions request,
        const std::vector<PLMD::Vector>& positions,
        const PLMD::Tensor& cell,
        bool periodic
    );

    // execute the model for the given system
    metatensor_torch::TensorBlock executeModel(metatomic_torch::System system);

    torch::jit::Module model_;

    metatomic_torch::ModelCapabilities capabilities_;

    // neighbor lists requests made by the model
    std::vector<metatomic_torch::NeighborListOptions> nl_requests_;

    // dtype/device to use to execute the model
    torch::ScalarType dtype_;
    torch::Device device_;

    torch::Tensor atomic_types_;
    // store the strain to be able to compute the virial with autograd
    torch::Tensor strain_;

    metatomic_torch::System system_;
    metatomic_torch::ModelEvaluationOptions evaluations_options_;
    bool check_consistency_;

    metatensor_torch::TensorMap output_;
    // shape of the output of this model
    unsigned n_samples_;
    unsigned n_properties_;
};


MetatomicPlumedAction::MetatomicPlumedAction(const ActionOptions& options):
    Action(options),
    ActionAtomistic(options),
    ActionWithValue(options),
    device_(torch::kCPU)
{
    if (metatomic_torch::version().find("0.1.") != 0) {
        this->error(
            "this code requires version 0.1.x of metatomic-torch, got version " +
            metatomic_torch::version()
        );
    }

    // first, load the model
    std::string extensions_directory_str;
    this->parse("EXTENSIONS_DIRECTORY", extensions_directory_str);

    torch::optional<std::string> extensions_directory = torch::nullopt;
    if (!extensions_directory_str.empty()) {
        extensions_directory = std::move(extensions_directory_str);
    }

    std::string model_path;
    this->parse("MODEL", model_path);

    try {
        this->model_ = metatomic_torch::load_atomistic_model(model_path, extensions_directory);
    } catch (const std::exception& e) {
        this->error("failed to load model at '" + model_path + "': " + e.what());
    }

    // extract information from the model
    auto metadata = this->model_.run_method("metadata").toCustomClass<metatomic_torch::ModelMetadataHolder>();
    this->capabilities_ = this->model_.run_method("capabilities").toCustomClass<metatomic_torch::ModelCapabilitiesHolder>();
    auto requests_ivalue = this->model_.run_method("requested_neighbor_lists");
    for (auto request_ivalue: requests_ivalue.toList()) {
        auto request = request_ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
        this->nl_requests_.push_back(request);
    }

    log.printf("\n%s\n", metadata->print().c_str());
    // add the model references to PLUMED citation handling mechanism
    for (const auto& it: metadata->references) {
        for (const auto& ref: it.value()) {
            this->cite(ref);
        }
    }

    // parse the atomic types from the input file
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
                log.printf("%d ", atoms[j]);
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

    this->atomic_types_ = torch::tensor(std::move(atomic_types));

    this->check_consistency_ = false;
    this->parseFlag("CHECK_CONSISTENCY", this->check_consistency_);
    if (this->check_consistency_) {
        log.printf("  checking for internal consistency of the model\n");
    }

    // create evaluation options for the model. These won't change during the
    // simulation, so we initialize them once here.
    evaluations_options_ = torch::make_intrusive<metatomic_torch::ModelEvaluationOptionsHolder>();
    evaluations_options_->set_length_unit(getUnits().getLengthString());

    auto outputs = this->capabilities_->outputs();
    if (!outputs.contains("features")) {
        auto existing_outputs = std::vector<std::string>();
        for (const auto& it: this->capabilities_->outputs()) {
            existing_outputs.push_back(it.key());
        }

        this->error(
            "expected a 'features' output in the capabilities of the model, "
            "could not find it. the following outputs exist: " + torch::str(existing_outputs)
        );
    }

    auto output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
    // this output has no quantity or unit to set

    output->per_atom = this->capabilities_->outputs().at("features")->per_atom;
    // we are using torch autograd system to compute gradients,
    // so we don't need any explicit gradients.
    output->explicit_gradients = {};
    evaluations_options_->outputs.insert("features", output);

    // Determine which device we should use based on user input, what the model
    // supports and what's available
    auto available_devices = std::vector<torch::Device>();
    for (const auto& device: this->capabilities_->supported_devices) {
        if (device == "cpu") {
            available_devices.push_back(torch::kCPU);
        } else if (device == "cuda") {
            if (torch::cuda::is_available()) {
                available_devices.push_back(torch::Device("cuda"));
            }
        } else if (device == "mps") {
            #if TORCH_VERSION_MAJOR >= 2
            if (torch::mps::is_available()) {
                available_devices.push_back(torch::Device("mps"));
            }
            #endif
        } else {
            this->warning(
                "the model declared support for unknown device '" + device +
                "', it will be ignored"
            );
        }
    }

    if (available_devices.empty()) {
        this->error(
            "failed to find a valid device for the model at '" + model_path + "': "
            "the model supports " + torch::str(this->capabilities_->supported_devices) +
            ", none of these where available"
        );
    }

    std::string requested_device;
    this->parse("DEVICE", requested_device);
    if (requested_device.empty()) {
        // no user request, pick the device the model prefers
        this->device_ = available_devices[0];
    } else {
        bool found_requested_device = false;
        for (const auto& device: available_devices) {
            if (device.is_cpu() && requested_device == "cpu") {
                this->device_ = device;
                found_requested_device = true;
                break;
            } else if (device.is_cuda() && requested_device == "cuda") {
                this->device_ = device;
                found_requested_device = true;
                break;
            } else if (device.is_mps() && requested_device == "mps") {
                this->device_ = device;
                found_requested_device = true;
                break;
            }
        }

        if (!found_requested_device) {
            this->error(
                "failed to find requested device (" + requested_device + "): it is either "
                "not supported by this model or not available on this machine"
            );
        }
    }

    this->model_.to(this->device_);
    this->atomic_types_ = this->atomic_types_.to(this->device_);

    log.printf(
        "  running model on %s device with %s data\n",
        this->device_.str().c_str(),
        this->capabilities_->dtype().c_str()
    );

    if (this->capabilities_->dtype() == "float64") {
        this->dtype_ = torch::kFloat64;
    } else if (this->capabilities_->dtype() == "float32") {
        this->dtype_ = torch::kFloat32;
    } else {
        this->error(
            "the model requested an unsupported dtype '" + this->capabilities_->dtype() + "'"
        );
    }

    auto tensor_options = torch::TensorOptions().dtype(this->dtype_).device(this->device_);
    this->strain_ = torch::eye(3, tensor_options.requires_grad(true));

    // determine how many properties there will be in the output by running the
    // model once on a dummy system
    auto dummy_system = torch::make_intrusive<metatomic_torch::SystemHolder>(
        /*types = */ torch::zeros({0}, tensor_options.dtype(torch::kInt32)),
        /*positions = */ torch::zeros({0, 3}, tensor_options),
        /*cell = */ torch::zeros({3, 3}, tensor_options),
        /*pbc = */ torch::zeros({3}, tensor_options.dtype(torch::kBool))
    );

    log.printf("  the following neighbor lists have been requested:\n");
    auto length_unit = this->getUnits().getLengthString();
    auto model_length_unit = this->capabilities_->length_unit();
    for (auto request: this->nl_requests_) {
        log.printf("    - %s list, %g %s cutoff (requested %g %s)\n",
            request->full_list() ? "full" : "half",
            request->engine_cutoff(length_unit),
            length_unit.c_str(),
            request->cutoff(),
            model_length_unit.c_str()
        );

        auto neighbors = this->computeNeighbors(
            request,
            {PLMD::Vector(0, 0, 0)},
            PLMD::Tensor(0, 0, 0, 0, 0, 0, 0, 0, 0),
            false
        );
        metatomic_torch::register_autograd_neighbors(dummy_system, neighbors, this->check_consistency_);
        dummy_system->add_neighbor_list(request, neighbors);
    }

    this->n_properties_ = static_cast<unsigned>(
        this->executeModel(dummy_system)->properties()->count()
    );

    // parse and handle atom sub-selection. This is done AFTER determining the
    // output size, since the selection might not be valid for the dummy system
    std::vector<int32_t> selected_atoms;
    this->parseVector("SELECTED_ATOMS", selected_atoms);
    if (!selected_atoms.empty()) {
        auto selection_value = torch::zeros(
            {static_cast<int64_t>(selected_atoms.size()), 2},
            torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
        );

        for (unsigned i=0; i<selected_atoms.size(); i++) {
            auto n_atoms = static_cast<int32_t>(this->atomic_types_.size(0));
            if (selected_atoms[i] <= 0 || selected_atoms[i] > n_atoms) {
                this->error(
                    "Values in metatomic's SELECTED_ATOMS should be between 1 "
                    "and the number of atoms (" + std::to_string(n_atoms) + "), "
                    "got " + std::to_string(selected_atoms[i]));
            }
            // PLUMED input uses 1-based indexes, but metatomic wants 0-based
            selection_value[i][1] = selected_atoms[i] - 1;
        }

        evaluations_options_->set_selected_atoms(
            torch::make_intrusive<metatensor_torch::LabelsHolder>(
                std::vector<std::string>{"system", "atom"}, selection_value
            )
        );
    }

    // Now that we now both n_samples and n_properties, we can setup the
    // PLUMED-side storage for the computed CV
    if (output->per_atom) {
        if (selected_atoms.empty()) {
            this->n_samples_ = static_cast<unsigned>(this->atomic_types_.size(0));
        } else {
            this->n_samples_ = static_cast<unsigned>(selected_atoms.size());
        }
    } else {
        this->n_samples_ = 1;
    }

    if (n_samples_ == 1 && n_properties_ == 1) {
        log.printf("  the output of this model is a scalar\n");

        this->addValue();
    } else if (n_samples_ == 1) {
        log.printf("  the output of this model is 1x%d vector\n", n_properties_);

        this->addValue({this->n_properties_});
    } else if (n_properties_ == 1) {
        log.printf("  the output of this model is %dx1 vector\n", n_samples_);

        this->addValue({this->n_samples_});
    } else {
        log.printf("  the output of this model is a %dx%d matrix\n", n_samples_, n_properties_);

        this->addValue({this->n_samples_, this->n_properties_});
        this->getPntrToComponent(0)->reshapeMatrixStore(n_properties_);
    }

    this->setNotPeriodic();
}

unsigned MetatomicPlumedAction::getNumberOfDerivatives() {
    // gradients w.r.t. positions (3 x N values) + gradients w.r.t. strain (9 values)
    return 3 * this->getNumberOfAtoms() + 9;
}


void MetatomicPlumedAction::createSystem() {
    if (this->getTotAtoms() != static_cast<unsigned>(this->atomic_types_.size(0))) {
        std::ostringstream oss;
        oss << "METATOMIC action needs to know about all atoms in the system. ";
        oss << "There are " << this->getTotAtoms() << " atoms overall, ";
        oss << "but we only have atomic types for " << this->atomic_types_.size(0) << " of them.";
        plumed_merror(oss.str());
    }

    // this->getTotAtoms()

    const auto& cell = this->getPbc().getBox();

    auto cpu_f64_tensor = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
    auto torch_cell = torch::zeros({3, 3}, cpu_f64_tensor);

    torch_cell[0][0] = cell(0, 0);
    torch_cell[0][1] = cell(0, 1);
    torch_cell[0][2] = cell(0, 2);

    torch_cell[1][0] = cell(1, 0);
    torch_cell[1][1] = cell(1, 1);
    torch_cell[1][2] = cell(1, 2);

    torch_cell[2][0] = cell(2, 0);
    torch_cell[2][1] = cell(2, 1);
    torch_cell[2][2] = cell(2, 2);

    using torch::indexing::Slice;

    auto pbc_a = torch_cell.index({0, Slice()}).norm().abs().item<double>() > 1e-9;
    auto pbc_b = torch_cell.index({1, Slice()}).norm().abs().item<double>() > 1e-9;
    auto pbc_c = torch_cell.index({2, Slice()}).norm().abs().item<double>() > 1e-9;
    auto torch_pbc = torch::tensor({pbc_a, pbc_b, pbc_c});

    const auto& positions = this->getPositions();

    auto torch_positions = torch::from_blob(
        const_cast<PLMD::Vector*>(positions.data()),
        {static_cast<int64_t>(positions.size()), 3},
        cpu_f64_tensor
    );

    torch_positions = torch_positions.to(this->dtype_).to(this->device_);
    torch_cell = torch_cell.to(this->dtype_).to(this->device_);
    torch_pbc = torch_pbc.to(this->device_);

    // setup torch's automatic gradient tracking
    if (!this->doNotCalculateDerivatives()) {
        torch_positions.requires_grad_(true);

        // pretend to scale positions/cell by the strain so that it enters the
        // computational graph.
        torch_positions = torch_positions.matmul(this->strain_);
        torch_positions.retain_grad();

        torch_cell = torch_cell.matmul(this->strain_);
    }

    this->system_ = torch::make_intrusive<metatomic_torch::SystemHolder>(
        this->atomic_types_,
        torch_positions,
        torch_cell,
        torch_pbc
    );

    auto periodic = torch::all(torch_pbc).item<bool>();
    if (!periodic && torch::any(torch_pbc).item<bool>()) {
        std::string periodic_directions;
        std::string non_periodic_directions;
        if (pbc_a) {
            periodic_directions += "A";
        } else {
            non_periodic_directions += "A";
        }

        if (pbc_b) {
            periodic_directions += "B";
        } else {
            non_periodic_directions += "B";
        }

        if (pbc_c) {
            periodic_directions += "C";
        } else {
            non_periodic_directions += "C";
        }

        plumed_merror(
            "mixed periodic boundary conditions are not supported, this system "
            "is periodic along the " + periodic_directions + " cell vector(s), "
            "but not along the " + non_periodic_directions + " cell vector(s)."
        );
    }

    // compute the neighbors list requested by the model, and register them with
    // the system
    for (auto request: this->nl_requests_) {
        auto neighbors = this->computeNeighbors(request, positions, cell, periodic);
        metatomic_torch::register_autograd_neighbors(this->system_, neighbors, this->check_consistency_);
        this->system_->add_neighbor_list(request, neighbors);
    }
}


metatensor_torch::TensorBlock MetatomicPlumedAction::computeNeighbors(
    metatomic_torch::NeighborListOptions request,
    const std::vector<PLMD::Vector>& positions,
    const PLMD::Tensor& cell,
    bool periodic
) {
    auto labels_options = torch::TensorOptions().dtype(torch::kInt32).device(this->device_);
    auto neighbor_component = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "xyz",
        torch::tensor({0, 1, 2}, labels_options).reshape({3, 1})
    );
    auto neighbor_properties = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "distance", torch::zeros({1, 1}, labels_options)
    );

    auto cutoff = request->engine_cutoff(this->getUnits().getLengthString());

    // use https://github.com/Luthaf/vesin to compute the requested neighbor
    // lists since we can not get these from PLUMED
    vesin::VesinOptions options;
    options.cutoff = cutoff;
    options.full = request->full_list();
    options.return_shifts = true;
    options.return_distances = false;
    options.return_vectors = true;

    vesin::VesinNeighborList* vesin_neighbor_list = new vesin::VesinNeighborList();
    memset(vesin_neighbor_list, 0, sizeof(vesin::VesinNeighborList));

    const char* error_message = NULL;
    int status = vesin_neighbors(
        reinterpret_cast<const double (*)[3]>(positions.data()),
        positions.size(),
        reinterpret_cast<const double (*)[3]>(&cell(0, 0)),
        periodic,
        vesin::VesinCPU,
        options,
        vesin_neighbor_list,
        &error_message
    );

    if (status != EXIT_SUCCESS) {
        plumed_merror(
            "failed to compute neighbor list (cutoff=" + std::to_string(cutoff) +
            ", full=" + (request->full_list() ? "true" : "false") + "): " + error_message
        );
    }

    // transform from vesin to metatomic format
    auto n_pairs = static_cast<int64_t>(vesin_neighbor_list->length);

    auto pair_vectors = torch::from_blob(
        vesin_neighbor_list->vectors,
        {n_pairs, 3, 1},
        /*deleter*/ [=](void*) {
            vesin_free(vesin_neighbor_list);
            delete vesin_neighbor_list;
        },
        torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
    );

    auto pair_samples_values = torch::empty({n_pairs, 5}, labels_options.device(torch::kCPU));
    auto pair_samples_values_ptr = pair_samples_values.accessor<int32_t, 2>();
    for (unsigned i=0; i<n_pairs; i++) {
        pair_samples_values_ptr[i][0] = static_cast<int32_t>(vesin_neighbor_list->pairs[i][0]);
        pair_samples_values_ptr[i][1] = static_cast<int32_t>(vesin_neighbor_list->pairs[i][1]);
        pair_samples_values_ptr[i][2] = vesin_neighbor_list->shifts[i][0];
        pair_samples_values_ptr[i][3] = vesin_neighbor_list->shifts[i][1];
        pair_samples_values_ptr[i][4] = vesin_neighbor_list->shifts[i][2];
    }

    auto neighbor_samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
        pair_samples_values.to(this->device_)
    );

    auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
        pair_vectors.to(this->dtype_).to(this->device_),
        neighbor_samples,
        std::vector<metatensor_torch::Labels>{neighbor_component},
        neighbor_properties
    );

    return neighbors;
}

metatensor_torch::TensorBlock MetatomicPlumedAction::executeModel(metatomic_torch::System system) {
    try {
        auto ivalue_output = this->model_.forward({
            std::vector<metatomic_torch::System>{system},
            evaluations_options_,
            this->check_consistency_,
        });

        auto dict_output = ivalue_output.toGenericDict();
        auto cv = dict_output.at("features");
        this->output_ = cv.toCustomClass<metatensor_torch::TensorMapHolder>();
    } catch (const std::exception& e) {
        plumed_merror("failed to evaluate the model: " + std::string(e.what()));
    }

    plumed_massert(this->output_->keys()->count() == 1, "output should have a single block");
    auto block = metatensor_torch::TensorMapHolder::block_by_id(this->output_, 0);
    plumed_massert(block->components().empty(), "components are not yet supported in the output");

    return block;
}


void MetatomicPlumedAction::calculate() {
    this->createSystem();

    auto block = this->executeModel(this->system_);
    auto torch_values = block->values().to(torch::kCPU).to(torch::kFloat64);

    if (static_cast<unsigned>(torch_values.size(0)) != this->n_samples_) {
        plumed_merror(
            "expected the model to return a TensorBlock with " +
            std::to_string(this->n_samples_) + " samples, got " +
            std::to_string(torch_values.size(0)) + " instead"
        );
    } else if (static_cast<unsigned>(torch_values.size(1)) != this->n_properties_) {
        plumed_merror(
            "expected the model to return a TensorBlock with " +
            std::to_string(this->n_properties_) + " properties, got " +
            std::to_string(torch_values.size(1)) + " instead"
        );
    }

    Value* value = this->getPntrToComponent(0);
    // reshape the plumed `Value` to hold the data returned by the model
    if (n_samples_ == 1) {
        if (n_properties_ == 1) {
            value->set(torch_values.item<double>());
        } else {
            // we have multiple CV describing a single thing (atom or full system)
            for (unsigned i=0; i<n_properties_; i++) {
                value->set(i, torch_values[0][i].item<double>());
            }
        }
    } else {
        auto samples = block->samples();
        plumed_assert((samples->names() == std::vector<std::string>{"system", "atom"}));

        auto samples_values = samples->values().to(torch::kCPU);
        auto selected_atoms = this->evaluations_options_->get_selected_atoms();

        // handle the possibility that samples are returned in
        // a non-sorted order.
        auto get_output_location = [&](unsigned i) {
            if (selected_atoms.has_value()) {
                // If the users picked some selected atoms, then we store the
                // output in the same order as the selection was given
                auto sample = samples_values.index({static_cast<int64_t>(i), torch::indexing::Slice()});
                auto position = selected_atoms.value()->position(sample);
                plumed_assert(position.has_value());
                return static_cast<unsigned>(position.value());
            } else {
                return static_cast<unsigned>(samples_values[i][1].item<int32_t>());
            }
        };

        if (n_properties_ == 1) {
            // we have a single CV describing multiple things (i.e. atoms)
            for (unsigned i=0; i<n_samples_; i++) {
                auto output_i = get_output_location(i);
                value->set(output_i, torch_values[i][0].item<double>());
            }
        } else {
            // the CV is a matrix
            for (unsigned i=0; i<n_samples_; i++) {
                auto output_i = get_output_location(i);
                for (unsigned j=0; j<n_properties_; j++) {
                    value->set(output_i * n_properties_ + j, torch_values[i][j].item<double>());
                }
            }
        }
    }
}


void MetatomicPlumedAction::apply() {
    const auto* value = this->getPntrToComponent(0);
    if (!value->forcesWereAdded()) {
        return;
    }

    auto block = metatensor_torch::TensorMapHolder::block_by_id(this->output_, 0);
    auto torch_values = block->values().to(torch::kCPU).to(torch::kFloat64);

    if (!torch_values.requires_grad()) {
        this->warning(
            "the output of the model does not requires gradients, this might "
            "indicate a problem"
        );
        return;
    }

    auto output_grad = torch::zeros_like(torch_values);
    if (n_samples_ == 1) {
        if (n_properties_ == 1) {
            output_grad[0][0] = value->getForce();
        } else {
            for (unsigned i=0; i<n_properties_; i++) {
                output_grad[0][i] = value->getForce(i);
            }
        }
    } else {
        auto samples = block->samples();
        plumed_assert((samples->names() == std::vector<std::string>{"system", "atom"}));

        auto samples_values = samples->values().to(torch::kCPU);
        auto selected_atoms = this->evaluations_options_->get_selected_atoms();

        // see above for an explanation of why we use this function
        auto get_output_location = [&](unsigned i) {
            if (selected_atoms.has_value()) {
                auto sample = samples_values.index({static_cast<int64_t>(i), torch::indexing::Slice()});
                auto position = selected_atoms.value()->position(sample);
                plumed_assert(position.has_value());
                return static_cast<unsigned>(position.value());
            } else {
                return static_cast<unsigned>(samples_values[i][1].item<int32_t>());
            }
        };

        if (n_properties_ == 1) {
            for (unsigned i=0; i<n_samples_; i++) {
                auto output_i = get_output_location(i);
                output_grad[i][0] = value->getForce(output_i);
            }
        } else {
            for (unsigned i=0; i<n_samples_; i++) {
                auto output_i = get_output_location(i);
                for (unsigned j=0; j<n_properties_; j++) {
                    output_grad[i][j] = value->getForce(output_i * n_properties_ + j);
                }
            }
        }
    }

    this->system_->positions().mutable_grad() = torch::Tensor();
    this->strain_.mutable_grad() = torch::Tensor();

    torch_values.backward(output_grad);
    auto positions_grad = this->system_->positions().grad();
    auto strain_grad = this->strain_.grad();

    positions_grad = positions_grad.to(torch::kCPU).to(torch::kFloat64);
    strain_grad = strain_grad.to(torch::kCPU).to(torch::kFloat64);

    plumed_assert(positions_grad.sizes().size() == 2);
    plumed_assert(positions_grad.is_contiguous());

    plumed_assert(strain_grad.sizes().size() == 2);
    plumed_assert(strain_grad.is_contiguous());

    auto derivatives = std::vector<double>(
        positions_grad.data_ptr<double>(),
        positions_grad.data_ptr<double>() + 3 * this->system_->size()
    );

    // add virials to the derivatives
    derivatives.push_back(-strain_grad[0][0].item<double>());
    derivatives.push_back(-strain_grad[0][1].item<double>());
    derivatives.push_back(-strain_grad[0][2].item<double>());

    derivatives.push_back(-strain_grad[1][0].item<double>());
    derivatives.push_back(-strain_grad[1][1].item<double>());
    derivatives.push_back(-strain_grad[1][2].item<double>());

    derivatives.push_back(-strain_grad[2][0].item<double>());
    derivatives.push_back(-strain_grad[2][1].item<double>());
    derivatives.push_back(-strain_grad[2][2].item<double>());

    unsigned index = 0;
    this->setForcesOnAtoms(derivatives, index);
}

} // namespace metatomic
} // namespace PLMD


#endif


namespace PLMD {
namespace metatomic {

// use the same implementation for both the actual action and the dummy one
// (when libtorch and libmetatomic could not be found).
void MetatomicPlumedAction::registerKeywords(Keywords& keys) {
    Action::registerKeywords(keys);
    ActionAtomistic::registerKeywords(keys);
    ActionWithValue::registerKeywords(keys);

    keys.add("compulsory", "MODEL", "path to the exported metatomic model");
    keys.add("optional", "EXTENSIONS_DIRECTORY", "path to the directory containing TorchScript extensions to load");
    keys.add("optional", "DEVICE", "Torch device to use for the calculations");

    keys.addFlag("CHECK_CONSISTENCY", false, "should we enable internal consistency checks when executing the model");

    keys.add("numbered", "SPECIES", "the indices of atoms in each PLUMED species");
    keys.reset_style("SPECIES", "atoms");

    keys.add("optional", "SELECTED_ATOMS", "subset of atoms that should be used for the calculation");
    keys.reset_style("SELECTED_ATOMS", "atoms");

    keys.add("optional", "SPECIES_TO_TYPES", "mapping from PLUMED SPECIES to metatomic's atom types");

    keys.addOutputComponent("outputs", "default", "scalar/vector/matrix", "collective variable created by the metatomic model");
    keys.setValueDescription("scalar/vector/matrix", "collective variable created by the metatomic model");
}

PLUMED_REGISTER_ACTION(MetatomicPlumedAction, "METATOMIC")

} // namespace metatomic
} // namespace PLMD
