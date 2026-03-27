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

//+PLUMEDOC METATOMICMOD_COLVAR METATOMIC
/*
Use arbitrary machine learning models as collective variables.

\note This action requires the metatomic-torch library. Check the
instructions in the \ref METATOMICMOD page to enable this module.

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

\par Examples

The following input shows how you can call metatomic and evaluate the model that
is described in the file `custom_cv.pt` from PLUMED.

\plumedfile metatomic_cv: METATOMIC ... MODEL=custom_cv.pt

    SPECIES1=1-26
    SPECIES2=27-62
    SPECIES3=63-76
    SPECIES_TO_TYPES=6,1,8
...
\endplumedfile

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

\plumedfile soap: METATOMIC ... MODEL=soap.pt EXTENSION_DIRECTORY=extensions
CHECK_CONSISTENCY

    SPECIES1=1-10
    SPECIES2=11-20
    SPECIES_TO_TYPES=8,13

    # only run the calculation for the Aluminium (type 13) atoms, but
    # include the Oxygen (type 8) as potential neighbors.
    SELECTED_ATOMS=11-20
...
\endplumedfile

\par Collective variables and metatomic  models

PLUMED can use the [`"features"` output][features_output] of metatomic models as
a collective variables.

*/ /*

[TorchScript]: https://pytorch.org/docs/stable/jit.html
[metatomic]: https://docs.metatensor.org/metatomic/
[mta_tutorials]: https://docs.metatensor.org/metatomic/latest/examples/
[features_output]: https://docs.metatensor.org/metatomic/latest/outputs/features.html
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

/// Small helper class to compute one neighbor list requested by the metatomc
/// model using vesin
class NeighborListCalculator {
public:
    NeighborListCalculator(
        metatomic_torch::NeighborListOptions options,
        const std::string& engine_length_unit
    );
    ~NeighborListCalculator();

    NeighborListCalculator(const NeighborListCalculator& other) = delete;
    NeighborListCalculator& operator=(const NeighborListCalculator& other) = delete;

    NeighborListCalculator(NeighborListCalculator&& other) noexcept;
    NeighborListCalculator& operator=(NeighborListCalculator&& other) noexcept;

    // compute the neighbor list following metatomic format, using data from PLUMED
    metatensor_torch::TensorBlock compute(
        const std::vector<PLMD::Vector>& positions,
        const PLMD::Tensor& cell,
        std::array<bool, 3> periodic,
        torch::ScalarType dtype,
        torch::Device device
    );

    metatomic_torch::NeighborListOptions options;
private:
    double engine_cutoff_;
    vesin::VesinNeighborList neighbors_;
};

NeighborListCalculator::NeighborListCalculator(
    metatomic_torch::NeighborListOptions options_,
    const std::string& engine_length_unit
):
    options(options_),
    engine_cutoff_(options_->engine_cutoff(engine_length_unit))
{
    memset(&this->neighbors_, 0, sizeof(vesin::VesinNeighborList));
}

NeighborListCalculator::NeighborListCalculator(NeighborListCalculator&& other) noexcept {
    this->options = other.options;
    this->engine_cutoff_ = other.engine_cutoff_;
    this->neighbors_ = other.neighbors_;

    memset(&other.neighbors_, 0, sizeof(vesin::VesinNeighborList));
}

NeighborListCalculator& NeighborListCalculator::operator=(NeighborListCalculator&& other) noexcept {
    if (this != &other) {
        vesin::vesin_free(&this->neighbors_);

        this->options = other.options;
        this->engine_cutoff_ = other.engine_cutoff_;
        this->neighbors_ = other.neighbors_;

        memset(&other.neighbors_, 0, sizeof(vesin::VesinNeighborList));
    }
    return *this;
}

NeighborListCalculator::~NeighborListCalculator() {
    vesin::vesin_free(&this->neighbors_);
}

metatensor_torch::TensorBlock NeighborListCalculator::compute(
    const std::vector<PLMD::Vector>& positions,
    const PLMD::Tensor& cell,
    std::array<bool, 3> periodic,
    torch::ScalarType dtype,
    torch::Device device
) {
    auto labels_options = torch::TensorOptions().dtype(torch::kInt32).device(device);
    auto neighbor_component = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "xyz",
        torch::tensor({0, 1, 2}, labels_options).reshape({3, 1})
    );
    auto neighbor_properties = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "distance", torch::zeros({1, 1}, labels_options)
    );

    // use https://github.com/Luthaf/vesin to compute the requested neighbor
    // lists since we can not get these from PLUMED
    vesin::VesinOptions vesin_options;
    vesin_options.cutoff = this->engine_cutoff_;
    vesin_options.full = this->options->full_list();
    vesin_options.return_shifts = true;
    vesin_options.return_distances = false;
    vesin_options.return_vectors = true;
    vesin_options.algorithm = vesin::VesinAutoAlgorithm;

    const char* error_message = nullptr;
    int status = vesin_neighbors(
        reinterpret_cast<const double (*)[3]>(positions.data()),
        positions.size(),
        reinterpret_cast<const double (*)[3]>(&cell(0, 0)),
        periodic.data(),
        {vesin::VesinCPU, 0},
        vesin_options,
        &this->neighbors_,
        &error_message
    );

    if (status != EXIT_SUCCESS) {
        plumed_merror(
            "failed to compute neighbor list (cutoff=" + std::to_string(this->engine_cutoff_) +
            ", full=" + (this->options->full_list() ? "true" : "false") + "): " + error_message
        );
    }

    // transform from vesin to metatomic format
    auto n_pairs = static_cast<int64_t>(this->neighbors_.length);

    auto pair_vectors = torch::from_blob(
        this->neighbors_.vectors,
        {n_pairs, 3, 1},
        torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU)
    );

    auto pair_samples_values = torch::empty({n_pairs, 5}, labels_options.device(torch::kCPU));
    auto pair_samples_values_ptr = pair_samples_values.accessor<int32_t, 2>();
    for (unsigned i=0; i<n_pairs; i++) {
        pair_samples_values_ptr[i][0] = static_cast<int32_t>(this->neighbors_.pairs[i][0]);
        pair_samples_values_ptr[i][1] = static_cast<int32_t>(this->neighbors_.pairs[i][1]);
        pair_samples_values_ptr[i][2] = this->neighbors_.shifts[i][0];
        pair_samples_values_ptr[i][3] = this->neighbors_.shifts[i][1];
        pair_samples_values_ptr[i][4] = this->neighbors_.shifts[i][2];
    }

    auto neighbor_samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
        pair_samples_values.to(device),
        // vesin should create unique pairs
        metatensor::assume_unique{}
    );

    auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
        pair_vectors.to(dtype).to(device),
        neighbor_samples,
        std::vector<metatensor_torch::Labels>{neighbor_component},
        neighbor_properties
    );

    return neighbors;
}

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


    // execute the model for the given system
    metatensor_torch::TensorBlock executeModel(metatomic_torch::System system);

    metatensor_torch::Module model_;

    metatomic_torch::ModelCapabilities capabilities_;
    // name of the output we request
    std::string features_key;

    // neighbor lists requests made by the model and the corresponding data
    std::vector<NeighborListCalculator> neighbor_lists_;

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
    model_(torch::jit::Module()),
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
    auto nl_requests_ivalue = this->model_.run_method("requested_neighbor_lists");
    for (auto nl_request_ivalue: nl_requests_ivalue.toList()) {
        auto nl_request = nl_request_ivalue.get().toCustomClass<metatomic_torch::NeighborListOptionsHolder>();
        this->neighbor_lists_.emplace_back(nl_request, this->getUnits().getLengthString());
    }

    auto extra_inputs = this->model_.run_method("requested_inputs").toGenericDict();
    auto standard_inputs = std::vector<std::string>{};
    auto custom_inputs = std::vector<std::string>{};
    for (const auto& item: extra_inputs) {
        auto key = item.key().toStringRef();
        if (key.find("::") != std::string::npos) {
            custom_inputs.push_back(key);
        } else {
            standard_inputs.push_back(key);
        }
    }

    if (!standard_inputs.empty()) {
        this->error(
            "The model requested extra inputs that are not yet supported in PLUMED. "
            "Please open an issue to request support for the following inputs: " +
            torch::str(standard_inputs)
        );
    }

    if (!custom_inputs.empty()) {
        this->error(
            "The model requested custom inputs (" + torch::str(custom_inputs) + ") "
            "that can not be provided by PLUMED. Please change your model to use "
            "standard inputs only."
        );
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

    this->atomic_types_ = torch::tensor(atomic_types);

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

    std::string requested_variant;
    torch::optional<std::string> requested_variant_opt = torch::nullopt;
    this->parse("VARIANT", requested_variant);
    if (!requested_variant.empty()) {
        requested_variant_opt = requested_variant;
    }

    this->features_key = metatomic_torch::pick_output("features", outputs, requested_variant_opt);

    auto output = torch::make_intrusive<metatomic_torch::ModelOutputHolder>();
    // this output has no quantity or unit to set

    output->per_atom = this->capabilities_->outputs().at(this->features_key)->per_atom;
    // we are using torch autograd system to compute gradients,
    // so we don't need any explicit gradients.
    output->explicit_gradients = {};
    evaluations_options_->outputs.insert(this->features_key, output);

    std::string requested_device;
    torch::optional<std::string> requested_device_opt = torch::nullopt;
    this->parse("DEVICE", requested_device);
    if (!requested_device.empty()) {
        requested_device_opt = requested_device;
    }

    this->device_ = torch::Device(
        metatomic_torch::pick_device(this->capabilities_->supported_devices, requested_device_opt),
        /*index=*/ 0
    );

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
    for (auto& nl: this->neighbor_lists_) {
        log.printf("    - %s list, %g %s cutoff (requested %g %s)\n",
            nl.options->full_list() ? "full" : "half",
            nl.options->engine_cutoff(length_unit),
            length_unit.c_str(),
            nl.options->cutoff(),
            model_length_unit.c_str()
        );

        auto neighbors = nl.compute(
            {PLMD::Vector(0, 0, 0)},
            PLMD::Tensor(0, 0, 0, 0, 0, 0, 0, 0, 0),
            {false, false, false},
            this->dtype_,
            this->device_
        );
        metatomic_torch::register_autograd_neighbors(dummy_system, neighbors, this->check_consistency_);
        dummy_system->add_neighbor_list(nl.options, neighbors);
    }

    this->n_properties_ = static_cast<unsigned>(
        this->executeModel(dummy_system)->properties()->count()
    );

    // parse and handle atom sub-selection. This is done AFTER determining the
    // output size, since the selection might not be valid for the dummy system
    std::vector<AtomNumber> selected_atoms;
    this->parseAtomList("SELECTED_ATOMS", selected_atoms);
    if (!selected_atoms.empty()) {
        auto selection_value = torch::zeros(
            {static_cast<int64_t>(selected_atoms.size()), 2},
            torch::TensorOptions().dtype(torch::kInt32).device(this->device_)
        );

        for (unsigned i=0; i<selected_atoms.size(); i++) {
            auto n_atoms = this->atomic_types_.size(0);
            if (selected_atoms[i].index() > n_atoms) {
                this->error(
                    "Values in metatomic's SELECTED_ATOMS should be between 1 "
                    "and the number of atoms (" + std::to_string(n_atoms) + "), "
                    "got " + std::to_string(selected_atoms[i].serial()));
            }
            selection_value[i][1] = static_cast<int32_t>(selected_atoms[i].index());
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
        this->getPntrToComponent(0)->buildDataStore();
    } else if (n_properties_ == 1) {
        log.printf("  the output of this model is %dx1 vector\n", n_samples_);

        this->addValue({this->n_samples_});
        this->getPntrToComponent(0)->buildDataStore();
    } else {
        log.printf("  the output of this model is a %dx%d matrix\n", n_samples_, n_properties_);

        this->addValue({this->n_samples_, this->n_properties_});
        this->getPntrToComponent(0)->buildDataStore();
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

    auto norm_a = torch_cell.index({0, Slice()}).norm().abs().item<double>();
    auto norm_b = torch_cell.index({1, Slice()}).norm().abs().item<double>();
    auto norm_c = torch_cell.index({2, Slice()}).norm().abs().item<double>();

    auto periodic = std::array<bool, 3>{true, true, true};

    // make sure the cell and pbc argument agree with each other
    if (norm_a < 1e-9) {
        if (norm_a > 1e-30) {
            this->warning(
                "the cell vector A has a very small norm (" + std::to_string(norm_a) + "), "
                "this direction will be treated as non periodic. If this is intentional, ensure "
                "the vector is exactly zero to silence this warning"
            );
        }
        torch_cell.index({0, Slice()}).fill_(0);
        periodic[0] = false;
    }

    if (norm_b < 1e-9) {
        if (norm_b > 1e-30) {
            this->warning(
                "the cell vector B has a very small norm (" + std::to_string(norm_b) + "), "
                "this direction will be treated as non periodic. If this is intentional, ensure "
                "the vector is exactly zero to silence this warning"
            );
        }
        torch_cell.index({1, Slice()}).fill_(0);
        periodic[1] = false;
    }

    if (norm_c < 1e-9) {
        if (norm_c > 1e-30) {
            this->warning(
                "the cell vector C has a very small norm (" + std::to_string(norm_c) + "), "
                "this direction will be treated as non periodic. If this is intentional, ensure "
                "the vector is exactly zero to silence this warning"
            );
        }
        torch_cell.index({2, Slice()}).fill_(0);
        periodic[2] = false;
    }

    auto torch_pbc = torch::zeros({3}, torch::TensorOptions().dtype(torch::kBool).device(torch::kCPU));
    for (unsigned i=0; i<3; i++) {
        torch_pbc[i] = periodic[i];
    }

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

    // compute the neighbors list requested by the model, and register them with
    // the system
    for (auto& nl: this->neighbor_lists_) {
        auto neighbors = nl.compute(positions, cell, periodic, this->dtype_, this->device_);
        metatomic_torch::register_autograd_neighbors(this->system_, neighbors, this->check_consistency_);
        this->system_->add_neighbor_list(nl.options, neighbors);
    }
}


metatensor_torch::TensorBlock MetatomicPlumedAction::executeModel(metatomic_torch::System system) {
    try {
        auto ivalue_output = this->model_.forward({
            std::vector<metatomic_torch::System>{system},
            evaluations_options_,
            this->check_consistency_,
        });

        auto dict_output = ivalue_output.toGenericDict();
        auto cv = dict_output.at(this->features_key);
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

    keys.add("optional", "VARIANT", "which variant of the 'features' output to pick");

    keys.addOutputComponent("outputs", "default", "collective variable created by the metatomic model");

    keys.setValueDescription("collective variable created by the metatomic model");
}

PLUMED_REGISTER_ACTION(MetatomicPlumedAction, "METATOMIC")

} // namespace metatomic
} // namespace PLMD
