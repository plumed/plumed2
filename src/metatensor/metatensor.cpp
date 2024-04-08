/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Copyright (c) 2024 Guillaume Fraux

    This module is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"


#if !defined(__PLUMED_HAS_LIBTORCH) || !defined(__PLUMED_HAS_METATENSOR)

namespace PLMD { namespace metatensor {
class MetatensorPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatensorPlumedAction(const ActionOptions&) {
        throw std::runtime_error(
            "Can not use metatensor action without the corresponding libraries. \n"
            "Make sure to configure with `--enable-metatensor --enable-libtorch` "
            "and that the corresponding libraries are found"
        );
    }
};

}} // namespace PLMD::metatensor

#else

#include <type_traits>

#include <torch/script.h>
#include "torch/csrc/autograd/autograd.h"

#include <metatensor/torch.hpp>
#include <metatensor/torch/atomistic.hpp>

#include "vesin.h"


// TEMPORARY HACK
#include <dlfcn.h>
// TEMPORARY HACK

namespace PLMD { namespace metatensor {

// We will cast Vector/Tensor to pointers to arrays and doubles, so let's make
// sure this is legal to do
static_assert(std::is_standard_layout<PLMD::Vector>::value);
static_assert(sizeof(PLMD::Vector) == sizeof(std::array<double, 3>));
static_assert(alignof(PLMD::Vector) == alignof(std::array<double, 3>));

static_assert(std::is_standard_layout<PLMD::Tensor>::value);
static_assert(sizeof(PLMD::Tensor) == sizeof(std::array<std::array<double, 3>, 3>));
static_assert(alignof(PLMD::Tensor) == alignof(std::array<std::array<double, 3>, 3>));

class MetatensorPlumedAction: public ActionAtomistic, public ActionWithValue {
public:
    static void registerKeywords(Keywords& keys);
    explicit MetatensorPlumedAction(const ActionOptions&);

    void calculate() override;
    void apply() override;
    unsigned getNumberOfDerivatives() override;

private:
    void createSystem();

    // compute a neighbor list following metatensor format, using data from PLUMED
    metatensor_torch::TorchTensorBlock computeNeighbors(
        metatensor_torch::NeighborsListOptions request,
        const std::vector<PLMD::Vector>& positions,
        const PLMD::Tensor& cell
    );

    torch::jit::Module model_;

    torch::Tensor atomic_types_;
    // store the strain to be able to compute the virial with autograd
    torch::Tensor strain_;

    metatensor_torch::System system_;
    metatensor_torch::ModelEvaluationOptions evaluations_options_;
    bool check_consistency_ = true;
    metatensor_torch::TorchTensorMap output_;
};


MetatensorPlumedAction::MetatensorPlumedAction(const ActionOptions& options):
    Action(options),
    ActionAtomistic(options),
    ActionWithValue(options)
{
    std::string extensions_directory;
    this->parse("EXTENSIONS_DIRECTORY", extensions_directory);

    // TEMPORARY BAD CODE, TO BE REMOVED
    dlopen(
        (extensions_directory + "/rascaline/lib/librascaline.dylib").c_str(),
        RTLD_LOCAL | RTLD_NOW
    );

    dlopen(
        (extensions_directory + "/rascaline/torch/lib/librascaline_torch.dylib").c_str(),
        RTLD_LOCAL | RTLD_NOW
    );
    // END OF TEMPORARY BAD CODE, TO BE REMOVED

    // load the model
    std::string model_path;
    this->parse("MODEL", model_path);

    try {
        this->model_ = metatensor_torch::load_atomistic_model(model_path);
    } catch (const std::exception& e) {
        error("failed to load model at '" + model_path + "': " + e.what());
    }


    // parse the atomic types from the input file
    std::vector<int32_t> atomic_types;
    std::vector<int32_t> species_to_metatensor_types;
    parseVector("SPECIES_TO_METATENSOR_TYPES", species_to_metatensor_types);
    bool has_custom_types = !species_to_metatensor_types.empty();

    std::vector<AtomNumber> all_atoms;
    parseAtomList("SPECIES", all_atoms);

    auto n_species = 0;
    if (all_atoms.empty()) {
        std::vector<AtomNumber> t;
        for (int i=1;;i++) {
            parseAtomList("SPECIES", i, t);
            if (t.empty()) {
                break;
            }

            int32_t type = i;
            if (has_custom_types) {
                if (species_to_metatensor_types.size() < i) {
                    error(
                        "SPECIES_TO_METATENSOR_TYPES is too small, "
                        "it should have one entry for each species (we have at least "
                        + std::to_string(i) + " species and " +
                        std::to_string(species_to_metatensor_types.size()) +
                        "entries in SPECIES_TO_METATENSOR_TYPES)"
                    );
                }

                type = species_to_metatensor_types[i - 1];
            }

            log.printf("  Species %d includes atoms : ", i);
            for(unsigned j=0; j<t.size(); j++) {
                all_atoms.push_back(t[j]);
                atomic_types.push_back(type);
            }
            log.printf("\n"); t.resize(0);

            n_species += 1;
        }
    } else {
        n_species = 1;

        int32_t type = 1;
        if (has_custom_types) {
            type = species_to_metatensor_types[0];
        }
        atomic_types.resize(all_atoms.size(), type);
    }

    if (has_custom_types && species_to_metatensor_types.size() != n_species) {
        this->warning(
            "SPECIES_TO_METATENSOR_TYPES contains more entries (" +
            std::to_string(species_to_metatensor_types.size()) +
            ") than there where species (" + std::to_string(n_species) + ")"
        );
    }

    this->atomic_types_ = torch::tensor(std::move(atomic_types));

    // Request the atoms and check we have read in everything
    requestAtoms(all_atoms);

    // TODO: selected_atoms
    // evaluations_options_->set_selected_atoms()

    // setup the output
    // TODO: define the size/type of output a bit better
    this->addValue({1, 1});
    this->setNotPeriodic();
    this->getPntrToComponent(0)->buildDataStore();

    // create evaluation options for the model. These won't change during the
    // simulation, so we initialize them once here.
    evaluations_options_ = torch::make_intrusive<metatensor_torch::ModelEvaluationOptionsHolder>();
    evaluations_options_->set_length_unit(getUnits().getLengthString());

    auto output = torch::make_intrusive<metatensor_torch::ModelOutputHolder>();
    // TODO: should this be configurable?
    output->per_atom = true;
    // we are using torch autograd system to compute gradients, so we don't need
    // any explicit gradients.
    output->explicit_gradients = {};
    evaluations_options_->outputs.insert("collective_variable", output);
}

unsigned MetatensorPlumedAction::getNumberOfDerivatives() {
    // gradients w.r.t. positions (3 x N values) + gradients w.r.t. strain (9 values)
    return 3 * this->getNumberOfAtoms() + 9;
}


void MetatensorPlumedAction::createSystem() {
    const auto& cell = this->getPbc().getBox();

    auto tensor_options = torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
    auto torch_cell = torch::zeros({3, 3}, tensor_options);

    // TODO: check if cell is stored in row or column major order
    // TODO: check if cell is zero for non-periodic systems
    torch_cell[0][0] = cell(0, 0);
    torch_cell[0][1] = cell(0, 1);
    torch_cell[0][2] = cell(0, 2);

    torch_cell[1][0] = cell(1, 0);
    torch_cell[1][1] = cell(1, 1);
    torch_cell[1][2] = cell(1, 2);

    torch_cell[2][0] = cell(2, 0);
    torch_cell[2][1] = cell(2, 1);
    torch_cell[2][2] = cell(2, 2);

    const auto& positions = this->getPositions();

    auto torch_positions = torch::from_blob(
        const_cast<PLMD::Vector*>(positions.data()),
        {static_cast<int64_t>(positions.size()), 3},
        tensor_options
    );

    // setup torch's automatic gradient tracking
    if (!this->doNotCalculateDerivatives()) {
        torch_positions.requires_grad_(true);

        this->strain_ = torch::eye(3, tensor_options.requires_grad(true));

        // pretend to scale positions/cell by the strain so that it enters the
        // computational graph.
        torch_positions = torch_positions.matmul(this->strain_);
        torch_positions.retain_grad();

        torch_cell = torch_cell.matmul(this->strain_);
    }


    // TODO: move data to another dtype/device as requested by the model or user
    this->system_ = torch::make_intrusive<metatensor_torch::SystemHolder>(
        this->atomic_types_,
        torch_positions,
        torch_cell
    );

    // compute the neighbors list requested by the model, and register them with
    // the system
    auto nl_requests = this->model_.run_method("requested_neighbors_lists");
    for (auto request_ivalue: nl_requests.toList()) {
        auto request = request_ivalue.get().toCustomClass<metatensor_torch::NeighborsListOptionsHolder>();

        auto neighbors = this->computeNeighbors(request, positions, cell);
        metatensor_torch::register_autograd_neighbors(this->system_, neighbors, this->check_consistency_);
        this->system_->add_neighbors_list(request, neighbors);
    }
}


metatensor_torch::TorchTensorBlock MetatensorPlumedAction::computeNeighbors(
    metatensor_torch::NeighborsListOptions request,
    const std::vector<PLMD::Vector>& positions,
    const PLMD::Tensor& cell
) {
    auto labels_options = torch::TensorOptions().dtype(torch::kInt32).device(torch::kCPU);
    auto neighbor_component = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "xyz",
        torch::tensor({0, 1, 2}, labels_options).reshape({3, 1})
    );
    auto neighbor_properties = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        "distance", torch::zeros({1, 1}, labels_options)
    );

    auto cutoff = request->engine_cutoff(this->getUnits().getLengthString());

    auto periodic = (
        cell(0, 0) == 0.0 && cell(0, 1) == 0.0 && cell(0, 2) == 0.0 &&
        cell(1, 0) == 0.0 && cell(1, 1) == 0.0 && cell(1, 2) == 0.0 &&
        cell(2, 0) == 0.0 && cell(2, 2) == 0.0 && cell(2, 2) == 0.0
    );

    // use https://github.com/Luthaf/vesin to compute the requested neighbor
    // lists since we can not get these from PLUMED
    VesinOptions options;
    options.cutoff = cutoff;
    options.full = request->full_list();
    options.return_shifts = true;
    options.return_distances = false;
    options.return_vectors = true;

    VesinNeighborsList* vesin_neighbor_list = new VesinNeighborsList();
    memset(vesin_neighbor_list, 0, sizeof(VesinNeighborsList));

    const char* error_message = NULL;
    int status = vesin_neighbors(
        reinterpret_cast<const double (*)[3]>(positions.data()),
        positions.size(),
        reinterpret_cast<const double (*)[3]>(&cell(0, 0)),
        periodic,
        VesinCPU,
        options,
        vesin_neighbor_list,
        &error_message
    );

    if (status != EXIT_SUCCESS) {
        this->error(
            "failed to compute neighbor list (cutoff=" + std::to_string(cutoff) +
            "full=" + (request->full_list() ? "true" : "false") + "): " + error_message
        );
    }

    // transform from vesin to metatensor format
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

    auto pair_samples_values = torch::zeros({n_pairs, 5}, labels_options);
    for (unsigned i=0; i<n_pairs; i++) {
        pair_samples_values[i][0] = static_cast<int32_t>(vesin_neighbor_list->pairs[i][0]);
        pair_samples_values[i][1] = static_cast<int32_t>(vesin_neighbor_list->pairs[i][1]);
        pair_samples_values[i][2] = vesin_neighbor_list->shifts[i][0];
        pair_samples_values[i][3] = vesin_neighbor_list->shifts[i][1];
        pair_samples_values[i][4] = vesin_neighbor_list->shifts[i][2];
    }

    auto neighbor_samples = torch::make_intrusive<metatensor_torch::LabelsHolder>(
        std::vector<std::string>{"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
        pair_samples_values
    );

    auto neighbors = torch::make_intrusive<metatensor_torch::TensorBlockHolder>(
        pair_vectors,
        neighbor_samples,
        std::vector<metatensor_torch::TorchLabels>{neighbor_component},
        neighbor_properties
    );

    return neighbors;
}


void MetatensorPlumedAction::calculate() {
    this->createSystem();

    try {
        auto ivalue_output = this->model_.forward({
            std::vector<metatensor_torch::System>{this->system_},
            evaluations_options_,
            this->check_consistency_,
        });

        auto dict_output = ivalue_output.toGenericDict();
        auto cv = dict_output.at("collective_variable");
        this->output_ = cv.toCustomClass<metatensor_torch::TensorMapHolder>();
    } catch (const std::exception& e) {
        error("failed to evaluate the model: " + std::string(e.what()));
    }

    // send the output back to plumed
    plumed_massert(this->output_->keys()->count() == 1, "output should have a single block");
    auto block = metatensor_torch::TensorMapHolder::block_by_id(this->output_, 0);
    plumed_massert(block->components().empty(), "components are not yet supported in the output");
    auto torch_values = block->values().to(torch::kCPU).to(torch::kFloat64);
    auto n_samples = torch_values.size(0);
    auto n_properties = torch_values.size(1);

    Value* value = this->getPntrToComponent(0);
    const auto& value_shape = value->getShape();
    // reshape the plumed `Value` to hold the data returned by the model
    if (n_samples == 1) {
        if (n_properties == 1) {
            // the CV is a single scalar
            if (value->getRank() != 0) {
                log.printf("  output of metatensor model is a scalar\n");
                value->setShape({});
            }

            value->set(torch_values.item<double>());
        } else {
            // we have multiple CV describing a single thing (atom or full system)
            if (value->getRank() != 1 || value_shape[0] != n_properties) {
                log.printf("  output of metatensor model is a 1x%d vector\n", n_properties);
                value->setShape({static_cast<unsigned>(n_properties)});
            }

            for (unsigned i=0; i<n_properties; i++) {
                value->set(i, torch_values[0][i].item<double>());
            }
        }
    } else {
        if (n_properties == 1) {
            // we have a single CV describing multiple things (i.e. atoms)
            if (value->getRank() != 1 || value_shape[0] != n_samples) {
                log.printf("  output of metatensor model is a %dx1 vector\n", n_samples);
                value->setShape({static_cast<unsigned>(n_samples)});
            }

            // TODO: check sample order?
            for (unsigned i=0; i<n_samples; i++) {
                value->set(i, torch_values[i][0].item<double>());
            }
        } else {
            // the CV is a matrix
            if (value->getRank() != 2 || value_shape[0] != n_samples || value_shape[1] != n_properties) {
                log.printf("  output of metatensor model is a %dx%d matrix\n", n_samples, n_properties);
                value->setShape({
                    static_cast<unsigned>(n_samples),
                    static_cast<unsigned>(n_properties),
                });
                value->reshapeMatrixStore(n_properties);
            }

            // TODO: check sample order?
            for (unsigned i=0; i<n_samples; i++) {
                for (unsigned j=0; j<n_properties; j++) {
                    value->set(i * n_properties + j, torch_values[i][j].item<double>());
                }
            }
        }
    }
}


void MetatensorPlumedAction::apply() {
    auto* value = this->getPntrToComponent(0);
    if (!value->forcesWereAdded()) {
        return;
    }

    auto block = metatensor_torch::TensorMapHolder::block_by_id(this->output_, 0);
    auto torch_values = block->values().to(torch::kCPU).to(torch::kFloat64);
    auto n_samples = torch_values.size(0);
    auto n_properties = torch_values.size(1);

    auto output_grad = torch::zeros_like(torch_values);
    if (n_samples == 1) {
        if (n_properties == 1) {
            output_grad[0][0] = value->getForce();
        } else {
            for (unsigned i=0; i<n_properties; i++) {
                output_grad[0][i] = value->getForce(i);
            }
        }
    } else {
        if (n_properties == 1) {
            // TODO: check sample order?
            for (unsigned i=0; i<n_samples; i++) {
                output_grad[i][0] = value->getForce(i);
            }
        } else {
            // TODO: check sample order?
            for (unsigned i=0; i<n_samples; i++) {
                for (unsigned j=0; j<n_properties; j++) {
                    output_grad[i][j] = value->getForce(i * n_properties + j);
                }
            }
        }
    }

    auto input_grad = torch::autograd::grad(
        {torch_values},
        {this->system_->positions(), this->strain_},
        {output_grad}
    );
    plumed_assert(input_grad[0].is_cpu());
    plumed_assert(input_grad[0].is_contiguous());

    plumed_assert(input_grad[1].is_cpu());
    plumed_assert(input_grad[1].is_contiguous());

    auto positions_grad = input_grad[0];
    auto strain_grad = input_grad[1];

    auto derivatives = std::vector<double>(
        positions_grad.data_ptr<double>(),
        positions_grad.data_ptr<double>() + 3 * this->system_->size()
    );

    // add virials to the derivatives
    derivatives.push_back(strain_grad[0][0].item<double>());
    derivatives.push_back(strain_grad[0][1].item<double>());
    derivatives.push_back(strain_grad[0][2].item<double>());

    derivatives.push_back(strain_grad[1][0].item<double>());
    derivatives.push_back(strain_grad[1][1].item<double>());
    derivatives.push_back(strain_grad[1][2].item<double>());

    derivatives.push_back(strain_grad[2][0].item<double>());
    derivatives.push_back(strain_grad[2][1].item<double>());
    derivatives.push_back(strain_grad[2][2].item<double>());


    unsigned index = 0;
    this->setForcesOnAtoms(derivatives, index);
}

}} // namespace PLMD::metatensor

#endif


namespace PLMD { namespace metatensor {
    // use the same implementation for both the actual action and the dummy one
    // (when libtorch and libmetatensor could not be found).
    void MetatensorPlumedAction::registerKeywords(Keywords& keys) {
        Action::registerKeywords(keys);
        ActionAtomistic::registerKeywords(keys);
        ActionWithValue::registerKeywords(keys);

        keys.add("compulsory", "MODEL", "path to the exported metatensor model");
        keys.add("optional", "EXTENSIONS_DIRECTORY", "path to the directory containing TorchScript extensions to load");

        keys.add("numbered", "SPECIES", "the atoms in each PLUMED species");
        keys.reset_style("SPECIES", "atoms");

        keys.add("optional", "SPECIES_TO_METATENSOR_TYPES", "mapping from PLUMED SPECIES to metatensor's atomic types");
    }

    PLUMED_REGISTER_ACTION(MetatensorPlumedAction, "METATENSOR")
}}
