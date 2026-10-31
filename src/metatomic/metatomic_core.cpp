/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2026 The METATOMIC-PLUMED team
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

/*INDENT-OFF*/
#ifdef __PLUMED_HAS_LIBMETATOMIC

#include "MetatomicAction.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <type_traits>

#include <unistd.h>
#ifdef __PLUMED_HAS_DLOPEN
#include <dlfcn.h>
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include <metatensor.hpp>
#include <metatomic.hpp>

#pragma GCC diagnostic pop

#include "vesin.h"

namespace PLMD {
namespace metatomic {

namespace mta = ::metatomic;

static_assert(std::is_standard_layout<PLMD::Vector>::value);
static_assert(sizeof(PLMD::Vector) == sizeof(std::array<double, 3>));
static_assert(sizeof(PLMD::Tensor) == sizeof(std::array<std::array<double, 3>, 3>));

namespace {

template<typename T>
constexpr DLDataType dlpackDType() {
    if constexpr (std::is_same_v<T, double>) {
        return {kDLFloat, 64, 1};
    } else if constexpr (std::is_same_v<T, float>) {
        return {kDLFloat, 32, 1};
    } else if constexpr (std::is_same_v<T, int32_t>) {
        return {kDLInt, 32, 1};
    } else {
        static_assert(std::is_same_v<T, uint8_t>, "unsupported DLPack type");
        return {kDLBool, 8, 1};
    }
}

/// DLPack tensor and the data it needs to stay valid. `owner` keeps copied
/// buffers alive; it is empty when `data` is a view of PLUMED's storage.
struct DLPackHolder {
    DLManagedTensorVersioned managed;
    std::vector<int64_t> shape;
    std::vector<int64_t> strides;
    std::shared_ptr<void> owner;
};

/// Row-major CPU DLPack tensor over `data`. The deleter frees the holder (and
/// `owner`), never `data` itself.
template<typename T>
mta::DLPackTensor dlpack(T* data, std::vector<int64_t> shape, std::shared_ptr<void> owner, bool read_only) {
    auto* holder = new DLPackHolder{};
    holder->shape = std::move(shape);
    holder->strides.resize(holder->shape.size(), 1);
    for (size_t i = holder->shape.size(); i > 1; i--) {
        holder->strides[i - 2] = holder->strides[i - 1] * holder->shape[i - 1];
    }
    holder->owner = std::move(owner);

    auto& managed = holder->managed;
    managed.version = {DLPACK_MAJOR_VERSION, DLPACK_MINOR_VERSION};
    managed.manager_ctx = holder;
    managed.deleter = [](DLManagedTensorVersioned* self) {
        delete static_cast<DLPackHolder*>(self->manager_ctx);
    };
    managed.flags = read_only ? DLPACK_FLAG_BITMASK_READ_ONLY : 0;
    managed.dl_tensor.data = const_cast<std::remove_const_t<T>*>(data);
    managed.dl_tensor.device = {kDLCPU, 0};
    managed.dl_tensor.ndim = static_cast<int32_t>(holder->shape.size());
    managed.dl_tensor.dtype = dlpackDType<std::remove_const_t<T>>();
    managed.dl_tensor.shape = holder->shape.data();
    managed.dl_tensor.strides = holder->strides.data();
    managed.dl_tensor.byte_offset = 0;

    return mta::DLPackTensor(&managed);
}

/// DLPack tensor owning a copy of `values`
template<typename T>
mta::DLPackTensor ownedDLPack(std::vector<T> values, std::vector<int64_t> shape) {
    auto storage = std::make_shared<std::vector<T>>(std::move(values));
    return dlpack(storage->data(), std::move(shape), storage, /*read_only=*/false);
}

/// Copy a strided DLPack array to a row-major vector of doubles
template<typename T>
std::vector<double> toDoubles(const metatensor::DLPackArray<T>& array) {
    const auto& shape = array.shape();
    const auto& strides = array.strides();
    size_t size = 1;
    for (auto dim: shape) {
        size *= dim;
    }

    auto result = std::vector<double>(size);
    for (size_t linear = 0; linear < size; linear++) {
        int64_t offset = 0;
        size_t remaining = linear;
        for (size_t dim = shape.size(); dim > 0; dim--) {
            offset += static_cast<int64_t>(remaining % shape[dim - 1]) * strides[dim - 1];
            remaining /= shape[dim - 1];
        }
        result[linear] = static_cast<double>(array.data()[offset]);
    }
    return result;
}

std::vector<double> valuesAsDoubles(metatensor::TensorBlock& block) {
    auto dtype = block.dtype();
    if (dtype.code == kDLFloat && dtype.bits == 64) {
        return toDoubles(block.values<double>());
    } else if (dtype.code == kDLFloat && dtype.bits == 32) {
        return toDoubles(block.values<float>());
    }
    plumed_merror("metatomic model returned data with an unsupported dtype, expected float32 or float64");
}

/// Pair list in metatomic format, built from a vesin neighbor list
template<typename T>
metatensor::TensorBlock pairsBlock(const vesin::VesinNeighborList& neighbors) {
    auto n_pairs = neighbors.length;

    auto samples = std::vector<int32_t>(5 * n_pairs);
    auto vectors = std::vector<T>(3 * n_pairs);
    for (size_t i = 0; i < n_pairs; i++) {
        samples[5 * i + 0] = static_cast<int32_t>(neighbors.pairs[i][0]);
        samples[5 * i + 1] = static_cast<int32_t>(neighbors.pairs[i][1]);
        samples[5 * i + 2] = neighbors.shifts[i][0];
        samples[5 * i + 3] = neighbors.shifts[i][1];
        samples[5 * i + 4] = neighbors.shifts[i][2];
        for (size_t k = 0; k < 3; k++) {
            vectors[3 * i + k] = static_cast<T>(neighbors.vectors[i][k]);
        }
    }

    return metatensor::TensorBlock(
        std::make_unique<metatensor::SimpleDataArray<T>>(
            std::vector<uintptr_t>{n_pairs, 3, 1}, std::move(vectors)
        ),
        // vesin should create unique pairs
        metatensor::Labels(
            {"first_atom", "second_atom", "cell_shift_a", "cell_shift_b", "cell_shift_c"},
            samples.data(), n_pairs, metatensor::assume_unique{}
        ),
        {metatensor::Labels({"xyz"}, {{0}, {1}, {2}})},
        metatensor::Labels({"distance"}, {{0}})
    );
}

/// The Python plugin installed next to the libmetatomic used by PLUMED (with
/// the same layout as `metatomic.utils.python_plugin_path`), or an empty string
std::string defaultPythonPlugin() {
#ifdef __PLUMED_HAS_DLOPEN
    Dl_info info;
    auto* symbol = dlsym(RTLD_DEFAULT, "mta_load_model");
    if (symbol != nullptr && dladdr(symbol, &info) != 0 && info.dli_fname != nullptr) {
        auto plugin = std::filesystem::path(info.dli_fname).parent_path().parent_path()
            / "libexec" / "metatomic" / "metatomic-python-plugin.so";
        if (std::filesystem::exists(plugin)) {
            return plugin.string();
        }
    }
#endif
    return "";
}

/// Temporary script for metatomic-core's Python plugin, loading the TorchScript
/// model at `path` with metatomic-torch's `LegacyModel`. The file is removed
/// when this goes out of scope.
class TorchScriptLoader {
public:
    TorchScriptLoader(const std::string& path, const std::string& extensions_directory) {
        // JSON strings are valid Python string literals
        auto script = std::string(
            "import metatensor\n"
            "from metatomic_torch.legacy import LegacyModel\n"
            "\n"
            "# PLUMED gives pair lists as metatensor's C++ SimpleDataArray\n"
            "metatensor.register_external_data_wrapper(\n"
            "    \"metatensor::SimpleDataArray\", metatensor.ExternalCpuArray\n"
            ")\n"
            "\n"
            "model = LegacyModel.load(\n"
            "    ") + nlohmann::json(path).dump() + ",\n"
            "    extensions_directory=" + (
                extensions_directory.empty() ? "None" : nlohmann::json(extensions_directory).dump()
            ) + ",\n"
            ")\n";

        file_ = (std::filesystem::temp_directory_path() / "plumed-metatomic-XXXXXX.py").string();
        auto fd = mkstemps(file_.data(), 3);
        if (fd == -1) {
            throw std::runtime_error("could not create a temporary file to load the TorchScript model");
        }
        close(fd);
        std::ofstream(file_) << script;
    }

    ~TorchScriptLoader() {
        std::error_code ignored;
        std::filesystem::remove(file_, ignored);
    }

    TorchScriptLoader(const TorchScriptLoader&) = delete;
    TorchScriptLoader& operator=(const TorchScriptLoader&) = delete;

    const std::string& path() const {
        return file_;
    }

private:
    std::string file_;
};

mta::ExternalModel loadModel(
    MetatomicPlumedAction& action,
    std::vector<std::string> plugin_libraries,
    std::string plugin,
    const std::string& model
) {
    std::string device;
    action.parse("DEVICE", device);
    if (!device.empty() && device != "cpu") {
        action.error("metatomic models only run on CPU in PLUMED, DEVICE=" + device + " is not supported");
    }

    // TorchScript models exported with metatomic-torch run in the Python plugin
    auto torchscript = plugin.empty() && model.size() > 3 && model.compare(model.size() - 3, 3, ".pt") == 0;

    std::string extensions_directory;
    action.parse("EXTENSIONS_DIRECTORY", extensions_directory);
    if (!extensions_directory.empty() && !torchscript) {
        action.error("EXTENSIONS_DIRECTORY can only be used with TorchScript (.pt) models");
    }

    if (plugin_libraries.empty()) {
        auto python = defaultPythonPlugin();
        if (!python.empty()) {
            plugin_libraries.push_back(python);
        } else if (torchscript) {
            action.error(
                "could not find metatomic's Python plugin, which is required to load TorchScript models. "
                "Give its path with PLUGIN_LIBRARY (see `metatomic.utils.python_plugin_path`)"
            );
        }
    }

    try {
        for (const auto& library: plugin_libraries) {
            action.log.printf("  loading metatomic plugin from %s\n", library.c_str());
            mta::load_plugin(library);
        }

        auto loader = std::optional<TorchScriptLoader>();
        if (torchscript) {
            action.log.printf("  loading TorchScript model with metatomic-torch's LegacyModel\n");
            loader.emplace(model, extensions_directory);
            plugin = "python";
        }

        auto plugin_name = plugin.empty() ? std::nullopt : std::optional<std::string>(plugin);
        return mta::load_model(loader ? loader->path() : model, /*options_json=*/ std::nullopt, plugin_name);
    } catch (const std::exception& e) {
        action.error("failed to load model '" + model + "' with metatomic-core: " + e.what());
    }
}

} // namespace

/// Run models loaded through metatomic-core plugins. Forces use the gradients
/// of the model output, computed by the model (explicitly or with autograd).
class CoreBackend final: public MetatomicPlumedAction::Backend {
public:
    CoreBackend(
        MetatomicPlumedAction& action,
        const std::vector<std::string>& plugin_libraries,
        const std::string& plugin
    );
    ~CoreBackend() override;

    CoreBackend(const CoreBackend&) = delete;
    CoreBackend& operator=(const CoreBackend&) = delete;

    void calculate() override;
    void apply() override;

private:
    struct PairList {
        mta::PairListOptions options;
        double engine_cutoff;
        vesin::VesinNeighborList neighbors;
    };

    template<typename T>
    mta::System createSystem(bool empty);

    metatensor::TensorMap execute(
        mta::System system,
        std::optional<metatensor::Labels> selected_atoms,
        const std::vector<mta::Quantity>& outputs
    );

    MetatomicPlumedAction& action_;
    std::string model_name_;
    mta::ExternalModel model_;

    // requested output, with and without the gradients needed for forces
    std::vector<mta::Quantity> outputs_;
    std::vector<mta::Quantity> values_only_;
    std::vector<PairList> pair_lists_;
    std::vector<int32_t> atomic_types_;
    std::optional<metatensor::Labels> selected_atoms_;
    bool float64_ = true;
    bool check_consistency_ = false;
    // multiply gradients w.r.t. positions by this to go from model to engine units
    double positions_gradient_factor_ = 1.0;

    std::optional<metatensor::TensorMap> output_;
    // row of the PLUMED value where each sample of `output_` is stored
    std::vector<unsigned> output_rows_;
};

CoreBackend::CoreBackend(
    MetatomicPlumedAction& action,
    const std::vector<std::string>& plugin_libraries,
    const std::string& plugin
):
    action_(action),
    model_name_([&]() { std::string model; action.parse("MODEL", model); return model; }()),
    model_(loadModel(action, plugin_libraries, plugin, model_name_))
{
    auto& log = action.log;

    try {
        auto metadata = model_.metadata();
        log.printf("\n%s\n", mta::format_metadata(nlohmann::json(metadata).dump()).c_str());
        const auto& references = metadata.references();
        for (const auto* refs: {&references.model(), &references.architecture(), &references.implementation()}) {
            for (const auto& ref: *refs) {
                action.cite(ref);
            }
        }

        auto capabilities = model_.capabilities();

        const auto& devices = capabilities.supported_devices();
        if (std::find(devices.begin(), devices.end(), mta::ModelCapabilities::Device::CPU) == devices.end()) {
            action.error("this model does not support running on CPU, which is required for metatomic-core models in PLUMED");
        }
        float64_ = capabilities.dtype() == mta::ModelCapabilities::DType::Float64;
        log.printf("  running model on cpu device with %s data\n", float64_ ? "float64" : "float32");

        this->atomic_types_ = action.parseAtomicTypes();

        action.parseFlag("CHECK_CONSISTENCY", check_consistency_);
        if (check_consistency_) {
            log.printf("  checking for internal consistency of the model\n");
        }

        std::string variant;
        action.parse("VARIANT", variant);
        auto feature_name = variant.empty() ? std::string("feature") : "feature/" + variant;

        const auto& available = capabilities.outputs();
        auto declared = std::find_if(available.begin(), available.end(), [&](const auto& q) {
            return q.name() == feature_name;
        });
        if (declared == available.end()) {
            std::string names;
            for (const auto& q: available) {
                names += (names.empty() ? "" : ", ") + q.name();
            }
            action.error("this model does not have a '" + feature_name + "' output, available outputs are [" + names + "]");
        }

        // only request the gradients the model knows how to compute
        auto gradients = std::vector<mta::Gradients>{};
        for (auto gradient: {mta::Gradients::Positions, mta::Gradients::Strain}) {
            const auto& declared_gradients = declared->gradients();
            if (std::find(declared_gradients.begin(), declared_gradients.end(), gradient) != declared_gradients.end()) {
                gradients.push_back(gradient);
            }
        }

        auto request = [&](std::vector<mta::Gradients> gradients) {
            return mta::Quantity::builder()
                .name(feature_name)
                .unit("")
                .sample_kind(declared->sample_kind())
                .gradients(std::move(gradients))
                .build();
        };
        outputs_.push_back(request(gradients));
        values_only_.push_back(request({}));

        auto inputs = model_.requested_inputs(outputs_);
        if (!inputs.empty()) {
            std::string names;
            for (const auto& q: inputs) {
                names += (names.empty() ? "" : ", ") + q.name();
            }
            action.error(
                "The model requested extra inputs (" + names + ") that are not "
                "yet supported in PLUMED. Please open an issue to request support."
            );
        }

        auto length_unit = action.getUnits().getLengthString();
        auto model_length_unit = capabilities.length_unit();
        auto to_engine = mta::unit_conversion_factor(model_length_unit, length_unit);
        positions_gradient_factor_ = mta::unit_conversion_factor(length_unit, model_length_unit);

        log.printf("  the following neighbor lists have been requested:\n");
        for (auto& options: model_.requested_pair_lists(outputs_)) {
            log.printf("    - %s list, %g %s cutoff (requested %g %s)\n",
                options.full_list() ? "full" : "half",
                options.cutoff() * to_engine,
                length_unit.c_str(),
                options.cutoff(),
                model_length_unit.c_str()
            );
            pair_lists_.push_back({options, options.cutoff() * to_engine, vesin::VesinNeighborList()});
        }

        // determine how many properties there will be in the output by running
        // the model once on an empty system
        auto empty = float64_ ? createSystem<double>(/*empty=*/true) : createSystem<float>(/*empty=*/true);
        auto dummy = this->execute(std::move(empty), std::nullopt, values_only_);
        auto n_properties = static_cast<unsigned>(dummy.block_by_id(0).properties().count());

        // parse and handle atom sub-selection. This is done AFTER determining
        // the output size, since the selection might not be valid for the
        // empty system
        auto selected_atoms = action.parseSelectedAtoms(atomic_types_.size());
        if (!selected_atoms.empty()) {
            auto values = std::vector<int32_t>();
            for (const auto& atom: selected_atoms) {
                values.push_back(0);
                values.push_back(static_cast<int32_t>(atom.index()));
            }
            selected_atoms_ = metatensor::Labels({"system", "atom"}, values.data(), selected_atoms.size());
        }

        unsigned n_samples = 1;
        if (declared->sample_kind() == mta::SampleKind::Atom) {
            n_samples = static_cast<unsigned>(selected_atoms.empty() ? atomic_types_.size() : selected_atoms.size());
        } else if (declared->sample_kind() != mta::SampleKind::System) {
            action.error("only 'system' and 'atom' sample kinds are supported for the '" + feature_name + "' output");
        }

        action.addOutputValue(n_samples, n_properties);
    } catch (const mta::Error& e) {
        action.error(std::string("metatomic-core error: ") + e.what());
    } catch (const metatensor::Error& e) {
        action.error(std::string("metatensor error: ") + e.what());
    }
}

CoreBackend::~CoreBackend() {
    for (auto& pairs: pair_lists_) {
        vesin::vesin_free(&pairs.neighbors);
    }
}

template<typename T>
mta::System CoreBackend::createSystem(bool empty) {
    auto [cell, periodic] = empty ?
        std::pair<Tensor, std::array<bool, 3>>{Tensor(), {false, false, false}} :
        action_.cellAndPeriodicity();

    const auto& positions = action_.getPositions();
    auto n_atoms = empty ? size_t(0) : positions.size();
    if (!empty && n_atoms != atomic_types_.size()) {
        plumed_merror(
            "METATOMIC action needs to know about all atoms in the system. There are " +
            std::to_string(n_atoms) + " atoms overall, but we only have atomic types for " +
            std::to_string(atomic_types_.size()) + " of them."
        );
    }
    auto n = static_cast<int64_t>(n_atoms);
    const auto* raw_positions = n_atoms == 0 ? nullptr : &positions[0][0];

    // positions are shared with PLUMED when possible, everything else is copied
    auto dl_positions = [&]() {
        if constexpr (std::is_same_v<T, double>) {
            return dlpack(raw_positions, {n, 3}, nullptr, /*read_only=*/true);
        } else {
            return ownedDLPack(std::vector<T>(raw_positions, raw_positions + 3 * n_atoms), {n, 3});
        }
    }();

    auto system = mta::System(
        action_.getUnits().getLengthString(),
        ownedDLPack(empty ? std::vector<int32_t>() : atomic_types_, {n}),
        std::move(dl_positions),
        ownedDLPack(std::vector<T>(&cell(0, 0), &cell(0, 0) + 9), {3, 3}),
        ownedDLPack(std::vector<uint8_t>(periodic.begin(), periodic.end()), {3})
    );

    for (auto& pairs: pair_lists_) {
        if (empty) {
            system.add_pairs(pairs.options, pairsBlock<T>(vesin::VesinNeighborList()));
            continue;
        }

        vesin::VesinOptions options{};
        options.cutoff = pairs.engine_cutoff;
        options.full = pairs.options.full_list();
        options.return_shifts = true;
        options.return_distances = false;
        options.return_vectors = true;
        options.algorithm = vesin::VesinAutoAlgorithm;

        const char* error_message = nullptr;
        int status = vesin::vesin_neighbors(
            reinterpret_cast<const double (*)[3]>(raw_positions),
            n_atoms,
            reinterpret_cast<const double (*)[3]>(&cell(0, 0)),
            periodic.data(),
            {vesin::VesinCPU, 0},
            options,
            &pairs.neighbors,
            &error_message
        );
        if (status != EXIT_SUCCESS) {
            plumed_merror(
                "failed to compute neighbor list (cutoff=" + std::to_string(pairs.engine_cutoff) +
                ", full=" + (options.full ? "true" : "false") + "): " + error_message
            );
        }

        system.add_pairs(pairs.options, pairsBlock<T>(pairs.neighbors));
    }

    return system;
}

metatensor::TensorMap CoreBackend::execute(
    mta::System system,
    std::optional<metatensor::Labels> selected_atoms,
    const std::vector<mta::Quantity>& outputs
) {
    auto systems = std::vector<mta::System>();
    systems.push_back(std::move(system));

    auto results = std::vector<metatensor::TensorMap>();
    try {
        results = mta::execute_model(model_, systems, std::move(selected_atoms), outputs, check_consistency_);
    } catch (const std::exception& e) {
        plumed_merror("failed to evaluate the model: " + std::string(e.what()));
    }

    auto& output = results.at(0);
    plumed_massert(output.keys().count() == 1, "output should have a single block");
    plumed_massert(output.block_by_id(0).components().empty(), "components are not yet supported in the output");
    return std::move(output);
}

void CoreBackend::calculate() {
    auto system = float64_ ? createSystem<double>(/*empty=*/false) : createSystem<float>(/*empty=*/false);
    // gradients can be expensive, only ask for them when forces might be applied
    const auto& outputs = action_.doNotCalculateDerivatives() ? values_only_ : outputs_;
    output_ = this->execute(std::move(system), selected_atoms_, outputs);

    auto block = output_->block_by_id(0);
    auto values = valuesAsDoubles(block);
    auto samples = block.samples();

    auto n_samples = action_.n_samples_;
    auto n_properties = action_.n_properties_;
    if (samples.count() != n_samples) {
        plumed_merror(
            "expected the model to return a TensorBlock with " + std::to_string(n_samples) +
            " samples, got " + std::to_string(samples.count()) + " instead"
        );
    } else if (block.properties().count() != n_properties) {
        plumed_merror(
            "expected the model to return a TensorBlock with " + std::to_string(n_properties) +
            " properties, got " + std::to_string(block.properties().count()) + " instead"
        );
    }

    output_rows_ = {0};
    if (n_samples != 1) {
        const auto sample_names = samples.names();
        plumed_assert((std::vector<std::string>(sample_names.begin(), sample_names.end()) == std::vector<std::string>{"system", "atom"}));
        const auto samples_values = samples.values_cpu();

        // samples might be returned in any order. If the users picked some
        // selected atoms, we store the output in the same order as the selection
        output_rows_.resize(n_samples);
        for (unsigned i=0; i<n_samples; i++) {
            auto atom = samples_values(i, 1);
            if (selected_atoms_.has_value()) {
                auto position = selected_atoms_->position({0, atom});
                plumed_assert(position.has_value());
                output_rows_[i] = static_cast<unsigned>(position.value());
            } else {
                output_rows_[i] = static_cast<unsigned>(atom);
            }
        }
    }

    action_.setOutputValues(output_rows_, values);
}

void CoreBackend::apply() {
    const auto* value = action_.getPntrToComponent(0);
    if (!value->forcesWereAdded()) {
        return;
    }

    auto block = output_->block_by_id(0);
    auto gradients = block.gradients_list();
    for (const auto* parameter: {"positions", "strain"}) {
        if (std::find(gradients.begin(), gradients.end(), parameter) == gradients.end()) {
            auto per_atom = outputs_[0].sample_kind() == mta::SampleKind::Atom;
            plumed_merror(
                std::string("forces were applied to this METATOMIC action, but the model did not compute "
                "gradients of '") + outputs_[0].name() + "' with respect to " + parameter +
                (per_atom ? ". Per-atom features can only be biased if the model computes these gradients "
                "explicitly, metatomic-core does not compute them with autograd" : "")
            );
        }
    }

    auto n_properties = action_.n_properties_;
    auto forces = action_.outputForces(output_rows_);
    // contract one gradient row (`n_properties` values) with the output forces
    auto contract = [&](const std::vector<double>& gradient, size_t row, size_t sample) {
        double sum = 0;
        for (unsigned p=0; p<n_properties; p++) {
            sum += gradient[row * n_properties + p] * forces[sample * n_properties + p];
        }
        return sum;
    };

    auto derivatives = std::vector<double>(3 * atomic_types_.size() + 9, 0.0);

    // positions gradients: samples (sample, system, atom), component xyz
    auto positions_gradient = block.gradient("positions");
    auto positions_values = valuesAsDoubles(positions_gradient);
    const auto positions_samples = positions_gradient.samples().values_cpu();
    for (size_t g=0; g<positions_gradient.samples().count(); g++) {
        auto sample = static_cast<size_t>(positions_samples(g, 0));
        auto atom = static_cast<size_t>(positions_samples(g, 2));
        for (size_t x=0; x<3; x++) {
            derivatives[3 * atom + x] += positions_gradient_factor_ * contract(positions_values, 3 * g + x, sample);
        }
    }

    // strain gradients: sample (sample), components (xyz_1, xyz_2). PLUMED
    // expects minus the strain derivative for the virial
    auto strain_gradient = block.gradient("strain");
    auto strain_values = valuesAsDoubles(strain_gradient);
    const auto strain_samples = strain_gradient.samples().values_cpu();
    auto* virial = &derivatives[3 * atomic_types_.size()];
    for (size_t g=0; g<strain_gradient.samples().count(); g++) {
        auto sample = static_cast<size_t>(strain_samples(g, 0));
        for (size_t ab=0; ab<9; ab++) {
            virial[ab] -= contract(strain_values, 9 * g + ab, sample);
        }
    }

    unsigned index = 0;
    action_.setForcesOnAtoms(derivatives, index);
}

std::unique_ptr<MetatomicPlumedAction::Backend> makeCoreBackend(
    MetatomicPlumedAction& action,
    const std::vector<std::string>& plugin_libraries,
    const std::string& plugin
) {
    return std::make_unique<CoreBackend>(action, plugin_libraries, plugin);
}

} // namespace metatomic
} // namespace PLMD

#endif
