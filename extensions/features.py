"""Compile/link probes for optional native features, independent of conda policy."""

import shlex
import subprocess

from symbols import module_symbol

PROBES = {
    "BOOST_SERIALIZATION": {
        "source": "#include <boost/archive/text_oarchive.hpp>\n#include <sstream>\nint main(){std::ostringstream s; boost::archive::text_oarchive a(s); int n=0; a << n;}",
        "libraries": ["boost_serialization"],
    },
    "FFTW": {
        "source": "#include <fftw3.h>\nint main(){return fftw_plan_dft_1d(1, nullptr, nullptr, FFTW_FORWARD, FFTW_ESTIMATE)==nullptr;}",
        "libraries": ["fftw3"],
    },
    "LIBTORCH": {
        "source": "#include <torch/torch.h>\n#include <torch/script.h>\nint main(){return torch::ones({1}).numel()!=1;}",
        "libraries": ["torch_cpu", "c10"],
        "includes": ["include/torch/csrc/api/include"],
    },
    "LIBMETATOMIC": {
        "source": """#include <torch/script.h>
#include <metatensor.hpp>
#include <metatensor/torch.hpp>
#include <metatomic/torch.hpp>
#if METATOMIC_TORCH_VERSION_MAJOR != 0 || METATOMIC_TORCH_VERSION_MINOR != 1 || METATOMIC_TORCH_VERSION_PATCH < 12
#error metatomic-torch >=0.1.12,<0.2 is required
#endif
int main(){metatomic_torch::version(); metatensor_torch::version(); mts_version();}
""",
        "libraries": [
            "metatomic_torch",
            "metatensor_torch",
            "metatensor",
            "torch_cpu",
            "c10",
        ],
        "includes": ["include/torch/csrc/api/include"],
    },
    "LIBHDF5": {
        "source": "#include <hdf5.h>\nint main(){unsigned a,b,c; return H5get_libversion(&a,&b,&c);}",
        "libraries": ["hdf5_hl", "hdf5"],
    },
    "LIBITENSOR": {
        "source": '#include <itensor/all.h>\nint main(){auto i=itensor::Index(2,"i"); auto t=itensor::ITensor(i); return 0;}',
        "libraries": ["itensor", "openblas"],
    },
    "KENREF": {
        "source": "#include <KEnRefBias.h>\nint main(){return 0;}",
        "libraries": ["kenref_core"],
        "pkg_config": ["kenref_plumed", "kenref_core"],
    },
}


class FeatureProbes:
    def __init__(self, compiler, cppflags, cxxflags, ldflags, prefix, directory):
        self.compiler = compiler
        self.cppflags = cppflags
        self.cxxflags = cxxflags
        self.ldflags = ldflags
        self.prefix = prefix
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=True)
        self.results = {}

    def probe(self, feature):
        if feature in self.results:
            return self.results[feature]
        definition = PROBES[feature]
        includes = ["-I" + str(self.prefix / "include")]
        includes += [
            "-I" + str(self.prefix / path) for path in definition.get("includes", [])
        ]
        native = ["-l" + name for name in definition["libraries"]]
        pkg = definition.get("pkg_config", [])
        if pkg:
            try:
                includes += shlex.split(
                    subprocess.check_output(
                        ["pkg-config", "--cflags", *pkg],
                        text=True,
                        stderr=subprocess.DEVNULL,
                    )
                )
                native += shlex.split(
                    subprocess.check_output(
                        ["pkg-config", "--libs", *pkg],
                        text=True,
                        stderr=subprocess.DEVNULL,
                    )
                )
            except (OSError, subprocess.CalledProcessError):
                # pkg-config is optional; the compile/link probe below still
                # checks the conventional include paths and library names.
                pass
        source = self.directory / (feature + ".cpp")
        output = self.directory / (feature + ".exe")
        source.write_text(definition["source"])
        command = [
            *self.compiler,
            *self.cppflags,
            *self.cxxflags,
            *includes,
            source,
            *self.ldflags,
            "-L" + str(self.prefix / "lib"),
            *native,
            "-o",
            output,
        ]
        result = subprocess.run(
            list(map(str, command)), check=False, text=True, capture_output=True
        )
        log = result.stdout + result.stderr
        # Recent LibTorch headers require C++20. Upgrade only the modules using
        # these headers, retaining the kernel's C++ ABI and compiler settings.
        if result.returncode and feature in ("LIBTORCH", "LIBMETATOMIC"):
            includes.append("-std=c++20")
            command.insert(
                len(self.compiler) + len(self.cppflags) + len(self.cxxflags),
                "-std=c++20",
            )
            result = subprocess.run(
                list(map(str, command)), check=False, text=True, capture_output=True
            )
            log += "\nRetry with C++20:\n" + result.stdout + result.stderr
        (self.directory / (feature + ".log")).write_text(log)
        enabled = result.returncode == 0
        self.results[feature] = (
            enabled,
            includes if enabled else [],
            native if enabled else [],
        )
        print(
            f"Optional feature {feature}: {'enabled' if enabled else 'unavailable; module remains enabled'}",
            flush=True,
        )
        return self.results[feature]

    def module_flags(self, features, required=(), optional=(), disabled=()):
        unknown = (set(required) | set(optional)) - PROBES.keys()
        if unknown:
            raise ValueError(
                "No native compile/link probe for " + ", ".join(sorted(unknown))
            )
        features = (set(features) | set(required) | set(optional)) - set(disabled)
        cpp = []
        libraries = []
        enabled = {}
        for feature in PROBES:
            if feature not in features:
                continue
            available, includes, native = self.probe(feature)
            cpp.append("-U__PLUMED_HAS_" + feature)
            if available:
                cpp += ["-D__PLUMED_HAS_" + feature + "=1", *includes]
                libraries += native
            enabled[feature] = available
            if feature in required and not available:
                raise ValueError(
                    "Required conda feature " + feature + " is unavailable"
                )
        cpp += ["-U__PLUMED_HAS_" + feature for feature in disabled]
        enabled.update({feature: False for feature in disabled})
        return cpp, list(dict.fromkeys(libraries)), enabled


def feature_identity(name, features):
    """Export actual compiled features instead of trusting installed metadata."""
    code = f'extern "C" const char* {module_symbol(name)}_features() {{ return ""\n'
    for feature in sorted(set(features)):
        code += f"#if defined(__PLUMED_HAS_{feature}) && __PLUMED_HAS_{feature}\n"
        names = {"MPI": "mpi real_mpi", "THREAD_MPI": "thread_mpi mpi"}.get(
            feature, feature.lower()
        )
        code += '  "' + names + ' "\n#endif\n'
    return code + "; }\n"
