"""Smoke checks run against the repaired wheel in cibuildwheel's test environment."""

import array
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import plumed
from plumed_runtime import kernel_path, package_root, resource_root


def check_kernel(kernel=None, lines=("d: DISTANCE ATOMS=1,2",)):
    p = plumed.Plumed(kernel=kernel)
    p.cmd("setMDEngine", "python")
    p.cmd("setNatoms", 2)
    p.cmd("init")
    for line in lines:
        p.cmd("readInputLine", line)
    p.cmd("readInputLine", "RESTRAINT ARG=d AT=0 KAPPA=1")
    positions = array.array("d", [0, 0, 0, 1, 2, 3])
    forces = array.array("d", [0, 0, 0, 0, 0, 0])
    # a bias on atomic positions also acts on the cell (virial), so PLUMED
    # needs the box and the virial array even when NOPBC is used
    box = array.array("d", [10, 0, 0, 0, 10, 0, 0, 0, 10])
    virial = array.array("d", [0] * 9)
    p.cmd("setStep", 0)
    p.cmd("setBox", box)
    p.cmd("setPositions", positions)
    p.cmd("setForces", forces)
    p.cmd("setVirial", virial)
    p.cmd("calc")
    bias = array.array("d", [0])
    p.cmd("getBias", bias)
    assert abs(bias[0] - 7.0) < 1.0e-12
    assert abs(forces[0] - 1.0) < 1.0e-12
    p.finalize()


def enabled_modules():
    info = (resource_root() / "build-info.txt").read_text()
    return re.search(r"^Enabled modules: (.*)$", info, re.M).group(1).split(", ")


def documented_names(source_tree, modules):
    """Actions and command line tools documented in the sources of the modules.

    Documentation inside #if blocks depending on optional features
    (__PLUMED_HAS_..., e.g. EMMIVOX with libtorch) is skipped.
    """
    actions, tools = set(), set()
    for module in modules:
        for path in sorted((source_tree / "src" / module).glob("*.cpp")):
            optional = []
            for line in path.read_text(errors="replace").splitlines():
                line = line.strip()
                if line.startswith("#if"):
                    optional.append(
                        "__PLUMED_HAS_" in line and not line.startswith("#ifndef")
                    )
                elif line.startswith("#endif") and optional:
                    optional.pop()
                match = re.match(r"//\+PLUMEDOC\s+(\S+)\s+(\S+)", line)
                # INTERNAL pages document concepts, not actions
                if match and not any(optional) and match.group(1) != "INTERNAL":
                    names = tools if match.group(1) == "TOOLS" else actions
                    names.add(match.group(2))
    return actions, tools


def check_actions(source_tree):
    """All the actions and tools of the enabled modules are available."""
    actions, tools = documented_names(source_tree, enabled_modules())
    assert len(actions) > 100, actions
    available = set(plumed.InputBuilder()._vimdict)
    missing = sorted(actions - available)
    assert not missing, f"actions missing from the kernel: {missing}"
    usage = subprocess.run(
        ["plumed", "--no-mpi", "--help"], check=True, capture_output=True, text=True
    ).stdout
    missing = sorted(
        tool
        for tool in tools
        if not re.search(rf"^\s*plumed {re.escape(tool)} ", usage, re.M)
    )
    assert not missing, f"tools missing from the plumed command: {missing}"
    print(f"{len(actions)} actions and {len(tools)} tools available")


def check_external_kernel(bundled):
    """InputBuilder reads the metadata from an explicitly selected new kernel
    located elsewhere (a copy of the bundled one)."""
    with tempfile.TemporaryDirectory() as temporary:
        kernel = str(Path(temporary) / ("external-" + bundled.name))
        shutil.copy2(bundled, kernel)
        assert plumed._readkernelmetadata(kernel=kernel) is not None
        builder = plumed.InputBuilder(kernel=kernel)
        assert builder._vimdict["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"


def kernel_json_copy(plumed_instance, length):
    """Copy the last JSON snapshot using the kernel's text-transfer command."""
    buffer = array.array("i", [0]) * length
    plumed_instance.cmd("getInputBuilderJSON", buffer)
    return array.array("B", buffer).tobytes().decode("utf-8")


def kernel_json(plumed_instance):
    """Query the current JSON length and copy its matching snapshot."""
    length = array.array("i", [0])
    plumed_instance.cmd("getInputBuilderJSONLength", length)
    return kernel_json_copy(plumed_instance, length[0])


def check_metadata_transfer_errors(bundled):
    """The transfer requires a length query and rejects undersized buffers."""
    with plumed.Plumed(kernel=str(bundled)) as instance:
        try:
            instance.cmd("getInputBuilderJSON", array.array("i", [0]) * 8)
        except RuntimeError as error:
            assert "getInputBuilderJSONLength" in str(error)
        else:
            raise AssertionError("metadata copy worked before a length query")

        length = array.array("i", [0])
        instance.cmd("getInputBuilderJSONLength", length)
        try:
            too_small = array.array("i", [0]) * (length[0] - 1)
            instance.cmd("getInputBuilderJSON", too_small)
        except RuntimeError as error:
            assert "too small" in str(error)
        else:
            raise AssertionError("metadata copy accepted an undersized buffer")


def check_vim_generation(executable):
    """Generated Vim help and dictionary retain the legacy manual format."""
    actions = sorted(plumed.InputBuilder()._vimdict)
    with tempfile.TemporaryDirectory() as temporary:
        output = Path(temporary) / "vim-output"
        subprocess.run(
            [str(executable), "--no-mpi", "gen_vim", "--outdir", str(output)],
            check=True,
            capture_output=True,
            cwd=temporary,
        )

        syntax = (output / "syntax" / "plumed.vim").read_text(encoding="utf-8")
        generated_dictionary = {
            match.group(1): line
            for line in syntax.splitlines()
            if (match := re.match(r'let b:plumedDictionary\["(.*)"\]=\[.*\]$', line))
        }
        assert set(generated_dictionary) == set(actions)

        for action in actions:
            manual = subprocess.run(
                [
                    str(executable),
                    "--no-mpi",
                    "manual",
                    "--action",
                    action,
                    "--vim",
                ],
                check=True,
                capture_output=True,
            ).stdout
            first_line, separator, help_text = manual.partition(b"\n")
            assert separator, f"manual output for {action} has no first line"
            tokens = first_line.decode("utf-8").split(",")
            assert tokens[0] == action, (action, tokens[0])

            rows = ['{"word":"LABEL=","menu":"(label)"}']
            for token in tokens[1:]:
                kind, name = token.split(":", 1)
                if name == "LABEL":
                    continue
                if kind == "flag":
                    word, menu = name, "(flag)"
                elif kind == "numbered":
                    word, menu = name, "(numbered)"
                else:
                    word, menu = name + "=", "(option)"
                rows.append('{"word":"%s","menu":"%s"}' % (word, menu))
            expected_dictionary = 'let b:plumedDictionary["%s"]=[%s]' % (
                action,
                ",".join(sorted(rows)),
            )
            assert generated_dictionary[action] == expected_dictionary, action

            expected_help = (
                b"****************************************\n"
                + f"Short helpfile for action {action}\n".encode("utf-8")
                + b"****************************************\n"
                + help_text
            )
            generated_help = (output / "help" / f"{action}.txt").read_bytes()
            assert generated_help == expected_help, action

        help_dir = output / "help"
        obsolete = help_dir / "OBSOLETE_TEST_ACTION.txt"
        obsolete.write_text("stale\n", encoding="utf-8")
        subprocess.run(
            [str(executable), "--no-mpi", "gen_vim", "--outdir", str(output)],
            check=True,
            capture_output=True,
            cwd=temporary,
        )
        assert not obsolete.exists()

        previous_files = {
            path.relative_to(output): path.read_bytes()
            for path in output.rglob("*")
            if path.is_file()
        }
        empty_templates = Path(temporary) / "missing-templates"
        empty_templates.mkdir()
        failed = subprocess.run(
            [
                str(executable),
                "--no-mpi",
                "gen_vim",
                "--outdir",
                str(output),
                "--template-dir",
                str(empty_templates),
            ],
            capture_output=True,
            cwd=temporary,
        )
        assert failed.returncode != 0
        current_files = {
            path.relative_to(output): path.read_bytes()
            for path in output.rglob("*")
            if path.is_file()
        }
        assert current_files == previous_files

        # Template validation happens before replacing either output directory.
        templates = resource_root() / "vim" / "templates"
        shutil.copy2(templates / "plumedf.vim", empty_templates / "plumedf.vim")
        original_template = (templates / "plumed.vim.in").read_text(encoding="utf-8")
        marker = "@PLUMED_DICTIONARY@"
        for replacement in ("", marker + "\n" + marker, "prefix" + marker):
            (empty_templates / "plumed.vim.in").write_text(
                original_template.replace(marker, replacement), encoding="utf-8"
            )
            failed = subprocess.run(
                [
                    str(executable),
                    "--no-mpi",
                    "gen_vim",
                    "--outdir",
                    str(output),
                    "--template-dir",
                    str(empty_templates),
                ],
                capture_output=True,
                cwd=temporary,
            )
            assert failed.returncode != 0
            assert {
                path.relative_to(output): path.read_bytes()
                for path in output.rglob("*")
                if path.is_file()
            } == previous_files


def check_json_matches_vim():
    """JSON-derived InputBuilder classifications match Vim, except LABEL."""
    from plumed_runtime import resource_root

    json_dict, json_doc = plumed._readkernelmetadata()
    vim_dict, vim_doc = plumed._readvimdict(plumedroot=str(resource_root()))
    assert json_doc == vim_doc
    # Identical generated source checks signatures/defaults, argument handling,
    # formatted-input construction and docstrings for every action and alias mode.
    for underscores in (False, True):
        assert plumed._create_functions(
            json_dict, doc=json_doc, append_underscores=underscores
        ) == plumed._create_functions(
            vim_dict, doc=vim_doc, append_underscores=underscores
        )
    assert set(json_dict) == set(vim_dict)
    for action, vim_keywords in vim_dict.items():
        expected = {
            name: kind for name, kind in vim_keywords.items() if name != "LABEL"
        }
        actual = {
            name: kind for name, kind in json_dict[action].items() if name != "LABEL"
        }
        assert actual == expected, (action, actual, expected)

    assert json_dict["ANGLE_SCALAR"]["ATOMS"] == "(option)"
    assert json_dict["ANGLE_VECTOR"]["ATOMS"] == "(numbered)"
    assert json_dict["DUMPPDB"]["FMT"] == vim_dict["DUMPPDB"]["FMT"]
    assert json_dict["DUMPMASSCHARGE"]["STRIDE"] == vim_dict["DUMPMASSCHARGE"]["STRIDE"]
    assert json_dict["DUMPMASSCHARGE"]["STRIDE"] == "(option)"


def check_cli_json_matches_kernel(bundled, executable):
    """The CLI and in-memory exports serialize the same registry identically."""
    with plumed.Plumed(kernel=str(bundled)) as instance:
        in_memory = kernel_json(instance)
    completed = subprocess.run(
        [str(executable), "--no-mpi", "gen_json"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout == in_memory

    # The legacy --actions override remains available to external callers.
    metadata = json.loads(in_memory)
    assert {"vimlink", "replicalink", "groups", "modules", "cltools"} <= set(metadata)
    for name in ("@mda:", "@mdt:", "@vmd:", "@vmdexec:"):
        description = metadata["groups"][name]["description"]
        assert 'href="https://' in description
        assert '\\"' not in description
    with tempfile.TemporaryDirectory() as temporary:
        actions = Path(temporary) / "action_list"
        actions.write_text(
            "".join(
                f"{name}:{entry.get('description', '')}\n"
                for name, entry in metadata.items()
                if isinstance(entry, dict) and "syntax" in entry
            ),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [str(executable), "--no-mpi", "gen_json", "--actions", str(actions)],
            check=True,
            capture_output=True,
            text=True,
        )
    assert completed.stdout == in_memory


def check_vim_fallback():
    """Emulate an older kernel without embedded metadata: InputBuilder must use
    the Vim files of the kernel installation, consistently with the metadata."""
    reference = plumed.InputBuilder()._vimdict
    original = plumed._readkernelmetadata
    plumed._readkernelmetadata = lambda kernel=None, instance=None: None
    try:
        legacy = plumed.InputBuilder()._vimdict
    finally:
        plumed._readkernelmetadata = original
    assert legacy is not reference
    assert set(legacy) == set(reference)
    for action, keywords in reference.items():
        # the Vim files also list the hidden keywords
        missing = set(keywords) - set(legacy[action])
        assert not missing, (action, missing)
    assert legacy["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"


def initialize_builder_test_kernel(instance):
    instance.cmd("setMDEngine", "python")
    instance.cmd("setNatoms", 2)
    instance.cmd("init")


def check_builder_load(toplumed, *, append_underscores=False):
    """A LOAD sent through an attached builder refreshes its action methods."""
    builder = plumed.InputBuilder(
        toplumed=toplumed, append_underscores=append_underscores
    )
    assert "WHEELTEST" not in builder._vimdict
    load_name = "LOAD__" if append_underscores else "LOAD"
    load_line = getattr(builder, load_name)(FILE="WheelTest.cpp")
    assert load_line.startswith("LOAD FILE=WheelTest.cpp")
    assert "WHEELTEST" in builder._vimdict
    method_name = "WHEELTEST__" if append_underscores else "WHEELTEST"
    method = getattr(builder, method_name)
    assert method.__doc__ and "Create action WHEELTEST" in method.__doc__
    assert "WHEELTEST" in method(ATOMS="1,2")
    assert json.loads(kernel_json(toplumed)).get("WHEELTEST") is not None


def check_inputbuilder_load():
    """Builders use their attached instance's live, isolated action registry."""
    with plumed.Plumed() as auto_loaded:
        initialize_builder_test_kernel(auto_loaded)
        check_builder_load(auto_loaded)
        # The builder borrows this object; metadata queries still work afterward.
        assert json.loads(kernel_json(auto_loaded)).get("WHEELTEST") is not None

    with plumed.Plumed() as loaded_first:
        initialize_builder_test_kernel(loaded_first)
        loaded_first.cmd("readInputLine", "LOAD FILE=WheelTest.cpp")
        builder = plumed.InputBuilder(toplumed=loaded_first)
        assert "WHEELTEST" in builder._vimdict
        assert "WHEELTEST" in builder.WHEELTEST(ATOMS="1,2")

    with plumed.Plumed() as auto_loaded_underscored:
        initialize_builder_test_kernel(auto_loaded_underscored)
        check_builder_load(auto_loaded_underscored, append_underscores=True)

    with plumed.Plumed() as failed_load:
        initialize_builder_test_kernel(failed_load)
        builder = plumed.InputBuilder(toplumed=failed_load)
        before_dictionary = {
            action: keywords.copy() for action, keywords in builder._vimdict.items()
        }
        before_doc = builder._doc.copy()
        before_functions = builder._functions.copy()
        failed_line = "LOAD FILE=definitely_missing_plugin_for_test.so"
        try:
            builder.LOAD(FILE="definitely_missing_plugin_for_test.so")
        except RuntimeError:
            pass
        else:
            raise AssertionError("LOAD accepted a missing plugin")
        assert builder.last_failure == failed_line
        assert builder._vimdict == before_dictionary
        assert builder._doc == before_doc
        assert builder._functions == before_functions

    with plumed.Plumed() as verbatim_loaded:
        initialize_builder_test_kernel(verbatim_loaded)
        builder = plumed.InputBuilder(toplumed=verbatim_loaded)
        builder.verbatim("LOAD FILE=WheelTest.cpp")
        assert "WHEELTEST" in builder._vimdict
        assert "WHEELTEST" in builder.WHEELTEST(ATOMS="1,2")

    with plumed.Plumed() as generic_loaded:
        initialize_builder_test_kernel(generic_loaded)
        builder = plumed.InputBuilder(toplumed=generic_loaded)
        builder("LOAD", FILE="WheelTest.cpp")
        assert "WHEELTEST" in builder._vimdict
        assert "WHEELTEST" in builder.WHEELTEST(ATOMS="1,2")

    with plumed.Plumed() as refreshed:
        initialize_builder_test_kernel(refreshed)
        builder = plumed.InputBuilder(toplumed=refreshed)
        assert "WHEELTEST" not in builder._vimdict
        refreshed.cmd("readInputLine", "LOAD FILE=WheelTest.cpp")
        assert "WHEELTEST" not in builder._vimdict
        assert builder.refresh() is True
        method = builder.WHEELTEST
        assert method.__doc__ and "Create action WHEELTEST" in method.__doc__
        assert "WHEELTEST" in method(ATOMS="1,2")
        assert json.loads(kernel_json(refreshed)).get("WHEELTEST") is not None

        with plumed.Plumed() as independent:
            initialize_builder_test_kernel(independent)
            independent_builder = plumed.InputBuilder(toplumed=independent)
            assert "WHEELTEST" not in independent_builder._vimdict

    with plumed.Plumed() as no_dictionary:
        initialize_builder_test_kernel(no_dictionary)
        builder = plumed.InputBuilder(toplumed=no_dictionary, load_dict=False)
        assert builder.refresh() is False
        no_dictionary.cmd("readInputLine", "LOAD FILE=WheelTest.cpp")
        assert builder.refresh() is False
        assert "WHEELTEST" not in builder.__dict__.get("_vimdict", {})
        assert builder.WHEELTEST("d", ATOMS="1,2") == "d: WHEELTEST ATOMS=1,2\n"


def check_load():
    """Compile an action at runtime with LOAD, as users do with their own CVs."""
    options = (resource_root() / "src" / "config" / "compile_options.sh").read_text()
    compiler = shlex.split(re.search(r"^compile='(.*?) -c ", options, re.M).group(1))
    if shutil.which(compiler[0]) is None:
        if os.environ.get("PLUMED_WHEEL_TEST_NO_COMPILER"):
            print(f"skipping the LOAD test: {compiler[0]} not found")
            return
        raise AssertionError(
            f"{compiler[0]} not found: set PLUMED_WHEEL_TEST_NO_COMPILER to skip"
        )
    template = (resource_root() / "src" / "colvar" / "Template.cpp").read_text()
    # out of the source tree, headers are included as "module/header.h"
    code = template.replace('"Colvar.h"', '"colvar/Colvar.h"')
    code = code.replace("TEMPLATE", "WHEELTEST")
    current = os.getcwd()
    with tempfile.TemporaryDirectory() as temporary:
        os.chdir(temporary)
        try:
            Path("WheelTest.cpp").write_text(code)
            check_inputbuilder_load()
            # Metadata is available before initialization and refreshes after
            # LOAD registers another action on the same PLUMED instance.
            with plumed.Plumed() as instance:
                length = array.array("i", [0])
                instance.cmd("getInputBuilderJSONLength", length)
                before_load = kernel_json_copy(instance, length[0])
                instance.cmd("setNatoms", 2)
                instance.cmd("init")
                instance.cmd("readInputLine", "LOAD FILE=WheelTest.cpp")
                stale_snapshot = kernel_json_copy(instance, length[0])
                assert stale_snapshot == before_load
                metadata = json.loads(kernel_json(instance))
                assert "WHEELTEST" in metadata
                with plumed.Plumed() as fresh_instance:
                    fresh_metadata = json.loads(kernel_json(fresh_instance))
                    assert "WHEELTEST" not in fresh_metadata
            check_kernel(lines=("LOAD FILE=WheelTest.cpp", "d: WHEELTEST ATOMS=1,2"))
        finally:
            os.chdir(current)


def check_dynamic_cltool_module(executable):
    """A dynamically loaded CLTool alone still has a defined module entry."""
    options = (resource_root() / "src" / "config" / "compile_options.sh").read_text()
    compiler = shlex.split(re.search(r"^compile='(.*?) -c ", options, re.M).group(1))
    if shutil.which(compiler[0]) is None:
        if os.environ.get("PLUMED_WHEEL_TEST_NO_COMPILER"):
            print(f"skipping the dynamic CLTool test: {compiler[0]} not found")
            return
        raise AssertionError(
            f"{compiler[0]} not found: set PLUMED_WHEEL_TEST_NO_COMPILER to skip"
        )

    source = r"""#include "core/CLTool.h"
#include "core/CLToolRegister.h"

namespace PLMD {
class WheelOnlyCLTool : public CLTool {
public:
  static void registerKeywords(Keywords& keys) {
    CLTool::registerKeywords(keys);
  }
  explicit WheelOnlyCLTool(const CLToolOptions& options): CLTool(options) {
    inputdata=inputType::commandline;
  }
  int main(FILE*,FILE*,Communicator&) override {
    return 0;
  }
  std::string description() const override {
    return "wheel test tool without a registered action";
  }
};

PLUMED_REGISTER_CLTOOL(WheelOnlyCLTool,"wheel_only_tool")
}
"""
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        plugin_source = directory / "WheelOnlyCLTool.cpp"
        plugin = directory / "wheel_only_tool.so"
        plugin_source.write_text(source, encoding="utf-8")
        subprocess.run(
            [
                str(executable),
                "--no-mpi",
                "mklib",
                "--out",
                plugin.name,
                plugin_source.name,
            ],
            check=True,
            capture_output=True,
            cwd=directory,
        )
        result = subprocess.run(
            [str(executable), "--no-mpi", "--load", str(plugin), "gen_json"],
            check=True,
            capture_output=True,
            text=True,
            cwd=directory,
        )
    metadata = json.loads(result.stdout)
    assert "wheel_only_tool" in metadata["cltools"]
    assert metadata["cltools"]["wheel_only_tool"]["module"] == "external"
    assert "external" in metadata["modules"]
    assert "WHEEL_ONLY_TOOL" not in metadata


def main():
    os.environ.pop("PLUMED_KERNEL", None)
    os.environ.pop("PLUMED_ROOT", None)
    bundled = kernel_path()
    assert bundled is not None and bundled.is_file()
    check_kernel()
    check_kernel(str(bundled))

    builder = plumed.InputBuilder()
    assert "RESTRAINT" in builder._vimdict
    assert "KAPPA" in builder._vimdict["RESTRAINT"]
    assert builder._vimdict["MOVINGRESTRAINT"]["KAPPA"] == "(numbered)"
    assert "DISTANCE" in builder._vimdict

    with tempfile.TemporaryDirectory() as temporary:
        override = Path(temporary) / bundled.name
        shutil.copy2(bundled, override)
        os.environ["PLUMED_KERNEL"] = str(override)
        check_kernel()
        other = Path(temporary) / ("second-" + bundled.name)
        shutil.copy2(bundled, other)
        os.environ["PLUMED_KERNEL"] = str(other)
        check_kernel()
        try:
            plumed.Plumed(kernel=str(Path(temporary) / "missing-kernel"))
        except RuntimeError as error:
            assert "Error loading PLUMED kernel" in str(error)
        else:
            raise AssertionError("an invalid explicit kernel path was accepted")

    os.environ.pop("PLUMED_KERNEL", None)
    command = subprocess.run(
        ["plumed", "--no-mpi", "info", "--root"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert Path(command.stdout.strip()).resolve() == resource_root().resolve()

    check_external_kernel(bundled)
    check_json_matches_vim()
    check_cli_json_matches_kernel(bundled, package_root() / "bin" / "plumed")
    check_metadata_transfer_errors(bundled)
    check_vim_generation(package_root() / "bin" / "plumed")
    check_vim_fallback()
    check_load()
    check_dynamic_cltool_module(package_root() / "bin" / "plumed")

    source_tree = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None
    hidden_source = (
        source_tree.with_name(source_tree.name + ".hidden") if source_tree else None
    )
    vim_tree = resource_root() / "vim"
    hidden_vim = vim_tree.with_name("vim.hidden")
    json_tree = resource_root() / "json"
    hidden_json = json_tree.with_name("json.hidden")
    executable = package_root() / "bin" / "plumed"
    hidden_executable = executable.with_name("plumed.hidden")
    if source_tree and source_tree.is_dir():
        check_actions(source_tree)
    with tempfile.TemporaryDirectory() as empty_path:
        original_path = os.environ.get("PATH", "")
        if source_tree and source_tree.is_dir():
            source_tree.rename(hidden_source)
        if vim_tree.exists():
            vim_tree.rename(hidden_vim)
        if json_tree.exists():
            json_tree.rename(hidden_json)
        if executable.exists():
            executable.rename(hidden_executable)
        os.environ["PATH"] = empty_path
        try:
            metadata = plumed._readkernelmetadata()
            assert metadata is not None
            isolated = plumed.InputBuilder()
            assert "RESTRAINT" in isolated._vimdict
        finally:
            os.environ["PATH"] = original_path
            if hidden_executable.exists():
                hidden_executable.rename(executable)
            if hidden_vim.exists():
                hidden_vim.rename(vim_tree)
            if hidden_json.exists():
                hidden_json.rename(json_tree)
            if hidden_source and hidden_source.exists():
                hidden_source.rename(source_tree)

    forbidden = ("mpi", "gomp", "openblas", "blas", "lapack", "fftw")
    if sys.platform.startswith("linux"):
        dependencies = subprocess.check_output(["ldd", str(bundled)], text=True)
    elif sys.platform == "darwin":
        dependencies = subprocess.check_output(["otool", "-L", str(bundled)], text=True)
    else:
        dependencies = ""
    lowered = dependencies.lower()
    assert not any(name in lowered for name in forbidden), dependencies
    print((resource_root() / "build-info.txt").read_text())


if __name__ == "__main__":
    main()
