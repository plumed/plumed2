#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-extension-test.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
cp "$kernel" "$task/kernel-before"
# Compatibility packages contain no plugin: the module moved into the kernel.
conda install -y -p "$prefix" --override-channels -c local -c conda-forge "plumed-ext-opes=$VERSION"
"$prefix/bin/plumed" config module opes -q
test ! -e "$prefix/lib/plumed/plugins/opes.json"
test ! -e "$prefix/lib/plumed/plugins/opes.plumed"
cmp "$kernel" "$task/kernel-before"
! "$prefix/bin/plumed" modules -q --has ves
conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
  "plumed-ext-ves=$VERSION" "plumed-ext-pamm=$VERSION"
cmp "$kernel" "$task/kernel-before"
# Output-specific host requirements must keep unrelated native stacks out of
# lightweight extensions, even though all modules share one build environment.
for record in "$prefix"/conda-meta/{libtorch,libboost,libmetatomic-torch}-*.json; do
  test ! -e "$record"
done
python extensions/tests/chemfiles.py --prefix "$prefix"
python extensions/tests/abi.py --prefix "$prefix"
python extensions/tests/kernel_requirements.py
python extensions/tests/requirements.py
python extensions/tests/module_features.py --base "$prefix"
python extensions/tests/installed.py --base "$prefix" --extensions "$prefix"
# The public C wrapper permits a host compiler to exercise several independent
# objects against the same installed kernel, including concurrent activation.
read -r -a host_compiler <<< "${CXX:-c++}"
"${host_compiler[@]}" -std=c++17 -pthread -I"$prefix/include" \
  extensions/tests/visibility.cpp "$prefix/lib/libplumedWrapper.a" \
  -L"$prefix/lib" -Wl,-rpath,"$prefix/lib" -ldl -o "$task/visibility"
PLUMED_AUTOLOAD=0 PLUMED_KERNEL="$kernel" \
  PLUMED_PLUGIN_PATH="$prefix/lib/plumed/plugins" "$task/visibility"

# Available conda libraries must actually enable their module features; merely
# loading an identity-only library would miss broken include/link probes.
read -r -a regression_packages <<< "$(python extensions/tests/regressions.py --list-packages)"
for index in "${!regression_packages[@]}"; do
  regression_packages[index]="${regression_packages[index]}=$VERSION"
done
conda install -y -p "$prefix" --override-channels -c local -c conda-forge \
  "plumed-ext-drr=$VERSION" "plumed-ext-fourier=$VERSION" \
  "plumed-ext-pytorch=$VERSION" "plumed-ext-metatomic=$VERSION" \
  "plumed-ext-ttsketch=$VERSION" "plumed-ext-kenref=$VERSION" \
  "${regression_packages[@]}"
python - "$prefix" <<'PYTHON'
import json
import sys
from pathlib import Path
reports = Path(sys.argv[1]) / "lib/plumed/plugins"
for module, feature in [
    ("drr", "BOOST_SERIALIZATION"), ("fourier", "FFTW"),
    ("isdb-pytorch", "LIBTORCH"), ("pytorch", "LIBTORCH"), ("metatomic", "LIBMETATOMIC"),
    ("metatomic", "LIBTORCH"), ("ttsketch", "LIBHDF5"),
]:
    report = json.loads((reports / (module + ".json")).read_text())
    assert report["optional_features"][feature], report
PYTHON
"$prefix/bin/plumed" modules -q --has drr fourier pytorch metatomic ttsketch kenref
cmp "$kernel" "$task/kernel-before"
# Runtime features follow the selected module closure; config stays static.
! "$prefix/bin/plumed" config -q has libtorch
! PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" modules -q --has-feature libtorch
PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=pytorch modules -q --has-feature libtorch
! PLUMED_AUTOLOAD=0 "$prefix/bin/plumed" --load-module=pytorch modules -q --has-feature libmetatomic
python extensions/tests/isdb_pytorch.py --prefix "$prefix"
python extensions/tests/regressions.py --prefix "$prefix"

conda remove -y -p "$prefix" --override-channels -c local -c conda-forge plumed-ext-ves
cmp "$kernel" "$task/kernel-before"
! "$prefix/bin/plumed" modules -q --has ves
"$prefix/bin/plumed" modules -q --has pamm adjmat
# Verify removal leaves all base functionality intact.
"$prefix/bin/plumed" --has-thread-mpi
"$prefix/bin/plumed" modules -q --has bias
