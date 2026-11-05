#!/usr/bin/env bash
# Build only the provider adapter: never link a second PLUMED kernel.
set -euo pipefail

MANUAL='Build an MPI provider for a serial-built PLUMED kernel.

Usage: plumed mkmpi --compiler /absolute/path/to/mpicc --out provider.so

Options:
  --compiler PATH        Absolute path to an MPI C or C++ compiler driver
  -o, --out PATH         Output shared library (required)
  --cflag ARG            Add one compiler argument; may be repeated
  --ldflag ARG           Add one linker argument; may be repeated
  -h, --help             Show this help

The adapter is C++98-compatible. The .cpp suffix selects the language; no
standard or language-selection flags are added. Configured shared-library
flags are reused. Extra flags are individual arguments without shell evaluation.
The generated provider uses the selected MPI runtime and no PLUMED library.'

compiler=
output=
cflags=()
ldflags=()
while (($#)); do
  case "$1" in
    --description) echo 'build an MPI provider using an MPI compiler'; exit 0 ;;
    --options) echo '--compiler -o --out --cflag --ldflag -h --help'; exit 0 ;;
    -h|--help) echo "$MANUAL"; exit 0 ;;
    --compiler|-o|--out|--cflag|--ldflag)
      option=$1
      if (($# < 2)); then
        echo "mkmpi: $option requires an argument" >&2
        exit 1
      fi
      case "$option" in
        --compiler) compiler=$2 ;;
        -o|--out) output=$2 ;;
        --cflag) cflags+=("$2") ;;
        --ldflag) ldflags+=("$2") ;;
      esac
      shift 2 ;;
    *) echo "mkmpi: unknown option $1 (see --help)" >&2; exit 1 ;;
  esac
done
if [[ "$compiler" != /* || ! -f "$compiler" || ! -x "$compiler" ]]; then
  echo 'mkmpi: --compiler must name an absolute executable MPI compiler path' >&2
  exit 1
fi
if [[ -z "$output" ]]; then
  echo 'mkmpi: --out is required' >&2
  exit 1
fi
plumed_mpi_root=${PLUMED_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
header=
header_candidates=()
if [[ -n "${PLUMED_INCLUDEDIR:-}" ]]; then
  header_candidates+=("$PLUMED_INCLUDEDIR/${PLUMED_PROGRAM_NAME:-plumed}/wrapper/Plumed.h")
fi
header_candidates+=("$plumed_mpi_root/src/wrapper/Plumed.h"
  "$plumed_mpi_root/src/include/${PLUMED_PROGRAM_NAME:-plumed}/Plumed.h")
for candidate in "${header_candidates[@]}"; do
  if [[ -r "$candidate" ]]; then header=$candidate; break; fi
done
if [[ -z "$header" ]]; then
  echo 'mkmpi: cannot find the public Plumed.h header; check PLUMED_INCLUDEDIR and PLUMED_PROGRAM_NAME' >&2
  exit 1
fi
# Reuse only flags checked by configure, as the mklib command does. The
# selected MPI driver replaces the configured compiler, not its MPI runtime.
options_file="$plumed_mpi_root/src/config/compile_options.sh"
if [[ ! -r "$options_file" ]]; then
  echo "mkmpi: missing configured build options $options_file" >&2
  exit 1
fi
source "$options_file"
if [[ -z "${mpi_provider_link_flags:-}" ]]; then
  echo 'mkmpi: shared-library compiler flags are unavailable; reconfigure PLUMED' >&2
  exit 1
fi
read -r -a configured_compile_flags <<< "${mpi_provider_compile_flags:-}"
read -r -a configured_link_flags <<< "$mpi_provider_link_flags"
# Use the destination filesystem so the final rename is atomic. Keeping the
# output basename also preserves the default library identity on macOS.
build_dir=$(mktemp -d "$(dirname "$output")/.plumed-mkmpi.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
temporary_output="$build_dir/$(basename "$output")"
adapter_source="$build_dir/mpiwrap_implementation.cpp"
cat > "$adapter_source" <<'ADAPTER'
#define PLUMED_MPIW_NATIVE_ADAPTER 1
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#define __PLUMED_WRAPPER_CXX 0
#include <Plumed.h>
#include <Plumed.h>
ADAPTER
platform_link_flags=()
if [[ $(uname -s) == Darwin ]]; then
  platform_link_flags+=("-Wl,-install_name,@rpath/$(basename "$output")")
fi
"$compiler" ${configured_compile_flags[@]+"${configured_compile_flags[@]}"} \
  -I"$(dirname "$header")" ${cflags[@]+"${cflags[@]}"} "$adapter_source" \
  "${configured_link_flags[@]}" "${platform_link_flags[@]}" \
  ${ldflags[@]+"${ldflags[@]}"} -o "$temporary_output" || {
    echo "mkmpi: MPI adapter compilation/link failed with $compiler" >&2
    exit 1
  }
if [[ ! -s "$temporary_output" ]]; then
  echo 'mkmpi: compiler did not produce a provider library' >&2
  exit 1
fi
mv -f "$temporary_output" "$output"
