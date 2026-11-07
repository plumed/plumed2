#!/usr/bin/env bash
# Build a provider adapter or recompile the standard main with host MPI.
set -euo pipefail

MANUAL='Build an MPI provider or MPI-enabled PLUMED executable.

Usage: plumed mkmpi [--compiler /absolute/path/to/mpicc] --out provider.so
       plumed mkmpi --executable [--link shared|runtime] --out plumed-mpi

Options:
  --compiler PATH        Absolute MPI compiler path (default: search below)
  --executable           Recompile the installed src/main/main.cpp with MPI
  --link MODE            Executable linking: shared (default) or runtime
  -o, --out PATH         Output library or executable (required)
  --cflag ARG            Add one compiler argument; may be repeated
  --ldflag ARG           Add one linker argument; may be repeated
  -h, --help             Show this help

Without --compiler, try $MPICC, $MPICXX, mpicc, then mpicxx, selecting the
first executable found in either mode. Environment values may be paths or
names on PATH.

The adapter and main are C99-compatible. The .c suffix selects the language; no
standard or language-selection flags are added. Configured shared-library
flags are reused. Extra flags are individual arguments without shell evaluation.
The generated provider uses the selected MPI runtime and no PLUMED library.

Executable mode uses the compiler default and the standard main, including
--mpi/--no-mpi. The same source remains usable by C++ compiler drivers.
It injects host MPI into provider-capable kernels; native-MPI kernels ignore
injection. Shared mode links the current libplumed. Runtime mode embeds the
public loader; PLUMED_KERNEL overrides the default kernel selected at build
time. Use PLUMED_MPI_WRAPPER_FORCE=1 together with PLUMED_MPI_WRAPPER to test
an environment-selected provider instead of injection.'

executable=no
link_mode=shared
link_given=no
compiler=
compiler_given=no
output=
cflags=()
ldflags=()
while (($#)); do
  case "$1" in
    --description) echo 'build an MPI provider or MPI-enabled executable'; exit 0 ;;
    --options) echo '--compiler --executable --link -o --out --cflag --ldflag -h --help'; exit 0 ;;
    -h|--help) echo "$MANUAL"; exit 0 ;;
    --executable) executable=yes; shift ;;
    --compiler|--link|-o|--out|--cflag|--ldflag)
      option=$1
      if (($# < 2)); then
        echo "mkmpi: $option requires an argument" >&2
        exit 1
      fi
      case "$option" in
        --compiler) compiler=$2; compiler_given=yes ;;
        --link) link_mode=$2; link_given=yes ;;
        -o|--out) output=$2 ;;
        --cflag) cflags+=("$2") ;;
        --ldflag) ldflags+=("$2") ;;
      esac
      shift 2 ;;
    *) echo "mkmpi: unknown option $1 (see --help)" >&2; exit 1 ;;
  esac
done
if [[ "$link_mode" != shared && "$link_mode" != runtime ]]; then
  echo 'mkmpi: --link must be shared or runtime' >&2
  exit 1
fi
if [[ "$link_given" = yes && "$executable" != yes ]]; then
  echo 'mkmpi: --link requires --executable' >&2
  exit 1
fi
if [[ "$compiler_given" = no ]]; then
  candidates=("${MPICC:-}" "${MPICXX:-}" mpicc mpicxx)
  for candidate in "${candidates[@]}"; do
    [[ -n "$candidate" ]] || continue
    candidate=$(type -P -- "$candidate") || continue
    if [[ -f "$candidate" && -x "$candidate" ]]; then
      # Keep the wrapper's basename (including symlinks), but resolve relative
      # paths before compiling into a temporary directory.
      compiler="$(cd "$(dirname "$candidate")" && pwd)/$(basename "$candidate")"
      break
    fi
  done
  if [[ -z "$compiler" ]]; then
    echo 'mkmpi: no MPI compiler found; set MPICC or MPICXX, or use --compiler' >&2
    exit 1
  fi
fi
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
if [[ "$executable" != yes && -z "${mpi_provider_link_flags:-}" ]]; then
  echo 'mkmpi: shared-library compiler flags are unavailable; reconfigure PLUMED' >&2
  exit 1
fi
read -r -a configured_compile_flags <<< "${mpi_provider_compile_flags:-}"
read -r -a configured_link_flags <<< "${mpi_provider_link_flags:-}"
# Use the destination filesystem so the final rename is atomic. Keeping the
# output basename also preserves the default library identity on macOS.
build_dir=$(mktemp -d "$(dirname "$output")/.plumed-mkmpi.XXXXXX")
trap 'rm -rf "$build_dir"' EXIT
temporary_output="$build_dir/$(basename "$output")"
if [[ "$executable" = yes ]]; then
  main_source="$plumed_mpi_root/src/main/main.cpp"
  if [[ ! -r "$main_source" ]]; then
    echo "mkmpi: missing standard main source $main_source" >&2
    exit 1
  fi
  # Compile identical source as C without requiring compiler-specific -x flags.
  cp "$main_source" "$build_dir/main.c"
  program_name=${PLUMED_PROGRAM_NAME:-plumed}
  library_dir=
  for candidate in "$plumed_mpi_root/src/lib" "$(dirname "$plumed_mpi_root")"; do
    if [[ -r "$candidate/lib${program_name}Kernel.${soext:-so}" ]]; then
      library_dir=$candidate
      break
    fi
  done
  if [[ -z "$library_dir" ]]; then
    echo 'mkmpi: cannot locate the current PLUMED kernel library' >&2
    exit 1
  fi
  main_flags=(-D__PLUMED_HAS_MPI=1)
  main_link_flags=()
  if [[ "$link_mode" = shared ]]; then
    if [[ ! -r "$library_dir/lib${program_name}.${soext:-so}" ]]; then
      echo 'mkmpi: libplumed is unavailable; use --link runtime' >&2
      exit 1
    fi
    main_link_flags=(-L"$library_dir" "-Wl,-rpath,$library_dir" "-l$program_name")
  else
    main_flags+=(-D__PLUMED_WRAPPER_IMPLEMENTATION=1
      -D__PLUMED_WRAPPER_LINK_RUNTIME=1 -D__PLUMED_HAS_DLOPEN=1
      "-D__PLUMED_DEFAULT_KERNEL=$library_dir/lib${program_name}Kernel.${soext:-so}")
    main_link_flags=(-ldl)
  fi
  "$compiler" "${main_flags[@]}" -I"$(dirname "$header")/.." \
    ${cflags[@]+"${cflags[@]}"} "$build_dir/main.c" "${main_link_flags[@]}" \
    ${ldflags[@]+"${ldflags[@]}"} -o "$temporary_output" || {
      echo "mkmpi: executable compilation/link failed with $compiler" >&2
      exit 1
    }
  if [[ ! -s "$temporary_output" ]]; then
    echo 'mkmpi: compiler did not produce an executable' >&2
    exit 1
  fi
  mv -f "$temporary_output" "$output"
  exit 0
fi
adapter_source="$build_dir/mpiwrap_implementation.c"
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
