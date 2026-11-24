#!/usr/bin/env bash
# TU-local host adapters may coexist for either wrapper linkage/namespace mode.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-host-linkage.XXXXXX")
trap 'rm -rf "$task"' EXIT
cc=${MPICC:-mpicc}
cxx=${MPICXX:-mpicxx}
cat > "$task/host.c" <<'C'
/* This probe uses the C MPI API; legacy C++ MPI bindings are unrelated. */
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#define __PLUMED_WRAPPER_MPI_HOST_ADAPTER 1
#define __PLUMED_WRAPPER_IMPLEMENTATION (!__PLUMED_WRAPPER_EXTERN)
#define __PLUMED_WRAPPER_GLOBAL 0
#include "Plumed.h"
#include "Plumed.h"
#ifdef __cplusplus
extern "C" {
#endif
const void* HOST_FUNCTION(void) { return plumed_mpi_implementation(); }
#ifdef __cplusplus
}
#endif
C
cat > "$task/main.c" <<'C'
#include <assert.h>
#ifdef __cplusplus
extern "C" {
#endif
const void* host_a(void);
const void* host_b(void);
#ifdef __cplusplus
}
#endif
int main(void) { assert(host_a()!=0); assert(host_b()!=0); return 0; }
C
for language in c99 c++98 c++11 c++17; do
  compiler=$cxx; source_language=c++
  if [[ $language = c99 ]]; then compiler=$cc; source_language=c; fi
  for anonymous in 0 1; do
    for external in 0 1; do
      for unit in a b; do
        "$compiler" -std="$language" -Wall -Wextra -Werror -O0 -I"$root/src/wrapper" \
          -D__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE=$anonymous -D__PLUMED_WRAPPER_EXTERN=$external \
          -D__PLUMED_HAS_DLOPEN=1 -DHOST_FUNCTION=host_$unit -x "$source_language" -c "$task/host.c" -o "$task/$unit.o"
        # Accessors and private adapter functions must never be exported.
        if nm -g "$task/$unit.o" | grep -E '[[:space:]]_?plumed_mpi_implementation$'; then
          echo 'host accessor unexpectedly exported' >&2; exit 1
        fi
      done
      "$compiler" -std="$language" -x "$source_language" -c "$task/main.c" -o "$task/main.o"
      "$compiler" "$task/a.o" "$task/b.o" "$task/main.o" -o "$task/check"
      "$task/check"
    done
  done
done
# Numeric zero disables both adapters and must not require MPI declarations.
printf '%s\n' '#define __PLUMED_WRAPPER_MPI_HOST_ADAPTER 0' \
  '#define __PLUMED_WRAPPER_MPI_PROVIDER_ADAPTER 0' '#include "Plumed.h"' > "$task/disabled.c"
"${CC:-cc}" -std=c99 -Wall -Wextra -Werror -I"$root/src/wrapper" -fsyntax-only "$task/disabled.c"
echo 'C99/C++98+ multi-TU host adapter linkage and numeric disable passed'
