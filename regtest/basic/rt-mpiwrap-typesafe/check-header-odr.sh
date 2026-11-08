#!/usr/bin/env bash
# MPI visibility must never alter shared public wrapper definitions.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-header-odr.XXXXXX")
trap 'rm -rf "$task"' EXIT
cxx=${CXX:-c++}
cat > "$task/a.cpp" <<'CPP'
#include "mpi.h"
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#include "Plumed.h"
extern "C" void probe_a(plumed p, double* x) { PLMD::Plumed::plumed_cmd_cxx(p,"a",x); }
CPP
cat > "$task/b.cpp" <<'CPP'
#ifdef ALL_MPI_AWARE
#include "mpi.h"
#endif
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#include "Plumed.h"
extern "C" void probe_b(plumed p, double* x) { PLMD::Plumed::plumed_cmd_cxx(p,"b",x); }
CPP
cat > "$task/main.cpp" <<'CPP'
#ifdef ALL_MPI_AWARE
#include "mpi.h"
#endif
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#include "Plumed.h"
#include <cassert>
#include <cstdlib>
extern "C" void probe_a(plumed,double*);
extern "C" void probe_b(plumed,double*);
extern "C" void* plumed_malloc(size_t n) { return std::malloc(n); }
extern "C" void plumed_free(void* p) { std::free(p); }
extern "C" void plumed_cmd_safe_nothrow(plumed,const char* key,plumed_safeptr safe,plumed_nothrow_handler) {
  const bool expected = __PLUMED_WRAPPER_CXX_MPI_COMM_CHECK && (key[0]=='a'
#ifdef ALL_MPI_AWARE
    || key[0]=='b'
#endif
  );
  assert(!!(safe.flags&0x20000000u)==expected);
}
int main() { plumed p={0}; double x=0; probe_a(p,&x); probe_b(p,&x); }
CPP
# Missing declarations must be an error for explicitly checked shared wrappers.
if "$cxx" -std=c++98 -D__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK=1 -I"$root/src/wrapper" \
    -fsyntax-only "$task/b.cpp" >"$task/no-mpi.log" 2>&1; then
  echo 'shared checked wrapper accepted missing MPI declarations' >&2; exit 1
fi
grep -q 'requires MPI declarations' "$task/no-mpi.log"
# Disabling the C++ wrapper must also disable its MPI prerequisites.
printf '%s\n' '#define __PLUMED_WRAPPER_IMPLEMENTATION 0' '#include "Plumed.h"' \
  'void probe(plumed p, double* x) { plumed_cmd(p,"x",x); }' > "$task/no-cxx.cpp"
for standard in c++98 c++17; do
  "$cxx" -std="$standard" -Wall -Wextra -Werror -D__PLUMED_WRAPPER_CXX=0 \
    -D__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK=1 -I"$root/src/wrapper" \
    -x c++ -fsyntax-only "$task/no-cxx.cpp"
done
# Both absent declarations and explicit disable must compile in C11.
printf '%s\n' '#define __PLUMED_WRAPPER_IMPLEMENTATION 0' '#include "Plumed.h"' \
  'void probe(plumed p, double* x) { plumed_cmd(p,"x",x); }' > "$task/no-mpi.c"
"${CC:-cc}" -std=c11 -Werror -I"$root/src/wrapper" -fsyntax-only "$task/no-mpi.c"
# Test LTO only if the compiler/linker support it; actual tests then must pass.
printf 'int main() { return 0; }\n' > "$task/lto.cpp"
optimizations=(-O0)
if "$cxx" -O2 -flto "$task/lto.cpp" -o "$task/lto" >/dev/null 2>&1; then
  optimizations+=("-O2 -flto")
fi
for representation in 'int' 'struct opaque_comm *' 'void *' 'struct WrapperComm { int value; }'; do
  printf '#define MPI_VERSION 5\ntypedef %s MPI_Comm;\n' "$representation" > "$task/mpi.h"
  for mode in shared anonymous checked; do
    options=
    case $mode in
      anonymous) options=-D__PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE=1 ;;
      checked) options='-D__PLUMED_WRAPPER_CXX_MPI_COMM_CHECK=1 -DALL_MPI_AWARE' ;;
    esac
    for optimize in "${optimizations[@]}"; do
      for unit in a b main; do
        "$cxx" -std=c++98 -Wall -Wextra -Werror $options $optimize -I"$task" -I"$root/src/wrapper" \
          -c "$task/$unit.cpp" -o "$task/$unit.o"
      done
      for order in 'a b' 'b a'; do
        objects=()
        for unit in $order main; do objects+=("$task/$unit.o"); done
        "$cxx" $optimize "${objects[@]}" -o "$task/check"
        "$task/check"
      done
    done
  done
done
echo 'mixed MPI-aware/unaware wrapper translation units passed, including supported LTO'
