#!/usr/bin/env bash
# Exercise the wrapper-include exception and state reset between input files.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-wrapper-lint.XXXXXX")
trap 'rm -rf "$task"' EXIT
mkdir -p "$task/probe"
printf 'show_used_modules:\n\t@echo used_modules\n' > "$task/probe/Makefile"
cat > "$task/probe/good.cpp" <<'CPP'
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_EXTERN 0
#define __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE 1
#include "../wrapper/Plumed.h"
CPP
cat > "$task/probe/c-only.cpp" <<'CPP'
#define __PLUMED_WRAPPER_IMPLEMENTATION 0
#define __PLUMED_WRAPPER_CXX 0
#include "../wrapper/Plumed.h"
CPP
cat > "$task/probe/shared.cpp" <<'CPP'
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_EXTERN 0
#include "../wrapper/Plumed.h"
CPP
cat > "$task/probe/zero.cpp" <<'CPP'
#define __PLUMED_WRAPPER_IMPLEMENTATION 1
#define __PLUMED_WRAPPER_EXTERN 0
#define __PLUMED_WRAPPER_CXX_ANONYMOUS_NAMESPACE 0
#include "../wrapper/Plumed.h"
CPP
printf '#include "../wrapper/Plumed.h"\n' > "$task/probe/no-options.cpp"
bash "$root/src/maketools/plumedcheck" --astyle="$root/astyle/astyle" \
  --astyle-options="$root/.astyle.options" "$task/probe/good.cpp" "$task/probe/c-only.cpp" \
  "$task/probe/shared.cpp" "$task/probe/zero.cpp" "$task/probe/no-options.cpp" >"$task/result" 2>&1
for file in good c-only; do
  if grep "$file.cpp:" "$task/result" | grep -q ':plmd_include_dots:'; then
    cat "$task/result" >&2; exit 1
  fi
done
for file in shared zero no-options; do
  if ! grep "$file.cpp:" "$task/result" | grep -q ':plmd_include_dots:'; then
    cat "$task/result" >&2; exit 1
  fi
done
echo 'wrapper include lint accepts private/C-only boundaries and rejects shared or inherited options'
