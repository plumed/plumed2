#!/bin/bash
#formatted with shfmt -ci -s -bn -i 2
#checked with shellcheck
#
# Compare two cython-generated C files (e.g. the committed python/plumed.c and a
# freshly generated one), ignoring the compressed copies of the string table.
#
# Since version 3.1, cython stores the module strings compressed with zlib, bz2
# and zstd, using the compression modules of the python interpreter that runs
# cython. These blocks thus depend on the python version (zstd needs >=3.14) and
# on the versions of the compression libraries (e.g. libzstd 1.4 and 1.5 give
# different bytes), not only on plumed.pyx and cplumed.pxd.
# The uncompressed copy of the same strings (the "#else /* compression: none"
# branch) is kept in the comparison, so changes in the strings are still detected.
#
# Exit status: 0 if the files are equivalent, 1 if they differ, 2 on errors.
# When the files differ, the first lines of the difference are printed on stderr.

if (($# != 2)); then
  echo "usage: $0 file1.c file2.c" >&2
  exit 2
fi

for file in "$1" "$2"; do
  if [[ ! -f "$file" ]]; then
    echo "$0: cannot find $file" >&2
    exit 2
  fi
done

strip_compressed() {
  awk '/^[[:space:]]*#if \(CYTHON_COMPRESS_STRINGS\)/ { skip = 1 }
       /^[[:space:]]*#else \/\* compression: none/ { skip = 0 }
       !skip' "$1"
}

if ! difference="$(diff <(strip_compressed "$1") <(strip_compressed "$2"))"; then
  echo "$difference" | head -n 40 >&2
  exit 1
fi
