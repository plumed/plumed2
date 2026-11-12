#! /usr/bin/env bash

if [ "$1" = --description ] ; then
  echo "compile one or more *.cpp files into a shared library"
  exit 0
fi

MANUAL='Compile one or more *.cpp files into a shared library.

Usage:

  plumed mklib [options] files1.cpp [file2.cpp ...]

Options:
  -h, --help
                         Print this help and exit
  -o LIBNAME, --out LIBNAME
                         Name of the output library. If missing, the name
                         of the first input file will me used, with its
                         suffix properly adjusted.
  -I DIR, --include DIR   Add a header search directory (repeatable).
                         Relative to the working directory; quote paths with spaces.
                         Cannot be combined with -n (include paths and headers
                         are not tracked by the timestamp cache).
  -n                     Reuse a library newer than all input files.
                         Otherwise build a replacement, keeping the old library
                         available until linking succeeds. If another builder
                         publishes an up-to-date library, keep that result.

Environment variables:
  PLUMED_MKLIB_CXX       Override the configured C++ compiler (then CXX)
  PLUMED_MKLIB_CFLAGS    Extra compile time flags to be used
  PLUMED_MKLIB_LDFLAGS   Extra link time flags (and libraries) to be used
'

OPTIONS="--description --options -h --help -o --out -I --include -n"

if [ $# == 0 ]
then
  echo "ERROR"
  echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
  exit 1
fi

prefix=""
lib=""
no_clobber=""
files=() # empty array
includes=()
for opt
do
  prefixopt="$prefix$opt"
  prefix=""
  case "$prefixopt" in
    (-o|--out) prefix="--out=";;
    (-I|--include) prefix="--include=";;
    (--include=*)
      dir="${prefixopt#--include=}"
      if [[ -z "$dir" ]]; then
        echo "ERROR: --include requires a nonempty directory"
        exit 1
      fi
      includes+=("-I$dir");;
    (-h|--help) echo "$MANUAL" ; exit ;;
    (--options) echo "$OPTIONS" ; exit ;;
    (--out=*) lib="${prefixopt#--out=}";;
    (-n) no_clobber=yes;;
    (-*)
      echo "ERROR: Unknown option $opt. Use --help for help."
      exit 1 ;;
    (*)
      files+=("${prefixopt}");;
  esac
done

if [[ -n "$prefix" ]]; then
  echo "ERROR: missing argument for ${prefix%=}"
  exit 1
fi
if [[ ${#includes[@]} -gt 0 && -n "$no_clobber" ]]; then
  echo "ERROR: --include cannot be combined with -n"
  exit 1
fi

if [ ${#files[@]} = 0 ] ; then
  echo ERROR
  echo "pass at least one file"
  exit 1
fi

# Completion generation queries this script before compiler settings exist.
# Only an actual library build needs the configured flags and suffix.
source "$PLUMED_ROOT"/src/config/compile_options.sh || exit 1

# Preserve the configured ABI flags while allowing an installed toolchain.
if [[ -n ${mklib_cxx:-} ]]; then
  selected_cxx=${PLUMED_MKLIB_CXX:-${CXX:-$mklib_cxx}}
  compile="$selected_cxx${compile#"$mklib_cxx"}"
  link_shared="$selected_cxx${link_shared#"$mklib_cxx"}"
  read -r -a compiler_command <<< "$selected_cxx"
  if ! command -v "${compiler_command[0]}" >/dev/null 2>&1; then
    echo "ERROR: C++ compiler ${compiler_command[0]} not found. For conda, install plumed-mklib and activate its environment." >&2
    exit 1
  fi
fi
signature=$({
  cat "$PLUMED_ROOT/src/config/module-build-id" 2>/dev/null
  printf '%s\n' "$compile" "$link_shared" "$PLUMED_MKLIB_CFLAGS" "$PLUMED_MKLIB_LDFLAGS"
} | cksum)

if [ -z "$lib" ]
then
  firstfile="${files[0]}"
  lib="${firstfile%%.cpp}".$soext
fi

recompile=no
[[ -f "$lib.plumed-build" && $(cat "$lib.plumed-build") == "$signature" ]] || recompile=yes
for file in "${files[@]}"
do
  if ! test "$lib" -nt "$file" ;
  then
    recompile=yes
  fi
done

if test -z "$no_clobber" ; then
  recompile=yes
fi

if test $recompile = no ; then
  echo "$lib is already up to date"
  exit 0
fi

toRemove=()
remover() {
  #if compilation fails or crashes the trap will delete the temporary files
  for f in "${toRemove[@]}"; do
    rm -rf "${f}"
  done
}
trap "remover" EXIT
objs=()
nobj=0

tmpdir=$(mktemp -d "$(dirname "$lib")/plumed_mklib.XXXXXX") || exit 1
tmp_lib="$tmpdir/$(basename "$lib")"

toRemove+=("$tmpdir")

for file in "${files[@]}"
do

  if [[ "$file" != *.cpp ]] ;
  then
    echo "ERROR"
    echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
    exit 1
  fi

  # objects are written in $tmpdir under their file name only (with a counter,
  # so that files with the same name in different directories do not collide):
  # $file can contain a directory (e.g. ../code/X.cpp), which does not exist in $tmpdir
  nobj=$((nobj+1))
  obj="${nobj}_$(basename "${file%%.cpp}").o"
  
  if [ ! -f "$file" ]
  then
    echo "ERROR: I cannot find file $file"
    exit 1
  fi
  #adding a simple tmpfile, to preprocess "in place" the input file,
  #this assumes the user has write permission in the current directory
  #which should be true since we are going to compile and output something here
  tmpfile=$(mktemp "${file%.cpp}.XXXXXX")
  mv "${tmpfile}" "${tmpfile}.cpp"
  tmpfile=${tmpfile}.cpp
  cp "${file}" "${tmpfile}"
  toRemove+=("$tmpfile" "${tmpfile}.bak")
  #this deprecates the shortcut to action register
  if grep -q '^#include "\(bias\|colvar\|function\|sasa\|vatom\)\/ActionRegister.h"' "${tmpfile}"; then
     >&2 echo 'WARNING: using a legacy ActionRegister.h include path, please use <<#include "core/ActionRegister.h">>'
     sed -i.bak 's%^#include ".*/ActionRegister.h"%#include "core/ActionRegister.h"%g' "${tmpfile}"
  fi
  #this deprecates the shortcut to CLtoolregister
  if grep -q '^#include "\(cltools\)\/CLToolRegister.h"' "${tmpfile}"; then
     >&2  echo 'WARNING: using a legacy  CLToolRegister.h include path, please use <<#include "core/CLToolRegister.h">>'
     sed -i.bak 's%^#include ".*/CLToolRegister.h"%#include "core/CLToolRegister.h"%g' "${tmpfile}"
  fi
  #this removes "#include "core/Atoms.h" and makes the compilation failing on the calls of Atoms
  #instead of on the include
  if grep -q '^#include\s*"core/Atoms.h"' "${tmpfile}"; then
     #\s match anly type of white space
     >&2  echo 'WARNING: "core/Atoms.h" does not exist anymore in  version >=2.10, you should change your code.'
     sed -i.bak '/^#include\s*"core\/Atoms.h"/d' "${tmpfile}"
  fi
  rm -f "$tmpdir/$obj"

  # Keep configured command/flag evaluation, but expand user paths only in
  # the evaluated command, as array arguments, never as shell program text.
  eval "$compile" "$PLUMED_MKLIB_CFLAGS" '"${includes[@]}" -o "$tmpdir/$obj" "$tmpfile"' || {
    echo "ERROR: compiling $file"
    printf '%s\n' "Compiler search PATH: $PATH" >&2
    exit 1
  }

  objs+=("$tmpdir/$obj")

done

link_libs="$link_libs_uninstalled"

if test "$PLUMED_IS_INSTALLED" = yes; then
  link_libs="${mklib_link_libs_installed:-$link_libs_installed}"
  # Inline matrix helpers can import BLAS/LAPACK directly. Preserve those actual
  # imports without requiring every private kernel dependency's linker files.
  if [[ -n ${mklib_link_libs_installed:-} ]]; then
    for header in blas/blas.h lapack/lapack.h; do
      sed -n -E 's/.*PLUMED_BLAS_F77_FUNC\( *([a-z0-9_]+) *,.*/\1/p' \
        "$PLUMED_INCLUDEDIR/$PLUMED_PROGRAM_NAME/$header"
    done | sort -u > "$tmpdir/numeric-symbols"
    nm -u $objs | awk '{n=$NF;sub(/@.*/,"",n);gsub(/^_+|_+$/,"",n);print n}' > "$tmpdir/imports"
    if grep -q -F -x -f "$tmpdir/numeric-symbols" "$tmpdir/imports"; then
      for flag in $link_libs_installed; do
        case "$flag" in
          *blas*|*lapack*|*mkl*|*blis*) link_libs+=" $flag" ;;
          Accelerate|vecLib) link_libs+=" -framework $flag" ;;
        esac
      done
    fi
  fi
fi

# The libraries (PLUMED_MKLIB_LDFLAGS, the kernel and its dependencies) follow the objects:
# linkers using --as-needed (the default on some distributions, e.g. Ubuntu) drop the
# libraries listed before the objects using them. The kernel would then not be recorded
# as a dependency, and LOAD would fail with undefined symbols when the kernel was loaded
# with RTLD_LOCAL (as done by the python interface).
eval "$link_shared" '"${objs[@]}"' "$PLUMED_MKLIB_LDFLAGS" "$link_libs" '-o "$tmp_lib"' || exit 1

if test -n "$no_clobber" ; then
  # Another caller may have finished compiling these inputs while we linked.
  recompile=no
  [[ -f "$lib.plumed-build" && $(cat "$lib.plumed-build") == "$signature" ]] || recompile=yes
  for file in "${files[@]}" ; do
    if ! test "$lib" -nt "$file" ; then
      recompile=yes
    fi
  done
  if test "$recompile" = no ; then
    exit 0
  fi
fi

# Keep the previous library available until its replacement is complete.
# The temporary directory is on the destination filesystem, so rename publishes
# the complete library atomically even when another rank is loading it.
printf '%s\n' "$signature" > "$tmpdir/build-id"
mv "$tmp_lib" "$lib"
mv "$tmpdir/build-id" "$lib.plumed-build"
