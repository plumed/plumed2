#! /usr/bin/env bash

if [ "$1" = --description ] ; then
  echo "compile one or more *.cpp files into a shared library"
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
  -n                     No-clobber mode, similar to `mv -n`.
                         Does not overwrite an existing library.
                         If the library exists when the command is started,
                         skip also the compilation phase.

Environment variables:
  PLUMED_MKLIB_CFLAGS    Extra compile time flags to be used
  PLUMED_MKLIB_LDFLAGS   Extra link time flags (and libraries) to be used
'

OPTIONS="--description --options -h --help -o --out -n"

if [ $# == 0 ]
then
  echo "ERROR"
  echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
  exit 1
fi

source "$PLUMED_ROOT"/src/config/compile_options.sh

prefix=""
lib=""
no_clobber=""
files=() # empty array
for opt
do
  prefixopt="$prefix$opt"
  prefix=""
  case "$prefixopt" in
    (-o|--out) prefix="--out=";;
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

if [ ${#files[@]} = 0 ] ; then
  echo ERROR
  echo "pass at least one file"
  exit 1
fi

if [ -z "$lib" ]
then
  firstfile="${files[0]}"
  lib="${firstfile%%.cpp}".$soext
fi

recompile=no
for file in "${files[@]}"
do
  if ! test $lib -nt $file ;
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

rm -f "$lib"

toRemove=""
remover() {
  #if compilation fails or crashes the trap will delete the temporary files
  for f in $toRemove; do
    rm -rf "${f}"
  done
}
trap "remover" EXIT
objs=""

tmpdir=$(mktemp -d "plumed_mklib.XXXXXX")

toRemove="${toRemove} $tmpdir"

for file in "${files[@]}"
do

  if [[ "$file" != *.cpp ]] ;
  then
    echo "ERROR"
    echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
    exit 1
  fi

  obj="${file%%.cpp}".o
  
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
  toRemove="${toRemove} ${tmpfile} ${tmpfile}.bak"
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

  eval "$compile" "$PLUMED_MKLIB_CFLAGS" -o "$tmpdir/$obj" "$tmpfile" || {
    echo "ERROR: compiling $file"
    exit 1
  }

  objs="$objs $tmpdir/$obj"

done

link_command="$link_uninstalled"

if test "$PLUMED_IS_INSTALLED" = yes; then
  link_command="$link_installed"
fi

eval "$link_command" "$PLUMED_MKLIB_LDFLAGS" $objs -o "$tmpdir/$lib"

if test -n "$no_clobber" ; then
  # || true is necessary with recent coreutils
  mv -n "$tmpdir/$lib" $lib || true
else
  mv "$tmpdir/$lib" $lib
fi
