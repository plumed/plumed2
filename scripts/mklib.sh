#! /usr/bin/env bash

if [ "$1" = --description ]; then
  echo "compile one or more *.cpp files into a shared library"
  echo " you can create and export the variable PLUMED_MKLIB_CFLAGS with some extra compile time flags to be used"
  echo " you can create and export the variable PLUMED_MKLIB_LDFLAGS with some extra link time flags (and libraries) to be used"
  exit 0
fi

if [ "$1" = --options ]; then
  echo "--description --options"
  exit 0
fi

if [ $# == 0 ]; then
  echo "ERROR"
  echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
  exit 1
fi

source "$PLUMED_ROOT"/src/config/compile_options.sh

lib="${1%%.cpp}".$soext
rm -f "$lib"

toRemove=""
remover() {
  #if compilation fails or crashes the trap will delete the temporary files
  for f in $toRemove; do
    rm -f "${f}"
  done
}
trap "remover" EXIT
objs=""

for file; do

  if [[ "$file" != *.cpp ]]; then
    echo "ERROR"
    echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
    exit 1
  fi

  obj="${file%%.cpp}".o

  if [ ! -f "$file" ]; then
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

  if grep -q '^#include "\(bias\|colvar\|function\|sasa\|vatom\)\/ActionRegister.h"' "${tmpfile}"; then
    echo >&2 'WARNING: using a legacy ActionRegister.h include path, please use <<#include "core/ActionRegister.h">>'
    sed -i.bak 's%^#include ".*/ActionRegister.h"%#include "core/ActionRegister.h"%g' "${tmpfile}"
  fi

  if grep -q '^#include "\(cltools\)\/CLToolRegister.h"' "${tmpfile}"; then
    echo >&2 'WARNING: using a legacy  CLToolRegister.h include path, please use <<#include "core/CLToolRegister.h">>'
    sed -i.bak 's%^#include ".*/CLToolRegister.h"%#include "core/CLToolRegister.h"%g' "${tmpfile}"
  fi

  rm -f "$obj"
  eval "$compile" "$PLUMED_MKLIB_CFLAGS" "$tmpfile" -o "$obj" || {
    echo "ERROR: compiling $file"
    exit 1
  }

  objs="$objs $obj"
done

link_command="$link_uninstalled"

if test "$PLUMED_IS_INSTALLED" = yes; then
  link_command="$link_installed"
fi

eval "$link_command" "$PLUMED_MKLIB_LDFLAGS" $objs -o "$lib"
