#! /usr/bin/env bash

if [ "$1" = --description ] ; then
  echo "compile a .cpp file into a shared library"
  exit 0
fi

if [ "$1" = --options ] ; then
  echo "--description --options"
  exit 0
fi

source "$PLUMED_ROOT"/src/config/compile_options.sh

if [ $# != 1 ] || [[ "$1" != *.cpp ]] ;
then
  echo "ERROR"
  echo "type 'plumed mklib file.cpp'"
  exit 1
fi


file="$1"
obj="${file%%.cpp}".o
lib="${file%%.cpp}".$soext

if [ ! -f "$file" ]
then
  echo "ERROR: I cannot find file $file"
  exit 1
fi

# Compile and link beside the destination, then publish only a complete library.
tmpdir=$(mktemp -d "$(dirname "$lib")/plumed_mklib.XXXXXX") || exit 1
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT

if test "$PLUMED_IS_INSTALLED" = yes ; then
  link_libs="$link_libs_installed"
else
  link_libs="$link_libs_uninstalled"
fi

tmpobj="$tmpdir/$(basename "$obj")"
tmplib="$tmpdir/$(basename "$lib")"
if eval "$compile" "$tmpobj" "$file" && eval "$link_shared" "$tmpobj" "$link_libs" -o "$tmplib" ; then
  mv "$tmplib" "$lib"
else
  echo "ERROR: compiling or linking $file" >&2
  exit 1
fi



