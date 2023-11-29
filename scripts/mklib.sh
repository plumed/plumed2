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

if [ $# == 0 ]
then
  echo "ERROR"
  echo "type 'plumed mklib file1.cpp [file2.cpp etc]'"
  exit 1
fi

lib="${1%%.cpp}".$soext
rm -f "$lib"

objs=""

for file
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
  tmpfile=$(mktemp ${file%.cpp}.XXXXXX).cpp
  cp "${file}" "${tmpfile}"
  
  if grep -q '^#include "\(bias\|colvar\|function\|sasa\|vatom\)\/ActionRegister.h"' "${tmpfile}"; then
     >&2 echo 'WARNING: using a legacy ActionRegister.h include path, please use <<#include "core/ActionRegister.h">>'
     sed -i.bak 's%^#include ".*/ActionRegister.h"%#include "core/ActionRegister.h"%g' "${tmpfile}"
  fi
  
  if grep -q '^#include "\(cltools\)\/CLToolRegister.h"' "${tmpfile}"; then
     >&2  echo 'WARNING: using a legacy  CLToolRegister.h include path, please use <<#include "core/CLToolRegister.h">>'
     sed -i.bak 's%^#include ".*/CLToolRegister.h"%#include "core/CLToolRegister.h"%g' "${tmpfile}"
  fi
  
  rm -f "$obj"

  eval "$compile" "$PLUMED_MKLIB_CFLAGS" -o "$obj" "$tmpfile" || {
    echo "ERROR: compiling $file"
    exit 1
  }

  rm -f ${tmpfile} ${tmpfile}.bak ${tmpfile%.cpp}
  objs="$objs $obj"

done

if test "$PLUMED_IS_INSTALLED" = yes ; then
  eval "$link_installed" "$PLUMED_MKLIB_LDFLAGS" -o "$lib" "$objs"
else
  eval "$link_uninstalled" "$PLUMED_MKLIB_LDFLAGS" -o "$lib" "$objs"
fi

