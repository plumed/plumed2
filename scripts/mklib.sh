#! /bin/bash

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

rm -f "$obj" "$lib"

eval "$compile" "$obj" "$file" && $link "$lib" "$obj"



