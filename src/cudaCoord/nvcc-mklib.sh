#! /usr/bin/env bash

if [ "$1" = --description ]; then
  echo "compile a .cu file into a shared library"
  exit 0
fi

if [ "$1" = --options ]; then
  echo "--description --options"
  exit 0
fi

source "$PLUMED_ROOT"/src/config/compile_options.sh

if [ $# != 1 ] || [[ "$1" != *.cu ]]; then
  echo "ERROR"
  echo "type 'mklib file.cu'"
  exit 1
fi

file="$1"
obj="${file%%.cu}".o
lib="${file%%.cu}".$soext

if [ ! -f "$file" ]; then
  echo "ERROR: I cannot find file $file"
  exit 1
fi

rm -f "$obj" "$lib"

#nvcc "$2" -Xcompiler -fPIC -c -o "$kernel"
compile="nvcc -ccbin ${compile}"
compile=${compile//-Wall/}
compile=${compile//-pedantic/}
for opt in -W -pedantic -f; do
  #echo $compile
  compile=${compile//${opt}/-Xcompiler ${opt}}
  #echo $compile
done
#compile=${compile/-fPIC/-Xcompiler -fPIC}

link_installed="nvcc -shared${link_installed#*-shared}"
link_installed=${link_installed/-rdynamic/-Xcompiler -rdynamic}
link_installed=${link_installed/-Wl,/-Xlinker }
link_installed=${link_installed/-fopenmp/-Xcompiler -fopenmp}

echo ${link_installed}
link_uninstalled="nvcc -shared${link_uninstalled#*-shared}"
link_uninstalled=${link_uninstalled/-rdynamic/-Xcompiler -rdynamic}
link_uninstalled=${link_uninstalled/-Wl,/-Xlinker }
link_uninstalled=${link_uninstalled/-fopenmp/-Xcompiler -fopenmp}
if test "$PLUMED_IS_INSTALLED" = yes; then
  eval "$compile" "$obj" "$file" && eval "$link_installed" "$lib" "$obj"
else
  eval "$compile" "$obj" "$file" && echo "$link_uninstalled" "$lib" "$obj"
fi