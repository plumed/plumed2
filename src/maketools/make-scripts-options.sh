#! /bin/bash

for file in ../../scripts/*.sh
do
  name=${file##*/}
  name=${name%.sh}
  name=${name//-/_}
  echo -n "local cmd_keys_${name}=\""
  echo -n $(PLUMED_ROOT=../../ $file --options)
  echo "\""
done
