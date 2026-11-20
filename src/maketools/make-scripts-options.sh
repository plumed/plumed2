#! /usr/bin/env bash

# Shell-neutral records: exact tool name, tab, whitespace-separated options.
for file in ../../scripts/*.sh; do
  name=${file##*/}
  name=${name%.sh}
  options=$(PLUMED_ROOT=../../ "$file" --options)
  printf '%s\t%s\n' "$name" "${options//$'\n'/ }"
done
