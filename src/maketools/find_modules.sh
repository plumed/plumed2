#!/bin/bash

cd ../
for dir in *
do
  
  if test -f "$dir/module.type"
  then
    case "$(cat "$dir/module.type")" in
    (always) echo $dir ;;
    (default-on) echo "$plumed_disabled_modules" | grep :$dir: > /dev/null 2> /dev/null || echo $dir ;;
    (default-off) echo "$plumed_enabled_modules" | grep :$dir: > /dev/null 2> /dev/null && echo $dir ;;
    esac
  fi
done
