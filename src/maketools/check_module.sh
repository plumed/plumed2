#!/bin/bash

dir=$1

msg="ERROR: this module depends on $dir, which is presently disabled"

cd ../

if test -f "$dir/module.type"
then
  case "$(cat "$dir/module.type")" in
  (always) exit 0 ;;
  (default-on)
    if echo "$plumed_disabled_modules" | grep :$dir: > /dev/null 2> /dev/null
    then
      echo "$msg"
      exit 1
    else 
      exit 0
    fi ;;
  (default-off)
    if echo "$plumed_enabled_modules" | grep -q :$dir: > /dev/null 2> /dev/null
    then
      exit 0
    else 
      echo "$msg"
      exit 1
    fi ;;
  esac
fi

