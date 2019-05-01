#!/bin/bash

cd ../
for dir in *
do
  
  if test -f "$dir/module.type"
  then
    case "$(cat "$dir/module.type")" in
    (always) echo $dir ;;
    (default-on) test -f $dir.off || echo $dir ;;
    (default-off) test -f $dir.on && echo $dir ;;
    esac
  fi
done
