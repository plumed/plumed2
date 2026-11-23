#! /usr/bin/env bash

cd ../
maketools/module-header check-obsolete . || exit 1
for dir in *
do
  
  if test -f "$dir/module.yml"
  then
    kind=$(maketools/module-header type "$dir/module.yml") || exit 1
    case "$kind" in
    (always) echo $dir ;;
    (default-on) test -f $dir.off || echo $dir ;;
    (default-off) test -f $dir.on && echo $dir ;;
    esac
  fi
done

exit 0
