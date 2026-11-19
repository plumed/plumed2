#! /usr/bin/env bash

module_header="$(cd "$(dirname "$0")" && pwd)/module-header"
cd ../
for dir in *
do
  
  if test -f "$dir/module.yml"
  then
    kind=$("$module_header" type "$dir/module.yml") || exit 1
    case "$kind" in
    (always) echo $dir ;;
    (default-on) test -f $dir.off || echo $dir ;;
    (default-off) test -f $dir.on && echo $dir ;;
    esac
  fi
done

exit 0
