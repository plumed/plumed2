#!/bin/bash

for dir in ../src/*
do
  if [ -e $dir/module.type ] ; then
    if [ "$(cat "$dir/module.type")" == always ] ; then
      continue
    fi
    short="$(echo $dir | sed -e 's/..\/src\///g')"
    if ! test -d $short ; then
      continue
    fi
    case "$(cat "$dir/module.type")" in
     (default-on) test -f $dir.off || echo $short ;;
     (default-off) test -f $dir.on && echo $short ;;
    esac

#    echo `echo $dir | sed -e 's/..\/src\///g'`
  fi
done
