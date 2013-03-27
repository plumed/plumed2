#!/bin/bash

for dir in ../src/*
do
  if [ -e $dir/module.type ] ; then
    if [ "$(cat "$dir/module.type")" == always ] ; then
      continue
    fi
#Check not colvar
    if `echo $dir | grep colvar | grep -v multi 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not vatom
    if `echo $dir | grep vatom 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not cltools
    if `echo $dir | grep cltools 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not vesselbase
    if `echo $dir | grep vesselbase 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not function
    if `echo $dir | grep function 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not bias
    if `echo $dir | grep bias 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not generic
    if `echo $dir | grep setup 1>/dev/null 2>& 1` ; then
       continue
    fi
# Check not setup
    if `echo $dir | grep generic 1>/dev/null 2>& 1` ; then
       continue
    fi
    case "$(cat "$dir/module.type")" in
     (default-on) test -f $dir.off || echo $dir | sed -e 's/..\/src\///g';;
     (default-off) test -f $dir.on && echo $dir | sed -e 's/..\/src\///g';;
    esac

#    echo `echo $dir | sed -e 's/..\/src\///g'`
  fi
done
