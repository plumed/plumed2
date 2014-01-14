#! /bin/bash

{
echo "#define PLUMED_VERSION_SHORT \"$(
  if test -d ../../.git && hash git 2> /dev/null ; then
    git rev-parse --abbrev-ref HEAD
  elif test -f ../../VERSION ; then
    grep SHORT: ../../VERSION | sed 's/SHORT: //'
  else
    echo "Unknown"
  fi
)\""

echo "#define PLUMED_VERSION_LONG \"$(
  if test -d ../../.git && hash git 2> /dev/null ; then
    git describe --tags
  elif test -f ../../VERSION ; then
    grep LONG: ../../VERSION | sed 's/LONG: //'
  else
    echo "Unknown"
  fi
)\""
} > $1~

cmp -s $1~ $1 || cp $1~ $1
rm $1~



