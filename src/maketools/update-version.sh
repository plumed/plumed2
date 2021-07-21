#! /bin/bash

{
echo "#ifndef __PLUMED_config_version_h"
echo "#define __PLUMED_config_version_h"
echo "#define PLUMED_VERSION_SHORT \"$(
  if test -f ../../VERSION ; then
    grep -v "#" ../../VERSION | sed  's/^\([0-9][0-9]*\.[0-9][0-9]*\).*/\1/'
  else
    echo "Unknown"
  fi
)\""

echo "#define PLUMED_VERSION_LONG \"$(
  if test -f ../../VERSION ; then
    grep -v "#" ../../VERSION
  else
    echo "Unknown"
  fi
)\""

echo "#define PLUMED_VERSION_GIT \"$(
  if test -d ../../.git && hash git 2> /dev/null ; then
# in case it does not work, fallback to normal hash (12 char long)
    git describe --long --dirty --always || git rev-parse  --short=12 HEAD
  else
    echo "Unknown"
  fi
)\""

echo "#define PLUMED_VERSION_MAJOR $(
  grep -v "#" ../../VERSION | sed  's/^\([0-9][0-9]*\).*/\1/'
)"

echo "#define PLUMED_VERSION_MINOR $(
  grep -v "#" ../../VERSION | sed  's/^[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'
)"

echo "#define PLUMED_VERSION_PATCH $(
  grep -v "#" ../../VERSION | sed  's/^[0-9][0-9]*\.[0-9][0-9]*\.\([0-9][0-9]*\).*/\1/'
)"

echo "#endif"

} > $1~

cmp -s $1~ $1 || cp $1~ $1
rm $1~



