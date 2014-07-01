#! /bin/bash

test -n "$1" || {
  echo "Usage: $0 outfile"
  exit 1
}

# if environment variable "prefix" is set, use it.
# otherwise defaults to /usr/local
prefix="${prefix:=/usr/local}"

# if environment variable PLUMED_PREFIX is set,
# override the present prefix
PLUMED_PREFIX="${PLUMED_PREFIX:=$prefix}"

PLUMED_LIBSUFFIX="${PLUMED_LIBSUFFIX:=}"
test -n "$PLUMED_LIBSUFFIX" && PLUMED_LIBSUFFIX="-${PLUMED_LIBSUFFIX}"
PLUMED_ROOT="${PLUMED_PREFIX}/lib/plumed${PLUMED_LIBSUFFIX}/"

{
echo "PLUMED_INSTALL_ROOT=${PLUMED_ROOT}"
echo "PLUMED_INSTALL_PREFIX=${PLUMED_PREFIX}"
echo "PLUMED_INSTALL_LIBSUFFIX=${PLUMED_LIBSUFFIX}"
} > install.conf

sed "s|@PLUMED_ROOT@|${PLUMED_ROOT}|g" > $1~

cmp -s $1~ $1 || cp $1~ $1
rm $1~



