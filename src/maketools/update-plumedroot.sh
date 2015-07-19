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
if test -n "$PLUMED_PREFIX" ; then
  echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "WARNING: using PLUMED_PREFIX variable is deprecated"
  echo "         please use one of the following choices:"
  echo "         (1) at configure time:"
  echo "           ./configure --prefix=$PLUMED_PREFIX"
  echo "         (2) or later, at install time:"
  echo "           make install prefix=$PLUMED_PREFIX"
  echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# this will result in an error in a later version
# now we fall back to the previous behavior
  prefix="${PLUMED_PREFIX}"
fi

# if environment variable PLUMED_LIBSUFFIX is set, complain
if test -n "$PLUMED_LIBSUFFIX" ; then
  echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  echo "WARNING: using PLUMED_LIBSUFFIX variable is deprecated, please use"
  echo "         ./configure --program-suffix='$PLUMED_LIBSUFFIX'"
  echo "         at configure time"
  echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
# this will result in an error in a later version
# now we fall back to the previous behavior
  PLUMED_LIBSUFFIX="${PLUMED_LIBSUFFIX:=}"
  test -n "$PLUMED_LIBSUFFIX" && PLUMED_LIBSUFFIX="_${PLUMED_LIBSUFFIX}"
fi

PLUMED_PROGRAM_NAME="$(echo plumed | sed "${program_transform_name}")${PLUMED_LIBSUFFIX}"
PLUMED_ROOT="$prefix/lib/${PLUMED_PROGRAM_NAME}"

{
echo "PLUMED_INSTALL_ROOT=${PLUMED_ROOT}"
echo "PLUMED_INSTALL_PREFIX=$prefix"
echo "PLUMED_PROGRAM_NAME=${PLUMED_PROGRAM_NAME}"
} > install.conf

sed "s|@PLUMED_ROOT@|${PLUMED_ROOT}|g" > $1~

cmp -s $1~ $1 || cp $1~ $1
rm $1~



