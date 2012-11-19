#! /bin/bash

if [ "$1" = --description ] ; then
  echo "create a new collective variable from a template"
  exit 0
fi

if [ $# != 2 ] ;
then
  echo "ERROR"
  echo "type 'plumed newcv directive cvname'"
  echo "E.g. 'plumed newcv TORSION Torsion'"
  exit 1
fi

directive=$1
cvname=$2

sed "s/TEMPLATE/$directive/g
     s/Template/$cvname/g" $PLUMED_ROOT/src/basic/ColvarTemplate.cpp > Colvar${cvname}.cpp

