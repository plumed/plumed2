#! /bin/bash

if [ "$1" = --description ] ; then
  echo "create a new collective variable from a template"
  exit 0
fi

if [ "$1" = --options ] ; then
  echo "--description --options"
  exit 0
fi

if [ $# != 2 ] ;
then
  echo "ERROR"
  echo "type 'plumed newcv directive classname'"
  echo "E.g. 'plumed newcv TORSION Torsion'"
  exit 1
fi

directive=$1
classname=$2

sed "s/TEMPLATE/$directive/g
     s/Template/$classname/g" "$PLUMED_ROOT"/src/colvar/Template.cpp > ${classname}.cpp

