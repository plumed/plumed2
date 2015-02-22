#! /bin/bash

file=$1

# only do fix for text files with "#! FIELDS " header
header="$(head -c 10 $file)"

if [ "$header" = "#! FIELDS " ] ; then
  sed "s/ -0.0000/  0.0000/g" $file > $file.$$.tmp
  mv $file.$$.tmp $file
fi

