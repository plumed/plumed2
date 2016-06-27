#! /bin/bash

file=$1

sed "s/-0.000/ 0.000/g" $file > $file.$$.tmp
mv $file.$$.tmp $file

