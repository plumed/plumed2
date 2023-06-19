#! /usr/bin/env bash

file=$1

grep -E '[^[:print:]]' $file 2> /dev/null && exit 0

sed "s/-\(0\.0*0 \)/ \1/g; 
     s/-\(0\.0*0$\)/ \1/g" $file > $file.$$.tmp
mv $file.$$.tmp $file

