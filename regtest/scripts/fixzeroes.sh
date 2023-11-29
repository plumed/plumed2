#! /usr/bin/env bash

file=$1

grep -E '[^[:print:]]' $file 2> /dev/null && exit 0

sed "s/-\(0\.0*0 \)/ \1/g; 
     s/-\(0\.0*0$\)/ \1/g" $file > $file.$$.tmp
sed "s/-\(nan \)/ \1/g; 
     s/-\(nan$\)/ \1/g" $file.$$.tmp > $file.$$.tmp2
mv $file.$$.tmp2 $file

