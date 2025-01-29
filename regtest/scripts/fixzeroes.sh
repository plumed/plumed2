#! /usr/bin/env bash

file=$1

grep -q -E '[^[:print:]]' $file 2> /dev/null && exit 0

sed "s/-\(0\.0*0 \)/ \1/g; 
     s/-\(0\.0*0$\)/ \1/g" $file |
sed "s/-nan/ nan/g" |
awk '{if($1=="#!" && $2=="FIELDS") gsub("@[1-9][0-9]*","@___")
      print($0)}' > $file.$$.tmp
mv $file.$$.tmp $file

