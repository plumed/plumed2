#! /bin/bash

file=$1

sed "s/-\(0\.0*0 \)/ \1/g; 
     s/-\(0\.0*0$\)/ \1/g" $file > $file.$$.tmp
mv $file.$$.tmp $file

