#! /usr/bin/env bash

cat ../*/*/*cpp |
awk '{
     if(inside && NF>0 && $1!="/*" && $1!="*/") { 
        printf "%s\n", $0 
        inside=0
     }
     if($1=="//+PLUMEDOC"){ 
        printf "%s:", $3 
        inside=1
     }
  }'

rm modules
for mod in `ls ../src/`; do
    if [ -d "../src/$mod" ] ; then 
       for file in `ls ../src/$mod`; do 
           if [ ! -d "../src/$mod/$file" ] && [ $(grep -c "PLUMED_REGISTER_ACTION" ../src/$mod/$file) -ge 1 ] ; then
              nreg=`grep -c "PLUMED_REGISTER_ACTION" ../src/$mod/$file`
              for (( i=1; i<=$nreg; i++ )) ; do 
                  modname=`grep "PLUMED_REGISTER_ACTION" ../src/$mod/$file | head -n $i | tail -n 1 | awk -F "\"" '{print $2}'` 
                  echo $modname:$mod >> modules
              done
           fi 
       done 
    fi
done
