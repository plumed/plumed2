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
