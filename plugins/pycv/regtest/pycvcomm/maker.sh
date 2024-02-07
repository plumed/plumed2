#!/bin/bash

failed=""
for a in rt-*; do
  echo $a
  (
  cd "$a" || exit
  make
  
  if make; then
    echo ok
  else
    failed="${failed} ${a}"
  fi
  )
done
echo $failed
