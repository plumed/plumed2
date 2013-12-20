#! /bin/bash

case "$1" in
(light)
  rm -fr *.on *.off
  touch $(grep default */module.type | awk 'BEGIN{FS="/"}{printf($1".off ")}') ;;
(heavy)
  rm -fr *.on *.off
  touch $(grep default */module.type | awk 'BEGIN{FS="/"}{printf($1".on ")}') ;;
(reset)
  rm -fr *.on *.off;;
(*)
  echo "ERROR"
  exit ;;
esac

