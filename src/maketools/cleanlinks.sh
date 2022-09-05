#! /usr/bin/env bash

for file in *
do
  test -d $file && test -e $file/.tmpdir && rm -fr $file/
  test -L $file && rm -fr $file
done
# this echo is to avoid triggering errors when last test fails
echo links removed
