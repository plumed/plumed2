#! /bin/bash

DIRS=$1

# keep only dirname
DIRS=${DIRS##*/}
echo "$DIRS"

test -z "$DIRS" && DIRS=*

for dir in $DIRS
do
test -d "$dir" || continue

test "$dir" = lapack && continue
test "$dir" = blas && continue
test "$dir" = molfile && continue
test "$dir" = lepton && continue
test "$dir" = asmjit && continue

cd $dir


for file in *.c *.cpp *.h *.inc.in
do

test -f "$file" || continue

echo -n "astyle $file"

../../astyle/astyle --options=../../.astyle.options < $file > $file.tmp && {
if cmp -s $file $file.tmp ; then
  echo 
else
  cp $file.tmp $file
  echo " +++ PATCHED"
  git add $file
fi
}

rm $file.tmp

done

cd -

done


