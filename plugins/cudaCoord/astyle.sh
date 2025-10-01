#! /usr/bin/env bash

astyle=$(which astyle || ../../astyle/astyle)

for file in *.c *.cpp *.h *.inc.in *.cu *.cuh; do

  test -f "$file" || continue

  echo -n "astyle $file"

  $astyle --options=../../.astyle.options <"$file" >"$file.tmp" && {
    if cmp -s "$file" "$file.tmp"; then
      echo
    else
      cp "$file.tmp" "$file"
      echo " +++ PATCHED"
      #git add $file
    fi
  }

  rm "$file.tmp"

done
