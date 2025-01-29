#! /usr/bin/env bash

#formatted with shfmt  https://github.com/mvdan/sh/releases
#checked with shellcheck

cd src || exit 1
for file in *.c *.cpp *.h *.inc.in; do

  test -f "$file" || continue

  echo -n "astyle $file"

  ../../../astyle/astyle --options=../../../.astyle.options <"$file" >"${file}.tmp" && {
    if cmp -s "$file" "${file}.tmp"; then
      echo
    else
      cp "${file}.tmp" "$file"
      echo " +++ PATCHED"
    fi
  }

  rm "${file}.tmp"

done
