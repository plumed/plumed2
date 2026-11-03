#! /usr/bin/env bash

file=$1

if [ -z "$file" ] || [ ! -f "$file" ]; then
  echo "usage: $0 FILE" >&2
  exit 2
fi

grep -q -E '[^[:print:][:space:]]' "$file" 2> /dev/null && exit 0

awk '
  function normalize_line(line, result, i, j, separator, token, start) {
    result=""
    i=1
    while(i<=length(line)) {
      separator=""
      while(i<=length(line) && substr(line,i,1) ~ /[[:space:]]/) {
        separator=separator substr(line,i,1)
        i++
      }
      if(i>length(line)) return result separator
      start=i
      while(i<=length(line) && substr(line,i,1) !~ /[[:space:]]/) i++
      token=substr(line,start,i-start)
      if(token ~ /^[-+]?0+(\.0*)?([eE][-+]?[0-9]+)?$/) {
        # Canonicalize spacing before zero fields so padded and unpadded
        # numeric formats compare equally without changing nonzero fields.
        if(result=="") separator=""
        else if(separator!="") separator=" "
        sub(/^[+-]/,"",token)
      }
      result=result separator token
    }
    return result
  }
  {
    print normalize_line($0)
  }
' "$file" > "$file.$$.tmp"
mv "$file.$$.tmp" "$file"

