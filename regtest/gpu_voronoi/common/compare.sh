#! /usr/bin/env bash
# Usage: compare.sh FILE FIRST
# Compares, in each row of FILE (written by PRINT or DUMPDERIVATIVES with ARG=c0,g0,c1,g1,...),
# the columns FIRST and FIRST+1, FIRST+2 and FIRST+3, ...: the CPU and the GPU values
# must agree within 2e-8 + 2e-10*|CPU|. Prints one line, which goes into the reference.
awk -v first="$2" -v name="$1" '
function abs(x) { return x<0 ? -x : x }
/^#/ { next }
{
  rows++
  for(i=first; i<NF; i+=2) {
    if($i !~ /^[-+]?[0-9]/ || $(i+1) !~ /^[-+]?[0-9]/ || !(abs($i-$(i+1)) <= 2e-8+2e-10*abs($i))) {
      if(++bad<=5) printf "%s: row %d column %d: CPU %s GPU %s\n", name, rows, i, $i, $(i+1)
    }
  }
}
END {
  if(rows==0) printf "%s: no rows\n", name
  else if(bad) printf "%s: %d values differ\n", name, bad
  else printf "%s: %d rows, CPU and GPU agree\n", name, rows
}' "$1"
