#!/usr/bin/env bash
# Test-only normalization for funnel/rt-sphere, whose default parameters are
# z_cc=1.8, R_cyl=0.1, and slope=1. The branch boundary is z=s-1.7; the
# potential and radial derivative are continuous there, while the longitudinal
# derivative is undefined. Mask only its five exact grid samples. Round
# both derivative columns to five decimal places to absorb smooth-point
# last-digit differences without choosing an individual coordinate.
set -euo pipefail

if (($# != 1)); then
  echo 'usage: normalize-funnel-bias.sh FILE' >&2
  exit 2
fi

file=$1
tmp="$file.$$.tmp"
awk '
  function replace_field(line, field_to_replace, replacement, result, i, separator, token, start, field) {
    result=""
    i=1
    field=0
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
      field++
      if(field==field_to_replace) {
        if(result!="") separator=" "
        token=replacement
      }
      result=result separator token
    }
    return result
  }
  !/^#/ && NF {
    rows++
    if(NF!=5 || $4 !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ ||
       $5 !~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
      print "unexpected funnel BIAS row " rows > "/dev/stderr"
      exit 1
    }
    s=$1+0
    z=$2+0
    difference=z-(s-1.7)
    if(difference<0) difference=-difference
    if(difference<1e-10) {
      $0=replace_field($0,4,"UNDEFINED____")
      kinks++
    } else {
      $0=replace_field($0,4,sprintf("%.5f",$4))
    }
    $0=replace_field($0,5,sprintf("%.5f",$5))
  }
  {print}
  END {
    if(kinks!=5) {
      print "expected 5 exact funnel kinks, found " kinks > "/dev/stderr"
      exit 1
    }
  }
' "$file" > "$tmp"
mv "$tmp" "$file"
