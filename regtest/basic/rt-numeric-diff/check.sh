#!/usr/bin/env bash

check() {
  local description="$1" expected="$2" actual
  shift 2
  if ../../../scripts/numeric-diff "$@" > diagnostics 2>&1; then actual=0; else actual=$?; fi
  if [[ "$actual" = "$expected" ]]; then
    echo "$description: OK"
  else
    echo "FAILURE: $description (expected $expected, got $actual)"
    cat diagnostics
  fi
}

printf '#! FIELDS x y\n0 1.000000\n1 label\n' > reference
printf '#! FIELDS x y\n0 1.000001\n1 label\n' > actual
check 'Last decimal place' 0 reference actual 1e-6
check 'Beyond tolerance' 1 reference actual 1e-7
printf '#! FIELDS x z\n0 1.000000\n1 label\n' > actual
check 'Header mismatch' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 1.000000\n1 other\n' > actual
check 'Label mismatch' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 nan\n1 label\n' > actual
check 'Unexpected NaN' 1 reference actual 1e-6
printf '1\n' > single-reference
printf 'nan\n' > single-actual
check 'Single-column NaN' 1 single-reference single-actual 1e-6
printf '#! FIELDS x y\n0 inf\n1 label\n' > actual
check 'Unexpected infinity' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 1e999\n1 label\n' > actual
check 'Overflowing number' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 1.000000 extra\n1 label\n' > actual
check 'Column count' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 1.000000\n' > actual
check 'Missing row' 1 reference actual 1e-6
printf '#! FIELDS x y\n0 1.000000\n1 label\n2 extra\n' > actual
check 'Extra row' 1 reference actual 1e-6
: > empty
check 'Empty actual file' 1 reference empty 1e-6
check 'Empty reference file' 1 empty reference 1e-6
check 'Two empty files' 0 empty empty 1e-6
check 'Invalid tolerance' 2 reference actual nan
check 'Overflowing tolerance' 2 reference actual 1e999
