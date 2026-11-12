#!/usr/bin/env bash
# A serial PLUMED kernel must coexist with an environment selecting MPI HDF5.
set -euo pipefail
: "${VERSION:?set VERSION to the built package version}"
task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-hdf5-variant.XXXXXX")
trap 'rm -rf "$task"' EXIT
prefix="$task/environment"
conda create -y -p "$prefix" --override-channels -c local -c conda-forge "plumed=$VERSION"
kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
cp "$kernel" "$task/kernel-before"
conda install -y -p "$prefix" --override-channels -c local -c conda-forge 'hdf5=*=mpi_openmpi_*' openmpi
cmp "$kernel" "$task/kernel-before"
"$prefix/bin/plumed" --has-thread-mpi
"$prefix/bin/plumed" config has libhdf5 -q
cat > "$task/input.dat" <<'INPUT'
d: DISTANCE ATOMS=1,2 NOPBC
PRINT ARG=d FILE=table.h5
INPUT
cat > "$task/positions.xyz" <<'XYZ'
2
10 10 10
H 0 0 0
H 1 0 0
2
10 10 10
H 0 0 0
H 2 0 0
XYZ
(cd "$task"; "$prefix/bin/plumed" driver --ixyz positions.xyz --length-units nm --plumed input.dat)
"$prefix/bin/plumed" h5convert --input "$task/table.h5" --output "$task/table.txt"
python - "$task/table.txt" <<'PY'
from pathlib import Path
import sys
rows=[line.split() for line in Path(sys.argv[1]).read_text().splitlines() if not line.startswith('#')]
assert [float(row[1]) for row in rows] == [1.0,2.0], rows
PY
echo 'MPI-enabled HDF5 solver substitution and numerical table I/O passed'
