#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 1 || ! -x $1 ]]; then
  echo "usage: $0 /path/to/plumed-cli" >&2
  exit 2
fi

plumed=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/mpiwrap-cli.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/plumed.dat" <<'EOF'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
EOF

"$plumed" --has-thread-mpi
for ranks in 1 2 4; do
  (cd "$tmp" && "$plumed" --tmpi-threads "$ranks" driver \
    --plumed plumed.dat --noatoms --multi "$ranks" --timestep 1 \
    --trajectory-stride 1 >"driver.$ranks.out" 2>"driver.$ranks.err")
done

# Benchmark uses the same generic launcher and exercises atom distribution.
cat > "$tmp/benchmark.dat" <<'EOF'
d: DISTANCE ATOMS=1,2
RESTRAINT ARG=d AT=1 KAPPA=1
PRINT ARG=d FILE=benchmark.colvar STRIDE=1
EOF
for ranks in 1 2 4; do
  mkdir "$tmp/benchmark.$ranks"
  cp "$tmp/benchmark.dat" "$tmp/benchmark.$ranks/plumed.dat"
  (cd "$tmp/benchmark.$ranks" && "$plumed" --tmpi-threads "$ranks" benchmark \
    --natoms 32 --nsteps 5 --domain-decomposition >benchmark.out 2>benchmark.err)
done
python3 - "$tmp" <<'PYTHON'
import math
from pathlib import Path
import sys

root = Path(sys.argv[1])
reference = None
for ranks in (1, 2, 4):
    lines = (root / f"benchmark.{ranks}" / "benchmark.colvar").read_text().splitlines()
    rows = [[float(x) for x in line.split()] for line in lines if not line.startswith("#")]
    if len(rows) != 5 or any(not math.isfinite(x) for row in rows for x in row):
        raise SystemExit(f"invalid benchmark output with {ranks} ranks")
    if reference is not None and rows != reference:
        raise SystemExit(f"benchmark output differs with {ranks} ranks")
    reference = rows
PYTHON

if (cd "$tmp" && "$plumed" --tmpi-threads 1 driver \
    --plumed plumed.dat --noatoms --multi 2 >mismatch.out 2>mismatch.err); then
  echo "launcher accepted a rank count inconsistent with driver --multi" >&2
  exit 1
fi
grep -q "invalid number of processes for multi environment" "$tmp/mismatch.err"

if (cd "$tmp" && "$plumed" --tmpi-threads 2 --no-mpi driver \
    >"conflict.out" 2>"conflict.err"); then
  echo "launcher accepted a contradictory --no-mpi option" >&2
  exit 1
fi
grep -q -- "--tmpi-threads cannot be combined with --mpi or --no-mpi" "$tmp/conflict.err"

cat > "$tmp/chdir.dat" <<'EOF'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
PLUMED FILE=guest.dat CHDIR=thread-mpi-directory
EOF
cat > "$tmp/guest.dat" <<'EOF'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
EOF
if (cd "$tmp" && "$plumed" --tmpi-threads 2 driver \
    --plumed chdir.dat --noatoms --multi 2 --timestep 1 \
    --trajectory-stride 1 >chdir.out 2>chdir.err); then
  echo "thread-MPI unexpectedly accepted PLUMED CHDIR" >&2
  exit 1
fi
grep -q "CHDIR is not supported while a thread-MPI runtime is active" "$tmp/chdir.err"
test ! -e "$tmp/thread-mpi-directory"

echo "thread-MPI CLI launch, benchmark distribution, rank handoff, and CHDIR guard passed"
