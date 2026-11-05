#!/usr/bin/env bash
# A real-MPI provider must allow CHDIR in a nested PLUMED instance.
set -euo pipefail
if [[ $# -ne 1 || -z ${PLUMED_MPI_WRAPPER:-} || -z ${PLUMED_KERNEL:-} ]]; then
  echo "usage: PLUMED_MPI_WRAPPER=... PLUMED_KERNEL=... $0 MPI_CLI_HELPER" >&2
  exit 2
fi
helper=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
tmp=$(mktemp -d "${TMPDIR:-/tmp}/plumed-mpi-chdir.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
mkdir "$tmp/guest"
cat > "$tmp/host.dat" <<'INPUT'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
PLUMED FILE=guest.dat CHDIR=guest
INPUT
cat > "$tmp/guest/guest.dat" <<'INPUT'
UNITS LENGTH=A ENERGY=kcal/mol TIME=ps
INPUT
(cd "$tmp" && "$helper" driver --plumed host.dat --parse-only --natoms 1 >out 2>err) || {
  cat "$tmp/out" "$tmp/err" >&2
  exit 1
}
grep -q 'running on separate directory guest' "$tmp/out"
echo 'real-MPI provider permits nested PLUMED CHDIR'
