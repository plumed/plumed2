#!/usr/bin/env bash
# Exercise capability reporting without selecting a provider or launching ranks.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo "usage: $0 THREAD_MPI_PLUMED REAL_MPI_PROVIDER" >&2
  exit 2
fi
program=$1
provider=$2
unset PLUMED_MPI_WRAPPER_FORCE
unset PLUMED_MPI_WRAPPER
"$program" --no-mpi --has-mpi
"$program" --no-mpi --has-thread-mpi
if "$program" --no-mpi --mpi-ranks-are-threads; then
  echo 'embedded thread-MPI was reported as selected without an explicit launch' >&2
  exit 1
fi
if "$program" --no-mpi --has-real-mpi; then
  echo 'internal thread-MPI was incorrectly reported as real MPI' >&2
  exit 1
fi
export PLUMED_MPI_WRAPPER="$provider"
"$program" --no-mpi --has-mpi
"$program" --no-mpi --has-real-mpi
# Both capabilities may coexist; rank kind describes only the selected one.
"$program" --no-mpi --has-thread-mpi
if "$program" --no-mpi --mpi-ranks-are-threads; then
  echo 'external real MPI was incorrectly reported as thread ranks' >&2
  exit 1
fi
echo 'runtime MPI capability and selected rank-kind queries passed'
