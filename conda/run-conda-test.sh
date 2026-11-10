#!/usr/bin/env bash
# Exercise solver selection and optional add-on installation using local builds.
set -euo pipefail
cd "$(dirname "$0")"
: "${VERSION:?set VERSION to the version built by run-conda-build.sh}"

bash run-extension-test.sh
bash run-pycv-test.sh
bash run-cuda-coordination-test.sh

task=$(mktemp -d "${TMPDIR:-/tmp}/plumed-conda-packaging.XXXXXX")
trap 'rm -rf "$task"' EXIT
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_rmaps_base_oversubscribe=1
# Keep MPI smoke tests on software transports in CI containers. MPICH's UCX
# netmod otherwise probes host RDMA devices (including Azure MANA interfaces).
export UCX_TLS=${UCX_TLS:-tcp,self}
# UCX settings do not constrain Open MPI's other PML/BTL choices. Use a
# portable software transport for these packaging checks on hosted runners.
export OMPI_MCA_pml=${OMPI_MCA_pml:-ob1}
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,tcp}
unset PLUMED_KERNEL PLUMED_LIB PLUMED_ROOT PLUMED_INCLUDEDIR PLUMED_PROGRAM_NAME
unset PLUMED_MPI_WRAPPER PLUMED_MPI_WRAPPER_FORCE

for mpi in mpich openmpi; do
  prefix="$task/$mpi"
  conda create --yes --prefix "$prefix" --override-channels -c local -c conda-forge \
    "plumed=$VERSION"
  if [[ $(uname) == Darwin ]]; then
    test "$(/usr/bin/stat -f %l "$prefix/bin/plumed")" = 1
  fi
  "$prefix/bin/plumed" --has-thread-mpi
  # Verify that installing PLUMED alone does not pull an external MPI runtime.
  for record in "$prefix"/conda-meta/{mpi,mpich,openmpi}-*.json; do
    test ! -e "$record"
  done
  kernel=$(find "$prefix/lib" -maxdepth 1 -name 'libplumedKernel.*' -type f | head -n 1)
  test -n "$kernel"
  cp "$kernel" "$task/kernel-before"

  # Installing a runtime chooses the matching variant of the single add-on name.
  conda install --yes --prefix "$prefix" --override-channels -c local -c conda-forge \
    "$mpi" "plumed-mpi=$VERSION"
  if [[ $(uname) == Darwin && "$mpi" == openmpi ]]; then
    # The add-on solve can replace an earlier Open MPI with a different build.
    # Copy the exact runtime selected by that transaction: a hardlinked prte
    # makes dyld search conda's package cache rather than this environment.
    # Install the add-on normally above so its no_link check remains meaningful.
    mpi_package=$(conda list --prefix "$prefix" --export | grep '^openmpi=')
    conda install --yes --copy --force-reinstall --freeze-installed --prefix "$prefix" \
      --override-channels -c local -c conda-forge "$mpi_package"
  fi
  if [[ $(uname) == Darwin ]]; then
    test "$(/usr/bin/stat -f %l "$prefix/bin/plumed-mpi")" = 1
  fi
  cmp "$kernel" "$task/kernel-before"
  for record in "$prefix"/conda-meta/plumed-mpi-*.json; do
    [[ "$record" == *"-mpi_${mpi}_"*.json ]]
  done
  "$prefix/bin/plumed" --has-thread-mpi
  "$prefix/bin/mpiexec" --version
  "$prefix/bin/mpiexec" -n 2 "$prefix/bin/plumed-mpi" --has-real-mpi
  "$prefix/bin/plumed-mpi" --no-mpi info --root

  conda remove --yes --prefix "$prefix" --override-channels -c local -c conda-forge plumed-mpi
  cmp "$kernel" "$task/kernel-before"
  test ! -e "$prefix/bin/plumed-mpi"
  test ! -e "$prefix/lib/libplumed-mpi-provider.so"
  test ! -e "$prefix/lib/libplumed-mpi-provider.dylib"
  "$prefix/bin/plumed" --has-thread-mpi
  "$prefix/bin/plumed" info --root
done
