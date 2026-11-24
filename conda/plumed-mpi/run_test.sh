#!/usr/bin/env bash
set -euo pipefail

# Test both routes with ordinary installed files, not a build-tree kernel or
# an environment-selected provider. This also exercises conda's relocation.
unset PLUMED_KERNEL PLUMED_LIB PLUMED_ROOT PLUMED_INCLUDEDIR PLUMED_PROGRAM_NAME
unset PLUMED_MPI_WRAPPER PLUMED_MPI_WRAPPER_FORCE
# Allow two ranks in small/root-owned CI containers.
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export OMPI_MCA_rmaps_base_oversubscribe=1
# Keep MPI smoke tests on software transports in CI containers. MPICH's UCX
# netmod otherwise probes host RDMA devices (including Azure MANA interfaces).
export UCX_TLS=${UCX_TLS:-tcp,self}
# Keep Open MPI on software transports too, independently of its PML default.
export OMPI_MCA_pml=${OMPI_MCA_pml:-ob1}
export OMPI_MCA_btl=${OMPI_MCA_btl:-self,tcp}

export MPICH_CC="$CC" OMPI_CC="$CC"
"$PREFIX/bin/mpicc" ${CPPFLAGS:-} ${CFLAGS:-} native-host.c -I"$PREFIX/include" \
  -L"$PREFIX/lib" "-Wl,-rpath,$PREFIX/lib" -lplumed ${LDFLAGS:-} -o native-host
# Any attempted runtime compilation must fail: the prebuilt provider is required.
export MPICC=$(command -v false) MPICXX=$(command -v false)

cat > plumed.dat <<'INPUT'
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar FMT=%8.3f
INPUT
cat > trajectory.xyz <<'INPUT'
2
10 10 10
H 0 0 0
H 1 0 0
INPUT

# New MPI hosts inject through the standard main; old hosts discover the file.
for host in "$PREFIX/bin/plumed-mpi" ./native-host; do
  "$PREFIX/bin/mpiexec" -n 2 "$host" --has-real-mpi
  "$PREFIX/bin/mpiexec" -n 2 "$host" driver \
    --plumed plumed.dat --ixyz trajectory.xyz --length-units nm
  awk '!/^#/ { if ($2 != 1) exit 1; found=1 } END { if (!found) exit 1 }' colvar
  rm colvar
done
# Installing a real-MPI add-on must not hijack the explicit internal launcher.
"$PREFIX/bin/plumed" --tmpi-threads 2 driver \
  --plumed plumed.dat --ixyz trajectory.xyz --length-units nm
awk '!/^#/ { if ($2 != 1) exit 1; found=1 } END { if (!found) exit 1 }' colvar
"$PREFIX/bin/plumed-mpi" --no-mpi info --root
