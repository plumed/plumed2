#!/usr/bin/env bash
# Install one host MPI per job and record its exact compiler/launcher paths.
set -euo pipefail

: "${BACKEND:?Set BACKEND to serial, tmpi, openmpi, or mpich}"
case "$BACKEND" in serial|tmpi|openmpi|mpich) ;; *) exit 2 ;; esac

if [[ "$RUNNER_OS" == Linux ]]; then
  sudo apt-get update -qq
  packages=(gawk ninja-build pkg-config)
  case "$BACKEND" in
    openmpi) packages+=(openmpi-bin libopenmpi-dev) ;;
    mpich) packages+=(mpich libmpich-dev openmpi-bin) ;;
  esac
  sudo apt-get install -y "${packages[@]}"
  echo 'CC=gcc' >> "$GITHUB_ENV"
  echo 'CXX=g++' >> "$GITHUB_ENV"
  if [[ "$BACKEND" == openmpi || "$BACKEND" == mpich ]]; then
    echo "MD_MPICC=/usr/bin/mpicc.$BACKEND" >> "$GITHUB_ENV"
    if [[ "$BACKEND" == openmpi ]]; then
      echo 'MD_MPICXX=/usr/bin/mpic++.openmpi' >> "$GITHUB_ENV"
    else
      echo 'MD_MPICXX=/usr/bin/mpicxx.mpich' >> "$GITHUB_ENV"
      # Ubuntu's MPICH uses external PMIx, incompatible with its Hydra launcher:
      # https://github.com/pmodels/mpich/issues/7064
      # Keep MPICH compilers/libraries, but use a PMIx-capable process manager.
      # This is the same workaround validated by the PLUMED injection tests.
      echo 'UCX_TLS=tcp,sm,self' >> "$GITHUB_ENV"
    fi
    echo 'MD_MPIEXEC=/usr/bin/mpiexec.openmpi' >> "$GITHUB_ENV"
    echo 'MD_LAUNCHER_KIND=openmpi' >> "$GITHUB_ENV"
  fi
else
  brew install gawk ninja libomp pkgconf
  echo 'CC=/usr/bin/clang' >> "$GITHUB_ENV"
  echo 'CXX=/usr/bin/clang++' >> "$GITHUB_ENV"
  echo "MD_LIBOMP=$(brew --prefix libomp)" >> "$GITHUB_ENV"
  if [[ "$BACKEND" == openmpi || "$BACKEND" == mpich ]]; then
    formula=open-mpi
    other=mpich
    if [[ "$BACKEND" == mpich ]]; then
      formula=mpich
      other=open-mpi
    fi
    # Hosted images can already contain the other implementation. Unlink it
    # before Homebrew links the requested one; never rely on its default MPI.
    if brew list --versions "$other" >/dev/null 2>&1; then
      brew unlink "$other"
    fi
    brew install "$formula"
    prefix=$(brew --prefix "$formula")
    echo "MD_MPICC=$prefix/bin/mpicc" >> "$GITHUB_ENV"
    echo "MD_MPICXX=$prefix/bin/mpicxx" >> "$GITHUB_ENV"
    echo "MD_MPIEXEC=$prefix/bin/mpiexec" >> "$GITHUB_ENV"
    echo "MD_LAUNCHER_KIND=$BACKEND" >> "$GITHUB_ENV"
    # Include helper binaries needed by the MPI process manager.
    echo "$prefix/bin" >> "$GITHUB_PATH"
  fi
fi
