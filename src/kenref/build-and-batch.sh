#!/usr/bin/env bash
# build-and-batch.sh — the old name of build-and-patch.sh, kept so that callers using it keep working
# (KEnRef's cmake/InvokePlumed.cmake, behind KENREF_WITH_PLUMED_GROMACS, invokes the script by this name).
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/build-and-patch.sh" "$@"
