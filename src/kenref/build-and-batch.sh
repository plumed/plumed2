#!/usr/bin/env bash
# =============================================================================
# build-and-batch.sh — build PLUMED with the KEnRef module AND batch a GROMACS
# 2025.x with it (native GMX_USE_PLUMED alone is not full integration until the
# GROMACS source is patched with `plumed patch`).
#
# Ensures kenref_core (delegated to KEnRef), builds PLUMED (its own autotools),
# then find-or-fetches a GROMACS 2025.x, `plumed patch`es it, and builds/installs
# it. Power-user one-liners:
#     ./build-and-batch.sh -y                                 # fetch everything
#     ./build-and-batch.sh --gromacs-src ~/gromacs-2025.4 -y  # use a provided gromacs
# =============================================================================
set -euo pipefail
KN_SCRIPT="build-and-batch"; KN_BATCHES=1
# shellcheck source=_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_common.sh"

kn_parse_args "$@"
kn_finalize_toolchain

say "PLUMED-side build (with gromacs batch): kenref -> plumed -> gromacs"
ensure_kenref
build_plumed
build_gromacs
say "DONE. PLUMED -> ${PLUMED_PREFIX} ; GROMACS (batched) -> ${GROMACS_PREFIX}."
