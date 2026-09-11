#!/usr/bin/env bash
# =============================================================================
# build-only.sh — build PLUMED with the KEnRef module, WITHOUT GROMACS batching.
# (For the GROMACS-batched build, use build-and-batch.sh.)
#
# Ensures kenref_core (delegated to the KEnRef repo — kenref builds kenref), then
# builds PLUMED with its own autotools. Power-user one-liners:
#     ./build-only.sh -y                         # find kenref (or clone), build+install plumed
#     ./build-only.sh --kenref-src ~/KEnRef -y
# =============================================================================
set -euo pipefail
KN_SCRIPT="build-only"; KN_BATCHES=0
# shellcheck source=_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_common.sh"

kn_parse_args "$@"
kn_finalize_toolchain

say "PLUMED-side build (no gromacs): kenref -> plumed"
ensure_kenref
build_plumed
say "DONE. PLUMED (+kenref) installed to ${PLUMED_PREFIX}."
