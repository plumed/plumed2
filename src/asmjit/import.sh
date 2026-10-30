#!/usr/bin/env bash
# Import the AsmJit 1.9 source layout (including OpenMM's bundled copy).
set -euo pipefail
if [[ $# != 1 ]]; then
  echo "usage: $0 /path/to/asmjit (or /path/to/openmm/libraries/asmjit)" >&2
  exit 1
fi
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/import.py" "$1"
