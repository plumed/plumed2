#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 1 ]]; then
  echo "usage: $0 /path/to/openmm (or /path/to/lepton)" >&2
  exit 1
fi
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
python3 "$script_dir/import.py" "$1"
