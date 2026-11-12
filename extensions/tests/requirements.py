#!/usr/bin/env python3
"""Check old/new regression feature query dispatch and false-result handling."""

import os
import subprocess
import tempfile
from pathlib import Path

source = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="plumed-requirements-") as scratch:
    root = Path(scratch)
    executable = root / "plumed"
    executable.write_text("""#!/usr/bin/env bash
if [[ $1 == modules && $2 == --help ]]; then
  [[ $TEST_MODE == modern ]] || exit 1
  echo 'Usage: plumed modules --module NAME --has-feature FEATURE'
  exit 0
fi
if [[ $1 == modules ]]; then
  [[ $* == *'--module enabled --has-feature boost_serialization'* ]]
  exit $?
fi
if [[ $1 == config && $2 == module ]]; then
  [[ $3 == enabled ]]
  exit $?
fi
if [[ $1 == config && $2 == has ]]; then
  [[ $3 == boost_serialization ]]
  exit $?
fi
exit 9
""")
    executable.chmod(0o755)
    for mode in ("legacy", "modern"):
        for requirement, expected in (
            ("enabled:boost_serialization", 0),
            ("disabled:boost_serialization", 1),
            ("enabled:missing", 1),
        ):
            command = [
                "bash",
                "-c",
                'plumed_nompi="$1"; source "$2"; plumed_check_feature "$3"',
                "test",
                str(executable),
                str(source / "regtest/scripts/requirements.sh"),
                requirement,
            ]
            result = subprocess.run(
                command, env=dict(os.environ, TEST_MODE=mode), check=False, timeout=10
            )
            assert result.returncode == expected, (mode, requirement, result.returncode)
print("Legacy and module-scoped regression prerequisites passed")
