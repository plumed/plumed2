#!/usr/bin/env bash
# Build the MPI host after the kernel: no MPI link or initialization in PLUMED.
set -euo pipefail
if (($# != 3 && $# != 4)); then
  echo 'usage: check-mpi-host-injection.sh MPI_CXX_COMPILER MPIEXEC SERIAL_KERNEL [PROVIDER_LIBRARY]' >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-host-injection.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
"$1" -std=c++98 -DPLUMED_MPI_HOST_LINK_RUNTIME=1 "$ROOT/regtest/scripts/mpi-cli.cpp" -ldl -o "$TASK_DIR/mpi-host"
# Fail deterministically if --no-mpi reaches a lifecycle call, including the
# exception path. These executable definitions interpose the shared MPI library.
cat > "$TASK_DIR/no-mpi-lifecycle.cpp" <<'SOURCE'
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
namespace PLMD {
extern "C" {
int MPI_Init(int *, char ***) {
  std::fputs("unexpected MPI_Init\n", stderr);
  std::exit(91);
}
int MPI_Finalize(void) {
  std::fputs("unexpected MPI_Finalize\n", stderr);
  std::exit(91);
}
int MPI_Abort(MPI_Comm, int) {
  std::fputs("unexpected MPI_Abort\n", stderr);
  std::exit(91);
}
}
}
SOURCE
"$1" -std=c++98 -DPLUMED_MPI_HOST_LINK_RUNTIME=1 \
  "$ROOT/regtest/scripts/mpi-cli.cpp" "$TASK_DIR/no-mpi-lifecycle.cpp" \
  -ldl -o "$TASK_DIR/no-init-host"
no_init_environment=(env -u PLUMED_MPI_WRAPPER -u PLUMED_MPI_WRAPPER_FORCE "PLUMED_KERNEL=$3")
"${no_init_environment[@]}" "$TASK_DIR/no-init-host" --no-mpi --has-real-mpi
"${no_init_environment[@]}" "$TASK_DIR/no-init-host" --mpi --no-mpi info --root
status=0
"${no_init_environment[@]}" PLUMED_MPI_HOST_EXPLICIT_PATH="$TASK_DIR/missing-kernel" \
  "$TASK_DIR/no-init-host" --no-mpi info --root > "$TASK_DIR/no-init-error.log" 2>&1 || status=$?
if ((status != 1)); then
  cat "$TASK_DIR/no-init-error.log" >&2
  echo "uninitialized host did not return a clean error" >&2
  exit 1
fi
# A later --mpi must enable initialization; also verify our interceptor works.
status=0
"${no_init_environment[@]}" "$TASK_DIR/no-init-host" --no-mpi --mpi info --root \
  > "$TASK_DIR/init-control.log" 2>&1 || status=$?
if ((status != 91)) || ! grep -q 'unexpected MPI_Init' "$TASK_DIR/init-control.log"; then
  cat "$TASK_DIR/init-control.log" >&2
  echo 'MPI option precedence or lifecycle interceptor failed' >&2
  exit 1
fi
echo 'MPI host honors --no-mpi without initialization, finalization, or abort'

# Separate launcher/runtime failures from any host-header or kernel behavior.
"$1" -std=c++98 "$ROOT/regtest/scripts/mpi-native-probe.cpp" -o "$TASK_DIR/native-probe"
if command -v ldd >/dev/null 2>&1; then
  ldd "$TASK_DIR/native-probe"
  ldd "$TASK_DIR/mpi-host"
fi
if ! "$2" -n 2 "$TASK_DIR/native-probe" 2; then
  echo 'Native MPI failed before loading PLUMED; retrying the probe without inherited include/library paths' >&2
  env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH -u LIBRARY_PATH -u LD_LIBRARY_PATH \
    "$1" -std=c++98 "$ROOT/regtest/scripts/mpi-native-probe.cpp" -o "$TASK_DIR/clean-native-probe" && \
  env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH -u LIBRARY_PATH -u LD_LIBRARY_PATH \
    "$2" -n 2 "$TASK_DIR/clean-native-probe" 2 || true
  exit 1
fi
python3 - "$ROOT" "$TASK_DIR" "$2" "$3" "${4:-}" <<'PYTHON'
import difflib
import os
import pathlib
import shutil
import subprocess
import sys

root, directory, launcher, kernel = map(pathlib.Path, sys.argv[1:5])
provider = pathlib.Path(sys.argv[5]).resolve() if sys.argv[5] else None
root, directory, launcher, kernel = [p.resolve() for p in (root, directory, launcher, kernel)]
fixture = root / "regtest/basic/rt-mpi0"
environment = os.environ.copy()
for name in ("PLUMED_MPI_WRAPPER", "PLUMED_MPI_WRAPPER_FORCE", "PLUMED_LOAD_DEBUG"):
    environment.pop(name, None)
environment["PLUMED_KERNEL"] = str(kernel)
for mode in (("injection", "environment") if provider else ("injection",)):
    if mode == "environment":
        environment["PLUMED_MPI_WRAPPER"] = str(provider)
    else:
        environment.pop("PLUMED_MPI_WRAPPER", None)
    for explicit in (False, True):
        for ranks in (1, 2, 4):
            work = directory / f"{mode}-{'path' if explicit else 'selected'}-{ranks}"
            work.mkdir()
            shutil.copyfile(fixture / "plumed.dat", work / "plumed.dat")
            if explicit:
                environment["PLUMED_MPI_HOST_EXPLICIT_PATH"] = str(kernel)
                # Explicit-path creation must not accidentally use the environment.
                environment["PLUMED_KERNEL"] = str(work / "missing-kernel")
            else:
                environment.pop("PLUMED_MPI_HOST_EXPLICIT_PATH", None)
                environment["PLUMED_KERNEL"] = str(kernel)
            environment["PLUMED_MPI_HOST_EXPECTED_RANKS"] = str(ranks)
            result = subprocess.run(
                [str(launcher), "-n", str(ranks), str(directory / "mpi-host"),
                 *(["--no-injection"] if mode == "environment" else []),
                 "driver", "--plumed=plumed.dat", "--timestep=0.05", "--ixyz",
                 str(root / "regtest/trajectories/trajectory.xyz"), "--dump-forces", "ff"],
                cwd=work, env=environment, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, timeout=120,
            )
            if result.returncode:
                print(result.stdout, file=sys.stderr)
                raise SystemExit(f"MPI host failed: mode={mode}, explicit={explicit}, ranks={ranks}")
            for reference in fixture.glob("*.reference"):
                name = reference.name[:-len(".reference")]
                if name.startswith("ff.") and int(name.split(".")[1]) >= ranks:
                    continue
                # Driver adds the rank suffix only when the communicator has
                # more than one rank. The reference fixture uses four ranks.
                actual = work / ("ff" if ranks == 1 and name == "ff.0" else name)
                expected = work / (name + ".expected")
                if not actual.is_file():
                    print(result.stdout, file=sys.stderr)
                    print("Produced files: " + ", ".join(sorted(p.name for p in work.iterdir())), file=sys.stderr)
                    raise SystemExit(f"MPI host did not produce {name}: mode={mode}, explicit={explicit}, ranks={ranks}")
                shutil.copyfile(reference, expected)
                # Use the same portable normalization as the ordinary regtests.
                subprocess.run([str(root / "regtest/scripts/fixzeroes.sh"), str(actual)], check=True)
                subprocess.run([str(root / "regtest/scripts/fixzeroes.sh"), str(expected)], check=True)
                if actual.read_bytes() != expected.read_bytes():
                    print("".join(difflib.unified_diff(
                        expected.read_text().splitlines(True), actual.read_text().splitlines(True),
                        fromfile=str(reference), tofile=str(actual))), file=sys.stderr)
                    raise SystemExit(f"host injection differs: {name}, explicit={explicit}, ranks={ranks}")
            print(f"MPI host passed: mode={mode}, explicit={explicit}, ranks={ranks}")
PYTHON
