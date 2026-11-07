#!/usr/bin/env bash
# Build the MPI host after the kernel: no MPI link or initialization in PLUMED.
set -euo pipefail
if (($# != 3 && $# != 4)); then
  echo 'usage: check-mpi-host-injection.sh MPI_CXX_COMPILER MPIEXEC SERIAL_KERNEL [PROVIDER_LIBRARY]' >&2
  exit 2
fi
ROOT=$(cd "$(dirname "$0")/../../.." && pwd)
TASK_DIR=$(mktemp -d "${TMPDIR:-/tmp}/plumed-host-injection.XXXXXX")
trap 'rm -rf "$TASK_DIR"' EXIT
kernel_dir=$(cd "$(dirname "$3")" && pwd)
program_name=$(basename "$3")
program_name=${program_name#lib}
program_name=${program_name%Kernel.*}
for candidate in "$kernel_dir/$program_name" "$kernel_dir/../.."; do
  if [[ -r "$candidate/src/main/main.cpp" && -r "$candidate/src/config/compile_options.sh" ]]; then
    mpi_root=$(cd "$candidate" && pwd)
    break
  fi
done
if [[ -z ${mpi_root:-} ]]; then
  echo 'cannot find the installed standard main beside the requested kernel' >&2
  exit 1
fi
include_dir="$mpi_root/src/include"
if [[ ! -d "$include_dir" ]]; then include_dir="$mpi_root/../../include"; fi
build_main() {
  env PLUMED_ROOT="$mpi_root" PLUMED_INCLUDEDIR="$include_dir" \
    PLUMED_PROGRAM_NAME="$program_name" \
    "$ROOT/scripts/mkmpi.sh" --executable --compiler "$1" "${@:2}"
}
for link in shared runtime; do
  build_main "$1" --link "$link" --out "$TASK_DIR/mpi-host-$link"
done
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
build_main "$1" --link runtime --cflag "$TASK_DIR/no-mpi-lifecycle.cpp" \
  --out "$TASK_DIR/no-init-host"
no_init_environment=(env -u PLUMED_MPI_WRAPPER -u PLUMED_MPI_WRAPPER_FORCE "PLUMED_KERNEL=$3")
"${no_init_environment[@]}" "$TASK_DIR/no-init-host" --no-mpi --has-real-mpi
"${no_init_environment[@]}" "$TASK_DIR/no-init-host" --mpi --no-mpi info --root
status=0
"${no_init_environment[@]}" PLUMED_KERNEL="$TASK_DIR/missing-kernel" \
  "$TASK_DIR/no-init-host" --no-mpi info --root > "$TASK_DIR/no-init-error.log" 2>&1 || status=$?
if ((status == 0 || status == 91)); then
  cat "$TASK_DIR/no-init-error.log" >&2
  echo "uninitialized host did not fail without entering MPI" >&2
  exit 1
fi
# Exercise a kernel-thrown exception as well as the invalid-loader error above.
status=0
"${no_init_environment[@]}" PLUMED_MPI_WRAPPER_FORCE=invalid \
  "$TASK_DIR/no-init-host" --no-mpi info --root > "$TASK_DIR/kernel-error.log" 2>&1 || status=$?
if ((status == 0 || status == 91)) || \
    ! grep -q 'PLUMED_MPI_WRAPPER_FORCE must be 0 or 1' "$TASK_DIR/kernel-error.log"; then
  cat "$TASK_DIR/kernel-error.log" >&2
  echo 'host did not report the kernel exception without entering MPI' >&2
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
"$1" -std=c++98 "$ROOT/regtest/basic/rt-mpiwrap-api/mpi-native-probe.cpp" -o "$TASK_DIR/native-probe"
if command -v ldd >/dev/null 2>&1; then
  ldd "$TASK_DIR/native-probe"
  ldd "$TASK_DIR/mpi-host-shared"
  ldd "$TASK_DIR/mpi-host-runtime"
fi
if ! "$2" -n 2 "$TASK_DIR/native-probe" 2; then
  echo 'Native MPI failed before loading PLUMED; retrying the probe without inherited include/library paths' >&2
  env -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH -u LIBRARY_PATH -u LD_LIBRARY_PATH \
    "$1" -std=c++98 "$ROOT/regtest/basic/rt-mpiwrap-api/mpi-native-probe.cpp" -o "$TASK_DIR/clean-native-probe" && \
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
        environment["PLUMED_MPI_WRAPPER_FORCE"] = "1"
    else:
        environment.pop("PLUMED_MPI_WRAPPER", None)
        environment.pop("PLUMED_MPI_WRAPPER_FORCE", None)
    for link in ("shared", "runtime"):
        for ranks in (1, 2, 4):
            work = directory / f"{mode}-{link}-{ranks}"
            work.mkdir()
            shutil.copyfile(fixture / "plumed.dat", work / "plumed.dat")
            environment["PLUMED_KERNEL"] = str(kernel)
            result = subprocess.run(
                [str(launcher), "-n", str(ranks), str(directory / f"mpi-host-{link}"),
                 "driver", "--plumed=plumed.dat", "--timestep=0.05", "--ixyz",
                 str(root / "regtest/trajectories/trajectory.xyz"), "--dump-forces", "ff"],
                cwd=work, env=environment, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, timeout=120,
            )
            if result.returncode:
                print(result.stdout, file=sys.stderr)
                raise SystemExit(f"MPI host failed: mode={mode}, link={link}, ranks={ranks}")
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
                    raise SystemExit(f"MPI host did not produce {name}: mode={mode}, link={link}, ranks={ranks}")
                shutil.copyfile(reference, expected)
                # Use the same portable normalization as the ordinary regtests.
                subprocess.run([str(root / "regtest/scripts/fixzeroes.sh"), str(actual)], check=True)
                subprocess.run([str(root / "regtest/scripts/fixzeroes.sh"), str(expected)], check=True)
                if actual.read_bytes() != expected.read_bytes():
                    print("".join(difflib.unified_diff(
                        expected.read_text().splitlines(True), actual.read_text().splitlines(True),
                        fromfile=str(reference), tofile=str(actual))), file=sys.stderr)
                    raise SystemExit(f"host injection differs: {name}, link={link}, ranks={ranks}")
            print(f"MPI host passed: mode={mode}, link={link}, ranks={ranks}")
PYTHON
