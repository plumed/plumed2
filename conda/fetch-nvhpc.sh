#!/usr/bin/env bash
# Install an isolated build toolchain, without root or system MPI/math packages.
set -euo pipefail
version=${1:?provide the NVHPC version}
destination=${2:?provide a private extraction directory}
[[ $(uname -s) == Linux ]] || { echo 'NVHPC build tools require Linux' >&2; exit 1; }
mkdir -p "$destination"
metadata=$(mktemp)
archive=$(mktemp --suffix=.deb)
trap 'rm -f "$metadata" "$archive"' EXIT
python - "$version" "$metadata" <<'PY'
import json
import platform
import re
import sys
import urllib.request
version, output = sys.argv[1:]
assert re.fullmatch(r"\d+\.\d+", version)
architecture = {"x86_64": "amd64", "aarch64": "arm64"}[platform.machine()]
root = f"https://developer.download.nvidia.com/hpc-sdk/ubuntu/{architecture}/"
text = urllib.request.urlopen(root + "Packages", timeout=60).read().decode()
name = "nvhpc-" + version.replace(".", "-")
for block in text.split("\n\n"):
    fields = dict(line.split(": ", 1) for line in block.splitlines() if ": " in line and not line.startswith(" "))
    if fields.get("Package") == name:
        with open(output, "w") as stream:
            json.dump({"url": root + fields["Filename"].removeprefix("./"), "sha256": fields["SHA256"]}, stream)
        break
else:
    raise RuntimeError(f"Official NVIDIA repository has no {name} for {architecture}")
PY
url=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["url"])' "$metadata")
checksum=$(python -c 'import json,sys; print(json.load(open(sys.argv[1]))["sha256"])' "$metadata")
curl -fL --retry 3 "$url" -o "$archive"
printf '%s  %s\n' "$checksum" "$archive" | sha256sum --check
ar p "$archive" data.tar.xz | tar -xJf - -C "$destination" --wildcards \
  '*/compilers/*' '*/cuda/*'
root="$destination/opt/nvidia/hpc_sdk/Linux_$(uname -m)/$version"
"$root/compilers/bin/makelocalrc" -x "$root/compilers/bin"
rm -f "$root/compilers/glibc_version"
cuda=$(find "$root/cuda" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | sort -V | tail -n 1)
for component in bin include lib64 nvvm; do
  if [[ -e "$root/cuda/$cuda/$component" ]]; then
    ln -sfn "$cuda/$component" "$root/cuda/$component"
  fi
done
"$root/compilers/bin/nvc++" -V
printf 'NVHPC_HOME=%s\n' "$destination/opt/nvidia/hpc_sdk"
