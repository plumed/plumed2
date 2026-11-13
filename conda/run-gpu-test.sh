#!/usr/bin/env bash
# Run already-built packages on NVIDIA hardware; never override driver metadata.
set -euo pipefail
cd "$(dirname "$0")/.."
: "${VERSION:?set VERSION to the package version}"
: "${GPU_CONDA_CHANNEL:?set GPU_CONDA_CHANNEL to the package channel}"
: "${GPU_TEST_ROOT:?set GPU_TEST_ROOT to a writable test directory}"
unset CONDA_OVERRIDE_CUDA
nvidia-smi
mkdir -p "$GPU_TEST_ROOT"
python conda/gpu-plan.py > "$GPU_TEST_ROOT/plan.json"
failed=0
kernel_spec=${GPU_KERNEL_SPEC:-plumed=$VERSION}
while IFS='|' read -r variant build_variant cuda_version; do
  prefix="$GPU_TEST_ROOT/$variant"
  if [[ $variant == cuda* ]]; then
    package=plumed-ext-cuda-coordination
  else
    package=plumed-ext-openacc-ptm
  fi
  spec="$package=$VERSION=$build_variant"
  if ! conda create -y -p "$prefix" --strict-channel-priority --override-channels -c "$GPU_CONDA_CHANNEL" -c conda-forge "$kernel_spec" "$spec" python=3.13 > "$GPU_TEST_ROOT/$variant.log" 2>&1; then
    echo "$variant: installation failed (see $GPU_TEST_ROOT/$variant.log)"
    echo 'INSTALLATION FAILED' > "$GPU_TEST_ROOT/$variant.status"
    failed=1
    continue
  fi
  # Pin the first selected exact kernel for all remaining alternatives too.
  kernel_spec=$(python - "$prefix" <<'PY'
import json,sys
from pathlib import Path
records = [json.loads(p.read_text()) for p in (Path(sys.argv[1])/"conda-meta").glob("plumed-*.json")]
record = next(r for r in records if r["name"] == "plumed")
print("plumed=" + record["version"] + "=" + record["build"])
PY
)
  if [[ $variant == cuda* ]]; then
    command=("$prefix/bin/python" extensions/tests/cuda_coordination.py --prefix "$prefix" --cuda-version "$cuda_version" --gpu)
  else
    command=("$prefix/bin/python" extensions/tests/openacc.py --prefix "$prefix" --device nvidia)
  fi
  if "${command[@]}" >> "$GPU_TEST_ROOT/$variant.log" 2>&1; then
    echo "$variant: GPU numerical tests passed"
    echo PASS > "$GPU_TEST_ROOT/$variant.status"
  else
    echo "$variant: GPU numerical tests failed (see $GPU_TEST_ROOT/$variant.log)"
    echo 'NUMERICAL TEST FAILED' > "$GPU_TEST_ROOT/$variant.status"
    failed=1
  fi
  cat "$GPU_TEST_ROOT/$variant.log"
done < <(python - "$GPU_TEST_ROOT/plan.json" <<'PY'
import json,sys
for item in json.load(open(sys.argv[1])):
    print("|".join(item[key] for key in ("name", "build", "cuda")))
PY
)
python - "$GPU_TEST_ROOT" "$VERSION" <<'PY'
import json,sys
from pathlib import Path
root, version = Path(sys.argv[1]), sys.argv[2]
plan = json.loads((root/"plan.json").read_text())
summary = {"version":version, "variants":{}, "passed":False}
for item in plan:
    name = item["name"]
    status = root/(name+".status")
    summary["variants"][name] = {"status":status.read_text().strip() if status.exists() else "NOT VALIDATED"}
    packages = [json.loads(p.read_text()) for p in (root/name/"conda-meta").glob("*.json")]
    (root/("packages-"+name+".json")).write_text(json.dumps([{key:p[key] for key in ("name","version","build")} for p in packages],indent=2)+"\n")
summary["passed"] = bool(plan) and all(item["status"]=="PASS" for item in summary["variants"].values())
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
PY
exit "$failed"
