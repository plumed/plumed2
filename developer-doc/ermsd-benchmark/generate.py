#!/usr/bin/env python3
"""Deterministic 1S72 base-triplet references and smooth perturbed XYZ trajectories (nm)."""

import json
import math
import os
import random
import urllib.request
from pathlib import Path

root = Path(os.environ.get("ERMSD_BENCHMARK_DIR", Path.cwd())).resolve()
root.mkdir(parents=True, exist_ok=True)
if not (root / "1S72.pdb").exists():
    urllib.request.urlretrieve(
        "https://files.rcsb.org/download/1S72.pdb", root / "1S72.pdb"
    )
res = {}
for line in (root / "1S72.pdb").read_text().splitlines():
    if not line.startswith("ATOM  ") or line[16] not in " A":
        continue
    name = line[17:20].strip()
    if name not in ("A", "C", "G", "U"):
        continue
    key = (line[21], line[22:27])
    res.setdefault(key, (name, {}))[1][line[12:16].strip()] = tuple(
        float(line[a : a + 8]) * 0.1 for a in (30, 38, 46)
    )
triplets = []
names = []
for name, atoms in res.values():
    order = ("C2", "C6", "C4") if name in ("A", "G") else ("C2", "C4", "C6")
    if all(a in atoms for a in order):
        triplets.append([atoms[a] for a in order])
        names.append(name)
print("complete RNA bases:", len(triplets))
# Nested prefixes cover progressively more of the ribosomal RNA, without changing density.
for n in (100, 300, 1000, 3000, 6000):
    folder = root / str(n)
    folder.mkdir(exist_ok=True)
    coords = []
    for i in range(n):
        copy, base = divmod(i, len(triplets))
        coords.extend(
            tuple(x + (25.0 * copy if k == 0 else 0) for k, x in enumerate(p))
            for p in triplets[base]
        )
    with (folder / "ref.pdb").open("w") as f:
        for i, p in enumerate(coords):
            residue = names[(i // 3) % len(triplets)]
            order = ("C2", "C6", "C4") if residue in ("A", "G") else ("C2", "C4", "C6")
            atom = order[i % 3]
            f.write(
                f"ATOM  {i + 1:5d} {atom:^4s} {residue:>3s} A{i // 3 + 1:4d}    {p[0] * 10:8.3f}{p[1] * 10:8.3f}{p[2] * 10:8.3f}  1.00  0.00           C\n"
            )
        f.write("END\n")
    rng = random.Random(1729)
    phases = [rng.uniform(0, 2 * math.pi) for _ in coords]
    # Smooth center excursions (0.08 nm) plus independent atomic distortions (0.004 nm).
    # 120 frames cross the automatic buffer-rebuild threshold multiple times.
    with (folder / "traj.xyz").open("w") as f:
        for frame in range(120):
            f.write(f"{len(coords)}\n0 0 0\n")
            t = 2 * math.pi * frame / 120
            for i, p in enumerate(coords):
                v = [
                    p[k]
                    + 0.08 * math.sin(t + phases[3 * (i // 3)] + k)
                    + 0.004 * math.sin(3 * t + phases[i] + k)
                    for k in range(3)
                ]
                f.write("C " + " ".join(f"{x:.8f}" for x in v) + "\n")
    (folder / "plumed.dat").write_text(
        f"e: ERMSD ATOMS=1-{3 * n} REFERENCE=ref.pdb NOPBC\nRESTRAINT ARG=e AT=0 KAPPA=1\n"
    )
(root / "metadata.json").write_text(
    json.dumps(
        {
            "source": "https://files.rcsb.org/download/1S72.pdb",
            "complete_bases": len(triplets),
            "seed": 1729,
            "frames": 120,
            "coordinate_units": "nm",
            "sizes": [100, 300, 1000, 3000, 6000],
        },
        indent=2,
    )
    + "\n"
)
