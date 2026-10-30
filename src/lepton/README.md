# Vendored scalar Lepton

The upstream source is `libraries/lepton` in OpenMM 8.6.1 (released 2026-09-11), commit
`b399af4725573963b46d6c1083fdcf7a37615857`.
Use it together with the AsmJit 1.9.0 copy in that same OpenMM checkout.
The separately cloned latest AsmJit has a different API and is not the target
of this import.

Check out upstream tag `8.6.1` before importing. From this directory:

```sh
./import.sh /path/to/openmm
# Equivalently:
./import.sh /path/to/openmm/libraries/lepton
```

Python 3 and `patch` are required. `import.py` flattens the upstream scalar
sources, wraps their namespace in `PLMD`, replaces the JIT feature macro, and
applies `plumed.patch` with no fuzz. It stages the complete import before
replacing any vendored sources. Neither importer stages changes in Git.
Upstream copyright notices are retained. Vector expressions are intentionally
excluded: PLUMED uses the scalar API and does not expose the separate upstream
`CompiledVectorExpression` API.

`plumed.patch` records PLUMED's changes relative to the transformed upstream:

- Additional functions and derivatives, including `nandelta` and the extended
  trigonometric functions.
- PLUMED's delta/step semantics and mathematical constants.
- Thread-safe initialization of the parser's operation map.
- The `PLUMED_USE_ASMJIT=yes/no` runtime switch (JIT is enabled by default when
  compiled in).
- An opaque runtime owner, keeping AsmJit headers out of Lepton's public headers.
- Interpreter fallback when x86 AVX is unavailable, and error reporting for JIT
  compilation failures.
- Releasing generated code before rebuilding it on variable rebinding.
- Floating-point ARM64 step/delta comparisons, including negative zero and NaN.

To inspect or update the patch, generate an unpatched baseline:

```sh
python3 import.py /path/to/openmm --unpatched --output-dir /tmp/lepton-baseline
```

Compare that baseline with the final `.h`/`.cpp` files and regenerate a unified
patch using paths `a/filename` and `b/filename`. For a reproduction check, import
into another directory with `--output-dir` and compare its sources byte for byte.
Future upstream updates require reviewing patch application and semantics.

Native ARM64 validation includes `rt-lepton`, `rt-lepton-asmjit`,
`rt-make-lepton-jit`, and the switching-function tests. `rt42` provides detailed
C++ versus Lepton timings. Configure accepts x86/x64 and AArch64; other
architectures continue to use the interpreter.

## Validation on Apple Silicon

With this port, the basic regression suite completed with 241 passing tests,
59 dependency-related skips, and no failures. The focused JIT regression also
passed as an x86_64 executable under Rosetta and in a standalone ARM64 build
with JIT support compiled out. `rt42` passed with JIT both enabled and disabled.
Both importers reproduced their generated sources byte for byte.

Seven paired `rt42` runs gave the following median calculation times (108
frames; these are local measurements, not a general performance guarantee):

| Action | Interpreter | ARM64 JIT |
|---|---:|---:|
| C++ `c` | 21.3 ms | 21.4 ms |
| Lepton `d` | 46.0 ms | 28.8 ms |
| Lepton `dfast` | 41.9 ms | 26.0 ms |

`dfast` became about 1.6 times faster and its time relative to C++ fell from
about 1.96 to 1.22. The native build remains configured with AsmJit enabled;
set `PLUMED_USE_ASMJIT=no` to compare the interpreter.
