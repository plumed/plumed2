The KENREF module implements *kinetic ensemble refinement*: it biases an ensemble of replicas so that
NMR observables computed **across** the ensemble reproduce experimental values. It provides a single
action, [KENREF](KENREF.md), which is a `Bias` and an `ActionAtomistic` at the same time — it acts on
atomic coordinates directly and takes no `ARG`.

The numerics live in an external library, **KEnRef** (<https://github.com/Smith-Group/KEnRef>), which
this module links against. Three energy models are available and are selected by name with the `MODEL`
keyword:

| `MODEL`    | restrains                                                        |
|:-----------|:-----------------------------------------------------------------|
| `SIGMA`    | cross-relaxation rates (interproton NOE build-up)                 |
| `PLATEAUS` | NOE plateau values                                                |
| `RELAX`    | longitudinal/transverse relaxation over the ensemble              |

Because the observables are ensemble averages, `RELAX` requires at least two replicas; run PLUMED with
`--multi` in the usual way. `SIGMA` and `PLATEAUS` also work with a single replica.

## Installation

This module is not compiled by default. It requires the KEnRef core library, which supplies both the
energy models and the source for this action's constructor. Configure PLUMED with:

```bash
./configure --enable-kenref
```

`--enable-kenref` and `--enable-modules=+kenref` are equivalent: either one enables the whole pathway.
If `kenref_core` is already installed, put its prefix on `PKG_CONFIG_PATH` and configure will find it
(the install ships an `env.sh` that does this for you). Otherwise configure clones KEnRef and delegates
the build to KEnRef's own CMake, into `kenref-deps/` inside the PLUMED build tree. Point it at a local
checkout instead with `--with-kenref-src=DIR`, or override the clone with `KENREF_GIT_URL` /
`KENREF_GIT_TAG`.

KEnRef stores Eigen objects inside its own containers, so **the module must be compiled with the same
Eigen alignment as `libkenref_core`** — that is, the same `-march`/SIMD width. A mismatch is caught at
compile time by a `static_assert` in KEnRef's headers rather than corrupting memory at run time. The
`.pc` files KEnRef installs carry the right `-march`, so this is handled for you; you only need to think
about it if you assemble the flags by hand.

`src/kenref/install.md` documents the build in more detail, and `src/kenref/build-only.sh` automates it.

## Where the source lives

This module is deliberately split. `KEnRefBias.cpp` here holds the **stable** parts of the action —
keyword registration, the PLUMED↔Eigen glue, `calculate()`, and the action registration — and is
maintained in PLUMED. The **one-time constructor**, which changes whenever KEnRef gains an energy model
or alters its parameter schema, is hosted in the KEnRef repository and compiled in through the
forwarding translation unit `KEnRefBias_setup.cpp`.

The point of the split is that adding a model to KEnRef, or changing its input format, does not require
a change to PLUMED. KEnRef remains the source of truth for that half; please do not "inline" the
forwarder.

For the same reason, the regression tests here are deliberately thin: they check that the action parses,
registers, runs and produces its declared components. The numerical correctness of each energy model is
validated in KEnRef's own test suite against reference values from the original R implementation, which
is where new models should be tested.

## Example

An ensemble refinement of GB3 against cross-relaxation rates. `EXP_DATA_FOLDER` holds the experimental
data and the atom-pair definitions; `REF` is the reference structure used for the Kabsch fit.

```plumed
#SETTINGS MOLFILE=regtest/kenref/rt-kenref-parse/gb3.pdb
kenref: KENREF ...
  MODEL=SIGMA
  K=1.0
  N=0.25
  PROTON_MHZ=700.0
  EXP_DATA_FOLDER=./exp_data/
  REF=./gb3.pdb
  ATOMNAME_MAPPING=./gb3.pdb
  GUIDE_ATOMS=39,60,82,101,117,234,241,256,270,284
  MAX_FORCE=999
  FIT_TO_REFERENCE
  SATURATE_FORCES
...

PRINT ARG=kenref.bias,kenref.energy,kenref.rmsd FILE=kenref.out STRIDE=1
```

Note that there is no `ARG` keyword: KENREF biases coordinates, not collective variables.
