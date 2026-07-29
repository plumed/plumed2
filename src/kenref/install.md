# Building PLUMED with the KEnRef module (from the PLUMED side)

The `kenref` module adds the `KENREF` bias. This directory carries only the **bridge** (`KEnRefBias` + the
bridge implementations); the numerical **core library lives in a separate repo**
(<https://github.com/Smith-Group/KEnRef.git>). So a PLUMED build with KEnRef always: (1) obtains
`kenref_core` — **delegated to the KEnRef repo, because kenref knows how to build kenref** — then (2) builds
PLUMED here with its own autotools.

> The `kenref` module is `default-off`: a plain `./configure && make` builds PLUMED **without** it and with
> no dependency on KEnRef. You opt in with `--enable-kenref`.

To start from the **KEnRef** side instead, see `INSTALL.md` in the KEnRef repo.

---

## Option A — the two scripts (recommended)

Both are thin orchestrators sharing `_common.sh`; power users need only a command or two.

```bash
# PLUMED with the kenref module, NO gromacs:
src/kenref/build-only.sh -y                          # find kenref (or clone), build+install PLUMED
src/kenref/build-only.sh --kenref-src ~/KEnRef -y

# PLUMED + a GROMACS 2025.x batched with it:
src/kenref/build-and-batch.sh -y                     # also fetch+patch+build GROMACS 2025.x
src/kenref/build-and-batch.sh --gromacs-src ~/gromacs-2025.4 -y   # use a provided gromacs
```

What they do: **kenref_core** — reused if found via `pkg-config`, else the KEnRef repo is located
(`--kenref-src`) or **auto-cloned** and **its** `build.sh` builds+installs it (with the plumedinterface
exported). Starting from PLUMED never assumes KEnRef is already present — an empty machine just clones it. **PLUMED** — built by its own autotools (`--enable-kenref
--enable-modules=+kenref`). **GROMACS** (batch script) — a provided source, else fetched to a
plumed-related dir; `plumed patch`ed (native `GMX_USE_PLUMED` alone is *not* full integration) then built.
Toolchain from the environment (`CXX=mpicxx CXXFLAGS="-stdlib=libc++" …`); `--accel` auto-detects.

## Option B — PLUMED's build system alone (no script)

`--enable-kenref` is **self-sufficient**: if `kenref_core` isn't installed, `configure` clones the KEnRef
repo and **delegates the build to KEnRef's CMake**, then links it.

```bash
autoreconf --force
./configure --enable-kenref CXXFLAGS="-stdlib=libc++ -O3 -std=c++17 -march=skylake-avx512"
make -j && make install
```

Controls (only if you don't want the default clone): `--with-kenref-src=DIR`, `KENREF_GIT_URL`/`KENREF_GIT_TAG`,
`KENREF_ACCEL`. Already have `kenref_core` installed? Put it on `PKG_CONFIG_PATH` and `configure` reuses it.
`--enable-kenref` and `--enable-modules=+kenref` converge (either enables the whole pathway); a
`kenref_core >= <min>` floor rejects a stale core at configure time.

## Where it installs

Both options above install to the **same** places — each component where its own build system puts it by
default, so nothing here is a surprise to someone who already knows that project:

| component | default prefix | set by |
|---|---|---|
| PLUMED | `/usr/local` | autotools' default (`--prefix`) |
| `kenref_core` | `/usr/local/kenref` | KEnRef's CMake (`KENREF_INSTALL_BASE/kenref`) |
| GROMACS (batched) | `/usr/local/gromacs` | GROMACS's CMake |

Override any of them: `--plumed-prefix DIR`, `--kenref-prefix DIR`, `--gromacs-prefix DIR` (Option A), or
`--prefix` (Option B). Build directories always stay inside the source trees — they are artefacts, not
installs, and are never written to a system path.

**Privilege:** these defaults are system paths, so the *install step alone* runs under `sudo` — and only
when the destination is not writable. Redirect the prefixes into your own space and nothing ever asks for a
password:

```bash
src/kenref/build-and-batch.sh --kenref-src ~/KEnRef --gromacs-src ~/gromacs-2025.4 \
    --plumed-prefix ~/.local --kenref-prefix ~/.local/kenref --gromacs-prefix ~/.local/gromacs -y
```

## Self-contained binaries (rpath)

With the `kenref` module enabled, `configure` turns **rpath ON by default** (PLUMED's global default is off).
This bakes the runtime library directories into the RUNPATH of **both** the `plumed` executable and
`libplumedKernel.so`, so the **installed** binaries run with no `LD_LIBRARY_PATH` or module setup — including
`sudo make install` itself, whose doc-generation step runs the freshly built `plumed` (with the environment
stripped by `sudo`).

Why it matters even for a stock toolchain: the kernel links **`libkenref_core`**, which is installed in its
own prefix (e.g. `/usr/local/kenref/lib`) that is **not** on the system loader path — so without rpath the
binaries can't find it at run time. If you *also* choose a **non-system C++ runtime** (e.g. build with
`-stdlib=libc++`), the same applies to `libc++.so.1`.

rpath records **`LIBRARY_PATH`**, so put your non-system runtime library dirs on it at configure time — the
`libkenref_core` lib dir, plus, for a libc++ build, the dir holding `libc++.so.1` (with LLVM that is the
`…/lib/x86_64-unknown-linux-gnu` **subdir**, not the parent `…/lib`):

```bash
export LIBRARY_PATH="/usr/local/kenref/lib:$LIBRARY_PATH"                 # libkenref_core
export LIBRARY_PATH="$LLVM/lib/x86_64-unknown-linux-gnu:$LIBRARY_PATH"    # only if building with libc++
autoreconf --force
./configure --enable-kenref --with-kenref-src=~/KEnRef CXXFLAGS="-O3 -std=c++17 -march=native"
make -j && sudo make install
plumed --version        # runs with no module / no LD_LIBRARY_PATH
```

- If every linked library (including `libkenref_core`) is already on the system loader path — e.g. installed
  under a prefix covered by `ldconfig` — rpath is harmless and simply records nothing new.
- **Opt out** with `--disable-rpath` (e.g. if you deploy runtime paths via environment modules instead); then
  the installed binaries — and `sudo make install` — need those non-system libs on `LD_LIBRARY_PATH`.
- At **build time**, `libc++` must also be on `LD_LIBRARY_PATH` (configure runs test programs that have no
  rpath yet) — normally provided by loading your compiler/toolchain module.

## Batching GROMACS (2025.x)

Running KEnRef-on-PLUMED-**on-GROMACS** needs GROMACS **2025.x (pref. 2025.4)**, batched with `plumed patch`
(native `-DGMX_USE_PLUMED=ON` compiles but isn't full integration on its own). `build-and-batch.sh` does this;
manually:

```bash
cd /path/to/gromacs-2025.x           # a source DEDICATED to the PLUMED path (not the kenref-gmx one)
plumed patch -p -e gromacs-2025.0 --mode runtime
cmake -S . -B build -DGMX_USE_PLUMED=ON -DGMX_MPI=ON && cmake --build build && cmake --install build
```

## Notes
- The `Smith-Group/KEnRef` fork/URL and this in-tree module are **temporary**, until `KEnRefBias` is merged
  upstream — after which the module ships with PLUMED and `kenref_core` is a plain `pkg-config` dependency.
- KEnRef installs its own env helpers (`env.sh`, TCL modulefiles); see the KEnRef `INSTALL.md`.
