This plugin computes the Debye sum of [SAXS](SAXS.md) and [SANS](SANS.md) on the GPU with [OpenACC](https://www.openacc.org).
It provides the `OPENACC` backend that is used when SAXS or SANS are given the `USEGPU` flag. Only the Debye sum, which is
the expensive part of the calculation (it scales with the square of the number of atoms or beads), runs on the GPU: form
factors, solvent accessible surface (ONEBEAD), resolution function and metainference are computed by PLUMED as usual.

## [Experimental] Compiling the openaccSAXS plugin

You need the [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) and the plumed library available in your path. Then:

```bash
cd plugins/openaccSAXS
./configure
make
```

This creates `plumedOpenACCSAXS.so` (`.dylib` on macOS). The notes in plugins/openaccPTM/README.md about MPI and about the
GCC toolchain used by `nvc++` apply also here. `make ACCTARGET=host` compiles the kernels for the CPU (serially, no GPU
needed): this is what the CI does to run the regression tests.

The tests in `regtest` are the SAXS/SANS tests of `regtest/isdb` with `USEGPU`, compared with the same (CPU) references:

```bash
cd plugins/openaccSAXS/regtest
make
```

## [Experimental] Using the plugin

```plumed
LOAD FILE=/path/to/plumedOpenACCSAXS.so

MOLINFO STRUCTURE=template.pdb
SAXS ...
  LABEL=saxs
  ATOMS=1-3708
  ATOMISTIC
  USEGPU
  QVALUE1=0.03 QVALUE2=0.06 QVALUE3=0.09
... SAXS
```

- `DEVICEID` selects the GPU; by default the one assigned to PLUMED by the MD code (if any) or the default OpenACC device.
- With MPI the rows of the Debye sum are split over the ranks, each one using its own GPU (or sharing one), and the
  results are summed. With `SERIAL` only the first rank computes on its GPU and the results are shared with the others.
- With `DOSCORE` the derivatives of the intensities stay on the GPU and only the forces are copied back. Without it the
  derivatives of every intensity with respect to every particle are copied back, which costs more for many q values.
- `GPU_BACKEND=OPENACC` is needed only if another backend (a different plugin) is loaded too.

## Implementation

The form factors and the q values are uploaded only when they change; the positions are uploaded at every step.
Each GPU thread computes a row (particle i, q value k) of the sum looping over all the other particles, so each pair is
computed twice, but there are no atomic operations and the derivative of particle i is complete in its thread.
The rows are summed on the GPU for each q value.

There is no guarantee that the GPU is faster than the CPU for small systems: each step pays for the transfers and the
kernel launches. Please let us know about your experience.
