This plugin allows you to use the USEGPU flag with certain actions in PLUMED. This flag turns on an experimental GPU parallized version of the 
command. GPU parallelism in PLUMED has been implemented using [openACC](https://www.openacc.org) and is currently experimental. We are actively working
on these features at the moment. __There is thus no guarantee that the GPU accelerated versions of actions are any faster than 
the CPU versions.__ If you have experimented with these features on your own calculations we would love to hear from you (even
if your experience was negative.)

## [Experimental] Compiling the openaccPTM plugin

The openaccPTM plugin is compiled separately to the main PLUMED code. To compile the plugin you need to have the 
[NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) and the plumed library avaiable in your path.  You then compile 
the plugin by issuing the following commands:

```bash
cd plugins/openaccPTM
./configure 
make clean
make
```  

If PLUMED has been compiled with MPI, `configure` sets `NVCXX=mpic++`. The wrapper must call `nvc++` and not the compiler
that was used to build MPI: the Makefile takes care of this by exporting `OMPI_CXX=nvc++`, `MPICH_CXX=nvc++` and `I_MPI_CXX=nvc++`
(for OpenMPI, MPICH and Intel MPI respectively) to the compile and link commands, so a plain `make` is enough. If these variables are
already set in your environment they are not overridden. With OpenMPI you can check what the wrapper calls with `OMPI_CXX=nvc++ mpic++ --showme`.

If PLUMED has been compiled with a GCC that is not the default one known by `nvc++` (for example one loaded via `module load`),
you may need to point `nvc++` to it, so that the same C++ standard library is used. After running `./configure` you can do:

```bash
GCCROOT=$(dirname $(dirname $(which g++)))
echo "CXXFLAGS += --gcc-toolchain=$GCCROOT" >> Makefile.conf
echo "LDSHARED += --gcc-toolchain=$GCCROOT" >> Makefile.conf
```

It is worth running these commands more than once if the compilation doesn't work the first time. If they complete correctly a dynamic libarary called `plumedOpenACC.so` 
will be created in the plumed2/plugins/openaccPTM directory.

We have tested this module with the following versions of the NVIDIA HPC SDK compiler and found that it works:

- [24.3](https://developer.nvidia.com/nvidia-hpc-sdk-243-downloads);
- 24.5, on GPUs (Leonardo Booster nodes, A100);
- 25.11, on the CPU (`ACCTARGET=host`), which is the version we used to run the regression tests in the GitHub Actions CI.
- 26.9, on the CPU (`ACCTARGET=host`), which is the version we use to run the regression tests in the GitHub Actions CI.

It currently does not work with NVHPC 26.5 compiler.

Note that we recommend running the test suite after compilation.  You can do this by issuing the following commands:

```bash
cd plugins/openaccPTM/regtest
make
```

## [Experimental] Using PLUMED with the openaccPTM plugin

To use openACC in PLUMED you need to load the dynamic libary that you compiled in the previous section.  You can do this as follows:

```plumed
LOAD FILE=/path/to/plumedOpenACC.so
```

where `/path/to` is the directory that contains the compiled plugins. If you choose not to move the library after compilation this will 
be the `plumed2/plugins/openaccPTM` directory. To use the plugin to calculate an action using the GPU you add the `USEGPU` keyword to the list of keywords
for that action.  For example, to calculate a contact matrix using the GPU you would use the following input:

```plumed
LOAD FILE=/path/to/plumedOpenACC.so
CONTACT_MATRIX GROUP=1-1000 SWITCH={EXP D_0=0.2 R_0=0.1 D_MAX=0.66} USEGPU
``` 

When PLUMED runs with MPI (e.g. an MD code with several MPI ranks), the actions with `USEGPU` are computed by rank 0 only,
on its GPU, and the results are sent to the other ranks. The other ranks do not use the GPU and wait for rank 0.

The next section contains a list of the actions that can take the USEGPU command in input and that can thus be run on the GPU using this plugin.

## List of actions that can be called with the USEGPU option:

 - module:
   - ACTION

 - colvar:
   - [ANGLE](ANGLE.md)
   - [DIPOLE](DIPOLE.md)
   - [DISTANCE](DISTANCE.md)
   - [PLANE](PLANE.md)
   - [POSITION](POSITION.md)
   - [TORSION](TORSION.md)
 - crystdistrib:
   - ~~[QUATERNION_BOND_PRODUCT_MATRIX](QUATERNION_BOND_PRODUCT_MATRIX.md)~~ setup, but deactivated
 - secondarystructure:
   - [SECONDARY_STRUCTURE_DRMSD](SECONDARY_STRUCTURE_DRMSD.md), and in particular:
     - [ALPHARMSD](ALPHARMSD.md) only with **TYPE=DRMSD**
     - [ANTIBETARMSD](ANTIBETARMSD.md) only with **TYPE=DRMSD**
     - [PARABETARMSD](PARABETARMSD.md) only with **TYPE=DRMSD**
 - volumes:
   - [AROUND](AROUND.md)
   - [INCYLINDER](INCYLINDER.md)
   - [INSPHERE](INSPHERE.md)
- adjmat
   - [CONTACT_MATRIX](CONTACT_MATRIX.md)
- function
   - [LESS_THAN](LESS_THAN.md)
   - [MORE_THAN](MORE_THAN.md)
   - [BETWEEN](BETWEEN.md)
   - [COMBINE](COMBINE.md)
- gridtools
   - [KDE](KDE.md)
   - [SPHERICAL_KDE](SPHERICAL_KDE.md)
- matrixtools
   - [MATRIX_PRODUCT](MATRIX_PRODUCT.md)
   - [MATRIX_VECTOR_PRODUCT](MATRIX_VECTOR_PRODUCT.md)
- isdb
   - [SAXS](SAXS.md) and [SANS](SANS.md), see below
- symfunc
   - [Q1](Q1.md)
   - [Q3](Q3.md)
   - [Q4](Q4.md)
   - [Q6](Q6.md)
   - [LOCAL_AVERAGE](LOCAL_AVERAGE.md)
   - [FCCUBIC](FCCUBIC.md)
   - [TETRAHEDRAL](TETRAHEDRAL.md)
   - [SIMPLECUBIC](SIMPLECUBIC.md)
   - [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)
   - [COORDINATION_SHELL_AVERAGE](COORDINATION_SHELL_AVERAGE.md)
   - [LOCAL_Q1](LOCAL_Q1.md)
   - [LOCAL_Q3](LOCAL_Q3.md)
   - [LOCAL_Q4](LOCAL_Q4.md)
   - [LOCAL_Q6](LOCAL_Q6.md)

## SAXS and SANS

[SAXS](SAXS.md) and [SANS](SANS.md) do not use the task manager: the plugin provides the `OPENACC` backend for their Debye
sum, which is the expensive part of the calculation (it scales with the square of the number of atoms or beads). Form
factors, solvent accessible surface (ONEBEAD), deuteration (SANS), resolution function and metainference are computed by
PLUMED as usual.

````
LOAD FILE=/path/to/plumedOpenACC.so

MOLINFO STRUCTURE=template.pdb
SAXS ...
  LABEL=saxs
  ATOMS=1-3708
  ATOMISTIC
  USEGPU
  QVALUE1=0.03 QVALUE2=0.06 QVALUE3=0.09
... SAXS
````

(This example is not highlighted only because the manual is built without the plugin.)

- `DEVICEID` selects the GPU; by default the one assigned to PLUMED by the MD code (if any) or the default OpenACC device.
- Differently from the other actions, with MPI the rows of the Debye sum are split over the ranks, each one using its own
  GPU (or sharing one), and the results are summed. With `SERIAL` only the first rank computes on its GPU and the results
  are shared with the others.
- With `DOSCORE` the derivatives of the intensities stay on the GPU and only the forces are copied back. Without it the
  derivatives of every intensity with respect to every particle are copied back, which costs more for many q values.
- `GPU_BACKEND=OPENACC` is needed only if another backend for SAXS/SANS (from a different plugin) is loaded too.

The form factors and the q values are uploaded only when they change; the positions are uploaded at every step.
Each GPU thread computes the rows of particle i for a chunk of 8 q values, looping over all the other particles, so
each pair is computed twice, but there are no atomic operations and the derivative of particle i is complete in its thread.
For each pair the distance is computed once per chunk and sin(qr) and cos(qr) are obtained by angle addition along the
q values, as on the CPU (exact sin and cos when the q values are far from a uniform grid or there are fewer than 4).
For small systems the loop over the other particles is also split in tiles, whose results are summed in a fixed order,
to have enough threads. The result does not depend on the scheduling, and the rows are summed on the GPU for each q value.
The chunk size and the number of threads below which tiles are used can be changed at compile time with
`-DPLUMED_SAXS_ACC_QCHUNK=...` (default 8) and `-DPLUMED_SAXS_ACC_TARGET_THREADS=...` (default 131072).

### Performance

Time per step of SAXS (ms) on a Leonardo Booster node (Xeon 8358 with 32 cores, 4 A100), 50 q values, double precision.
The CPU reference uses the 32 cores (OpenMP); "score" is with `DOSCORE`, "comp" without it (the intensities are used by
another action, e.g. a bias):

| system | particles | CPU 32 cores | 1 GPU score | 1 GPU comp | 4 GPUs score | 4 GPUs comp |
|---|---|---|---|---|---|---|
| protein, `ATOMISTIC` | 3439 | 179 | 6.2 | 7.6 | 2.3 | 9.2 |
| 16 copies, `ATOMISTIC` | 55024 | 38400 | 1360 | 1390 | 383 | 508 |
| 48 lysozymes, `MARTINI` | 14592 | 2810 | 98 | 105 | 33 | 68 |
| 48 copies, `ONEBEAD` | 9888 beads, 165072 atoms | 1190 | 52 | 75 | 20 | 56 |

The results agree with the CPU to 1e-12 (relative).

- One GPU is 20-40 times faster than the 32 cores, more for small numbers of q values.
- Several ranks, each with its own GPU, pay off with `DOSCORE` for large systems (more than about 10000 particles).
  Without `DOSCORE` every rank needs the derivatives of all the intensities with respect to all the particles, which
  are summed over the ranks: in this case a single GPU per simulation (or per replica) is the best choice.
- `ONEBEAD` without `DOSCORE` also stores the derivatives with respect to all the atoms, which is done on the CPU.

The tests in `regtest/targetSAXS` are the SAXS/SANS tests of `regtest/isdb` run with `USEGPU`: they take their input files
and their (CPU) references from `regtest/isdb`.
