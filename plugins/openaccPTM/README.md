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

When PLUMED runs with MPI (e.g. an MD code with several MPI ranks, or `mpirun -np 4 plumed driver`), the actions with `USEGPU`
are computed by the first rank, on its GPU, and the results are sent to the other ranks. With `MULTIGPU` (accepted by the actions
and by the shortcuts that take `USEGPU`) the tasks are split in blocks between the ranks, each rank computes its block on its GPU
(see below), and the results are summed, as for the CPU version; `SERIAL` restores the default. `MULTIGPU` is not the default
because, for most actions, a call is dominated by work on the CPU (copies of the input, of the results and of the forces) that
does not decrease with the number of ranks, while the sums over the ranks add more: on Leonardo (A100) splitting
[COORDINATIONNUMBER](COORDINATIONNUMBER.md) or [Q6](Q6.md) of 20000 atoms or a [DISTANCES](DISTANCES.md) of 100000 between
2 or 4 GPUs was slower than using one (measured with the contact matrices, before the sums described below). It pays off when
the calculation on the GPU dominates.

The forces of the matrices are gathered in one pass, adding with atomics the forces that several rows touch. The low memory
version (`AVOID_THREAD_GATHER`), which computes the matrix again column by column, was the default on the GPU, but it makes the
backward pass about two times slower ([CONTACT_MATRIX](CONTACT_MATRIX.md) of 20000 atoms on an A100).

### Coordination numbers and Steinhardt parameters

With `USEGPU`, [COORDINATIONNUMBER](COORDINATIONNUMBER.md) and [Q1](Q1.md), [Q3](Q3.md), [Q4](Q4.md), [Q6](Q6.md) (and so
[LOCAL_Q1](LOCAL_Q1.md), ..., [LOCAL_AVERAGE](LOCAL_AVERAGE.md) that use their output) compute the sums over the neighbours
of each atom on the GPU in a single pass, without building the contact matrix (`COORDINATIONNUMBER_SUMSACC`, `Q6_SUMSACC`, ...,
see `src/NeighborSums.cpp`): only the positions and the linked cells go to the GPU, and only the sums and the forces on the atoms
come back. These are the actions that the CPU uses without `USEGPU` (`src/symfunc/NeighborSums.h` in PLUMED), with the same
kernels: the linked cells are built on the CPU and searched on the GPU, in any box, one GPU thread per atom. With up to 512
atoms in `SPECIESA` (e.g. one ion against all the water molecules) the cells are not built: the GPU threads share the atoms of
`SPECIESB` in blocks, and the partial sums of each atom are added in order by a second kernel. With `SPECIES` the force on each
atom is gathered from its own neighbours (no atomics); with `SPECIESA` and `SPECIESB` the forces on the second group are added
with atomics.
For Q1, ..., Q6 the coordination numbers (the denominators) are computed in the same pass as the sums of the spherical harmonics,
and the forces of the spherical harmonics are contracted pair by pair, without storing the derivatives of the $2(2l+1)$ components.
The results are those of the CPU version (in double precision). On one A100 (8 CPU threads), ms per step:
COORDINATIONNUMBER 3.6 (20000 atoms) and 16 (100000), Q6 12.6 and 58; on 32 CPU cores 4.0, 22, 15.7 and 90.
This is not used with `MASK`, and for [COORDINATION_MOMENTS](COORDINATION_MOMENTS.md): the contact matrix is then computed
on the GPU as before. The actions of the contact matrix (e.g. `q6_mat`, `q6_sh`) are still created, and are computed (on the GPU)
only when another action uses them.
For Q1, ..., Q6 the norm of the vector of the $2l+1$ sums (`Q6_NORM`, as [COMBINE](COMBINE.md) and the square root) is always
computed on the CPU, which is faster than moving the $2(2l+1)$ vectors to the GPU and back.

### Secondary structure

With `USEGPU` and `TYPE=DRMSD`, the number of segments that are similar to a helix or a sheet computed by [ALPHARMSD](ALPHARMSD.md),
[ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md) (the value, or the `lessthan` component) comes from
`ALPHARMSD_SCALAR_ACC`, `ANTIBETARMSD_SCALAR_ACC` and `PARABETARMSD_SCALAR_ACC` (`src/SecondaryStructureCount.cpp`), the GPU
version of the single action that PLUMED uses on the CPU
(`src/secondarystructure/SecondaryStructureCount.h`), with the same functions for a segment: a first kernel lists the pairs of strands
within `STRANDS_CUTOFF` (1 nm by default) and whose switching function is not certainly zero (see the documentation of
ANTIBETARMSD), one GPU thread per listed segment computes its distance from the references, the switching function
and the derivatives, which are added to the atoms with atomics, and a last kernel sums the values and the virial. The forces
need no second pass. The vectors of the expanded input are computed (on the GPU) only when another action uses them.

On one A100 (8 CPU threads), for 144 peptides of 8 residues (20160 atoms, 370800 pairs of strands, about 5000 of them
within `STRANDS_CUTOFF=1`), [ALPHARMSD](ALPHARMSD.md), ANTIBETARMSD and PARABETARMSD together take 1.8 ms per step
(6.1 ms on 32 CPU cores, 12.5 ms on 8). For a single protein,
with a few hundred residues, the CPU is faster: the time on the GPU is then set by the cost of starting the kernels and of
the transfers (about 0.3 ms per step for ALPHARMSD).

### Which GPU is used

All the actions of the plugin (SAXS and SANS included) use the same GPU as the MD code, so that they work on a node with
several GPUs and do not interfere with an MD code that keeps its data on the GPU:

1. the GPU assigned to PLUMED by the MD code (with the `setGpuDeviceId` command, as the GROMACS and NAMD patches do when
   GROMACS or NAMD use CUDA), otherwise
2. the GPU the MD code is using (the one whose CUDA context is current when PLUMED is called), otherwise
3. the one set with `ACC_DEVICE_NUM`, otherwise
4. with several GPUs on the node, GPU number (local rank % number of GPUs), where the local rank is the rank of the process within
   its node, as set by the MPI launcher (`OMPI_COMM_WORLD_LOCAL_RANK`, `MPI_LOCALRANKID`, `SLURM_LOCALID`, ...).
   This is what happens with [driver](driver.md), where there is no MD code, otherwise
5. the first GPU. With `CUDA_VISIBLE_DEVICES` (or `srun --gpus-per-task=1`) each rank can be given its own GPU.

The choice is made by PLUMED itself (`PLMD::GpuDevice`, in `src/core/GpuDevice.h`), so it is the same for all the plugins.
The GPU of an action is chosen the first time it is computed, and it does not change afterwards. It is written in the log,
in the same way for all the actions, with the rule that was applied, e.g. `d runs on GPU 1 (assigned by the MD code)` for an
action labelled `d` (or `d runs on the host` when OpenACC sees no GPU). After each call the MD code finds its own GPU
(its CUDA context) as it left it.

The cudaCoord plugin chooses its GPU with the same rule. `make check-gpu` (on a node with GPUs) checks that the actions
run on the expected GPU and that the MD code keeps its own; `make -C gpucheck PLUGINS="/path/plumedOpenACC.so /path/CudaCoordination.so"`
checks the two plugins loaded together.

For developers: new actions should follow the same behaviour, as described in `src/ACCDevice.h`. The actions that use the
OpenACC task manager get it for free. Code that uses OpenACC outside of it (as the SAXS/SANS backend) keeps an `ACC::Device`,
chooses the GPU with `resolve(action)` before using it, and runs every use of OpenACC (kernels, transfers, allocations,
and the destructor that frees its device memory) inside an `ACC::DeviceScope`.

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
     - [ALPHARMSD](ALPHARMSD.md) only with **TYPE=DRMSD** (the number of segments: `ALPHARMSD_SCALAR_ACC`, see below)
     - [ANTIBETARMSD](ANTIBETARMSD.md) only with **TYPE=DRMSD** (the number of segments: `ANTIBETARMSD_SCALAR_ACC`, see below)
     - [PARABETARMSD](PARABETARMSD.md) only with **TYPE=DRMSD** (the number of segments: `PARABETARMSD_SCALAR_ACC`, see below)
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
   - [COORDINATIONNUMBER](COORDINATIONNUMBER.md)
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

- The GPU is chosen as for the other actions (see above): there is no `DEVICEID` keyword.
- Differently from the other actions, which need `MULTIGPU` for this, with MPI the rows of the Debye sum are split over the ranks, each one using its own
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
