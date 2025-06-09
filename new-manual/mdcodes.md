## Codes interfaced with PLUMED

PLUMED can be incorporated into an MD code and used to analyze or bias a molecular dynamics run on the fly.
Some MD codes already include calls to the PLUMED library
and be PLUMED-ready in their original distribution.
As far as we know, the following MD codes can be used with PLUMED out of the box:

- [Amber](http://ambermd.org/), pmemd module, since version 20.
- [AmberTools](http://ambermd.org/), sander module, since version 15.
- [CP2K](http://www.cp2k.org), since Feb 2015.
- [ESPResSo](http://espressomd.org), in a version that has been patched with PLUMED can be found
  [here](http://davidebr.github.io/espresso/).
- [PINY-MD](http://github.com/TuckermanGroup/PINY), in its plumed branch.
- [IPHIGENIE](http://sourceforge.net/projects/iphigenie/).
- [AceMD](http://www.multiscalelab.org/acemd/), see [this link](https://github.com/tonigi/ACEMD-PLUMED).
- [OpenMM](http://openmm.org), using the [openmm-plumed plugin](http://github.com/peastman/openmm-plumed).
- [DL_POLY4](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx).
- [VNL-ATK](https://quantumwise.com), see [this link](https://docs.quantumwise.com/tutorials/metadynamics_with_plumed/metadynamics_with_plumed.html).
- [ABIN](https://github.com/PHOTOX/ABIN).
- [i-pi](https://github.com/i-pi/i-pi).
- [LAMMPS](https://lammps.sandia.gov/) since Nov 2018.
- [Yaff](https://github.com/molmod/yaff), since Jul 2019.
- [DFTB+](https://www.dftbplus.org/), since release 20.1.
- [Metalwalls](https://gitlab.com/ampere2/metalwalls)
- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [GPUMD](https://gpumd.org/)
- [GROMACS](https://www.gromacs.org/) as of version 2025, with limited support (no replica exchange, no [ENERGY](ENERGY.md) collective variable, no lambda dynamics).

We also provide patches for [namd](https://www.ks.uiuc.edu/Research/namd/) and [quantum espresso](https://www.quantum-espresso.org) as well as a 
patch for [GROMACS](https://www.gromacs.org/) that you should use in place of the native implementation if you want to do replica exchange, 
use the [ENERGY](ENERGY.md) collective variable or do lambda dynamics.

If you maintain another MD code that is PLUMED-ready let us know and we will add it to this list.

The status of the interface between some of these MD codes and PLUMED are tested [here](https://plumed-testcenter.github.io).  You can find examples of how we build
the interfaces between PLUMED and these codes on that site.  However, you will in general need to refer to the documentation of the MD code to know how to use it with the latest PLUMED release.

PLUMED can also be used as a [tool](module_cltools.md) for post processing the results from molecular dynamics
or enhanced sampling calculations.  Notice that PLUMED can be used as an analysis tool from the following packages:

- [PLUMED-GUI](http://github.com/tonigi/vmd_plumed) is a [VMD](http://www.ks.uiuc.edu/Research/vmd/) plugin that computes PLUMED collective variables.
- [HTMD](http://www.htmd.org/) can use PLUMED collective variables for analysis.
- [OpenPathSampling](http://openpathsampling.org/), using the [PLUMED Wrapper for OpenPathSampling](https://e-cam.readthedocs.io/en/latest/Classical-MD-Modules/modules/OpenPathSampling/ops_plumed_wrapper/readme.html).

## GROMACS and PLUMED with GPU

Since version 4.6.x GROMACS can run in an hybrid mode making use of both
your CPU and your GPU (either using CUDA or OpenCL for newer versions of
GROMACS). The calculation of the short-range non-bonded interactions is
performed on the GPU while long-range and bonded interactions are at the
same time calculated on the CPU. By varying the cut-off for short-range
interactions GROMACS can optimize the balance between GPU/CPU loading
and obtain amazing performances.

GROMACS patched with PLUMED takes into account PLUMED in its load-balancing,
adding the PLUMED timings to the one resulting from bonded interactions and long-
range interactions. This means that the CPU/GPU balance will be optimized
automatically to take into account PLUMED!

It is important to notice that the optimal setup to use GROMACS alone
on the GPU or GROMACS + PLUMED can be different, try to change the number
of MPI/OpenMP processes (\ref Openmp) used by GROMACS and PLUMED to find
optimal performances. Remember that in GROMACS multiple MPI threads
can use the same GPU:

i.e. if you have 4 cores and 2 GPU you can:

- use 2 MPI/2GPU/2OPENMP:

```bash
export PLUMED_NUM_THREADS=2
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
```

- use 4 MPI/2GPU:

```bash
export PLUMED_NUM_THREADS=1
mpiexec -np 4 gmx_mpi mdrun -nb gpu -ntomp 1 -pin on -gpu_id 0011
```

Of notice that since plumed 2.5 and gromacs 2018.3 the number of openMP threads can automatically set by gromacs (so `PLUMED_NUM_THREADS` is not needed, and the number of OpenMP threads used by plumed is set by `-ntomp`)

```bash
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
```

