\page Performances Performances 

In this page we collect hints on how to use the features available in PLUMED to speed
up your calculations. Please note that PLUMED performs many different tasks, it
can calculate a number of different collective variables, functions of collective 
variables, bias, on-the-fly analysis, etc in a way that is compatible with a number of
different molecular dynamics codes. This means that there cannot be a single 
strategy to speed up all the possible calculations. 

\ref performance-optimization "Here" 
you can find a step-by-step tutorial on optimizing PLUMED performances, discussing some of the topics
below in more detail and using practical examples.

PLUMED makes use of MPI and OpenMP to parallelize some of its functions, try to always
compile it with these features enabled. Furthermore, newer compilers with proper optimization 
flags can provide a dramatic boost to performances.

PLUMED collects atoms from an external code and sends back forces, so it is key to minimize
the effect of PLUMED on highly parallel calculations to keep to the minimum the number of atoms 
used by PLUMED at every calculation step. The less is the number of atoms you need to send 
to PLUMED the less will be the overhead in the communication between PLUMED and the code.

In the following you can find specific strategies for specific calculations, these could
help in taking the most by using PLUMED for your simulations.

- \subpage GMXGPU 
- \subpage Metadyn
- \subpage MTS
- \subpage Multicolvar 
- \subpage Neighbour 
- \subpage Openmp
- \subpage Secondary
- \subpage Time
- \subpage Lepton

\page GMXGPU GROMACS and PLUMED with GPU

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

\verbatim
export PLUMED_NUM_THREADS=2
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
\endverbatim

- use 4 MPI/2GPU:

\verbatim
export PLUMED_NUM_THREADS=1
mpiexec -np 4 gmx_mpi mdrun -nb gpu -ntomp 1 -pin on -gpu_id 0011
\endverbatim

Of notice that since plumed 2.5 and gromacs 2018.3 the number of openMP threads can automatically set by gromacs (so PLUMED_NUM_THREADS is not needed, and the number of OpenMP threads used by plumed is set by -ntomp)

\verbatim
mpiexec -np 2 gmx_mpi mdrun -nb gpu -ntomp 2 -pin on -gpu_id 01
\endverbatim


\page Metadyn Metadynamics

Metadynamics can be sped up significantly using grids,
which are activated setting the GRID_MIN and GRID_MAX keywords of \ref METAD.
This makes addition of a hill to the list a bit slower (since
the Gaussian has to be evaluated for many grid points)
but the evaluation of the potential very fast. Since
the former is usually done every few hundred steps, whereas the latter 
typically at every step, using grids will make the simulation
 faster in particular for long runs.

Notice that when restarting a simulation the history is read  by default
from a file and hills are added again to the grid.
This allows one to change the grid boundaries upon restart. However,
the first step after restart is usually very slow.
Since PLUMED 2.3 you can also store the grid on a file
and read it upon restart. This can be particularly
useful if you perform many restarts and if your hills are large.

For the precise syntax, see \ref METAD

\page MTS Multiple time stepping

By setting a STRIDE different from 1, you change how frequently
an action is calculated. In the case of actions such as \ref PRINT, this just
means how frequently you dump some quantity on the disk.
Notice that variables are only computed when necessary. Thus,
if a variable is only appearing as the argument of a \ref PRINT statement with
STRIDE=10, it will be computed every 10 steps.

In a similar fashion, the STRIDE keyword can be used in a bias potential
so as to apply the bias potential every few steps.
In this case, forces from this bias potential are scaled up by
a factor equal to STRIDE.

This technique can allow your simulation to run faster if you need
the apply a bias potential on some very expensive collective variable.
Consider the following input:
\plumedfile
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500
\endplumedfile
This performs a \ref METAD simulation biasing the distance between two
centers of mass. Since computing these centers requires a lot of atoms
to be imported from the MD engine, it could slow down significantly the
simulation. Notice that whereas the bias is changed every PACE=500 steps,
it is applied every STRIDE step, where STRIDE=1 by default.
The following input could lead to a significantly faster simulation at the price
of a negligible systematic error
\plumedfile
c1: COM ATOMS=1-1000
c2: COM ATOMS=1001-2000
d:  DISTANCE ATOMS=c1,c2
METAD ARG=d HEIGHT=1 SIGMA=0.1 BIASFACTOR=5 PACE=500 STRIDE=2
\endplumedfile
Similarly, the STRIDE keyword can be used with other biases (e.g. \ref RESTRAINT).

The technique is discussed in details here \cite Ferrarotti2015.
See also \subpage EFFECTIVE_ENERGY_DRIFT.

\page Multicolvar Multicolvar

Whenever you have a multicolvar action such as:

\plumedfile
COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1. D_MAX=3.0} MORE_THAN={RATIONAL R_0=6.0 NN=6 MM=12 D_0=0}
\endplumedfile

You will get a colossal speedup by specifying the D_MAX keyword in all switching functions that act on distances.
D_MAX tells PLUMED that the switching function is strictly zero if the distance is greater than this value.  As a result
PLUMED knows that it does not need to calculate these zero terms in what are essentially sums with a very large number of terms.
In fact when D_MAX is set PLUMED uses linked lists when calculating these coordination numbers, which is what 
gives you such a dramatic increase in performance.

\page Neighbour Neighbor Lists

Collective variables that can be speed up making us of neighbor lists:
- \ref COORDINATION
- \ref DHENERGY
- \ref PATHMSD

By tuning the cut-off for the neighbor list and the frequency for the recalculation of the list it is
possible to balance between accuracy and performances.

Notice that for \ref COORDINATION and \ref DHENERGY using a neighbor list could imply that a smaller
number of atoms are requested to the host MD engine. This is typically true when considering
\ref COORDINATION of a small number of atoms (e.g. a ligand) again many atoms (e.g. water).
When the neighbor list is used, only the water atoms close to the ligand will be requested at each step.

\warning
Notice that the calculation of the neighbor list is not not parallelized for \ref COORDINATION and \ref DHENERGY.
As a consequence, if you run
with many processors and/or OpenMP threads, the neighbor list might even make the calculation slower.


\page Openmp OpenMP

PLUMED is partly parallelized using OpenMP.
This should be enabled by default if your compiler supports it,
and can be disabled with `--disable-openmp`..
At runtime, you should set the environment variable
PLUMED_NUM_THREADS to the number of threads you wish to use with PLUMED.
The number of OpenMP threads can be set either by the MD code, if implemented in the patch, or generally by setting PLUMED_NUM_THREADS.
If they are not set openmp will be disabled at runtime. 

E.g., to run with gromacs you can do:
\verbatim
export PLUMED_NUM_THREADS=8
mdrun -plumed
\endverbatim

or as well

\verbatim
mdrun -plumed -ntomp 8
\endverbatim

In the first case the number of OpenMP threads used by plumed is 8 while the one used by gromacs can be 1 or something else, this is usually sub optimal.
In the second case GROMACS and plumed will use the same number of OpenMP threads.

Notice that:
- This option is likely to improve the performance, but could also slow down
  the code in some case.
- Results could be slightly different because of numerical round off and
  different order in summations. This should be harmless.
- The optimum number of threads is not necessary "all of them", nor should be
  equal to the number of threads used to parallelize MD.
- Only a few CVs are parallelized with openMP (currently, \ref COORDINATION and
  \ref DHENERGY).
- You might want to tune also the environmental variable PLUMED_CACHELINE_SIZE,
  by default 512, to set the size of cache lines on your machine. This is used
  by PLUMED to decrease the number of threads to be used in each loop so as to
  avoid clashes in memory access. This variable is expected to affect
  performance only, not results.


\page Secondary Secondary Structure

Secondary Structure collective variables (\ref ALPHARMSD, \ref PARABETARMSD and \ref ANTIBETARMSD)
can be particularly demanding if you want to calculate them for all the residues of a protein. 
This is particularly true for the calculation of beta structures.

The FIRST thing to speed up \ref PARABETARMSD and \ref ANTIBETARMSD is to use the keyword
STRANDS_CUTOFF (i.e. STRANDS_CUTOFF=1), in this way only a subset of possible fragments, the one
less than 1. nm apart, are used in the calculation.

The metric used to calculate the distance from ideal secondary structure elements can also influence 
the performances, try to use TYPE=OPTIMAL or TYPE=OPTIMAL-FAST instead of TYPE=DRMSD.

At last, try to reduce the number of residues in the calculation.

\page Lepton Making lepton library faster

In case you are using a lot of \ref CUSTOM functions or \ref switchingfunction "switching functions",
notice that these functionalities depend on the lepton library included in PLUMED.
This library replaces libmatheval since PLUMED 2.5, and by itself it is significantly faster than libmatheval.
However, you can make it even faster using a [just-in-time compiler](https://github.com/asmjit/asmjit.git).
As of PLUMED 2.6, the correct version of ASMJIT is embedded in PLUMED. In order to enable it
it is sufficient to use a specific flag in configure:
\verbatim
./configure --enable-asmjit
make
make install
\endverbatim

You are done!

Once ASMJIT has been configured, you can disable it at runtime setting the environment variable
`PLUMED_USE_ASMJIT`:
\verbatim
export PLUMED_USE_ASMJIT=no
\endverbatim


In some case using a custom expression is almost as fast as using a hard-coded
function. For instance, with an input like this one:
\plumedfile
...
c: COORDINATION GROUPA=1-108 GROUPB=1-108 R_0=1
d_fast: COORDINATION GROUPA=1-108 GROUPB=1-108 SWITCH={CUSTOM FUNC=1/(1+x2^3) R_0=1}
...
\endplumedfile
I (GB) obtained the following timings (on a Macbook laptop):
\verbatim
...
PLUMED: 4A  1 c                                          108     0.126592     0.001172     0.000701     0.002532
PLUMED: 4A  2 d_fast                                      108     0.135210     0.001252     0.000755     0.002623
...
\endverbatim

Notice the usage of `x2` as a variable for the switching function (see \ref switchingfunction), which
avoids an unnecessary square root calculation (this is done automatically by the hard-coded switching functions
when you use only even powers). The asmjit calculation (`d_fast`) takes less than 10% more than the hard-coded
one (`c`).

\page Time Time your Input

Once you have prepared your plumed input file you can run a test simulation, or use driver, 
to see which collective variable, function, bias or analysis is consuming more time and can 
thus be the target for a different definition (use less atoms, change relevant parameters,
or just use something else)

To have an accurate timing of your input you can use the \ref DEBUG DETAILED_TIMERS.
  
