# Reactive Voronoi assignment backends

VORONOI_COORDINATION, VORONOI_DISTANCE and VORONOI_POSITION accept explicit
atom selections and preserve their mathematical definitions across the
assignment backends. The optional CUDA library implements FP64 soft
assignments and their analytic reverse derivatives.

## Input and backend selection

The default is shared CPU evaluation. Compatible actions in one PLUMED
instance reuse geometry, normalized weights and occupancies. Different
KAPPA values share CPU geometry and normalize separately. Each action keeps
its own REFERENCE, SELECT, coefficients, sign and reduction semantics.

~~~plumed
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2
d: VORONOI_DISTANCE CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 GROUP1=1,2
PRINT ARG=q,d FILE=COLVAR
~~~

Add NO_SHARE to select independent CPU evaluation for a reference comparison.
BACKEND=CPU also selects shared CPU evaluation.

BACKEND=CUDA explicitly requests CUDA and reports an error when the plugin has not been
loaded or the device or requested calculation is unsupported.

~~~plumed
LOAD FILE=./libReactiveVoronoiCuda.so
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 BACKEND=CUDA
~~~

BACKEND=AUTO opts into capability detection. It selects CUDA when the plugin has
been loaded and the device and input are supported. Otherwise it selects
CPU and logs the reason.

~~~plumed
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 BACKEND=AUTO
~~~

Place LOAD before the Voronoi actions. Backend names are case insensitive.
The log identifies the requested backend, selected backend and GPU device
when available.

AUTO checks capability and input compatibility. Actual speed depends on system
size, reduction groups, host transfers and the MD workload. The CPU default
keeps existing inputs reproducible, and small groups can favor CPU execution.
An initialization fallback persists for the action's lifetime. Changes between
orthorhombic and triclinic boxes select a compatible backend on each request.
Calculation errors after successful CUDA initialization stop the calculation.
For a rejected center–assigned pair, the message includes the MD step and
the atom numbers from the input selections.

The currently supported CUDA domain is orthorhombic PBC or NOPBC and a
single rank, with either complete assignments or NLIST. AUTO falls back to
CPU for triclinic PBC and multiple MPI ranks.

NLIST uses the same CPU-built candidate list and cutoff in both backends.
Each CUDA context uploads the list after a rebuild and reuses the uploaded
indices between evaluations. The center-distance reduction for
VORONOI_DISTANCE runs on the host over the configured GROUP1 and GROUP2
pairs. NLIST renormalizes the assignment over the retained centers, so check
both CVs and forces against complete assignments when choosing NL_CUTOFF.
NL_STRIDE=1 rebuilds the candidates each step. For longer strides use a
positive NL_SKIN to enable displacement- and box-triggered rebuilds.

~~~plumed
UNITS LENGTH=A
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 NLIST NL_CUTOFF=6 NL_SKIN=0.5 NL_STRIDE=10 BACKEND=AUTO
~~~

The distances in this example are in angstroms. Choose the cutoff and buffer
for the actual atom selections and length units.

## Using CUDA Voronoi with other actions and force fields

The existing cudaCoord plugin and the CUDA Voronoi backend can run in the
same PLUMED input. Load both plugins and select the Voronoi backend on
the Voronoi action.

~~~plumed
LOAD FILE=./CudaCoordination.so
LOAD FILE=./libReactiveVoronoiCuda.so
c: CUDACOORDINATION GROUPA=1,2 GROUPB=3,4 R_0=0.32 D_MAX=1.25
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 BACKEND=CUDA
PRINT ARG=c,q FILE=COLVAR
~~~

The example uses PLUMED's default length unit (nm). The two plugins keep
their own GPU data and can use the same device. The MD engine can also use
that device for the potential, as in the single-GPU DeePMD examples. The CV's backend is chosen independently of
the potential. For a CPU force field, the default shared CPU CV requires no
CUDA library. BACKEND=CUDA can instead offload the Voronoi assignment if
the run and GPU meet the requirements above. Coordinates and bias
derivatives pass through the regular host PLUMED interface.

## Building the optional library

The PLUMED core builds with a normal C++ compiler and keeps CUDA headers
confined to the plugin. Build the CUDA library on a Linux CUDA host.

~~~bash
make PLUMED=/path/to/plumed CUDA_ARCH=sm_70
~~~

Select the architecture of the target GPU, such as sm_70 for V100 or sm_80
for A100. The default CUDA_ARCH=native detects the GPU on the compilation
host. Set NVCC and NVCCFLAGS when the toolchain requires them. Install the
resulting library in a versioned directory and pass its path to LOAD. Use
the same PLUMED build at compile time and at runtime. The Makefile uses the
shared-library extension reported by PLUMED.

The build links the shared CUDA runtime and marks the Linux library NODELETE.
This follows the lifetime protection used by the official CUDA plugin and
preserves an MD engine's live device allocations when PLUMED is destroyed.

## Device ownership and data flow

Device selection follows the MD-provided setGpuDeviceId when supplied.
Otherwise it captures the current CUDA device on first evaluation. Every
allocation, transfer, kernel, synchronization and release uses that captured
device, and restores the caller's device afterward. CUDA_VISIBLE_DEVICES
provides the usual logical-device mapping.

Each compatible KAPPA context retains O(Ncenters+Nassigned) device buffers
for complete assignments. NLIST adds offsets by assigned atom and one center
index per buffered candidate. FP64 kernels stream over the selected centers
for each assigned atom and use shifted-softmax normalization. CUDA FP64
atomic additions can sum in a different order between runs, so the final
digits of values and derivatives can vary. Each scalar
action still contracts its own analytic derivatives. PLUMED combines bias
sensitivities using its existing force graph.

OPES, scalar reductions, direct distance/position terms and their direct
virial contributions execute on CPU. Full-pool VORONOI_DISTANCE therefore
includes a CPU center-pair reduction. The current host MD interface supplies
host coordinates and receives host forces.

The official cudaCoord plugin registers CUDACOORDINATION and
CUDACOORDINATIONFLOAT after LOAD. The openaccPTM plugin adds USEGPU selection
to its supported actions after LOAD. The Voronoi plugin registers its assignment
backend after LOAD. Its BACKEND keyword belongs to the three Voronoi actions.

## Validation

The core regtests rt-reactive-voronoi-shared and
rt-reactive-voronoi-backends run with CPU-only PLUMED. The latter covers
coordination powers/signs/selections, cross/within-group distances, position
axes/signs/normalization/absolute values, reordered atoms, mixed KAPPA,
PBC/NOPBC, strided lists, repeated-step changes and finite-difference forces
and cell derivatives. It includes a triclinic position-image regression.

Run the driver matrix from the plugin directory with the selected PLUMED
installation. Run the CUDA cases on a CUDA-capable device.

~~~bash
make CUDA_ARCH=sm_70 check PLUMED=/absolute/path/bin/plumed
~~~

The driver matrix includes a fixed-box trajectory that forces a rebuild
before the scheduled step. It compares values and derivatives to an
independent CPU list rebuilt every step. It also checks invalid stride/skin
inputs and pinpoints failed comparisons by file, row and column.
Run the CPU and missing-plugin fallback cases on a CPU-only host with:

~~~bash
python3 regtest/backend_checks.py --plumed /absolute/path/bin/plumed --cpu-only
~~~

With a CUDA toolkit, `make CUDA_ARCH=sm_80 check-nogpu` also tests loading
the real plugin with no visible GPU. It checks the numerical CPU fallback
for AUTO and the initialization error for an explicit CUDA request.

In a configured PLUMED source build, the native CPU regression is available
through PLUMED's test harness.

~~~bash
make -C regtest/basic/rt-reactive-voronoi-backends
~~~

The native test executable also accepts a backend suffix for a manual CUDA
or AUTO run. Use an empty working directory for each run because the test
writes logs and an output file in the current directory.

~~~bash
PLUMED_VORONOI_TEST_PLUGIN=/absolute/path/libReactiveVoronoiCuda.so \
  /absolute/path/to/regtest/basic/rt-reactive-voronoi-backends/tmp/exe " BACKEND=CUDA"
~~~

Strict CUDA tests compare every supported force/virial component to CPU and
finite-difference diagonal strains for PBC. AUTO also checks sheared boxes
through its CPU fallback. All cell strains are checked for NOPBC.

Application benchmarks and large trajectories belong in separate case records.
