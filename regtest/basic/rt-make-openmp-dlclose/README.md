This Linux regression links only the runtime wrapper, so the test executable
cannot keep the kernel or OpenMP runtime loaded through a link dependency.

First it creates and finalizes a kernel without starting OpenMP workers. It
checks that the kernel unloads and its dynamically linked OpenMP runtime stays
mapped. This gives a deterministic failure with the runtime-retention fix
removed, without depending on the timing of a worker-thread segmentation fault.
Then it runs 100 loading/calculation/unloading cycles with a four-thread
COORDINATION action. Other platforms simply produce the reference output.

Run in a normal environment without preloading OpenMP or retaining the kernel.
An OpenMP-disabled build still tests repeated kernel unloading.

The kernel retains a shared OpenMP runtime by default on Linux. To manage the
runtime lifetime explicitly, set PLUMED_OPENMP_RETAIN_RUNTIME=no before creating
any PLUMED instance. The application must then ensure all runtime resources are
released before unloading it. Changing the variable after retention cannot undo
RTLD_NODELETE. The regression clears this variable to test the default behavior.
It also starts a fresh subprocess with the opt-out enabled and checks that the
runtime unloads. That subprocess does not start workers, so no resource cleanup
is required before unloading.

Some ELF loader/compiler configurations retain the kernel even with runtime
retention explicitly disabled (for example through GNU-unique symbols). The
isolated opt-out child detects this baseline restriction with status 77. In
that case the parent still exercises the repeated calculations, but does not
require the loader to unmap the kernel or its dependencies. When the baseline
kernel can unload, all kernel and runtime mapping assertions remain enforced.
