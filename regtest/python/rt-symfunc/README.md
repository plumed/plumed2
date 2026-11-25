This regression compares the Python interface with the existing symmetry-function
reference. Keep the dynamically loaded kernel resident until process exit because
OpenMP worker threads can still be waiting inside its runtime at dlclose. The
PLUMED_LOAD_DLCLOSE setting is scoped to this regression and does not change its
thread count, calculations, numerical tolerance or reference output.
