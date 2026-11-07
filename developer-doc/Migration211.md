\page Migration211 Running inputs written for PLUMED 2.10 with PLUMED 2.11

Most inputs written for PLUMED 2.10 run with PLUMED 2.11 without changes. This page lists the cases in which an input
stops with an error, or runs but gives different results, and what to write instead.
The changes to the C++ interface used by the code compiled with LOAD are listed on the manual page of LOAD.

\section Migration211Renamed Keywords that were renamed

These keywords are still read, as the keyword that replaces them, with a warning in the log. Use the new name in new inputs.

| Action | Old keyword | New keyword |
|---|---|---|
| DUMPGRID, CONVERT_TO_FES | `GRID` | `ARG` |
| DFSCLUSTERING and the other actions that take a matrix | `MATRIX` | `ARG` |
| HISTOGRAM | `DATA` | `ARG` |
| RDF, CONTACT_MATRIX and the other adjacency matrices | `ATOMS` | `GROUP` |
| PAIRENTROPY | `SIGMA` | `BANDWIDTH` |
| DISTANCE_MATRIX | `CUTOFF` | `LINKCELL_CUTOFF` |
| TORSION | `VECTOR1`, `VECTOR2` | `VECTORA`, `VECTORB` |
| INSPHERE | `ATOM` | `CENTER` |
| AROUND | `ATOM` | `ORIGIN` |

\section Migration211Removed Keywords that were removed

- `CUTOFF` of RDF and PAIRENTROPY: the kernels are now truncated with the `D_MAX` of `KERNEL`, in units of `BANDWIDTH`.
  `CUTOFF=x` becomes `KERNEL={GAUSSIAN D_MAX=y}` with y=sqrt(2x); the error message gives the number.
- `NLIST_PACE_RESET` of OPES_METAD and OPES_METAD_EXPLORE has no effect: the neighbor list is now exact and on by default,
  and with WALKERS_MPI it is reset after each deposition.

\section Migration211Results Inputs that give different results

- The default kernel of KDE, HISTOGRAM, RDF and PAIRENTROPY is now a Gaussian that is stretched and scaled to go to zero at
  its cutoff, as in METAD and sum_hills. The kernel can be set with the switching function syntax, e.g. `KERNEL={GAUSSIAN D_MAX=3}`.
- The kernels of OPES_METAD and OPES_METAD_EXPLORE are force-shifted at `KERNEL_CUTOFF` (their derivative also goes to zero there).
- UPPER_WALLS and LOWER_WALLS acting on a periodic argument do not apply the periodic boundary conditions to the difference
  between the argument and `AT` any more.
- `ARG=*` leaves out the values that the action cannot use (e.g. in PRINT).
- JCOUPLING with `TYPE=CUSTOM` needs `A`, `B` and `C`, and `SHIFT` is zero by default.
- A negative number where a positive integer is expected (e.g. `STRIDE=-1`) is an error; it was read as a very large number.
- Many inputs that crashed or gave `nan` now stop with an error that names the keyword to change: see the file CHANGES/v2.11.md.

\section Migration211Older Inputs written for PLUMED 2.9 or older

These inputs stop with an error since PLUMED 2.10:

- HISTOGRAM with `LOGWEIGHTS` and without `NORMALIZATION`: set `NORMALIZATION=true` or `NORMALIZATION=false`
  (otherwise the weights were ignored), or use KDE and ACCUMULATE.
- The keywords of the old multicolvar syntax on actions that do not compute many values (e.g. `SUM` on HBOND_MATRIX):
  use the actions that compute the function of a vector or of a matrix (e.g. SUM, or MATRIX_VECTOR_PRODUCT for the rows of a matrix).
