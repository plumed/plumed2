<!--
description: Tensor train metadynamics (TT-MetaD)
authors: Nils E. Strand, Siyao Yang, Yuehaw Khoo, Aaron R. Dinner
reference: https://doi.org/10.1021/acs.jctc.6c00487
-->

# Overview

This module implements tensor train metadynamics (TT-MetaD), a method that
represents the metadynamics bias potential as a tensor train (TT) rather than
as a sum of Gaussian hills. The accumulated Gaussian hills are periodically
compressed into a low-rank TT approximation using the TT-Sketch algorithm,
keeping memory and evaluation cost bounded throughout long simulations.

For the full algorithmic details and for guidance on the choice of the parameters see the paper cited below
and the [TTMETAD](TTMETAD.md) page.


# Installation

This module requires the [ITensor] library (v3) with HDF5 support.

## Building ITensor

Follow the [ITensor installation instructions]. Copy `options.mk` to the
ITensor source directory and edit it to match your compiler and LAPACK/BLAS
setup. Then build:

```bash
cd /path/to/ITensor
cp options.mk.sample options.mk
# edit options.mk as needed
make
```

## Building PLUMED with ttsketch

Once ITensor is built, configure PLUMED from its source directory. Adjust
library paths and flags to match your `options.mk`:

```bash
./configure --enable-modules=ttsketch --enable-libitensor \
    LDFLAGS="-L/path/to/ITensor/lib" \
    CPPFLAGS="-I/path/to/ITensor -DITENSOR_USE_HDF5 -DITENSOR_USE_OMP -D__ITENSOR_LAPACK_WRAP_h" \
    LIBS="-litensor -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lhdf5 -lhdf5_hl" \
    CXXFLAGS="-O3 -fconcepts"
```

The flags above assume Intel MKL for LAPACK/BLAS. If you use OpenBLAS or
another provider, replace the `-lmkl_*` flags accordingly.


[ITensor]: https://github.com/ITensor/ITensor
[ITensor installation instructions]: https://itensor.org/docs.cgi?page=install
