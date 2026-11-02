Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

WARNING: You will need to also build and install the plumed library (see http://www.plumed.org) and make sure the file `libplumedKernel.so` (or `libplumedKernel.dylib`) is available on your system.

The PyPI package contains the Python interface only. It does not include or install the PLUMED runtime kernel. Install a compatible PLUMED kernel separately and set ``PLUMED_KERNEL`` to its library path as shown below.

You should then make sure the library is found setting the environment variable `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple plumed versions on your system using tcl environment modules, this should be taken care automatically
by the plumed module.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This will allow you to mix different plumed versions in the same python script.

Maintainer release notes
------------------------

The ``python-package.yml`` workflow builds and validates the wheel and source distribution on pull requests and publishes both when a ``v*`` release tag is pushed. Before the first automated release, configure PyPI Trusted Publishing for the ``plumed`` project with GitHub owner ``plumed``, repository ``plumed2``, workflow ``python-package.yml``, and the GitHub environment ``pypi``. No PyPI API token is required by the workflow.

CHANGES: See the PLUMED documentation.
