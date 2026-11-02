Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

The Linux and macOS wheels include a basic PLUMED kernel and the ``plumed``
command. Linux wheels target the manylinux2014 baseline (glibc 2.17), and macOS
wheels target macOS 10.15 or newer. The bundled kernel is used by default. You
can select another kernel for a process with ``PLUMED_KERNEL`` or for one
object with ``Plumed(kernel=...)``.

To override the bundled kernel for all new PLUMED objects, set `PLUMED_KERNEL`::

     export PLUMED_KERNEL=/path/to/libplumedKernel.so
     python
     >>> import plumed
     >>> p=plumed.Plumed()

If you manage multiple PLUMED versions on your system with environment modules, setting ``PLUMED_KERNEL`` selects that external version instead of the bundled one.

Alternatively, a pure python solution is::

    >>> import plumed
    >>> os.environ["PLUMED_KERNEL"]="/path/to/libplumedKernel.so"
    >>> p=plumed.Plumed()

Finally, notice that you can set the path to the plumed library directly when declaring a Plumed object::

    >>> import plumed
    >>> p=plumed.Plumed(kernel="/path/to/libplumedKernel.so")

This allows you to mix different PLUMED versions in the same Python script.

Source builds bundle a kernel by default. Set ``PLUMED_BUILD_KERNEL=0`` when building the package to produce a wrapper-only distribution.

Maintainer release notes
------------------------

The ``python-package.yml`` workflow builds and validates the source distribution on pull requests and publishes it when a ``v*`` release tag is pushed. Before the first automated release, configure PyPI Trusted Publishing for the ``plumed`` project with GitHub owner ``plumed``, repository ``plumed2``, workflow ``python-package.yml``, and the GitHub environment ``pypi``. For a TestPyPI rehearsal, register the same workflow and repository with the ``testpypi`` environment on TestPyPI, then run the workflow manually with ``publish_to_testpypi`` enabled. No upload token is stored in GitHub Actions.

CHANGES: See the PLUMED documentation.
