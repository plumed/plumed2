Python wrappers for plumed
==========================

Install using the following command::

     python -m pip install plumed

The package supports CPython 3.11 through 3.15. The Linux and macOS wheels
include a basic PLUMED kernel and the ``plumed`` command. Linux wheels target
the manylinux2014 baseline (glibc 2.17). The macOS deployment target is 10.15
for Intel and at least 11.0 for Apple Silicon; the wheel tag records the actual
minimum, which also depends on the Python version. The bundled kernel is used
by default. You
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

InputBuilder metadata
---------------------

With kernel API version 11 or newer, ``InputBuilder`` reads metadata directly
from the selected kernel. It requires no executable, source checkout, or Vim
files. The kernel generates JSON from its live action registry, including
hidden keywords and the exact numbered-keyword classification. The same
serializer serves ``plumed gen_json``; its optional ``--actions`` argument
still accepts a legacy action-summary file.

Older kernels continue to use their installation's Vim dictionary and help
files. Keep those files installed when selecting an older kernel.

With ``InputBuilder(toplumed=p)``, metadata is read from that existing PLUMED
object, including actions it has already loaded. A successful ``ib.LOAD(...)``
refreshes the builder's methods and help automatically. If an action is loaded
directly with ``p.cmd(...)``, call ``ib.refresh()`` afterward. The supplied
object remains owned by the caller. Refresh is unavailable with older kernels
using static Vim files, or when ``load_dict=False``. For an attached older kernel,
select its Vim installation through ``kernel=...`` or ``PLUMED_ROOT`` as before.

Vim files
---------

``plumed gen_vim --outdir DIR`` generates both Vim syntax plugins and all
per-action help files in one invocation. It reads two static templates from
``<PLUMED_ROOT>/vim/templates``; use ``--template-dir DIR`` to select another
location. The templates are installed resources, while the action help comes
from the kernel's keyword registrations. An invalid template leaves existing
syntax and help files in place; successful regeneration removes stale help.

Maintainer release notes
------------------------

The ``python-package.yml`` workflow builds and validates the source distribution on pull requests and publishes it when a ``v*`` release tag is pushed. Before the first automated release, configure PyPI Trusted Publishing for the ``plumed`` project with GitHub owner ``plumed``, repository ``plumed2``, workflow ``python-package.yml``, and the GitHub environment ``pypi``. For a TestPyPI rehearsal, register the same workflow and repository with the ``testpypi`` environment on TestPyPI, then run the workflow manually with ``publish_to_testpypi`` enabled. No upload token is stored in GitHub Actions.

CHANGES: See the PLUMED documentation.
