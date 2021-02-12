\page InstallationLayout Installation Layout

I write here some notes related to how plumed is installed.

As a first notice, plumed package is mostly designed to make these tools available:
- a `plumed` executable, that can be used to launch command line tools.
- a `libplumed.so` library, that can be linked to an externally supplied MD code.

These are the main entry points to plumed, but they require several other resources to be
properly located so as to work. Moreover, plumed is designed to be usable both
when "just compiled" (so as to allow for fast development) and when properly installed.
This results in the non trivial problem of knowing where the required resources are located.

Additionally, we provide shell-only alternatives to command line tools. E.g.,
to use `plumed patch` in a cross compiled environment, one can call `plumed-patch` which
does the same but bypass the `plumed` executable and directly executes as a bash script.

As a result, plumed routines and scripts can be entered in three ways:
- calling `plumed` from the command line.
- calling `plumed-patch` from the command line.
- entering the shared library from another code (typically an MD code).

This is achieved in the following way:
- `plumed-*` scripts contains hardcoded environment variables pointing at the correct paths when they are launched.
- `plumed` executable and `libplumed.so` have access to methods (in the config namespace)
  to access to the paths and locate resources.

\warning
Since paths are hardcoded, plumed executable and library are de facto not relocatable.
It is however possible to use them after relocation provided that some environment
variables are set as discussed below.

As an example, the `PLUMED_ROOT` variable is defined to tell to the plumed scripts where to find
most of the plumed-related files. Similarly, from C++ you can use \ref config::getPlumedRoot() to retrieve
the same path.

When a plumed command line tool implemented as script is invoked by the plumed executable,
thus transferring the control from C++ to
an external script, the environment should be consistently set. This is done in method \ref config::getEnvCommand()
which builds a string in the form `env PLUMED_ROOT=/path env PLUMED_INCLUDEDIR=/path ` etc.
In this ways, the scripts are run in an environment with the correct settings.

The following paths need to be set for plumed to work properly. Here they are listed
together with the name of the corresponding environment variables (that can be used in plumed
scripts) and the method that can be used to retrieve them from C++ code.
- Root of plumed: `$PLUMED_ROOT` \ref config::getPlumedRoot()
- Path to include files: `$PLUMED_INCLUDEDIR` \ref config::getPlumedIncludedir()
- Path to html files: `$PLUMED_HTMLDIR` \ref config::getPlumedHtmldir()
- Name of plumed program: `$PLUMED_PROGRAM_NAME` \ref config::getPlumedProgramName()

When using plumed from its build directory (without installing it) these paths will be set to the
value reported below:
- `PLUMED_ROOT=/build/directory`
- `PLUMED_INCLUDEDIR=$PLUMED_ROOT/src/include` (this works thanks to a symlink of `/build/directory/src` to `/build/directory/src/include/plumed`)
- `PLUMED_HTMLDIR=$PLUMED_ROOT`
- `PLUMED_PROGRAM_NAME=plumed`

These paths are hardcoded in `plumed` executable and in `plumed-*` scripts when they are compiled.
Notice that it is possible to set the `PLUMED_ROOT` variable before calling plumed overriding the hard code values.
E.g., you can compile plumed in directory `/build/directory1`, move it to `/build/directory2`, and launch it
with `PLUMED_ROOT=/build/directory2 /build/directory2/src/lib/plumed`. Notice however that although plumed will find all the
required resources in this way, it might not be possible to perform some task such as patching MD code. Also notice that
since the structure of the build directory is fixed the `PLUMED_ROOT` variable is sufficient to reconstruct the other paths.

When using plumed after it has been installed, these paths will be set to the value reported below:
- `PLUMED_ROOT=/usr/local/lib/plumed`
- `PLUMED_INCLUDEDIR=/usr/local/include`
- `PLUMED_HTMLDIR=/usr/local/share/doc/plumed`
- `PLUMED_PROGRAM_NAME=plumed`

These paths are hardcoded in `plumed` executable and in `plumed-*` scripts when they are installed.
Notice that these value can be customized at configure step using standard arguments to ./configure.
When using an installed copy of plumed one can override the hard code values by setting the variables
`PLUMED_ROOT`, `PLUMED_INCLUDEDIR` ,`PLUMED_HTMLDIR`, and `PLUMED_PROGRAM_NAME` before launching plumed.

Notice that to enforce a consistent behavior of scripts and plumed executable the same logic needed to
be implemented twice. One implementation is found in the `src/config/Config.inc.in` file, another implementation
is prependend to the installed scripts by the `src/lib/Makefile`.

Also consider that environment is inherited by subprocesses. That means that if you want to
launch another plumed version from a plumed script (crazy idea, perhaps nobody will ever do it)
you should unexport the relevant environment variables so that the second plumed executable
will find its paths correctly.

\section InstallationLayout-files Installed files

I here describe what's the content of the most important files installed by plumed.

`/usr/local/bin/plumed`: this is a static executable that can be used to launch plumed.
It is typically used to launch a command line tool (e.g. `plumed sum_hills`).
Notice that some command line tools are actually implemented as bash scripts (e.g. `plumed patch`).
Those scripts are located in `$PLUMED_ROOT/scripts/` with an extra `.sh`  suffix. E.g.
the `plumed patch` command will set properly the environment then call `$PLUMED_ROOT/scripts/patch.sh`.

`/usr/local/lib/libplumed.so`: this is a library containing all the plumed routines.
Notice that `/usr/local/bin/plumed` described above is equal to the combination of
`/usr/local/lib/libplumed.so` with a single object file compiled from `buildroot/src/main/main.cpp`.

`/usr/local/lib/libplumedKernel.so`: this is a library containing almost all the plumed routines,
with the exception of those called from MD engines.
Notice that `/usr/local/lib/libplumed.so` described above is equal to the combination of
`/usr/local/lib/libplumedKernel.so` with a single object file compiled from `buildroot/src/wrapper/PlumedStatic.cpp`

`/usr/local/lib/libplumedWrapper.a`: this is a static library containing exclusively the
object file compiled from `buildroot/src/wrapper/Plumed.cpp`

To summarize:
- `bin/plumed` = `buildroot/src/main/main.cpp` + `lib/libplumed.so`
- `lib/libplumed.so` = `buildroot/src/wrapper/PlumedStatic.cpp` + `lib/libplumedKernel.so`
- `lib/libplumedWrapper.a` = `buildroot/src/wrapper/Plumed.cpp`

The logic of this subdivision is that it is possible to either link the MD code to `/usr/local/lib/libplumed.so`
or to link it to a single object file (the one compiled from `buildroot/src/wrapper/Plumed.c` or the installed `libplumedWrapper.a`)
so as to allow linking at run time an a posteriori chosen plumed library. This is the trick behind the `--runtime` patching procedure.

Notice that the only differences between `buildroot/src/wrapper/PlumedStatic.cpp` and `buildroot/src/wrapper/Plumed.c` are that
runtime binding is disabled for the first one and that the second one is compiled as plain C.
This makes it less likely to do mistakes when linking lib/libplumed.so (by unintentionally using a different version
of plumed), and makes C++ library unnecessary if an external code is only interesting in linking the PLUMED
wrappers in `buildroot/src/wrapper/Plumed.c` or in `libplumedWrapper.a`.

We can then dissect more the content of `/usr/local/lib/libplumedKernel.so`.
This library puts together a large list of object files. The same object files will be located after install
in `/usr/local/lib/plumed/obj/k*.o`. I use a wildcard here because these might be many files (named `k0.o`, 'k1.o', etc) or
a single `kernel.o` file (when `ld -r -o` can be used to merge them together). The reason why we 
store object files in the installed directory is that this is the most portable way to link statically C++
objects to another executable. Indeed, merging them in a single .a file (such as libplumed.a) 
would require this library to be linked with special flags so as to allow dropping all the static constructors.
Whereas the special flags could be found by autoconf, it seems simpler to directly link `/usr/local/lib/plumed/obj/k*.o`.


Also notice that this library changes slighlty in the installed version (`/usr/local/lib/libplumedKernel.so`)
and in the pre-install version (`buildroot/src/lib/libplumedKernel.so`). Indeed, whereas the former
include the object file from `buildroot/src/config/ConfigInstall.cpp` the latter includes the object file from
`buildroot/src/config/Config.cpp`. This object file is the one containing the hardcoded paths discussed above,
and thus should include different strings in the installed and pre-install versions.

\note
New in PLUMED v2.5, the `./configure` script will check if it is possible to build a `/usr/local/lib/libplumed.a` library.
This library contains basically `buildroot/src/wrapper/PlumedStatic.cpp` and the single object obtained
merging all the objects in the kernel. When this library is linked, if at least one of the functions in the wrappers
is called (e.g. `plumed_cmd`) then all the objects are pulled in. In principle, this should solve the problem
with C++ static constructors. This feature can be disabled with `--disable-static-archive`.

\section InstallationLayout-installation Installation procedure

When `make` is invoked, several things are performed. First, all the source files are compiled.
The `plumed` executable and the library files are put in `buildroot/src/lib`.
Then, the "to be installed" versions
of the executable and library files are produced and located in `buildroot/src/lib/install`. These are different from
those located in `buildroot/src/lib` in that they include the `buildroot/src/config/ConfigInstall.o` object so as
to hardcode the proper paths.

When `make install` is invoked, the makefile checks if the objects in `buildroot/src/lib/install` should be updated
and, if necessary, recompiles them. If not, it just copies all the material in place. Notice that all the resulting
files are real files (no symlinks). This is a novelty with respect to PLUMED 2.1 and allows for a proper implementation
of the DESTDIR feature required by unix distributions.

Using the standard behavior explained in the autoconf documentation, it is possible to change the paths
for plumed install either during configure (with `--prefix`) or by setting `prefix` during `make install`.


