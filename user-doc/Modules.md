\page mymodules List of modules

The functionality in PLUMED 2 is divided into a small number of modules.  Some
users may only wish to use a subset of the functionality available within the 
code while others may wish to use some of the more complicated features that are available.
For this reason the plumed source code is divided into modules, which users can
activate or deactivate to their hearts content.  

You can activate a module at configure time using the keyword `--enable-modules`.
For example:
\verbatim
./configure --enable-modules=modulename
\endverbatim
will enable module called modulename. A module that is on by default can be disabled
using the following syntax
\verbatim
./configure --enable-modules=-modulename
\endverbatim
To enable or disable multiple modules one should provide them as a : separated
list. Notice that `+modulename` and `modulename` both activate the module, whereas
`-modulename` deactivates it. E.g.
\verbatim
./configure --enable-modules=+crystallization:-colvar
\endverbatim
will disable the colvar module and enable the crystallization module.
Also notice that `:` can be omitted when using `+` or `-`. Thus, the same can be obtained
with
\verbatim
./configure --enable-modules=+crystallization-colvar
\endverbatim

If you repeat the `--enable-modules` keyword only the last instance will be used. Thus
`./configure --enable-modules=crystallization --enable-modules=-colvar` will _not_ do what you expect!

There are also some shortcuts available:
- `./configure --enable-modules=all` to enable all optional modules. This includes the maximum number of features in PLUMED,
including modules that might not be properly functional.
- `./configure --enable-modules=none` or `./configure --disable-modules` to disable all optional modules. This produces a minimal
PLUMED which can be used as a library but which has no command line tools and no collective variables or biasing methods.
- `./configure --enable-modules=reset` or `./configure --enable-modules` to enable the default modules.

The two kinds of syntax can be combined and, for example, `./configure --enable-modules=none:colvar` will result
in a PLUMED with all the modules disabled with the exception of the colvar module.

Some modules are active by default in the version of PLUMED 2 that you download from 
the website while others are inactive.  The following lists all of the modules that
are available in plumed and tells you whether or not they are active by default.

@MODULES@

Until PLUMED 2.2, it was also possible to switch on or off modules by adding files
in the plumed2/src directory. Since PLUMED 2.3 this is discouraged, since any choice made
in this manner will be overwritten next time `./configure` is used.

