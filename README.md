
Content
-------

    CHANGES/         : change log
    COPYING.LESSER   : license
    Makefile         : makefile
    Makefile.conf.in : template configuration makefile
    PEOPLE           : list of authors
    README           : this file
    VERSION          : version file
    configurations/  : template configuration files
    configure        : configuration script
    configure.ac     : configuration script (autoconf)
    developer-doc    : developer documentation
    include          : symbolic link for include files
    macports         : directory where Portfiles are generated
    patches          : patch scripts
    release.sh       : developer utility to publish releases
    regtest          : regression tests, including reference results
    scripts          : shell tools
    src              : source code
    sourceme.sh      : template configuration script
    test             : examples
    user-doc         : user documentation
    vim              : directory where vim syntax is generated

Install
-------

Extensive installation instructions are in the user documentation.
Links to precompiled versions of the documentation can be found [here](http://www.plumed.org/documentation).

Required software:

* GNU make
* C/c++ compiler (c++11 support is required as of version 2.4).
* A modern version of the `patch` command line tool
* Support for POSIX library `dirent.h`
* `xxd` (present in most unix distributions)

Suggested software (libraries are checked by `./configure` and enabled if available):

* Optimized blas and lapack libraries. Can be replaced by an internal version if not available.
* MPI library to run parallel simulations. It should be the same library used by your MD code.
* [VMD molfile plugins](http://www.ks.uiuc.edu/Research/vmd/plugins) to read arbitrary file formats. Can be replaced by an internal version supporting a few formats if not available.
* [Matheval library](http://www.gnu.org/software/libmatheval) to use algebraic collective variables.
* [Zlib library](http://zlib.net/) to use compressed data files.
* [Xdrfile library](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library) to have read/write access to gromacs
  trajectory files.
* [Doxygen](http:://www.doxygen.org) to build user manual. Doxygen might need the following packages:
  * Latex to build the pdf user manual.
  * [Graphviz](http://www.graphviz.org) to show class hierarchy in
    developer manual.

Quick compilation instructions
------------------------------

Configure for your system

    ./configure --prefix=$HOME/opt
    
    
If necessary, edit `Makefile.conf`. 
Configure your environment

    source ./sourceme.sh
    
Compile plumed

    make
    
The `plumed` executable should be now in your execution path

    plumed help
    
Compile the manuals (pre-compiled manual is available online):

    make doc

User documentation can be found at `user-doc/html/index.html`.
Developer documentation can be found at `developer-doc/html/index.html`.

Install PLUMED in `$HOME/opt` (directory should be set during `./configure`):

    umask 022
    make install
    
A sample modulefile with environment variable will be placed in
`$HOME/opt/lib/plumed/src/lib/modulefile`. Path to the installed documentation can be found with command `plumed info --user-doc`.


Branches and releases
---------------------

Several branches and tags are stored on the git repository.

Branches named `v2.X` correspond to release branches.

Master branch may contain non tested features and is not expected to be used by non-developers.
It typically contains features that will be available on the next release.

Tags named `v2.XbY` correspond to beta releases, use it with care.
Tags named `v2.X.Y` correspond to official releases, use the latest available.

To report problems found on beta or official releases, use the normal `plumed-users@googlegroups.com`
mailing list. Just state exactly which version you are using.

To report problems found on `master` branch, use the `plumed2-git@googlegroups.com` mailing list.
This is also the correct place for discussions about new features etc.
When reporting please provide the git hash (you can obtain it with `git rev-parse HEAD`).

Status
------

| Branch   |      Status   | Supported |
|:--------:|:-------------:|:--------:|
| master   | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=master)](https://travis-ci.org/plumed/plumed2) | yes |
| v2.3     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2) | yes |
| v2.2     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.3)](https://travis-ci.org/plumed/plumed2) | no |
| v2.1     | [![Build Status](https://travis-ci.org/plumed/plumed2.svg?branch=v2.2)](https://travis-ci.org/plumed/plumed2) | no |
| v2.0     | Not available | no |


