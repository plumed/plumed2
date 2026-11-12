# /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   Copyright (c) 2011-2016 The plumed team
#   (see the PEOPLE file at the root of the distribution for a list of names)
#
#   See http://www.plumed.org for more information.
#
#   This file is part of plumed, version 2.
#
#   plumed is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   plumed is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
# This python routine builds an interface between plumed and python using cython
#
from setuptools import setup
from setuptools import Extension
from setuptools import find_namespace_packages
from setuptools.command.build_py import build_py
from pathlib import Path
import shutil
import os
import os.path
import sys
import platform
import sysconfig
import re

if sys.version_info < (3,):
    raise ImportError("PLUMED 2.6 only supports Python 3")

is_platform_mac = sys.platform == "darwin"

print("""You can manipulate the behaviour of the compilation
by exporting the following variables [default values for pypi]:

 - plumed_program_name          ["plumed"]
 - plumed_version               [-greps from ../Version.txt-]
 - plumed_default_kernel        [None] - set the full path
 - plumed_force_cython          [None] - active if set to "yes"
 - plumed_disable_rtld_deepbind [None] - active if exported
""")
# parameters that can be changed with command line
plumedname = os.getenv("plumed_program_name")
if plumedname is None:
    plumedname = "plumed"

plumedversion = os.getenv("plumed_version")
if plumedversion is None:
    # TODO: if PLUMED_VERSION is not present copy the file from the
    # base directory (or use Unknown ?)
    # as grep -v "#" PLUMED_VERSION, which is not available on Windows
    with open("./PLUMED_VERSION") as f:
        plumedversion = "".join(line for line in f if "#" not in line).rstrip()

print("Module name " + plumedname)
print("Version number " + plumedversion)

extra_compile_args = [
    "-D__PLUMED_HAS_DLOPEN",
    "-D__PLUMED_WRAPPER_LINK_RUNTIME=1",
    "-D__PLUMED_WRAPPER_IMPLEMENTATION=1",
    "-D__PLUMED_WRAPPER_EXTERN=0",
]

defaultkernel = os.getenv("plumed_default_kernel")
if defaultkernel is not None:
    defaultkernel = os.path.abspath(defaultkernel)
    if sys.platform == "win32":
        # Plumed.h turns the macro into a string with the # operator, which
        # keeps backslashes as they are: "C:\Users\..." would then contain
        # escape sequences. Windows accepts forward slashes as well.
        defaultkernel = defaultkernel.replace("\\", "/")
    extra_compile_args.append("-D__PLUMED_DEFAULT_KERNEL=" + defaultkernel)
    print("Hardcoded PLUMED_KERNEL " + defaultkernel)

plumed_disable_rtld_deepbind = os.getenv("plumed_disable_rtld_deepbind")
if plumed_disable_rtld_deepbind is not None:
    extra_compile_args.append("-D__PLUMED_WRAPPER_ENABLE_RTLD_DEEPBIND=0")
    print("Disabling RTLD_DEEPBIND")

# Fixes problem with compiling the PYTHON interface in Mac OS 10.14 and higher.
# Sets the deployment target to 10.9 when compiling on version 10.9 and above,
# overriding distutils behaviour to target the Mac OS version python was built for.
# This can be overridden by setting MACOSX_DEPLOYMENT_TARGET before compiling the
# python interface.
# This fix is taken from https://github.com/pandas-dev/pandas/pull/24274/files
if is_platform_mac:
    if "MACOSX_DEPLOYMENT_TARGET" not in os.environ:

        def version_tuple(value):
            return tuple(int(part) for part in re.findall(r"\d+", value))

        current_system = version_tuple(platform.mac_ver()[0])
        python_target = version_tuple(
            sysconfig.get_config_var("MACOSX_DEPLOYMENT_TARGET") or "0"
        )
        if python_target < (10, 9) and current_system >= (10, 9):
            os.environ["MACOSX_DEPLOYMENT_TARGET"] = "10.9"


def include_dir():
    """Return the directory containing Plumed.h.

    include/Plumed.h is a symbolic link to ../../src/wrapper/Plumed.h. Git
    without symbolic link support (the default on Windows) checks it out as a
    small text file containing the path of the target: follow it by hand.
    """
    header = Path("include") / "Plumed.h"
    if header.is_file() and header.stat().st_size < 1024:
        target = Path("include") / header.read_text().strip()
        if target.is_file():
            return str(target.parent)
    return "./include"


def readme():
    with open("README.rst") as f:
        return f.read()


# allow one to force using cython with env var plumed_force_cython=yes
USE_CYTHON = False
if os.environ.get("plumed_force_cython") == "yes":
    print("plumed_force_cython=yes")
    USE_CYTHON = True

# if plumed.c is available, do not need cython
if not USE_CYTHON:
    if not os.path.isfile("plumed.c"):
        print("plumed.c not found, cython is needed")
        USE_CYTHON = True

bundle_default = "0" if os.environ.get("plumed_default_kernel") else "1"
bundle_kernel = os.environ.get("PLUMED_BUILD_KERNEL", bundle_default).lower() not in {
    "0",
    "no",
    "false",
}

# try to import cython
if not USE_CYTHON:
    print("using available plumed.c file")

ext_modules = [
    Extension(
        name=plumedname + "._plumed",
        # without cython, compile the pregenerated plumed.c explicitly
        # (so that it is not regenerated even if cython is installed)
        sources=["plumed.pyx"] if USE_CYTHON else ["plumed.c"],
        include_dirs=[include_dir()],
        extra_compile_args=extra_compile_args,
    )
]


class BuildPy(build_py):
    """Replace staged runtime files, including resources from older layouts."""

    def run(self):
        # setuptools otherwise keeps removed package data across builds. Clear
        # only its output copy, also when switching to a wrapper-only wheel.
        shutil.rmtree(Path(self.build_lib) / plumedname / "runtime", ignore_errors=True)
        super().run()


setup(
    name=plumedname,
    version=plumedversion,
    description="Python interface to PLUMED",
    long_description=readme(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Programming Language :: Python :: 3.15",
    ],
    license="LGPL-3.0-or-later",
    author="Gareth A. Tribello",
    author_email="plumed-users@googlegroups.com",
    url="http://www.plumed.org",
    ext_modules=ext_modules,
    cmdclass={"build_py": BuildPy},
    # the directories of the private runtime (bin, lib, include/...) are declared as
    # namespace packages, otherwise setuptools warns for each of them;
    # their files are selected by MANIFEST.in (include_package_data)
    packages=[
        plumedname + name[len("plumed") :]
        for name in find_namespace_packages(include=["plumed", "plumed.*"])
    ],
    package_dir={plumedname: "plumed"},
    include_package_data=True,
    entry_points={"console_scripts": ["plumed=" + plumedname + ".runtime.cli:main"]}
    if bundle_kernel
    else {},
    license_files=["LICENSE"],
    setup_requires=["Cython"] if USE_CYTHON else None,
    extras_require={"hdf5": ["pandas", "h5py"]},
    zip_safe=False,
    # this will be deprecated shortly
    # test_suite='nose.collector',
    python_requires=">=3.11",
)
