#/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#
# This python routine builds an interface between plumed and python using cython
#
from setuptools import setup
from setuptools import Extension
import subprocess
import os
import os.path
import sys
from shutil import copyfile
import platform
from distutils.sysconfig import get_config_var
from distutils.version import LooseVersion

if sys.version_info < (3,):
    raise ImportError("PLUMED 2.6 only supports Python 3")

def is_platform_mac():
    return sys.platform == 'darwin'

if os.getenv("plumed_macports") is not None:
    copyfile("../VERSION.txt","PLUMED_VERSION")
    try:
        os.mkdir("include")
    except OSError:
        pass
    copyfile("../src/wrapper/Plumed.h","include/Plumed.h")

plumedname = os.getenv("plumed_program_name")
if plumedname is None:
    plumedname = "plumed"

plumedversion = os.getenv("plumed_version")
if plumedversion is None:
    plumedversion = subprocess.check_output(['grep','-v','#','./PLUMED_VERSION']).decode("utf-8").rstrip()

print( "Module name " + plumedname )
print( "Version number " + plumedversion )

extra_compile_args=['-D__PLUMED_HAS_DLOPEN','-D__PLUMED_WRAPPER_LINK_RUNTIME=1','-D__PLUMED_WRAPPER_IMPLEMENTATION=1','-D__PLUMED_WRAPPER_EXTERN=0']

defaultkernel=os.getenv("plumed_default_kernel")
if defaultkernel is not None:
    extra_compile_args.append("-D__PLUMED_DEFAULT_KERNEL=" + os.path.abspath(defaultkernel))
    print( "Hardcoded PLUMED_KERNEL " + os.path.abspath(defaultkernel))

# Fixes problem with compiling the PYTHON interface in Mac OS 10.14 and higher.
# Sets the deployment target to 10.9 when compiling on version 10.9 and above,
# overriding distutils behaviour to target the Mac OS version python was built for.
# This can be overridden by setting MACOSX_DEPLOYMENT_TARGET before compiling the
# python interface.
# This fix is taken from https://github.com/pandas-dev/pandas/pull/24274/files
if is_platform_mac():
    if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ:
        current_system = LooseVersion(platform.mac_ver()[0])
        python_target = LooseVersion(get_config_var('MACOSX_DEPLOYMENT_TARGET'))
        if python_target < '10.9' and current_system >= '10.9':
            os.environ['MACOSX_DEPLOYMENT_TARGET'] = '10.9'

def readme():
    with open('README.rst') as f:
        return f.read()


try:
    include_dirs=[os.environ["plumed_include_dir"]]
except KeyError:
    include_dirs=["./include"]

# allow one to force using cython with env var plumed_force_cython=yes
USE_CYTHON = False
try:
    if(os.environ["plumed_force_cython"]=="yes"):
        print('plumed_force_cython=yes')
        USE_CYTHON = True
except KeyError:
    pass

# if plumed.c is available, do not need cython
if not USE_CYTHON:
    if not os.path.isfile("plumed.c"):
        print('plumed.c not found, cython is needed')
        USE_CYTHON = True

# try to import cython
if USE_CYTHON:
    try:
        print('importing cython')
        from Cython.Build import cythonize
        extension="pyx"
    except ImportError:
        print('Error: building ' + plumedname + ' requires cython. Please install it first with pip install cython')
        sys.exit(-1)
else:
    print('using available plumed.c file')
    extension="c"

ext_modules=[Extension(
     name=plumedname,
     sources=["plumed." + extension],
     language="c",
     include_dirs=include_dirs,
     extra_compile_args=extra_compile_args
  )]

if USE_CYTHON:
    ext_modules=cythonize(ext_modules,language_level=3)

setup(
  name=plumedname,
  version=plumedversion,
  description='Python interface to PLUMED',
  long_description=readme(),
  classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics',
          'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
  ],
  author='Gareth A. Tribello',
  author_email='plumed-users@googlegroups.com',
  url='http://www.plumed.org',
  ext_modules = ext_modules,
  zip_safe=False,
  test_suite='nose.collector',
  python_requires='>=3'
)
