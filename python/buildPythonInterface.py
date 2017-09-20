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
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy
import subprocess
import os

# Function for checking if PLUMED is in path
def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

# Check if PLUMED is in PATH
plumedexe = '../src/lib/plumed'
for path in os.environ["PATH"].split(os.pathsep):
    path = path.strip('"')
    exe_file = os.path.join(path, 'plumed') 
    if is_exe(exe_file) : 
       plumedexe=exe_file
       break

# Get information on where plumed headers and libraries are installed and the version number
print( "Plumedexe is " + plumedexe )
plumedroot = subprocess.check_output([plumedexe, 'info', '--root']).decode("utf-8").strip("\n")
print( "Creating interface for plumed version in " + plumedroot )
plumedhead = subprocess.check_output([plumedexe, 'info', '--include-dir']).decode("utf-8").strip("\n") + "/plumed/wrapper/"
print( "Using headers in " + plumedhead )
plumedversion = subprocess.check_output([plumedexe, 'info', '--version']).decode("utf-8")
print( "Version number for this plumed is " + plumedversion )
# Get list containing all config variables so we can extract information on compilers to use during cython build
plumedconfig = subprocess.check_output([plumedexe, 'info', '--configuration']).decode("utf-8").split("\n")

for line in plumedconfig :
   if "CC=" in line : os.environ["CC"] = line.replace("CC=","").replace("\n","")
   if "CXX=" in line : os.environ["CXX"] = line.replace("CXX=","").replace("\n","")
   if "LDSHARED=" in line : os.environ["LDSHARED"] = line.replace("LDSHARED=","").replace("\n","")

print( "Building interface using CC=" + os.environ["CC"] + " , CXX=" + os.environ["CXX"] + " and LDSHARED=" + os.environ["LDSHARED"] )

setup(
  name='plumed',
  version=plumedversion,
  description='Python interface to PLUMED',
  author='Gareth A. Tribello',
  author_email='plumed-users@googlegroups.com',
  url='http://www.plumed.org',
  ext_modules = cythonize([
                  Extension( name="plumed", 
                             sources=["plumed.pyx"],
                             library_dirs=[plumedroot, plumedroot + "/src/lib/"],
                             libraries=["plumed"],
                             language="c++",
                             include_dirs=[plumedhead, numpy.get_include()]
                           )
                          ])
)
