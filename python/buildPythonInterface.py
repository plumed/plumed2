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

# Read makefile.conf and get CC and CXX flags
f = open('../Makefile.conf', 'r')
for line in f :
   if "CC=" in line : os.environ["CC"] = line.replace("CC=","").replace("\n","")
   if "CXX=" in line : os.environ["CXX"] = line.replace("CXX=","").replace("\n","")
   if "LDSO=" in line : os.environ["LDSHARED"] = line.replace("LDSO=","").replace("\n","")
f.close()
print( "Building interface using CC=" + os.environ["CC"] + " , CXX=" + os.environ["CXX"] + " and LDSHARED=" + os.environ["LDSHARED"] )
print( "Environment variables are:")
for key in os.environ.keys():
    print( key, os.environ[key] )

plumedversion = subprocess.check_output(['../src/lib/plumed', 'info', '--version']).decode("utf-8")
print( "Creating interface for plumed version " + plumedversion )

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
                             library_dirs=["../src/lib/"],
                             libraries=["plumed"],
                             language="c++",
                             include_dirs=["../src/wrapper/", numpy.get_include()]
                           )
                          ])
)
