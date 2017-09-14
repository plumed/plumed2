Welcome to the plumed2 wiki!

This fork of PLUMED has eABF/DRR implementation.

**Requirements**

Boost::serialization and C++11 compiler

**Compiling instruction:**

After clone this repository, please cd to the plumed2 directory and run:

1. `autoconf`
2. `./configure --enable-boost_serialization --enable-modules=drr`
3. Modify your Makefile.conf and add `-lboost_serialization` in `DYNAMIC_LIBS=`
4. `make` and `sudo make install`

**Usage**

Run `make doc` and the usage is in `user-doc/html/_d_r_r.html`

**Authors**

Chen Haochuan (All files in drr module except colvar_UIestimator.h)

Fu Haohao (colvar_UIestimator.h)

**COPYRIGHT**

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


