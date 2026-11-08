/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "AtomDistributionFiles.h"
#include "View.h"

namespace PLMD {

void fileTraj::rewind() {
  auto errormessage=parser.rewind();
  if (errormessage) {
    // A workarounf for not implemented rewind is to dump the trajectory in an xyz and then read that
    plumed_error()<<*errormessage;
  }
  //the extra false prevents an infinite loop in case of unexpected consecutice EOFs after a rewind
  step(false);
}

//read the next step
void fileTraj::step(bool doRewind) {
  read=false;
  TrajectoryFrame<double> trajFrame;
  trajFrame.readBox=!dont_read_pbc;
  std::optional<std::string> errormessage;
  // the first time the header is read here, to size the buffers
  const bool firstFrame=masses.empty();
  if (firstFrame) {
    errormessage=parser.readHeader(trajFrame);
    if (errormessage) {
      plumed_error()<<*errormessage;
    }
    const size_t natoms = parser.nOfAtoms();

    masses.assign(natoms,0.0);
    charges.assign(natoms,0.0);
    coordinates.assign(natoms,Vector(0.0,0.0,0.0));
    cell.assign(9,0.0);
  }
  trajFrame.positions=&coordinates[0][0];
  trajFrame.box=cell.data();
  trajFrame.masses=masses.data();
  trajFrame.charges=charges.data();
  if (firstFrame) {
    errormessage=parser.readAtoms(trajFrame);
  } else {
    errormessage=parser.readFrame(trajFrame);
  }

  if (errormessage) {
    if (*errormessage =="EOF" && doRewind) {
      rewind();
    } else {
      plumed_error()<<*errormessage;
    }
  }
}

void fileTraj::frame(View<Vector> posToUpdate,
                     View<double,9> box,
                     unsigned /*step*/,
                     Random& /*rng*/) {
  if (read) {
    step();
  }
  read=true;
  std::copy(coordinates.begin(),coordinates.end(),posToUpdate.begin());
  std::copy(cell.begin(),cell.end(),box.begin());
}

fileTraj::fileTraj(std::string_view fmt,
                   std::string_view fname,
                   std::string_view format) {
  auto errormessage=parser.init(fmt,
                                fname,
                                format);
  if(errormessage) {
    plumed_error()<<*errormessage;
  }
  step();
}

bool fileTraj::overrideNat(unsigned& natoms) {
  natoms = masses.size();
  return true;
}

} //namespace PLMD
