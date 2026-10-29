/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2026 The plumed team
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
#ifndef __PLUMED_bias_MetaDHills_h
#define __PLUMED_bias_MetaDHills_h

#include <string>

namespace PLMD {

class Communicator;

namespace bias {

/// shortest representation of x that is converted back to the same double
std::string exactString(double x);

/// broadcast a string from process 0: the length is sent first, as on the other processes the string
/// can have a different length (Communicator::Bcast of a string uses the local length)
void bcastString(Communicator& c, std::string& s);

/// remove the last line of a hills file if it is incomplete (left by a run that was interrupted while
/// writing it); returns the message for the log (empty if nothing was done)
std::string removeIncompleteLastLine(const std::string& path);

}
}

#endif
