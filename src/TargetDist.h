/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_TargetDist_h
#define __PLUMED_TargetDist_h

#include "Value.h"
#include "ActionWithValue.h"
#include "PDB.h"
#include <vector>
#include <string>

namespace PLMD{

class Log;
class PDB;

class TargetDist {
private:
  std::vector<Value*> args; 
  std::vector<double> target;
  Log &log;
public:
  TargetDist(Log& log) : log(log) {};
  void read( const PDB& pdb, std::vector<Value*> args ); 
  void read( const std::vector<double>& targ, std::vector<Value*> ar );
  double calculate( std::vector<double>& derivs );
};

}

#endif
