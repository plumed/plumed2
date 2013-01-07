/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_tools_RMSD_h
#define __PLUMED_tools_RMSD_h

#include "Vector.h"
#include <vector>
#include <string>

namespace PLMD{

class Log;
class PDB;
class OptimalAlignment;

/// \ingroup TOOLBOX
/// A class that implements RMSD calculations
class RMSD
{
  enum AlignmentMethod {SIMPLE, OPTIMAL};
  AlignmentMethod alignmentMethod;
  std::vector<Vector> reference;
  std::vector<double> align;
  std::vector<double> displace;
  OptimalAlignment *myoptimalalignment;
  Log &log;
public:
/// initialize the log in the constructor
  RMSD(Log & log ): myoptimalalignment(NULL),log(log){};
/// a copy constructor
  RMSD(const RMSD &);
/// assignment
  RMSD& operator=(const RMSD& );
/// the destructor needs to delete the myalignment object eventually
  ~RMSD();
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void set(const PDB&, std::string mytype);
/// set the type of alignment we are doing
  void setType(std::string mytype);
/// set reference coordinates
  void setReference(const std::vector<Vector> & reference);
/// set weights
  void setAlign(const std::vector<double> & align);
/// set align
  void setDisplace(const std::vector<double> & displace);
/// 
  std::string getMethod();	
///
  double simpleAlignment(const  std::vector<double>  & align,
  		                     const  std::vector<double>  & displace,
  		                     const std::vector<Vector> & positions,
  		                     const std::vector<Vector> & reference ,
  		                     Log &log,
  		                     std::vector<Vector>  & derivatives, bool squared=false);
/// Compute rmsd
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared=false);
};

}

#endif

