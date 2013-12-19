/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#ifndef __PLUMED_tools_RMSD_h
#define __PLUMED_tools_RMSD_h

#include "Vector.h"
#include <vector>
#include <string>

namespace PLMD{

class Log;
class PDB;

/** \ingroup TOOLBOX
A class that implements RMSD calculations
This is a class that implements the various infrastructure to calculate the 
RMSD or MSD respect a given frame. It can be done through an optimal alignment scheme
as Kearsley or, more simply, by resetting the center of mass. 
This is the class that decides this. A very simple use is  
\verbatim
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "tools/Vector.h"
using namespace PLMD;
RMSD rmsd;
PDB pdb;
// get the pdb (see PDB documentation)
pdb.read("file.pdb",true,1.0);
string type;
type.assign("OPTIMAL");
// set the reference and the type 
rmsd.set(pdb,type);
// this calculates the rmsd and the derivatives
vector<Vector> derivs;
double val; 
val=rmsd.calculate(getPositions(),derivs,true);
\endverbatim

**/

class RMSD
{
  enum AlignmentMethod {SIMPLE, OPTIMAL, OPTIMAL_FAST};
  AlignmentMethod alignmentMethod;
// Reference coordinates
  std::vector<Vector> reference;
// Weights for alignment
  std::vector<double> align;
// Weights for deviation
  std::vector<double> displace;
// Logfile
  Log *log;
public:
/// Constructor
  RMSD(Log & log );
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
  		                     std::vector<Vector>  & derivatives, bool squared=false)const;
template <bool safe,bool alEqDis>
  double optimalAlignment(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,
                          const std::vector<Vector> & reference ,
                          std::vector<Vector>  & derivatives, bool squared=false)const;
/// Compute rmsd
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared=false)const;
};

}

#endif

