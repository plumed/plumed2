/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#ifndef __PLUMED_tools_OptimalAlignment_h
#define __PLUMED_tools_OptimalAlignment_h

#include "Vector.h"
#include "Tensor.h"
#include <vector>

namespace PLMD{

class Log;
class Kearsley;

/// A class that is intended to include or combine various optimal alignment algorithms

class	OptimalAlignment 
{
private:
	/// a pointer to the object that performs the optimal alignment via quaternions
	Kearsley *mykearsley;
	/// displacement vector : a double that says if the coordinate should be used in calculating the RMSD/MSD
	std::vector<double> displace;
	/// alignment vector: a double that says if the atom has to be used in reset COM and makeing the alignment
	std::vector<double> align;
	/// position of one frame (generally the MD)
	std::vector<Vector> p0;
	/// position of the reference frames
	std::vector<Vector> p1;
	/// derivatives of the error  respect to the p0 (MD running frame)
	std::vector<Vector> derrdp0;
	/// derivatives of the error respect to the p1 (static frame, do not remove: useful for SM)
	std::vector<Vector> derrdp1;
	/// the pointer to the logfile
	Log* log;
	/// a bool that decides to make the fast version (alignment vec= displacement vec) or the slower case
	bool fast;

public:
	/// the contructor
	OptimalAlignment( const  std::vector<double>  & align,  const std::vector<double>   & displace, const std::vector<Vector> & p0, const std::vector<Vector> & p1 , Log* &log );
	/// the destructor: delete kearsley
	~OptimalAlignment();
	/// assignment of the running frame p0
	void assignP0(  const std::vector<Vector> & p0 );
	/// assignment to the reference frame p1
	void assignP1(  const std::vector<Vector> & p1 );
	// this updates align runtime
	void assignAlign(  const std::vector<double> & align );
	// this updates displace runtime
	void assignDisplace(  const std::vector<double> & displace );
	/// this does the real calculation
	double calculate( bool rmsd, std::vector<Vector> & derivatives);
	/// this should perform the weighted alignment
	double weightedAlignment( bool rmsd);
	// a finite difference test
	double weightedFindiffTest( bool rmsd);
};

}

#endif

