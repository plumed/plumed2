/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "Matrix.h"
#include "Tensor.h"
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
// Center for reference and flag for its calculation 
  Vector reference_center;
  bool reference_center_is_calculated;
  bool reference_center_is_removed; 
// Center for running position (not used in principle but here to reflect reference/positio symmetry
  Vector positions_center;
  bool positions_center_is_calculated; 
  bool positions_center_is_removed; 
// calculates the center from the position provided
  Vector calculateCenter(std::vector<Vector> &p,std::vector<double> &w){
	plumed_massert(p.size()==w.size(),"mismatch in dimension of position/align arrays while calculating the center");
	unsigned n; n=p.size();
	Vector c; c.zero();
	for(unsigned i=0;i<n;i++)c+=p[i]*w[i];
	return c;
  };
// removes the center for the position provided
  void removeCenter(std::vector<Vector> &p, Vector &c){
	unsigned n; n=p.size();
	for(unsigned i=0;i<n;i++)p[i]-=c;
  };

public:
/// Constructor
  RMSD();
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure
  void set(const PDB&, std::string mytype);
/// set the type of alignment we are doing
  void setType(std::string mytype);
/// set reference coordinates
  void setReference(const std::vector<Vector> & reference);
/// set weights
  void setAlign(const std::vector<double> & align, bool normalize=true, bool remove_center=true);
/// set align
  void setDisplace(const std::vector<double> & displace, bool normalize=true);
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

/// this is a class which is needed to share information across the various non-threadsafe routines
/// so that the public function of rmsd are threadsafe while the inner core can safely share information 
class RMSDCoreData
{
	private:
		bool alEqDis;
		bool distanceIsMSD; // default is RMSD but can deliver the MSD 
		bool hasDistance;  // distance is already calculated
		bool isInitialized;
		bool safe;

		// use reference assignment to speed up instead of copying
                const std::vector<Vector> &positions;
                const std::vector<Vector> &reference;
                const std::vector<double> &align; 
                const std::vector<double> &displace; 

		// the needed stuff for distance and more (one could use eigenvecs components and eigenvals for some reason)
		double dist;
		std::vector<double> eigenvals;
		Matrix<double> eigenvecs;
		double rr00; //  sum of positions squared (needed for dist calc)
		double rr11; //  sum of reference squared (needed for dist calc)
		Tensor rotation; // rotation derived from the eigenvector having the smallest eigenvalue
		Tensor drotation_drr01[3][3]; // derivative of the rotation only available when align!=displace
		Tensor ddist_drr01;
	        Tensor ddist_drotation;
		std::vector<Vector> d; // difference of components
 		Vector cpositions,creference; // geometric center of the running position and reference
	public:
		// the constructor (note: only references are passed, therefore is rather fast)
		// note:: this aligns the reference onto the positions
		RMSDCoreData(const std::vector<double> &a ,const std::vector<double> &d,const std::vector<Vector> &p, const std::vector<Vector> &r ):alEqDis(false),distanceIsMSD(false),hasDistance(false),isInitialized(false),safe(false),positions(p),reference(r),align(a),displace(d){};
		//  does the core calc : first thing to call after the constructor	
		void doCoreCalc(bool safe,bool alEqDis);
		// retrieve the distance if required after doCoreCalc 
		double getDistance(bool squared);
		// retrieve the derivative of the distance respect to the position
		std::vector<Vector> getDDistanceDPositions();
		// retrieve the derivative of the distance respect to the reference
		std::vector<Vector> getDDistanceDReference();
		// get aligned reference onto position
                std::vector<Vector> getAlignedReferenceToPositions();	
		// get aligned position onto reference
                std::vector<Vector> getAlignedPositionsToReference();	
		// get centered positions
                std::vector<Vector> getCenteredPositions();	
		// get centered reference
                std::vector<Vector> getCenteredReference();	
		// get rotation matrix (reference ->positions) 
		Tensor getRotationMatrixReferenceToPositions();
		// get rotation matrix (positions -> reference) 
		Tensor getRotationMatrixPositionsToReference();
		// get the derivative of the rotation matrix respect to positions
		// note that the this transformation overlap the  reference onto position
		// if inverseTransform=true then aligns the positions onto reference
		Matrix<std::vector<Vector> > getDRotationDPosition( bool inverseTransform=false );
		// get the derivative of the rotation matrix respect to reference 
		// note that the this transformation overlap the  reference onto position
		// if inverseTransform=true then aligns the positions onto reference
		Matrix<std::vector<Vector> >  getDRotationDReference(bool inverseTransform=false );
};

}

#endif

