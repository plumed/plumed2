/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_tools_RMSD_h
#define __PLUMED_tools_RMSD_h

#include "Vector.h"
#include "Matrix.h"
#include "Tensor.h"
#include <vector>
#include <string>
#include <array>

namespace PLMD {

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
  Vector calculateCenter(std::vector<Vector> &p,std::vector<double> &w) {
    plumed_massert(p.size()==w.size(),"mismatch in dimension of position/align arrays while calculating the center");
    unsigned n; n=p.size();
    Vector c; c.zero();
    for(unsigned i=0; i<n; i++)c+=p[i]*w[i];
    return c;
  };
// removes the center for the position provided
  void removeCenter(std::vector<Vector> &p, Vector &c) {
    unsigned n; n=p.size();
    for(unsigned i=0; i<n; i++)p[i]-=c;
  };
// add center
  void addCenter(std::vector<Vector> &p, Vector &c) {Vector cc=c*-1.; removeCenter(p,cc);};

public:
/// Constructor
  RMSD();
/// clear the structure
  void clear();
/// set reference, align and displace from input pdb structure: evtl remove com from the initial structure and normalize the input weights from the pdb
  void set(const PDB&,const std::string & mytype, bool remove_center=true, bool normalize_weights=true);
/// set align displace reference and type from input vectors
  void set(const std::vector<double> & align, const std::vector<double> & displace, const std::vector<Vector> & reference,const std::string & mytype, bool remove_center=true, bool normalize_weights=true );
/// set the type of alignment we are doing
  void setType(const std::string & mytype);
/// set reference coordinates, remove the com by using uniform weights
  void setReference(const std::vector<Vector> & reference);
  std::vector<Vector> getReference();
/// set weights and remove the center from reference with normalized weights. If the com has been removed, it resets to the new value
  void setAlign(const std::vector<double> & align, bool normalize_weights=true, bool remove_center=true);
  std::vector<double> getAlign();
/// set align
  void setDisplace(const std::vector<double> & displace, bool normalize_weights=true);
  std::vector<double> getDisplace();
///
  std::string getMethod();
/// workhorses
  double simpleAlignment(const  std::vector<double>  & align,
                         const  std::vector<double>  & displace,
                         const std::vector<Vector> & positions,
                         const std::vector<Vector> & reference,
                         std::vector<Vector>  & derivatives,
                         std::vector<Vector>  & displacement,
                         bool squared=false)const;
  template <bool safe,bool alEqDis>
  double optimalAlignment(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,
                          const std::vector<Vector> & reference,
                          std::vector<Vector>  & DDistDPos, bool squared=false)const;

  template <bool safe, bool alEqDis>
  double optimalAlignmentWithCloseStructure(const  std::vector<double>  & align,
      const  std::vector<double>  & displace,
      const std::vector<Vector> & positions,
      const std::vector<Vector> & reference,
      std::vector<Vector>  & derivatives,
      Tensor & rotationPosClose,
      Tensor & rotationRefClose,
      std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01,
      bool squared=false)const;

  template <bool safe, bool alEqDis>
  double optimalAlignment_Rot_DRotDRr01(const  std::vector<double>  & align,
                                        const  std::vector<double>  & displace,
                                        const std::vector<Vector> & positions,
                                        const std::vector<Vector> & reference,
                                        Tensor & Rotation,
                                        std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01,
                                        bool squared=false)const;

  template <bool safe, bool alEqDis>
  double optimalAlignment_Rot(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              const std::vector<Vector> & reference,
                              std::vector<Vector>  & derivatives,
                              Tensor & Rotation,
                              bool squared=false)const;

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef(const  std::vector<double>  & align,
                                    const  std::vector<double>  & displace,
                                    const std::vector<Vector> & positions,
                                    const std::vector<Vector> & reference,
                                    std::vector<Vector>  & DDistDPos,
                                    std::vector<Vector> &  DDistDRef,
                                    bool squared=false) const;

  template <bool safe,bool alEqDis>
  double optimalAlignment_SOMA(const  std::vector<double>  & align,
                               const  std::vector<double>  & displace,
                               const std::vector<Vector> & positions,
                               const std::vector<Vector> & reference,
                               std::vector<Vector>  & DDistDPos,
                               std::vector<Vector> &  DDistDRef,
                               bool squared=false) const;

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef_Rot_DRotDPos(const  std::vector<double>  & align,
      const  std::vector<double>  & displace,
      const std::vector<Vector> & positions,
      const std::vector<Vector> & reference,
      std::vector<Vector>  & DDistDPos,
      std::vector<Vector> &  DDistDRef,
      Tensor & Rotation,
      Matrix<std::vector<Vector> > &DRotDPos,
      bool squared=false) const;

  template <bool safe,bool alEqDis>
  double optimalAlignment_DDistDRef_Rot_DRotDPos_DRotDRef(const  std::vector<double>  & align,
      const  std::vector<double>  & displace,
      const std::vector<Vector> & positions,
      const std::vector<Vector> & reference,
      std::vector<Vector>  & DDistDPos,
      std::vector<Vector> &  DDistDRef,
      Tensor & Rotation,
      Matrix<std::vector<Vector> > &DRotDPos,
      Matrix<std::vector<Vector> > &DRotDRef,
      bool squared=false) const;

  template <bool safe,bool alEqDis>
  double optimalAlignment_PCA(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              const std::vector<Vector> & reference,
                              std::vector<Vector> & alignedpositions,
                              std::vector<Vector> & centeredpositions,
                              std::vector<Vector> & centeredreference,
                              Tensor & Rotation,
                              std::vector<Vector> & DDistDPos,
                              Matrix<std::vector<Vector> > & DRotDPos,
                              bool squared=false) const ;

  template <bool safe,bool alEqDis>
  double optimalAlignment_Fit(const  std::vector<double>  & align,
                              const  std::vector<double>  & displace,
                              const std::vector<Vector> & positions,
                              const std::vector<Vector> & reference,
                              Tensor & Rotation,
                              Matrix<std::vector<Vector> > & DRotDPos,
                              std::vector<Vector> & centeredpositions,
                              Vector & center_positions,
                              bool squared=false);


/// Compute rmsd: note that this is an intermediate layer which is kept in order to evtl expand with more alignment types/user options to be called while keeping the workhorses separated
  double calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, bool squared=false)const;
/// Other convenience methods:
/// calculate the derivative of distance respect to position(DDistDPos) and reference (DDistDPos)
  double calc_DDistDRef( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, std::vector<Vector>& DDistDRef, const bool squared=false   );
/// calculate the derivative for SOMA (i.e. derivative respect to reference frame without the matrix derivative)
  double calc_SOMA( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, std::vector<Vector>& DDistDRef, const bool squared=false   );
///
  double calc_DDistDRef_Rot_DRotDPos( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, std::vector<Vector>& DDistDRef, Tensor & Rotation,Matrix<std::vector<Vector> > &DRotDPos, const bool squared=false   );
  double calc_DDistDRef_Rot_DRotDPos_DRotDRef( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, std::vector<Vector>& DDistDRef, Tensor & Rotation,Matrix<std::vector<Vector> > &DRotDPos,Matrix<std::vector<Vector> > &DRotDRef, const bool squared=false   );
/// convenience method to retrieve all the bits and pieces for PCA
  double calc_PCAelements( const std::vector<Vector>& pos, std::vector<Vector> &DDistDPos, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,std::vector<Vector>  & alignedpositions, std::vector<Vector> & centeredpositions, std::vector<Vector> &centeredreference, const bool& squared=false) const ;
/// convenience method to retrieve all the bits and pieces needed by FitToTemplate
  double calc_FitElements( const std::vector<Vector>& pos, Tensor & Rotation, Matrix<std::vector<Vector> > & DRotDPos,std::vector<Vector> & centeredpositions,Vector & center_positions, const bool& squared=false );
///calculate rotation matrix, derivative of rotation matrix w.r.t. positions, derivative of rotation matrix w.r.t. rr01
  double calc_Rot_DRotDRr01( const std::vector<Vector>& positions, Tensor & Rotation, std::array<std::array<Tensor,3>,3> & DRotDRr01, const bool squared=false   );
///calculate rotation matrix, derivative of rotation matrix w.r.t. positions
  double calc_Rot( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, Tensor & Rotation, const bool squared=false   );
///calculate with close structure, i.e. approximate the RMSD without expensive computation of rotation matrix by reusing saved rotation matrices from previous iterations
  double calculateWithCloseStructure( const std::vector<Vector>& positions, std::vector<Vector> &DDistDPos, Tensor & rotationPosClose, Tensor & rotationRefClose, std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01, const bool squared=false   );
/// static convenience method to get the matrix i,a from drotdpos (which is a bit tricky)
  static  Tensor getMatrixFromDRot(Matrix< std::vector<Vector> > &drotdpos, const unsigned &i, const unsigned &a) {
    Tensor t;
    t[0][0]=drotdpos[0][0][i][a]; t[0][1]=drotdpos[0][1][i][a]; t[0][2]=drotdpos[0][2][i][a];
    t[1][0]=drotdpos[1][0][i][a]; t[1][1]=drotdpos[1][1][i][a]; t[1][2]=drotdpos[1][2][i][a];
    t[2][0]=drotdpos[2][0][i][a]; t[2][1]=drotdpos[2][1][i][a]; t[2][2]=drotdpos[2][2][i][a];
    return t;
  };
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

  // this need to be copied and they are small, should not affect the performance
  Vector creference;
  bool creference_is_calculated;
  bool creference_is_removed;
  Vector cpositions;
  bool cpositions_is_calculated;
  bool cpositions_is_removed;
  bool retrieve_only_rotation;

  // use reference assignment to speed up instead of copying
  const std::vector<Vector> &positions;
  const std::vector<Vector> &reference;
  const std::vector<double> &align;
  const std::vector<double> &displace;

  // the needed stuff for distance and more (one could use eigenvecs components and eigenvals for some reason)
  double dist;
  Vector4d eigenvals;
  Tensor4d eigenvecs;
  double rr00; //  sum of positions squared (needed for dist calc)
  double rr11; //  sum of reference squared (needed for dist calc)
  Tensor rotation; // rotation derived from the eigenvector having the smallest eigenvalue
  std::array<std::array<Tensor,3>,3> drotation_drr01; // derivative of the rotation only available when align!=displace
  Tensor ddist_drr01;
  Tensor ddist_drotation;
  std::vector<Vector> d; // difference of components
public:
  /// the constructor (note: only references are passed, therefore is rather fast)
  /// note: this aligns the reference onto the positions
  ///
  /// this method assumes that the centers are already calculated and subtracted
  RMSDCoreData(const std::vector<double> &a,const std::vector<double> &d,const std::vector<Vector> &p, const std::vector<Vector> &r, Vector &cp, Vector &cr ):
    alEqDis(false),distanceIsMSD(false),hasDistance(false),isInitialized(false),safe(false),
    creference(cr),creference_is_calculated(true),creference_is_removed(true),
    cpositions(cp),cpositions_is_calculated(true),cpositions_is_removed(true),retrieve_only_rotation(false),positions(p),reference(r),align(a),displace(d),dist(0.0),rr00(0.0),rr11(0.0) {};

  // this constructor does not assume that the positions and reference have the center subtracted
  RMSDCoreData(const std::vector<double> &a,const std::vector<double> &d,const std::vector<Vector> &p, const std::vector<Vector> &r):
    alEqDis(false),distanceIsMSD(false),hasDistance(false),isInitialized(false),safe(false),
    creference_is_calculated(false),creference_is_removed(false),
    cpositions_is_calculated(false),cpositions_is_removed(false),retrieve_only_rotation(false),positions(p),reference(r),align(a),displace(d),dist(0.0),rr00(0.0),rr11(0.0)
  {cpositions.zero(); creference.zero();};

  // set the center on the fly without subtracting
  void calcPositionsCenter() {
    plumed_massert(!cpositions_is_calculated,"the center was already calculated");
    cpositions.zero(); for(unsigned i=0; i<positions.size(); i++) {cpositions+=positions[i]*align[i];} cpositions_is_calculated=true;
  }
  void calcReferenceCenter() {
    plumed_massert(!creference_is_calculated,"the center was already calculated");
    creference.zero(); for(unsigned i=0; i<reference.size(); i++) {creference+=reference[i]*align[i];} creference_is_calculated=true;
  };
  // assume the center is given externally
  void setPositionsCenter(Vector v) {plumed_massert(!cpositions_is_calculated,"You are setting the center two times!"); cpositions=v; cpositions_is_calculated=true;};
  void setReferenceCenter(Vector v) {plumed_massert(!creference_is_calculated,"You are setting the center two times!"); creference=v; creference_is_calculated=true;};
  // the center is already removed
  void setPositionsCenterIsRemoved(bool t) {cpositions_is_removed=t;};
  void setReferenceCenterIsRemoved(bool t) {creference_is_removed=t;};
  bool getPositionsCenterIsRemoved() {return cpositions_is_removed;};
  bool getReferenceCenterIsRemoved() {return creference_is_removed;};
  //  does the core calc : first thing to call after the constructor:
  // only_rotation=true does not retrieve the derivatives, just retrieve the optimal rotation (the same calc cannot be exploit further)
  void doCoreCalc(bool safe,bool alEqDis, bool only_rotation=false);
  // do calculation with close structure data structures
  void doCoreCalcWithCloseStructure(bool safe,bool alEqDis, Tensor & rotationPosClose, Tensor & rotationRefClose, std::array<std::array<Tensor,3>,3> & drotationPosCloseDrr01);
  // retrieve the distance if required after doCoreCalc
  double getDistance(bool squared);
  // retrieve the derivative of the distance respect to the position
  std::vector<Vector> getDDistanceDPositions();
  // retrieve the derivative of the distance respect to the reference
  std::vector<Vector> getDDistanceDReference();
  // specific version for SOMA calculation (i.e. does not need derivative respect to rotation matrix)
  std::vector<Vector> getDDistanceDReferenceSOMA();
  // get aligned reference onto position
  std::vector<Vector> getAlignedReferenceToPositions();
  // get aligned position onto reference
  std::vector<Vector> getAlignedPositionsToReference();
  // get centered positions
  std::vector<Vector> getCenteredPositions();
  // get centered reference
  std::vector<Vector> getCenteredReference();
  // get center of positions
  Vector getPositionsCenter();
  // get center of reference
  Vector getReferenceCenter();
  // get rotation matrix (reference ->positions)
  Tensor getRotationMatrixReferenceToPositions();
  // get rotation matrix (positions -> reference)
  Tensor getRotationMatrixPositionsToReference();
  // get the derivative of the rotation matrix respect to positions
  // note that the this transformation overlap the  reference onto position
  // if inverseTransform=true then aligns the positions onto reference
  Matrix<std::vector<Vector> > getDRotationDPositions( bool inverseTransform=false );
  // get the derivative of the rotation matrix respect to reference
  // note that the this transformation overlap the  reference onto position
  // if inverseTransform=true then aligns the positions onto reference
  Matrix<std::vector<Vector> >  getDRotationDReference(bool inverseTransform=false );
  const std::array<std::array<Tensor,3>,3> & getDRotationDRr01() const;
};

}

#endif

