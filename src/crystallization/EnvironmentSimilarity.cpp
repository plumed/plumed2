/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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

/* ----------------------------------------------------------------------
   Contributing author: Pablo Piaggi (Princeton University)
------------------------------------------------------------------------- */

#include "multicolvar/MultiColvarBase.h"
#include "multicolvar/AtomValuePack.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

#include <string>
#include <cmath>
#include <cfloat>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR ENVIRONMENTSIMILARITY
/*

*/
//+ENDPLUMEDOC


class EnvironmentSimilarity : public MultiColvarBase {
private:
  // All global variables end with underscore
  // square of cutoff, square of broadening parameter
  double rcut2_, sigmaSqr_;
  // lambda parameter for softmax function
  double lambda_;
  // Array of Vectors to store the reference environments, i.e. the templates
  std::vector<std::vector<Vector>> environments_;
public:
  static void registerKeywords( Keywords& keys );
  explicit EnvironmentSimilarity(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
// Calculates maximum distance in an environment
  double maxDistance(std::vector<Vector> environment);
  // Parse everything connected to the definition of the reference environments
  // First argument is the array of Vectors that stores the reference environments
  // Second argument is the maximum distance in the ref environments and sets the
  // cutoff for the cell lists
  void parseReferenceEnvironments( std::vector<std::vector<Vector>>& environments, double& max_dist);
};

PLUMED_REGISTER_ACTION(EnvironmentSimilarity,"ENVIRONMENTSIMILARITY")

void EnvironmentSimilarity::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","SIGMA","0.1","Broadening parameter");
  keys.add("compulsory","CRYSTAL_STRUCTURE","FCC","Targeted crystal structure");
  keys.add("optional","LATTICE_CONSTANTS","Lattice constants");
  keys.add("compulsory","LAMBDA","100","Lambda parameter");
  keys.add("optional","REFERENCE","PDB file with relative distances from central atom.");
  keys.add("numbered","REFERENCE_","PDB files with relative distances from central atom."
		                   "Each file corresponds to one template");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

EnvironmentSimilarity::EnvironmentSimilarity(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Parse everything connected to the definition of the reference environments
  double max_dist_ref_vector;
  parseReferenceEnvironments(environments_, max_dist_ref_vector);

  double sigma;
  parse("SIGMA", sigma);
  log.printf("  representing local density as a sum of Gaussians with standard deviation %f\n",sigma);
  sigmaSqr_=sigma*sigma;

  lambda_=100;
  parse("LAMBDA", lambda_);
  if (environments_.size()>1) log.printf("  using a soft max function with lambda %f\n",lambda_);

  // Set the link cell cutoff
  double rcut = max_dist_ref_vector + 3*sigma;
  setLinkCellCutoff( rcut );
  rcut2_ = rcut * rcut;

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
}

double EnvironmentSimilarity::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  if (environments_.size()==1) {
    // One reference environment case
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over atoms in the reference environment
        for(unsigned k=0; k<environments_[0].size(); ++k) {
          Vector distanceFromRef=distance-environments_[0][k];
          double value = std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/environments_[0].size() ;
          // CAREFUL! Off-diagonal virial is incorrect. Do not perform NPT simulations with flexible box angles.
          accumulateSymmetryFunction( 1, i, value, (value/(2*sigmaSqr_))*(-distance+environments_[0][k]) , (value/(2*sigmaSqr_))*Tensor(distance-environments_[0][k],distance) , myatoms );
        }
      }
    }
    return myatoms.getValue(1);
  } else {
    // More than one reference environment case
    std::vector<double> values(environments_.size()); //value for each template
    // First time calculate sums
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over templates
        for(unsigned j=0; j<environments_.size(); ++j) {
          // Iterate over atoms in the template
          for(unsigned k=0; k<environments_[j].size(); ++k) {
            Vector distanceFromRef=distance-environments_[j][k];
            values[j] += std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/environments_[j].size() ;
          }
        }
      }
    }
    double sum=0;
    for(unsigned j=0; j<environments_.size(); ++j) {
       values[j] = std::exp(lambda_*values[j]);
       sum += values[j];
    }
    // Second time find derivatives
    for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
      Vector& distance=myatoms.getPosition(i);
      double d2;
      if ( (d2=distance[0]*distance[0])<rcut2_ &&
           (d2+=distance[1]*distance[1])<rcut2_ &&
           (d2+=distance[2]*distance[2])<rcut2_ &&
           d2>epsilon ) {
        // Iterate over reference environment
        for(unsigned j=0; j<environments_.size(); ++j) {
          // Iterate over atoms in the reference environment
          for(unsigned k=0; k<environments_[j].size(); ++k) {
            Vector distanceFromRef=distance-environments_[j][k];
            double value = std::exp(-distanceFromRef.modulo2()/(4*sigmaSqr_) )/environments_[j].size() ;
            accumulateSymmetryFunction( 1, i, value, -(values[j]/sum)*(value/(2*sigmaSqr_))*distanceFromRef  , (values[j]/sum)*(value/(2*sigmaSqr_))*Tensor(distanceFromRef,distance) , myatoms );
          }
        }
      }
    }
    return std::log(sum)/lambda_;
  }
}

double EnvironmentSimilarity::maxDistance( std::vector<Vector> environment ) {
  double max_dist = 0.0;
  for(unsigned i=0;i<environment.size(); ++i) {
    double norm=environment[i].modulo();
    if (norm>max_dist) max_dist=norm;
  }
  return max_dist;
}

void EnvironmentSimilarity::parseReferenceEnvironments( std::vector<std::vector<Vector>>& environments, double& max_dist) {
  std::vector<double> lattice_constants;
  parseVector("LATTICE_CONSTANTS", lattice_constants);
  std::string crystal_structure;
  parse("CRYSTAL_STRUCTURE", crystal_structure);
  // find crystal structure
  if (crystal_structure == "FCC") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for FCC");
    environments.resize(1);
    environments[0].resize(12);
    environments[0][0]  = Vector(+0.5,+0.5,+0.0)*lattice_constants[0];
    environments[0][1]  = Vector(-0.5,-0.5,+0.0)*lattice_constants[0];
    environments[0][2]  = Vector(+0.5,-0.5,+0.0)*lattice_constants[0];
    environments[0][3]  = Vector(-0.5,+0.5,+0.0)*lattice_constants[0];
    environments[0][4]  = Vector(+0.5,+0.0,+0.5)*lattice_constants[0];
    environments[0][5]  = Vector(-0.5,+0.0,-0.5)*lattice_constants[0];
    environments[0][6]  = Vector(-0.5,+0.0,+0.5)*lattice_constants[0];
    environments[0][7]  = Vector(+0.5,+0.0,-0.5)*lattice_constants[0];
    environments[0][8]  = Vector(+0.0,+0.5,+0.5)*lattice_constants[0];
    environments[0][9]  = Vector(+0.0,-0.5,-0.5)*lattice_constants[0];
    environments[0][10] = Vector(+0.0,-0.5,+0.5)*lattice_constants[0];
    environments[0][11] = Vector(+0.0,+0.5,-0.5)*lattice_constants[0];
    max_dist = std::sqrt(2)*lattice_constants[0]/2.;
  } else if (crystal_structure == "SC") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for SC");
    environments.resize(1);
    environments[0].resize(6);
    environments[0][0]  = Vector(+1.0,+0.0,+0.0)*lattice_constants[0];
    environments[0][1]  = Vector(-1.0,+0.0,+0.0)*lattice_constants[0];
    environments[0][2]  = Vector(+0.0,+1.0,+0.0)*lattice_constants[0];
    environments[0][3]  = Vector(+0.0,-1.0,+0.0)*lattice_constants[0];
    environments[0][4]  = Vector(+0.0,+0.0,+1.0)*lattice_constants[0];
    environments[0][5]  = Vector(+0.0,+0.0,-1.0)*lattice_constants[0];
    max_dist = lattice_constants[0];
  } else if (crystal_structure == "BCC") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for BCC");
    environments.resize(1);
    environments[0].resize(14);
    environments[0][0]  = Vector(+0.5,+0.5,+0.5)*lattice_constants[0];
    environments[0][1]  = Vector(-0.5,-0.5,-0.5)*lattice_constants[0];
    environments[0][2]  = Vector(-0.5,+0.5,+0.5)*lattice_constants[0];
    environments[0][3]  = Vector(+0.5,-0.5,+0.5)*lattice_constants[0];
    environments[0][4]  = Vector(+0.5,+0.5,-0.5)*lattice_constants[0];
    environments[0][5]  = Vector(-0.5,-0.5,+0.5)*lattice_constants[0];
    environments[0][6]  = Vector(+0.5,-0.5,-0.5)*lattice_constants[0];
    environments[0][7]  = Vector(-0.5,+0.5,-0.5)*lattice_constants[0];
    environments[0][8]  = Vector(+1.0,+0.0,+0.0)*lattice_constants[0];
    environments[0][9]  = Vector(+0.0,+1.0,+0.0)*lattice_constants[0];
    environments[0][10] = Vector(+0.0,+0.0,+1.0)*lattice_constants[0];
    environments[0][11] = Vector(-1.0,+0.0,+0.0)*lattice_constants[0];
    environments[0][12] = Vector(+0.0,-1.0,+0.0)*lattice_constants[0];
    environments[0][13] = Vector(+0.0,+0.0,-1.0)*lattice_constants[0];
    max_dist = lattice_constants[0];
  } else if (crystal_structure == "HCP") {
    if (lattice_constants.size() != 2) error("Number of LATTICE_CONSTANTS arguments must be two for HCP");
    environments.resize(2);
    environments[0].resize(12);
    environments[1].resize(12);
    double sqrt3=std::sqrt(3);
    environments[0][0]  = Vector(+0.5,+sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[0][1]  = Vector(-0.5,+sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[0][2]  = Vector(+0.5,-sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[0][3]  = Vector(-0.5,-sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[0][4]  = Vector(+1.0,+0.0,+0.0)      *lattice_constants[0];
    environments[0][5]  = Vector(-1.0,+0.0,+0.0)      *lattice_constants[0];
    environments[0][6]  = Vector(+0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[0][7]  = Vector(-0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[0][8]  = Vector(+0.0,-sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[0][9]  = Vector(+0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    environments[0][10] = Vector(-0.5,+sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    environments[0][11] = Vector(+0.0,-sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    environments[1][0]  = Vector(+0.5,+sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[1][1]  = Vector(-0.5,+sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[1][2]  = Vector(+0.5,-sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[1][3]  = Vector(-0.5,-sqrt3/2.0,+0.0)*lattice_constants[0];
    environments[1][4]  = Vector(+1.0,+0.0,+0.0)      *lattice_constants[0];
    environments[1][5]  = Vector(-1.0,+0.0,+0.0)      *lattice_constants[0];
    environments[1][6]  = Vector(+0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[1][7]  = Vector(-0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[1][8]  = Vector(+0.0,+sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,+0.5)*lattice_constants[1];
    environments[1][9]  = Vector(+0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    environments[1][10] = Vector(-0.5,-sqrt3/6.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    environments[1][11] = Vector(+0.0,+sqrt3/3.0,+0.0)*lattice_constants[0] + Vector(+0.0,+0.0,-0.5)*lattice_constants[1];
    max_dist = lattice_constants[0];
 } else if (crystal_structure == "DIAMOND") {
    if (lattice_constants.size() != 1) error("Number of LATTICE_CONSTANTS arguments must be one for DIAMOND");
    environments.resize(2);
    environments[0].resize(4); environments[1].resize(4);
    environments[0][0]  = Vector(+1.0,+1.0,+1.0)*lattice_constants[0]/4.0;
    environments[0][1]  = Vector(-1.0,-1.0,+1.0)*lattice_constants[0]/4.0;
    environments[0][2]  = Vector(+1.0,-1.0,-1.0)*lattice_constants[0]/4.0;
    environments[0][3]  = Vector(-1.0,+1.0,-1.0)*lattice_constants[0]/4.0;
    environments[1][0]  = Vector(+1.0,-1.0,+1.0)*lattice_constants[0]/4.0;
    environments[1][1]  = Vector(-1.0,+1.0,+1.0)*lattice_constants[0]/4.0;
    environments[1][2]  = Vector(+1.0,+1.0,-1.0)*lattice_constants[0]/4.0;
    environments[1][3]  = Vector(-1.0,-1.0,-1.0)*lattice_constants[0]/4.0;
    max_dist = std::sqrt(3)*lattice_constants[0]/4.0;
 } else if (crystal_structure == "CUSTOM") {
    std::string reffile;
    parse("REFERENCE",reffile);
    if (!reffile.empty()) {
      // Case with one reference environment
      environments.resize(1);
      PDB pdb; pdb.read(reffile,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength());
      unsigned natoms=pdb.getPositions().size(); environments[0].resize( natoms );
      for(unsigned i=0;i<natoms;++i) environments[0][i]=pdb.getPositions()[i];
      max_dist=maxDistance(environments[0]);
      log.printf("  reading %d reference vectors from %s \n", natoms, reffile.c_str() );
    } else {
      // Case with several reference environments
      max_dist=0;
      for(unsigned int i=1;; i++) {
        if(!parseNumbered("REFERENCE_",i,reffile) ) {break;}
        PDB pdb; pdb.read(reffile,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength());
        unsigned natoms=pdb.getPositions().size();   std::vector<Vector> environment; environment.resize( natoms );
        for(unsigned i=0;i<natoms;++i) environment[i]=pdb.getPositions()[i];
        environments.push_back(environment);
        double norm = maxDistance(environment);
        if (norm>max_dist) max_dist=norm;
        log.printf("  Reference environment %d : reading %d reference vectors from %s \n", i, natoms, reffile.c_str() );
      }
    }
    if (environments.size()==0) error("No environments have been found! Please specify a PDB file in the REFERENCE "
                                      "or in the REFERENCE_1, REFERENCE_2, etc keywords");
    log.printf("  Number of reference environments is %d\n",environments.size() );
    log.printf("  Number of vectors per reference environment is %d\n",environments[0].size() );
  } else {
    error("CRYSTAL_STRUCTURE=" + crystal_structure + " does not match any structures in the database");
  }

  log.printf("  targeting the %s crystal structure",crystal_structure.c_str());
  if (lattice_constants.size()>0) log.printf(" with lattice constants %f\n",lattice_constants[0]);
  else log.printf("\n");

  log.printf("  maximum distance in the reference environment is %f\n",max_dist);
}

}
}
